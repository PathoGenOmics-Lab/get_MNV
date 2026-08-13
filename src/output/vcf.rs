//! VCF/BCF output writer: generates VCF (plain or BGZF-compressed) with
//! INFO fields, FILTER tags, Tabix indexing, and BCF conversion.
//!
//! The metric/filter helpers, INFO-entry builders, and per-record writers live
//! in submodules as additional `impl VcfWriter` blocks; this file owns the
//! writer types, the lifecycle methods, and the standalone index/BCF helpers.

use crate::error::AppResult;
use crate::io::ReferenceMap;
use crate::variants::{VariantInfo, VariantType};
use std::fs::File;
use std::io::Write;
use std::path::Path;

use super::common::*;

mod entry;
mod metrics;
mod records;

#[cfg(test)]
mod tests;

/// Configuration for constructing a VCF writer.
pub struct VcfWriterConfig<'a> {
    pub filename: &'a str,
    pub bam_provided: bool,
    pub min_snp_reads: usize,
    pub min_snp_frequency: f64,
    pub min_mnv_reads: usize,
    pub min_mnv_frequency: f64,
    pub min_quality: u8,
    pub min_mapq: u8,
    pub command_line: &'a str,
    pub contigs: &'a [String],
    pub bgzf_output: bool,
    pub min_snp_strand_reads: usize,
    pub min_mnv_strand_reads: usize,
    pub min_strand_bias_p: f64,
    pub emit_filtered: bool,
    pub include_strand_bias_info: bool,
    pub original_info_headers: &'a [String],
}

pub struct VcfWriter {
    writer: Box<dyn Write>,
    /// Records for the contig being written, held until they can go out in
    /// coordinate order.
    ///
    /// A multi-position variant emits one record per constituent position, so a
    /// variant sitting between those positions was written after them: the file
    /// came out unsorted, `tabix` refused it, and because indexing runs after
    /// the VCF is already on disk the run then died before writing the BCF, the
    /// summary JSON, the manifest or the report. Buffering one contig costs the
    /// same order of memory as the variants it came from, which the caller is
    /// already holding.
    pending: Vec<(usize, String)>,
    bam_provided: bool,
    min_snp_reads: usize,
    min_snp_frequency: f64,
    min_mnv_reads: usize,
    min_mnv_frequency: f64,
    min_snp_strand_reads: usize,
    min_mnv_strand_reads: usize,
    min_strand_bias_p: f64,
    emit_filtered: bool,
    include_strand_bias_info: bool,
}

struct SupportFilterInput {
    support_reads: usize,
    min_reads: usize,
    depth: usize,
    min_frequency: f64,
    forward_reads: usize,
    reverse_reads: usize,
    min_strand_reads: usize,
    strand_bias_p: Option<f64>,
}

impl VcfWriter {
    pub fn new(cfg: VcfWriterConfig<'_>) -> AppResult<Self> {
        let filename = cfg.filename;
        let bam_provided = cfg.bam_provided;
        let min_snp_reads = cfg.min_snp_reads;
        let min_snp_frequency = cfg.min_snp_frequency;
        let min_mnv_reads = cfg.min_mnv_reads;
        let min_mnv_frequency = cfg.min_mnv_frequency;
        let min_quality = cfg.min_quality;
        let min_mapq = cfg.min_mapq;
        let command_line = cfg.command_line;
        let contigs = cfg.contigs;
        let bgzf_output = cfg.bgzf_output;
        let min_snp_strand_reads = cfg.min_snp_strand_reads;
        let min_mnv_strand_reads = cfg.min_mnv_strand_reads;
        let min_strand_bias_p = cfg.min_strand_bias_p;
        let emit_filtered = cfg.emit_filtered;
        let include_strand_bias_info = cfg.include_strand_bias_info;
        let original_info_headers = cfg.original_info_headers;
        let out_file = if bgzf_output {
            format!("{filename}.MNV.vcf.gz")
        } else {
            format!("{filename}.MNV.vcf")
        };
        let mut writer: Box<dyn Write> = if bgzf_output {
            Box::new(noodles::bgzf::io::Writer::new(File::create(&out_file)?))
        } else {
            // Buffer plain-VCF output so each record is not a separate write(2).
            Box::new(std::io::BufWriter::new(File::create(&out_file)?))
        };

        writeln!(writer, "##fileformat=VCFv4.2")?;
        writeln!(writer, "##source=get_mnv")?;
        writeln!(writer, "##get_mnv_version={}", env!("CARGO_PKG_VERSION"))?;
        writeln!(writer, "##get_mnv_command={command_line}")?;
        writeln!(writer, "##get_mnv_min_quality={min_quality}")?;
        writeln!(writer, "##get_mnv_min_mapq={min_mapq}")?;
        writeln!(writer, "##get_mnv_min_snp_reads={min_snp_reads}")?;
        writeln!(
            writer,
            "##get_mnv_min_snp_frequency={}",
            format_freq(min_snp_frequency)
        )?;
        writeln!(writer, "##get_mnv_min_mnv_reads={min_mnv_reads}")?;
        writeln!(
            writer,
            "##get_mnv_min_mnv_frequency={}",
            format_freq(min_mnv_frequency)
        )?;
        if min_snp_strand_reads > 0 {
            writeln!(
                writer,
                "##get_mnv_min_snp_strand_reads={min_snp_strand_reads}"
            )?;
        }
        if min_mnv_strand_reads > 0 {
            writeln!(
                writer,
                "##get_mnv_min_mnv_strand_reads={min_mnv_strand_reads}"
            )?;
        }
        if min_strand_bias_p > 0.0 {
            writeln!(writer, "##get_mnv_min_strand_bias_p={min_strand_bias_p}")?;
        }
        for contig in contigs {
            // Sanitize contig names to prevent VCF header corruption from
            // control characters in FASTA sequence identifiers
            let safe_id: String = contig
                .chars()
                .map(|c| if c.is_control() { '_' } else { c })
                .collect();
            writeln!(writer, "##contig=<ID={safe_id}>")?;
        }
        write_info_header(
            &mut writer,
            bam_provided,
            include_strand_bias_info,
            original_info_headers,
        )?;
        if emit_filtered {
            writeln!(
                writer,
                "##FILTER=<ID=LowSupport,Description=\"Insufficient supporting reads\">"
            )?;
            writeln!(
                writer,
                "##FILTER=<ID=LowFrequency,Description=\"Allele frequency below threshold\">"
            )?;
            writeln!(
                writer,
                "##FILTER=<ID=StrandSupport,Description=\"Insufficient support in one or both strands\">"
            )?;
            writeln!(
                writer,
                "##FILTER=<ID=StrandBias,Description=\"Strand-bias p-value below threshold\">"
            )?;
        }
        writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

        Ok(Self {
            pending: Vec::new(),
            writer,
            bam_provided,
            min_snp_reads,
            min_snp_frequency,
            min_mnv_reads,
            min_mnv_frequency,
            min_snp_strand_reads,
            min_mnv_strand_reads,
            min_strand_bias_p,
            emit_filtered,
            include_strand_bias_info,
        })
    }

    /// Flush buffered output. Must be called before the writer is dropped so a
    /// flush error is reported instead of being swallowed by `BufWriter::drop`.
    pub fn flush(&mut self) -> AppResult<()> {
        self.writer.flush()?;
        Ok(())
    }

    pub(super) fn write_variant_line(
        &mut self,
        chrom: &str,
        pos: usize,
        ref_allele: &str,
        alt_allele: &str,
        filter: &str,
        info: &str,
    ) -> AppResult<()> {
        self.pending.push((
            pos,
            format!("{chrom}\t{pos}\t.\t{ref_allele}\t{alt_allele}\t.\t{filter}\t{info}"),
        ));
        Ok(())
    }

    /// Write the buffered records in coordinate order and clear the buffer.
    fn flush_pending(&mut self) -> AppResult<()> {
        let entries = std::mem::take(&mut self.pending);
        write_sorted_vcf_entries(&mut self.writer, entries)
    }

    pub(super) fn reference_sequence_for_variant<'a>(
        &self,
        variant: &VariantInfo,
        references: &'a ReferenceMap,
    ) -> AppResult<&'a str> {
        references
            .get(&variant.chrom)
            .map(std::string::String::as_str)
            .ok_or_else(|| {
                format!(
                    "Missing reference sequence for contig '{}' ({})",
                    variant.chrom,
                    variant_context(variant)
                )
                .into()
            })
    }

    pub fn write_variants(
        &mut self,
        variants: &[VariantInfo],
        references: &ReferenceMap,
    ) -> AppResult<()> {
        self.write_variants_buffered(variants, references)?;
        self.flush_pending()
    }

    fn write_variants_buffered(
        &mut self,
        variants: &[VariantInfo],
        references: &ReferenceMap,
    ) -> AppResult<()> {
        for variant in variants {
            if variant.gene == "intergenic" {
                self.write_intergenic(variant)?;
                continue;
            }
            match variant.variant_type {
                VariantType::Indel => self.write_indel(variant)?,
                VariantType::Snp => self.write_snp(variant)?,
                VariantType::Mnv => {
                    let reference_sequence =
                        self.reference_sequence_for_variant(variant, references)?;
                    self.write_mnv(variant, reference_sequence)?;
                }
                VariantType::SnpMnv => {
                    let reference_sequence =
                        self.reference_sequence_for_variant(variant, references)?;
                    self.write_snp_mnv(variant, reference_sequence)?;
                }
            }
        }
        Ok(())
    }
}

pub fn build_tabix_index(vcf_gz_path: &str) -> AppResult<()> {
    if !vcf_gz_path.ends_with(".vcf.gz") {
        return Err(format!("Tabix indexing requires a .vcf.gz file, got '{vcf_gz_path}'").into());
    }
    if !Path::new(vcf_gz_path).exists() {
        return Err(
            format!("Cannot build Tabix index: file '{vcf_gz_path}' does not exist").into(),
        );
    }
    // Use external tabix command (from htslib/samtools)
    let status = std::process::Command::new("tabix")
        .args(["-p", "vcf", vcf_gz_path])
        .status();
    match status {
        Ok(s) if s.success() => Ok(()),
        Ok(s) => Err(format!("tabix exited with status {s} for '{vcf_gz_path}'").into()),
        Err(_) => {
            log::warn!("tabix not found in PATH. Skipping .tbi index for '{vcf_gz_path}'. Install samtools/htslib for automatic indexing.");
            Ok(())
        }
    }
}

pub fn convert_vcf_to_bcf(input_vcf_path: &str, output_bcf_path: &str) -> AppResult<()> {
    if !Path::new(input_vcf_path).exists() {
        return Err(
            format!("Cannot convert to BCF: input VCF '{input_vcf_path}' does not exist").into(),
        );
    }
    // Use external bcftools for VCF→BCF conversion
    let status = std::process::Command::new("bcftools")
        .args(["view", "-O", "b", "-o", output_bcf_path, input_vcf_path])
        .status();
    match status {
        Ok(s) if s.success() => Ok(()),
        Ok(s) => Err(format!("bcftools exited with status {s} converting to BCF").into()),
        Err(_) => {
            log::warn!("bcftools not found in PATH. Skipping BCF output. Install bcftools for BCF conversion.");
            Ok(())
        }
    }
}
