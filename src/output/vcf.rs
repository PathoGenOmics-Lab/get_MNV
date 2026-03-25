//! VCF/BCF output writer: generates VCF (plain or BGZF-compressed) with
//! INFO fields, FILTER tags, Tabix indexing, and BCF conversion.

use crate::error::AppResult;
use crate::io::ReferenceMap;
use crate::variants::{VariantInfo, VariantType};
use rust_htslib::bcf;
use rust_htslib::bcf::Read as BcfReadTrait;
use rust_htslib::bgzf;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use super::common::{
    build_alt_region, build_info_string, filter_value, get_mnv_depth_from_variant, get_required,
    mnv_strand_bias_p_value, reference_subsequence, snp_aa_for_index, snp_bam_vectors,
    snp_strand_bias_p_value, validate_variant_shape, variant_context, vcf_entry_line,
    write_info_header, write_sorted_vcf_entries, MnvCallMetrics, SnpBamVectors, SnpCallMetrics,
    VcfEntry,
};

pub struct VcfWriter {
    writer: Box<dyn Write>,
    bam_provided: bool,
    min_snp_reads: usize,
    min_mnv_reads: usize,
    min_snp_strand_reads: usize,
    min_mnv_strand_reads: usize,
    min_strand_bias_p: f64,
    emit_filtered: bool,
    include_strand_bias_info: bool,
}

impl VcfWriter {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        filename: &str,
        bam_provided: bool,
        min_snp_reads: usize,
        min_mnv_reads: usize,
        min_quality: u8,
        min_mapq: u8,
        command_line: &str,
        contigs: &[String],
        bgzf_output: bool,
        min_snp_strand_reads: usize,
        min_mnv_strand_reads: usize,
        min_strand_bias_p: f64,
        emit_filtered: bool,
        include_strand_bias_info: bool,
        original_info_headers: &[String],
    ) -> AppResult<Self> {
        let out_file = if bgzf_output {
            format!("{}.MNV.vcf.gz", filename)
        } else {
            format!("{}.MNV.vcf", filename)
        };
        let mut writer: Box<dyn Write> = if bgzf_output {
            Box::new(bgzf::Writer::from_path(&out_file)?)
        } else {
            Box::new(File::create(&out_file)?)
        };

        writeln!(writer, "##fileformat=VCFv4.2")?;
        writeln!(writer, "##source=get_mnv")?;
        writeln!(writer, "##get_mnv_version={}", env!("CARGO_PKG_VERSION"))?;
        writeln!(writer, "##get_mnv_command={}", command_line)?;
        writeln!(writer, "##get_mnv_min_quality={}", min_quality)?;
        writeln!(writer, "##get_mnv_min_mapq={}", min_mapq)?;
        writeln!(writer, "##get_mnv_min_snp_reads={}", min_snp_reads)?;
        writeln!(writer, "##get_mnv_min_mnv_reads={}", min_mnv_reads)?;
        if min_snp_strand_reads > 0 {
            writeln!(
                writer,
                "##get_mnv_min_snp_strand_reads={}",
                min_snp_strand_reads
            )?;
        }
        if min_mnv_strand_reads > 0 {
            writeln!(
                writer,
                "##get_mnv_min_mnv_strand_reads={}",
                min_mnv_strand_reads
            )?;
        }
        if min_strand_bias_p > 0.0 {
            writeln!(writer, "##get_mnv_min_strand_bias_p={}", min_strand_bias_p)?;
        }
        for contig in contigs {
            writeln!(writer, "##contig=<ID={}>", contig)?;
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
                "##FILTER=<ID=StrandSupport,Description=\"Insufficient support in one or both strands\">"
            )?;
            writeln!(
                writer,
                "##FILTER=<ID=StrandBias,Description=\"Strand-bias p-value below threshold\">"
            )?;
        }
        writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

        Ok(Self {
            writer,
            bam_provided,
            min_snp_reads,
            min_mnv_reads,
            min_snp_strand_reads,
            min_mnv_strand_reads,
            min_strand_bias_p,
            emit_filtered,
            include_strand_bias_info,
        })
    }

    fn write_variant_line(
        &mut self,
        chrom: &str,
        pos: usize,
        ref_allele: &str,
        alt_allele: &str,
        filter: &str,
        info: &str,
    ) -> AppResult<()> {
        writeln!(
            self.writer,
            "{}\t{}\t.\t{}\t{}\t.\t{}\t{}",
            chrom, pos, ref_allele, alt_allele, filter, info
        )?;
        Ok(())
    }

    fn reference_sequence_for_variant<'a>(
        &self,
        variant: &VariantInfo,
        references: &'a ReferenceMap,
    ) -> AppResult<&'a str> {
        references
            .get(&variant.chrom)
            .map(|s| s.as_str())
            .ok_or_else(|| {
                format!(
                    "Missing reference sequence for contig '{}' ({})",
                    variant.chrom,
                    variant_context(variant)
                )
                .into()
            })
    }

    fn should_evaluate_strand_bias(&self) -> bool {
        self.include_strand_bias_info || self.min_strand_bias_p > 0.0
    }

    fn snp_strand_bias(&self, variant: &VariantInfo, index: usize) -> Option<f64> {
        if self.should_evaluate_strand_bias() {
            snp_strand_bias_p_value(variant, index)
        } else {
            None
        }
    }

    fn mnv_strand_bias(&self, variant: &VariantInfo) -> Option<f64> {
        if self.should_evaluate_strand_bias() {
            mnv_strand_bias_p_value(variant)
        } else {
            None
        }
    }

    fn build_support_filters(
        &self,
        support_reads: usize,
        min_reads: usize,
        forward_reads: usize,
        reverse_reads: usize,
        min_strand_reads: usize,
        strand_bias_p: Option<f64>,
    ) -> Vec<&'static str> {
        let mut filters = Vec::new();
        if support_reads < min_reads {
            filters.push("LowSupport");
        }
        if forward_reads < min_strand_reads || reverse_reads < min_strand_reads {
            filters.push("StrandSupport");
        }
        if let Some(p_value) = strand_bias_p {
            if self.min_strand_bias_p > 0.0 && p_value < self.min_strand_bias_p {
                filters.push("StrandBias");
            }
        }
        filters
    }

    fn should_emit_record(&self, filters: &[&str]) -> bool {
        filters.is_empty() || self.emit_filtered
    }

    fn snp_metrics_at(
        &self,
        variant: &VariantInfo,
        bam_vectors: &SnpBamVectors<'_>,
        index: usize,
    ) -> AppResult<SnpCallMetrics> {
        Ok(SnpCallMetrics {
            support_reads: *get_required(bam_vectors.reads, index, "snp_reads", variant)?,
            forward_reads: *get_required(
                bam_vectors.forward_reads,
                index,
                "snp_forward_reads",
                variant,
            )?,
            reverse_reads: *get_required(
                bam_vectors.reverse_reads,
                index,
                "snp_reverse_reads",
                variant,
            )?,
            depth: *get_required(bam_vectors.total_reads, index, "total_reads", variant)?,
            strand_bias_p: self.snp_strand_bias(variant, index),
        })
    }

    fn mnv_metrics(
        &self,
        variant: &VariantInfo,
        total_reads: &[usize],
    ) -> AppResult<MnvCallMetrics> {
        let support_reads = variant
            .mnv_reads
            .ok_or_else(|| format!("Missing MNV read counts for {}", variant_context(variant)))?;
        let forward_reads = variant.mnv_forward_reads.ok_or_else(|| {
            format!(
                "Missing MNV forward read counts for {}",
                variant_context(variant)
            )
        })?;
        let reverse_reads = variant.mnv_reverse_reads.ok_or_else(|| {
            format!(
                "Missing MNV reverse read counts for {}",
                variant_context(variant)
            )
        })?;

        Ok(MnvCallMetrics {
            support_reads,
            forward_reads,
            reverse_reads,
            depth: get_mnv_depth_from_variant(variant, total_reads),
            strand_bias_p: self.mnv_strand_bias(variant),
        })
    }

    fn build_snp_entry(
        &self,
        variant: &VariantInfo,
        index: usize,
        aa: &str,
        metrics: SnpCallMetrics,
    ) -> AppResult<Option<VcfEntry>> {
        let filters = self.build_support_filters(
            metrics.support_reads,
            self.min_snp_reads,
            metrics.forward_reads,
            metrics.reverse_reads,
            self.min_snp_strand_reads,
            metrics.strand_bias_p,
        );
        if !self.should_emit_record(&filters) {
            return Ok(None);
        }

        let pos = *get_required(&variant.positions, index, "positions", variant)?;
        let ref_base = get_required(&variant.ref_bases, index, "ref_bases", variant)?;
        let alt_base = get_required(&variant.base_changes, index, "base_changes", variant)?;
        let info = build_info_string(
            variant,
            Some(aa),
            VariantType::Snp.as_str(),
            Some((
                metrics.support_reads,
                metrics.forward_reads,
                metrics.reverse_reads,
            )),
            None,
            Some(index),
            Some(metrics.depth),
            Some(metrics.support_reads),
            if self.include_strand_bias_info {
                metrics.strand_bias_p
            } else {
                None
            },
            None,
            variant.original_info.as_deref(),
        );
        let filter = filter_value(&filters);
        Ok(Some((
            pos,
            vcf_entry_line(&variant.chrom, pos, ref_base, alt_base, &filter, &info),
        )))
    }

    fn build_mnv_entry(
        &self,
        variant: &VariantInfo,
        reference_sequence: &str,
        aa: &str,
        metrics: MnvCallMetrics,
    ) -> AppResult<Option<VcfEntry>> {
        let filters = self.build_support_filters(
            metrics.support_reads,
            self.min_mnv_reads,
            metrics.forward_reads,
            metrics.reverse_reads,
            self.min_mnv_strand_reads,
            metrics.strand_bias_p,
        );
        if !self.should_emit_record(&filters) {
            return Ok(None);
        }

        let min_pos = *variant
            .positions
            .iter()
            .min()
            .ok_or_else(|| format!("Missing positions for {}", variant_context(variant)))?;
        let max_pos = *variant
            .positions
            .iter()
            .max()
            .ok_or_else(|| format!("Missing positions for {}", variant_context(variant)))?;
        let ref_region = reference_subsequence(reference_sequence, min_pos, max_pos, variant)?;
        let alt_region = build_alt_region(
            reference_sequence,
            &variant.positions,
            &variant.base_changes,
        )?;
        let info = build_info_string(
            variant,
            Some(aa),
            VariantType::Mnv.as_str(),
            None,
            Some((
                metrics.support_reads,
                metrics.forward_reads,
                metrics.reverse_reads,
            )),
            None,
            Some(metrics.depth),
            Some(metrics.support_reads),
            None,
            if self.include_strand_bias_info {
                metrics.strand_bias_p
            } else {
                None
            },
            variant.original_info.as_deref(),
        );
        let filter = filter_value(&filters);
        Ok(Some((
            min_pos,
            vcf_entry_line(
                &variant.chrom,
                min_pos,
                ref_region,
                &alt_region,
                &filter,
                &info,
            ),
        )))
    }

    fn write_indel(&mut self, variant: &VariantInfo) -> AppResult<()> {
        validate_variant_shape(variant)?;
        let ref_base = get_required(&variant.ref_bases, 0, "ref_bases", variant)?;
        let alt_base = get_required(&variant.base_changes, 0, "base_changes", variant)?;
        let pos = *get_required(&variant.positions, 0, "positions", variant)?;
        let info = build_info_string(
            variant,
            None,
            VariantType::Indel.as_str(),
            None,
            None,
            Some(0),
            None,
            None,
            None,
            None,
            variant.original_info.as_deref(),
        );
        self.write_variant_line(&variant.chrom, pos, ref_base, alt_base, "PASS", &info)
    }

    fn write_snp(&mut self, variant: &VariantInfo) -> AppResult<()> {
        validate_variant_shape(variant)?;
        if self.bam_provided {
            let bam_vectors = snp_bam_vectors(variant)?;
            for i in 0..variant.positions.len() {
                let metrics = self.snp_metrics_at(variant, &bam_vectors, i)?;
                let aa = variant.aa_changes.join(",");
                if let Some((_, line)) = self.build_snp_entry(variant, i, &aa, metrics)? {
                    writeln!(self.writer, "{}", line)?;
                }
            }
            Ok(())
        } else {
            for (i, &pos) in variant.positions.iter().enumerate() {
                let ref_base = get_required(&variant.ref_bases, i, "ref_bases", variant)?;
                let alt_base = get_required(&variant.base_changes, i, "base_changes", variant)?;
                let aa = variant.aa_changes.join(",");
                let info = build_info_string(
                    variant,
                    Some(&aa),
                    VariantType::Snp.as_str(),
                    None,
                    None,
                    Some(i),
                    None,
                    None,
                    None,
                    None,
                    variant.original_info.as_deref(),
                );
                self.write_variant_line(&variant.chrom, pos, ref_base, alt_base, "PASS", &info)?;
            }
            Ok(())
        }
    }

    fn write_mnv(&mut self, variant: &VariantInfo, reference_sequence: &str) -> AppResult<()> {
        validate_variant_shape(variant)?;
        if self.bam_provided {
            let total_reads = variant.total_reads.as_ref().ok_or_else(|| {
                format!("Missing total read depth for {}", variant_context(variant))
            })?;
            let metrics = self.mnv_metrics(variant, total_reads)?;
            let aa = variant.aa_changes.join(",");
            if let Some((_, line)) =
                self.build_mnv_entry(variant, reference_sequence, &aa, metrics)?
            {
                writeln!(self.writer, "{}", line)?;
            }
            Ok(())
        } else {
            let min_pos = *variant
                .positions
                .iter()
                .min()
                .ok_or_else(|| format!("Missing positions for {}", variant_context(variant)))?;
            let max_pos = *variant
                .positions
                .iter()
                .max()
                .ok_or_else(|| format!("Missing positions for {}", variant_context(variant)))?;
            let ref_region = reference_subsequence(reference_sequence, min_pos, max_pos, variant)?;
            let alt_region = build_alt_region(
                reference_sequence,
                &variant.positions,
                &variant.base_changes,
            )?;
            let aa = variant.aa_changes.join(",");
            let info = build_info_string(
                variant,
                Some(&aa),
                VariantType::Mnv.as_str(),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                variant.original_info.as_deref(),
            );
            self.write_variant_line(
                &variant.chrom,
                min_pos,
                ref_region,
                &alt_region,
                "PASS",
                &info,
            )
        }
    }

    fn write_snp_mnv(&mut self, variant: &VariantInfo, reference_sequence: &str) -> AppResult<()> {
        validate_variant_shape(variant)?;
        let mut entries: Vec<VcfEntry> = Vec::new();

        if self.bam_provided {
            let bam_vectors = snp_bam_vectors(variant)?;
            for i in 0..variant.positions.len() {
                let metrics = self.snp_metrics_at(variant, &bam_vectors, i)?;
                if metrics.support_reads == 0 {
                    continue;
                }
                let aa = snp_aa_for_index(variant, i);
                if let Some(entry) = self.build_snp_entry(variant, i, &aa, metrics)? {
                    entries.push(entry);
                }
            }

            let mnv_metrics = self.mnv_metrics(variant, bam_vectors.total_reads)?;
            if mnv_metrics.support_reads > 0 {
                let aa = variant.aa_changes.join(",");
                if let Some(entry) =
                    self.build_mnv_entry(variant, reference_sequence, &aa, mnv_metrics)?
                {
                    entries.push(entry);
                }
            }
        } else {
            for (i, &pos) in variant.positions.iter().enumerate() {
                let ref_base = get_required(&variant.ref_bases, i, "ref_bases", variant)?;
                let alt_base = get_required(&variant.base_changes, i, "base_changes", variant)?;
                let aa = snp_aa_for_index(variant, i);
                let info = build_info_string(
                    variant,
                    Some(&aa),
                    VariantType::Snp.as_str(),
                    None,
                    None,
                    Some(i),
                    None,
                    None,
                    None,
                    None,
                    variant.original_info.as_deref(),
                );
                entries.push((
                    pos,
                    vcf_entry_line(&variant.chrom, pos, ref_base, alt_base, "PASS", &info),
                ));
            }

            let min_pos = *variant
                .positions
                .iter()
                .min()
                .ok_or_else(|| format!("Missing positions for {}", variant_context(variant)))?;
            let max_pos = *variant
                .positions
                .iter()
                .max()
                .ok_or_else(|| format!("Missing positions for {}", variant_context(variant)))?;
            let ref_region = reference_subsequence(reference_sequence, min_pos, max_pos, variant)?;
            let alt_region = build_alt_region(
                reference_sequence,
                &variant.positions,
                &variant.base_changes,
            )?;

            let aa = variant.aa_changes.join(",");
            let info = build_info_string(
                variant,
                Some(&aa),
                VariantType::Mnv.as_str(),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                variant.original_info.as_deref(),
            );
            entries.push((
                min_pos,
                vcf_entry_line(
                    &variant.chrom,
                    min_pos,
                    ref_region,
                    &alt_region,
                    "PASS",
                    &info,
                ),
            ));
        }

        write_sorted_vcf_entries(&mut self.writer, entries)
    }

    fn write_intergenic(&mut self, variant: &VariantInfo) -> AppResult<()> {
        validate_variant_shape(variant)?;
        // When BAM-based filters are active, intergenic variants have no read
        // support so they are subject to the same thresholds (support=0).
        let filters = if self.bam_provided {
            self.build_support_filters(0, self.min_snp_reads, 0, 0, self.min_snp_strand_reads, None)
        } else {
            Vec::new()
        };
        if !self.should_emit_record(&filters) {
            return Ok(());
        }
        let ref_base = get_required(&variant.ref_bases, 0, "ref_bases", variant)?;
        let alt_base = get_required(&variant.base_changes, 0, "base_changes", variant)?;
        let pos = *get_required(&variant.positions, 0, "positions", variant)?;
        let info = build_info_string(
            variant,
            None,
            variant.variant_type.as_str(),
            None,
            None,
            Some(0),
            None,
            None,
            None,
            None,
            variant.original_info.as_deref(),
        );
        let filter = filter_value(&filters);
        self.write_variant_line(&variant.chrom, pos, ref_base, alt_base, &filter, &info)
    }

    pub fn write_variants(
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
        return Err(format!(
            "Tabix indexing requires a .vcf.gz file, got '{}'",
            vcf_gz_path
        )
        .into());
    }
    if !Path::new(vcf_gz_path).exists() {
        return Err(format!(
            "Cannot build Tabix index: file '{}' does not exist",
            vcf_gz_path
        )
        .into());
    }
    bcf::index::build(vcf_gz_path, None, 1, bcf::index::Type::Tbx)
        .map_err(|e| format!("Failed to build Tabix index for '{}': {}", vcf_gz_path, e).into())
}

pub fn convert_vcf_to_bcf(input_vcf_path: &str, output_bcf_path: &str) -> AppResult<()> {
    if !Path::new(input_vcf_path).exists() {
        return Err(format!(
            "Cannot convert to BCF: input VCF '{}' does not exist",
            input_vcf_path
        )
        .into());
    }
    let mut reader = bcf::Reader::from_path(input_vcf_path)?;
    let header = bcf::Header::from_template(reader.header());
    let mut writer = bcf::Writer::from_path(output_bcf_path, &header, false, bcf::Format::Bcf)?;
    for result in reader.records() {
        let record = result?;
        writer.write(&record)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::write_sorted_vcf_entries;
    use std::fs::{self, File};
    use std::io::Read;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_temp_path(prefix: &str, suffix: &str) -> std::path::PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time before UNIX_EPOCH")
            .as_nanos();
        std::env::temp_dir().join(format!("{}_{}_{}", prefix, nanos, suffix))
    }

    #[test]
    fn test_write_sorted_vcf_entries_orders_by_position() {
        let path = unique_temp_path("get_mnv_vcf_sort", "vcf");
        let mut file = File::create(&path).expect("failed to create temp file");
        write_sorted_vcf_entries(
            &mut file,
            vec![
                (20, "chr\t20\t.\tA\tT\t.\tPASS\tTYPE=SNP".to_string()),
                (10, "chr\t10\t.\tA\tC\t.\tPASS\tTYPE=SNP".to_string()),
            ],
        )
        .expect("failed to write sorted entries");
        drop(file);

        let mut contents = String::new();
        File::open(&path)
            .expect("failed to open temp file")
            .read_to_string(&mut contents)
            .expect("failed to read temp file");
        let lines: Vec<&str> = contents.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].contains("\t10\t"));
        assert!(lines[1].contains("\t20\t"));
        let _ = fs::remove_file(path);
    }
}
