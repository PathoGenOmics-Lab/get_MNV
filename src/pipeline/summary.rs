//! Run summary data types and aggregation logic.

use crate::error::{AppResult, IoResultExt};
use crate::variants::{VariantInfo, VariantType};
use serde::Serialize;
use std::fs::File;
use std::io::Write;

#[derive(Debug, Clone, Default, Serialize)]
pub struct ContigSummary {
    pub contig: String,
    pub snp_records_in_vcf: usize,
    pub mapped_genes: usize,
    pub produced_variants: usize,
    pub snp_variants: usize,
    pub mnv_variants: usize,
    pub snp_mnv_variants: usize,
    pub indel_variants: usize,
    pub intergenic_variants: usize,
    pub region_cache_hits: usize,
    pub region_cache_misses: usize,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct GlobalSummary {
    pub contig_count: usize,
    pub snp_records_in_vcf: usize,
    pub mapped_genes: usize,
    pub produced_variants: usize,
    pub snp_variants: usize,
    pub mnv_variants: usize,
    pub snp_mnv_variants: usize,
    pub indel_variants: usize,
    pub intergenic_variants: usize,
    pub region_cache_hits: usize,
    pub region_cache_misses: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunSummary {
    pub schema_version: String,
    pub sample: Option<String>,
    pub dry_run: bool,
    pub bam_provided: bool,
    pub translation_table: u8,
    pub inputs: RunInputs,
    pub output_tsv: Option<String>,
    pub output_vcf: Option<String>,
    pub output_bcf: Option<String>,
    pub contigs: Vec<ContigSummary>,
    pub timings: RunTimings,
    pub global: GlobalSummary,
}

impl Default for RunSummary {
    fn default() -> Self {
        Self {
            schema_version: String::new(),
            sample: None,
            dry_run: false,
            bam_provided: false,
            translation_table: 11,
            inputs: RunInputs::default(),
            output_tsv: None,
            output_vcf: None,
            output_bcf: None,
            contigs: Vec::new(),
            timings: RunTimings::default(),
            global: GlobalSummary::default(),
        }
    }
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct RunTimings {
    pub parse_inputs_ms: f64,
    pub process_ms: f64,
    pub emit_ms: f64,
    pub total_ms: f64,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct RunInputs {
    pub vcf: String,
    pub fasta: String,
    pub annotation: String,
    pub bam: Option<String>,
    pub checksums: InputChecksums,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct InputChecksums {
    pub vcf_sha256: String,
    pub fasta_sha256: String,
    pub annotation_sha256: String,
    pub bam_sha256: Option<String>,
}

/// Progress event emitted during pipeline execution (for desktop GUI).
#[derive(Debug, Clone, Serialize)]
pub struct ProgressEvent {
    pub phase: String,
    pub contig: Option<String>,
    pub current: usize,
    pub total: usize,
}

pub fn summarize_contig_variants(
    contig: &str,
    snp_records_in_vcf: usize,
    mapped_genes: usize,
    variants: &[VariantInfo],
    region_cache_hits: usize,
    region_cache_misses: usize,
) -> ContigSummary {
    let mut summary = ContigSummary {
        contig: contig.to_string(),
        snp_records_in_vcf,
        mapped_genes,
        produced_variants: variants.len(),
        region_cache_hits,
        region_cache_misses,
        ..ContigSummary::default()
    };
    for variant in variants {
        if variant.gene == "intergenic" {
            summary.intergenic_variants += 1;
        } else {
            match variant.variant_type {
                VariantType::Snp => summary.snp_variants += 1,
                VariantType::Mnv => summary.mnv_variants += 1,
                VariantType::SnpMnv => summary.snp_mnv_variants += 1,
                VariantType::Indel => summary.indel_variants += 1,
            }
        }
    }
    summary
}

pub fn update_global_summary(global: &mut GlobalSummary, contig_summary: &ContigSummary) {
    global.contig_count += 1;
    global.snp_records_in_vcf += contig_summary.snp_records_in_vcf;
    global.mapped_genes += contig_summary.mapped_genes;
    global.produced_variants += contig_summary.produced_variants;
    global.snp_variants += contig_summary.snp_variants;
    global.mnv_variants += contig_summary.mnv_variants;
    global.snp_mnv_variants += contig_summary.snp_mnv_variants;
    global.indel_variants += contig_summary.indel_variants;
    global.intergenic_variants += contig_summary.intergenic_variants;
    global.region_cache_hits += contig_summary.region_cache_hits;
    global.region_cache_misses += contig_summary.region_cache_misses;
}

#[cfg(test)]
pub(crate) fn summary_to_json(summary: &RunSummary) -> String {
    serde_json::to_string_pretty(summary).unwrap_or_else(|e| {
        format!("{{\"schema_version\":\"1.0.0\",\"error\":\"Failed to serialize summary: {e}\"}}")
    }) + "\n"
}

pub fn write_summary_json(path: &str, summary: &RunSummary) -> AppResult<()> {
    let mut file = File::create(path).with_path(path)?;
    serde_json::to_writer_pretty(&mut file, summary)?;
    file.write_all(b"\n").with_path(path)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::VariantType;

    fn make_variant(gene: &str, vtype: VariantType) -> VariantInfo {
        VariantInfo {
            chrom: "chr1".to_string(),
            gene: gene.to_string(),
            positions: vec![100],
            ref_bases: vec!["A".to_string()],
            base_changes: vec!["T".to_string()],
            aa_changes: vec!["V100A".to_string()],
            snp_aa_changes: vec!["V100A".to_string()],
            variant_type: vtype,
            change_type: crate::variants::ChangeType::NonSynonymous,
            snp_reads: None,
            snp_forward_reads: None,
            snp_reverse_reads: None,
            mnv_reads: None,
            mnv_forward_reads: None,
            mnv_reverse_reads: None,
            mnv_total_reads: None,
            total_reads: None,
            total_forward_reads: None,
            total_reverse_reads: None,
            mnv_total_forward_reads: None,
            mnv_total_reverse_reads: None,
            ref_codon: None,
            snp_codon: None,
            mnv_codon: None,
            original_dp: None,
            original_freq: None,
            original_info: None,
        }
    }

    #[test]
    fn test_summarize_contig_variants_counts() {
        let variants = vec![
            make_variant("geneA", VariantType::Snp),
            make_variant("geneA", VariantType::Mnv),
            make_variant("intergenic", VariantType::Snp),
            make_variant("geneB", VariantType::Indel),
        ];
        let summary = summarize_contig_variants("chr1", 10, 2, &variants, 5, 3);
        assert_eq!(summary.snp_variants, 1);
        assert_eq!(summary.mnv_variants, 1);
        assert_eq!(summary.indel_variants, 1);
        assert_eq!(summary.intergenic_variants, 1);
        assert_eq!(summary.produced_variants, 4);
        assert_eq!(summary.region_cache_hits, 5);
        assert_eq!(summary.region_cache_misses, 3);
    }

    #[test]
    fn test_update_global_summary_accumulates() {
        let mut global = GlobalSummary::default();
        let contig1 = ContigSummary {
            contig: "chr1".to_string(),
            snp_records_in_vcf: 10,
            mapped_genes: 5,
            produced_variants: 8,
            snp_variants: 4,
            mnv_variants: 2,
            snp_mnv_variants: 1,
            indel_variants: 0,
            intergenic_variants: 1,
            region_cache_hits: 3,
            region_cache_misses: 2,
        };
        let contig2 = ContigSummary {
            contig: "chr2".to_string(),
            snp_records_in_vcf: 5,
            mapped_genes: 3,
            produced_variants: 4,
            snp_variants: 2,
            mnv_variants: 1,
            snp_mnv_variants: 0,
            indel_variants: 1,
            intergenic_variants: 0,
            region_cache_hits: 1,
            region_cache_misses: 1,
        };
        update_global_summary(&mut global, &contig1);
        update_global_summary(&mut global, &contig2);
        assert_eq!(global.contig_count, 2);
        assert_eq!(global.snp_records_in_vcf, 15);
        assert_eq!(global.produced_variants, 12);
        assert_eq!(global.region_cache_hits, 4);
    }
}
