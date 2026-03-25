//! Run summary data types and aggregation logic.

use crate::error::AppResult;
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

#[derive(Debug, Clone, Default, Serialize)]
pub struct RunSummary {
    pub schema_version: String,
    pub sample: Option<String>,
    pub dry_run: bool,
    pub bam_provided: bool,
    pub inputs: RunInputs,
    pub output_tsv: Option<String>,
    pub output_vcf: Option<String>,
    pub output_bcf: Option<String>,
    pub contigs: Vec<ContigSummary>,
    pub timings: RunTimings,
    pub global: GlobalSummary,
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
        format!(
            "{{\"schema_version\":\"1.0.0\",\"error\":\"Failed to serialize summary: {}\"}}",
            e
        )
    }) + "\n"
}

pub fn write_summary_json(path: &str, summary: &RunSummary) -> AppResult<()> {
    let mut file = File::create(path)?;
    serde_json::to_writer_pretty(&mut file, summary)
        .map_err(|e| crate::error::AppError::msg(format!("Failed to serialize summary JSON: {}", e)))?;
    file.write_all(b"\n")?;
    Ok(())
}
