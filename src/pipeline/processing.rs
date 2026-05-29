//! Per-contig parallel processing: variant annotation, BAM read counting,
//! LRU cache, and output emission.
//!
//! Input parsing and per-gene read support live in submodules; this file owns
//! the worker state types and the per-contig orchestration (`process_contig`)
//! plus output emission.

use super::summary::{summarize_contig_variants, ContigSummary};
use crate::cli::Args;
use crate::error::{AppError, AppResult, ErrorCode};
use crate::io::{self, ReferenceMap, VcfPosition};
use crate::output;
use crate::read_count::RegionObservationCache;
use crate::variants::{self, Gene, VariantInfo};
use log::info;
use lru::LruCache;
use noodles::bam;
use noodles::sam::Header;
use rayon::prelude::*;
use std::collections::HashMap;
use std::num::NonZeroUsize;

mod inputs;
mod read_support;

pub(crate) use inputs::{parse_inputs, resolve_variant_input_format};
use read_support::{append_supported_phased_indel_haplotypes, count_gene_variant_reads};

#[derive(Debug)]
pub(crate) struct ParsedInputs {
    pub(crate) base_name: String,
    pub(crate) references: ReferenceMap,
    pub(crate) snp_by_contig: HashMap<String, Vec<VcfPosition>>,
    pub(crate) contigs: Vec<String>,
    pub(crate) command_line: String,
    pub(crate) preloaded_gff: Option<HashMap<String, Vec<Gene>>>,
    pub(crate) original_info_headers: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct RegionCacheKey {
    contig: String,
    start: usize,
    end: usize,
    positions: Vec<usize>,
    min_mapq: u8,
}

struct WorkerState {
    bam: Option<bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>>,
    bam_header: Option<Header>,
    region_cache: LruCache<RegionCacheKey, RegionObservationCache>,
}

#[derive(Default)]
struct WorkerResult {
    variants: Vec<VariantInfo>,
    region_cache_hits: usize,
    region_cache_misses: usize,
}

struct PhasedIndelInputs<'a, 'r> {
    contig: &'a str,
    gene: &'a Gene,
    reference: &'a io::Reference<'r>,
    snp_list: &'a [VcfPosition],
    genetic_code: crate::genetic_code::GeneticCode,
}

pub(crate) fn sort_variants(variants: &mut [VariantInfo]) {
    variants.sort_by(|a, b| {
        let a_min = a.positions.iter().copied().min().unwrap_or(usize::MAX);
        let b_min = b.positions.iter().copied().min().unwrap_or(usize::MAX);
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| a_min.cmp(&b_min))
            .then_with(|| a.positions.cmp(&b.positions))
            .then_with(|| a.variant_type.cmp(&b.variant_type))
            .then_with(|| a.gene.cmp(&b.gene))
    });
}

pub(crate) fn reclassify_generic_as_validation(error: AppError) -> AppError {
    if error.code == ErrorCode::Generic {
        AppError::validation(error.message)
    } else {
        error
    }
}

fn annotate_variants_for_gene(
    gene: &Gene,
    snp_list: &[VcfPosition],
    reference: &io::Reference,
    contig: &str,
    genetic_code: crate::genetic_code::GeneticCode,
    indel_config: &variants::IndelAnnotationConfig,
) -> Vec<VariantInfo> {
    variants::get_mnv_variants_for_gene_with_config(
        gene,
        snp_list,
        reference,
        contig,
        genetic_code,
        indel_config,
    )
}

pub(crate) fn process_contig(
    args: &Args,
    contig: &str,
    references: &ReferenceMap,
    snp_by_contig: &HashMap<String, Vec<VcfPosition>>,
    preloaded_gff: Option<&HashMap<String, Vec<Gene>>>,
    genetic_code: crate::genetic_code::GeneticCode,
) -> AppResult<(Vec<VariantInfo>, ContigSummary)> {
    let reference =
        io::reference_for_chrom(references, contig).map_err(reclassify_generic_as_validation)?;
    let snp_list = snp_by_contig
        .get(contig)
        .ok_or_else(|| AppError::validation(format!("Missing VCF data for contig '{contig}'")))?;
    io::validate_vcf_reference_alleles(contig, snp_list, &reference)
        .map_err(reclassify_generic_as_validation)?;

    let genes = if let Some(gff_genes) = preloaded_gff {
        let all_contig_genes = gff_genes.get(contig).cloned().unwrap_or_default();
        let filtered = io::filter_genes_with_snps(&all_contig_genes, snp_list);
        log::info!(
            "GFF/GFF3 contig '{}': {} gene entries, {} mapped to SNPs, {} without SNPs",
            contig,
            all_contig_genes.len(),
            filtered.len(),
            all_contig_genes.len() - filtered.len()
        );
        filtered
    } else {
        io::load_genes(
            args.genes_file(),
            snp_list,
            Some(contig),
            &args.gff_features(),
        )
        .map_err(reclassify_generic_as_validation)?
    };
    info!(
        "Contig '{}' -> {} SNP/variant records in VCF, {} mapped genes",
        contig,
        snp_list.len(),
        genes.len()
    );

    let bam_path = args.bam_file.clone();
    let should_count_reads = args.bam_file.is_some() && !args.dry_run;
    let worker_results: Result<Vec<WorkerResult>, AppError> = genes
        .par_iter()
        .map_init(
            || -> AppResult<WorkerState> {
                let (bam, bam_header) = if should_count_reads {
                    if let Some(path) = bam_path.as_ref() {
                        let mut reader = bam::io::indexed_reader::Builder::default()
                            .build_from_path(path)
                            .map_err(|e| {
                                AppError::validation(format!(
                                    "Failed to open BAM '{path}' in worker thread: {e}"
                                ))
                            })?;
                        let header = reader.read_header().map_err(|e| {
                            AppError::validation(format!(
                                "Failed to read BAM header from '{path}': {e}"
                            ))
                        })?;
                        (Some(reader), Some(header))
                    } else {
                        (None, None)
                    }
                } else {
                    (None, None)
                };
                Ok(WorkerState {
                    bam,
                    bam_header,
                    region_cache: LruCache::new(NonZeroUsize::new(64).expect("64 > 0")),
                })
            },
            |state_result, gene| -> AppResult<WorkerResult> {
                let state = state_result
                    .as_mut()
                    .map_err(|err| AppError::validation(err.to_string()))?;
                let indel_config = variants::IndelAnnotationConfig {
                    frameshift_min_freq: args.frameshift_min_freq,
                };
                let mut variants = annotate_variants_for_gene(
                    gene,
                    snp_list,
                    &reference,
                    contig,
                    genetic_code,
                    &indel_config,
                );
                if variants.is_empty() {
                    return Ok(WorkerResult::default());
                }
                let (cache_hits, cache_misses) =
                    count_gene_variant_reads(state, args, contig, gene, &mut variants)?;
                append_supported_phased_indel_haplotypes(
                    state,
                    args,
                    PhasedIndelInputs {
                        contig,
                        gene,
                        reference: &reference,
                        snp_list,
                        genetic_code,
                    },
                    &mut variants,
                )?;
                Ok(WorkerResult {
                    variants,
                    region_cache_hits: cache_hits,
                    region_cache_misses: cache_misses,
                })
            },
        )
        .collect();

    let mut all_variants = Vec::new();
    let mut cache_hits = 0usize;
    let mut cache_misses = 0usize;
    for result in worker_results.map_err(|e| {
        AppError::validation(format!("Error while processing contig '{contig}': {e}"))
    })? {
        cache_hits += result.region_cache_hits;
        cache_misses += result.region_cache_misses;
        all_variants.extend(result.variants);
    }

    // Collect intergenic variants (positions not covered by any gene).
    if !args.exclude_intergenic {
        let mut covered = vec![false; snp_list.len()];
        for gene in &genes {
            for (idx, snp) in snp_list.iter().enumerate() {
                if io::gene_overlaps_variant(gene, snp) {
                    covered[idx] = true;
                }
            }
        }
        let mut intergenic_count = 0usize;
        for (idx, snp) in snp_list.iter().enumerate() {
            if !covered[idx] {
                all_variants.push(variants::build_intergenic_variant(contig, snp));
                intergenic_count += 1;
            }
        }
        if intergenic_count > 0 {
            info!("Contig '{contig}' -> {intergenic_count} intergenic variant(s) added");
        }
    }

    sort_variants(&mut all_variants);
    let summary = summarize_contig_variants(
        contig,
        snp_list.len(),
        genes.len(),
        &all_variants,
        cache_hits,
        cache_misses,
    );
    Ok((all_variants, summary))
}

pub(crate) fn emit_contig_variants(
    tsv_writer: Option<&mut output::TsvWriter>,
    vcf_writer: Option<&mut output::VcfWriter>,
    variants: &[VariantInfo],
    references: &ReferenceMap,
) -> AppResult<()> {
    if let Some(writer) = tsv_writer {
        writer.write_variants(variants)?;
    }
    if let Some(writer) = vcf_writer {
        writer.write_variants(variants, references)?;
    }
    Ok(())
}
