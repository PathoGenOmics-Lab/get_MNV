//! Per-contig parallel processing: variant annotation, BAM read counting,
//! LRU cache, and output emission.

use super::config::{
    sanitized_command_line, selected_contigs, validate_contig_inputs,
    validate_strict_original_metrics,
};
// build_input_metadata used by mod.rs, not here
use super::summary::{summarize_contig_variants, ContigSummary};
use crate::cli::Args;
use crate::error::ErrorCode;
use crate::error::{AppError, AppResult};
use crate::io::{self, AnnotationFormat, ReferenceMap, VcfPosition};
use crate::output;
use crate::read_count::{self, ReadCountSummary, RegionObservationCache};
use crate::variants::{self, Gene, VariantInfo, VariantType};
use log::info;
use lru::LruCache;
use rayon::prelude::*;
use rust_htslib::bam::IndexedReader;
use std::collections::{HashMap, HashSet};
use std::num::NonZeroUsize;

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
    bam: Option<IndexedReader>,
    region_cache: LruCache<RegionCacheKey, RegionObservationCache>,
}

#[derive(Default)]
struct WorkerResult {
    variants: Vec<VariantInfo>,
    region_cache_hits: usize,
    region_cache_misses: usize,
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

pub(crate) fn parse_inputs(args: &Args, sample_override: Option<&str>) -> AppResult<ParsedInputs> {
    let base_name = io::get_base_name(&args.vcf_file).map_err(reclassify_generic_as_validation)?;
    let references =
        io::load_references(&args.fasta_file).map_err(reclassify_generic_as_validation)?;
    // Use fast text parser for plain .vcf files, htslib for .bcf/.vcf.gz
    let snp_by_contig = if io::vcf_fast::use_fast_parser(&args.vcf_file) {
        io::vcf_fast::load_vcf_text(
            &args.vcf_file,
            sample_override,
            args.split_multiallelic,
            args.normalize_alleles,
            args.keep_original_info,
        )
        .map_err(reclassify_generic_as_validation)?
    } else {
        io::load_vcf_positions_by_contig(
            &args.vcf_file,
            sample_override,
            args.split_multiallelic,
            args.normalize_alleles,
            args.keep_original_info,
        )
        .map_err(reclassify_generic_as_validation)?
    };
    let annotation_format =
        io::detect_annotation_format(args.genes_file()).map_err(reclassify_generic_as_validation)?;
    let contigs = selected_contigs(args, &snp_by_contig)?;
    validate_strict_original_metrics(&contigs, &snp_by_contig, args.strict)?;
    validate_contig_inputs(&contigs, &references, &snp_by_contig, annotation_format)?;

    let preloaded_gff = if annotation_format == AnnotationFormat::Gff {
        Some(
            io::preload_gff_genes(args.genes_file(), &args.gff_features())
                .map_err(reclassify_generic_as_validation)?,
        )
    } else {
        None
    };

    let original_info_headers = if args.keep_original_info {
        if io::vcf_fast::use_fast_parser(&args.vcf_file) {
            io::vcf_fast::extract_text_info_headers(&args.vcf_file)
                .map_err(reclassify_generic_as_validation)?
        } else {
            io::extract_original_info_headers(&args.vcf_file)
                .map_err(reclassify_generic_as_validation)?
        }
    } else {
        Vec::new()
    };

    Ok(ParsedInputs {
        base_name,
        references,
        snp_by_contig,
        contigs,
        command_line: sanitized_command_line(),
        preloaded_gff,
        original_info_headers,
    })
}

pub(crate) fn reclassify_generic_as_validation(error: AppError) -> AppError {
    if error.code == ErrorCode::Generic {
        AppError::validation(error.message)
    } else {
        error
    }
}

fn apply_read_summary(variant: &mut VariantInfo, summary: ReadCountSummary) {
    variant.snp_reads = Some(summary.snp_counts);
    variant.snp_forward_reads = Some(summary.snp_forward_counts);
    variant.snp_reverse_reads = Some(summary.snp_reverse_counts);
    variant.mnv_reads = Some(summary.mnv_count);
    variant.mnv_forward_reads = Some(summary.mnv_forward_count);
    variant.mnv_reverse_reads = Some(summary.mnv_reverse_count);
    variant.mnv_total_reads = Some(summary.mnv_total_reads);
    variant.total_reads = Some(summary.total_reads);
    variant.total_forward_reads = Some(summary.total_forward_reads);
    variant.total_reverse_reads = Some(summary.total_reverse_reads);
    variant.mnv_total_forward_reads = Some(summary.mnv_total_forward_reads);
    variant.mnv_total_reverse_reads = Some(summary.mnv_total_reverse_reads);
}

fn annotate_variants_for_gene(
    gene: &Gene,
    snp_list: &[VcfPosition],
    reference: &io::Reference,
    contig: &str,
) -> Vec<VariantInfo> {
    variants::get_mnv_variants_for_gene(gene, snp_list, reference, contig)
}

fn count_gene_variant_reads(
    state: &mut WorkerState,
    args: &Args,
    contig: &str,
    gene: &Gene,
    variants: &mut [VariantInfo],
) -> AppResult<(usize, usize)> {
    if args.dry_run || state.bam.is_none() {
        return Ok((0, 0));
    }

    let mut target_positions = variants
        .iter()
        .filter(|variant| variant.variant_type != VariantType::Indel)
        .flat_map(|variant| variant.positions.iter().copied())
        .collect::<Vec<_>>();
    if target_positions.is_empty() {
        return Ok((0, 0));
    }

    target_positions.sort_unstable();
    target_positions.dedup();

    let cache_key = RegionCacheKey {
        contig: contig.to_string(),
        start: gene.start,
        end: gene.end,
        positions: target_positions.clone(),
        min_mapq: args.min_mapq,
    };

    let (cache, cache_hits, cache_misses) =
        if let Some(cached) = state.region_cache.get(&cache_key) {
            (cached.clone(), 1, 0)
        } else {
            let bam = match state.bam.as_mut() {
                Some(b) => b,
                None => {
                    return Err(AppError::validation(
                        "BAM reader unavailable in worker thread",
                    ))
                }
            };
            let built = read_count::build_region_observation_cache(
                bam,
                contig,
                gene.start,
                gene.end,
                &target_positions,
                args.min_mapq,
            )
            .map_err(|e| {
                AppError::validation(format!(
                    "Failed building read cache for contig '{}' gene '{}' at interval {}-{}: {}",
                    contig, gene.name, gene.start, gene.end, e
                ))
            })?;

            let result = built.clone();
            state.region_cache.put(cache_key, built);
            (result, 0, 1)
        };

    for variant in variants {
        if variant.variant_type == VariantType::Indel {
            continue;
        }
        let summary = read_count::count_reads_from_cache(
            &cache,
            &variant.positions,
            &variant.base_changes,
            args.min_quality,
        )
        .map_err(|e| {
            AppError::validation(format!(
                "Failed counting reads from cache for contig '{}' gene '{}' at positions {:?}: {}",
                contig, gene.name, variant.positions, e
            ))
        })?;
        apply_read_summary(variant, summary);
    }

    Ok((cache_hits, cache_misses))
}

pub(crate) fn process_contig(
    args: &Args,
    contig: &str,
    references: &ReferenceMap,
    snp_by_contig: &HashMap<String, Vec<VcfPosition>>,
    preloaded_gff: Option<&HashMap<String, Vec<Gene>>>,
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
        io::load_genes(args.genes_file(), snp_list, Some(contig), &args.gff_features())
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
                let bam = if should_count_reads {
                    if let Some(path) = bam_path.as_ref() {
                        Some(IndexedReader::from_path(path).map_err(|e| {
                            AppError::validation(format!(
                                "Failed to open BAM '{path}' in worker thread: {e}"
                            ))
                        })?)
                    } else {
                        None
                    }
                } else {
                    None
                };
                Ok(WorkerState {
                    bam,
                    region_cache: LruCache::new(NonZeroUsize::new(64).expect("64 > 0")),
                })
            },
            |state_result, gene| -> AppResult<WorkerResult> {
                let state = state_result
                    .as_mut()
                    .map_err(|err| AppError::validation(err.to_string()))?;
                let mut variants = annotate_variants_for_gene(gene, snp_list, &reference, contig);
                if variants.is_empty() {
                    return Ok(WorkerResult::default());
                }
                let (cache_hits, cache_misses) =
                    count_gene_variant_reads(state, args, contig, gene, &mut variants)?;
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
        let mut covered: HashSet<usize> = HashSet::with_capacity(snp_list.len());
        for gene in &genes {
            for snp in snp_list {
                if snp.position >= gene.start && snp.position <= gene.end {
                    covered.insert(snp.position);
                }
            }
        }
        let mut intergenic_count = 0usize;
        for snp in snp_list {
            if !covered.contains(&snp.position) {
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
