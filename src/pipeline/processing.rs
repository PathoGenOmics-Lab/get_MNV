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
use crate::read_count::{self, RegionObservationCache};
use crate::variants::{self, Gene, VariantInfo, VariantType};
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
use read_support::{
    append_supported_phased_indel_haplotypes, apply_read_summary, compute_frameshift_phasing,
    count_gene_variant_reads,
};

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
    phasing: &variants::FrameshiftPhasing,
) -> Vec<VariantInfo> {
    variants::get_mnv_variants_for_gene_with_config(
        gene,
        snp_list,
        reference,
        contig,
        genetic_code,
        indel_config,
        phasing,
    )
}

/// Count BAM read support for intergenic variants. These fall outside every
/// gene, so the per-gene counting loop never reaches them. Intergenic SNPs are
/// counted so they can be filtered by their real support like any other SNP;
/// intergenic indels are left uncounted (as elsewhere).
///
/// Nearby intergenic SNPs are grouped into bounded windows and counted from a
/// single region observation cache per window, so clustered sites share one BAM
/// query instead of one query each. The window span and position count are
/// capped so memory stays bounded when sites are dense; distant sites simply
/// fall into separate windows.
/// First gene whose intron places this variant in a splice region, with the
/// splice consequence. The variant's anchor position is tested; only multi-exon
/// transcript models have internal junctions, so single-feature (prokaryotic)
/// annotations never match.
fn splice_site_for_variant(
    genes: &[Gene],
    snp: &VcfPosition,
) -> Option<(String, crate::variants::SpliceConsequence)> {
    genes.iter().find_map(|gene| {
        variants::splice::splice_consequence_for_position(gene, snp.position)
            .map(|consequence| (gene.name.clone(), consequence))
    })
}

fn count_intergenic_variant_reads(
    args: &Args,
    contig: &str,
    variants: &mut [VariantInfo],
) -> AppResult<()> {
    const MAX_WINDOW_SPAN: usize = 50_000;
    const MAX_WINDOW_POSITIONS: usize = 256;

    let bam_path = match args.bam_file.as_ref() {
        Some(path) => path,
        None => return Ok(()),
    };

    // Indices of the intergenic SNPs, ordered by position so neighbours batch.
    let mut snp_indices: Vec<usize> = variants
        .iter()
        .enumerate()
        .filter(|(_, v)| v.variant_type == VariantType::Snp && !v.positions.is_empty())
        .map(|(i, _)| i)
        .collect();
    if snp_indices.is_empty() {
        return Ok(());
    }
    snp_indices.sort_by_key(|&i| variants[i].positions[0]);

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| {
            AppError::validation(format!(
                "Failed to open BAM '{bam_path}' for intergenic read counting: {e}"
            ))
        })?;
    let header = reader.read_header().map_err(|e| {
        AppError::validation(format!("Failed to read BAM header from '{bam_path}': {e}"))
    })?;

    let mut start = 0usize;
    while start < snp_indices.len() {
        let window_start_pos = variants[snp_indices[start]].positions[0];
        let mut end = start + 1;
        while end < snp_indices.len()
            && end - start < MAX_WINDOW_POSITIONS
            && variants[snp_indices[end]].positions[0].saturating_sub(window_start_pos)
                <= MAX_WINDOW_SPAN
        {
            end += 1;
        }

        let window = &snp_indices[start..end];
        let region_end_pos = variants[window[window.len() - 1]].positions[0];
        let positions: Vec<usize> = window.iter().map(|&i| variants[i].positions[0]).collect();

        let cache = read_count::build_region_observation_cache(
            &mut reader,
            &header,
            contig,
            window_start_pos,
            region_end_pos,
            &positions,
            args.min_mapq,
        )
        .map_err(|e| {
            AppError::validation(format!(
                "Failed building read cache for intergenic window {contig}:{window_start_pos}-{region_end_pos}: {e}"
            ))
        })?;

        for &i in window {
            let summary = read_count::count_reads_from_cache(
                &cache,
                &variants[i].positions,
                &variants[i].base_changes,
                args.min_quality,
            )
            .map_err(|e| {
                AppError::validation(format!(
                    "Failed counting reads for intergenic position {contig}:{}: {e}",
                    variants[i].positions[0]
                ))
            })?;
            apply_read_summary(&mut variants[i], summary);
        }

        start = end;
    }
    Ok(())
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
                // BAM-derived cis/trans phasing so an upstream indel does not
                // frameshift a downstream codon it is not on the same molecule as.
                let phasing = compute_frameshift_phasing(state, args, contig, gene, snp_list)?;
                let mut variants = annotate_variants_for_gene(
                    gene,
                    snp_list,
                    &reference,
                    contig,
                    genetic_code,
                    &indel_config,
                    &phasing,
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
        let mut intergenic: Vec<VariantInfo> = Vec::new();
        for (idx, snp) in snp_list.iter().enumerate() {
            if !covered[idx] {
                // A variant not in any CDS may still fall in the splice region of
                // a gene's intron; annotate it as a splice variant (gene-named)
                // instead of intergenic. It stays in this vector so its read
                // support is counted with the intergenic SNPs.
                match splice_site_for_variant(&genes, snp) {
                    Some((gene_name, splice)) => intergenic.push(variants::build_splice_variant(
                        contig, snp, &gene_name, splice,
                    )),
                    None => intergenic.push(variants::build_intergenic_variant(contig, snp)),
                }
            }
        }
        if !intergenic.is_empty() {
            let intergenic_count = intergenic.len();
            if should_count_reads {
                count_intergenic_variant_reads(args, contig, &mut intergenic)?;
            }
            all_variants.extend(intergenic);
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
