//! Per-contig parallel processing: variant annotation, BAM read counting,
//! LRU cache, and output emission.

use super::config::{
    sanitized_command_line, selected_contigs, validate_contig_inputs,
    validate_strict_original_metrics,
};
// build_input_metadata used by mod.rs, not here
use super::summary::{summarize_contig_variants, ContigSummary};
use crate::cli::Args;
use crate::cli::VariantInputFormat;
use crate::error::ErrorCode;
use crate::error::{AppError, AppResult};
use crate::io::{self, AnnotationFormat, ReferenceMap, VcfPosition};
use crate::output;
use crate::read_count::{self, ReadCountSummary, RegionObservationCache};
use crate::variants::{self, Gene, VariantInfo, VariantType};
use log::info;
use lru::LruCache;
use noodles::bam;
use noodles::sam::Header;
use rayon::prelude::*;
use std::collections::HashMap;
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
    let variant_file = args.variant_file();
    let base_name = io::get_base_name(variant_file).map_err(reclassify_generic_as_validation)?;
    let references =
        io::load_references(&args.fasta_file, args.chrom.as_deref())
            .map_err(reclassify_generic_as_validation)?;
    let input_format = resolve_variant_input_format(args)?;
    if input_format == VariantInputFormat::Tsv && sample_override.is_some() {
        return Err(AppError::config(
            "--sample is only supported for VCF input; TSV input is treated as a single sample",
        ));
    }

    let snp_by_contig = match input_format {
        VariantInputFormat::Tsv => io::ivar::load_ivar_tsv(variant_file, &references)
            .map_err(reclassify_generic_as_validation)?,
        VariantInputFormat::Vcf | VariantInputFormat::Auto => {
            // Use the fast text parser for plain .vcf files and the
            // BGZF-aware parser for .vcf.gz. BCF input is rejected with a
            // conversion hint.
            if io::vcf_fast::use_fast_parser(variant_file) {
                io::vcf_fast::load_vcf_text(
                    variant_file,
                    sample_override,
                    args.split_multiallelic,
                    args.normalize_alleles,
                    args.keep_original_info,
                )
                .map_err(reclassify_generic_as_validation)?
            } else {
                io::load_vcf_positions_by_contig(
                    variant_file,
                    sample_override,
                    args.split_multiallelic,
                    args.normalize_alleles,
                    args.keep_original_info,
                )
                .map_err(reclassify_generic_as_validation)?
            }
        }
    };
    let annotation_format = io::detect_annotation_format(args.genes_file())
        .map_err(reclassify_generic_as_validation)?;
    let contigs = selected_contigs(args, &snp_by_contig)?;
    validate_strict_original_metrics(&contigs, &snp_by_contig, args.strict)?;
    validate_contig_inputs(&contigs, &references, &snp_by_contig, annotation_format)?;

    // Auto-select the phase/splice-aware CDS model when the user did not set
    // --gff-features and the GFF actually contains CDS features; otherwise a
    // eukaryotic GFF would be analysed over whole-gene spans (introns included),
    // silently mis-numbering amino acids. Explicitly passing --gff-features
    // (e.g. `gene`) keeps the legacy behaviour and the existing phase warning.
    let gff_features = if args.gff_features_raw.is_none()
        && annotation_format == AnnotationFormat::Gff
        && io::annotation::gff_has_cds_features(args.genes_file())
            .map_err(reclassify_generic_as_validation)?
    {
        log::info!(
            "GFF contains CDS features and --gff-features was not set; analysing 'CDS' \
             (phase- and splice-aware codon model). Pass --gff-features gene to use whole-gene spans."
        );
        vec!["CDS".to_string()]
    } else {
        args.gff_features()
    };

    let preloaded_gff = if annotation_format == AnnotationFormat::Gff {
        Some(
            io::preload_gff_genes(args.genes_file(), &gff_features)
                .map_err(reclassify_generic_as_validation)?,
        )
    } else {
        None
    };

    let original_info_headers =
        if args.keep_original_info && input_format == VariantInputFormat::Vcf {
            if io::vcf_fast::use_fast_parser(variant_file) {
                io::vcf_fast::extract_text_info_headers(variant_file)
                    .map_err(reclassify_generic_as_validation)?
            } else {
                io::extract_original_info_headers(variant_file)
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

pub(crate) fn resolve_variant_input_format(args: &Args) -> AppResult<VariantInputFormat> {
    if args.tsv_file.is_some() && args.input_format == VariantInputFormat::Vcf {
        return Err(AppError::config(
            "--tsv cannot be combined with --input-format vcf",
        ));
    }

    let variant_file = args.variant_file();
    match args.effective_input_format() {
        VariantInputFormat::Auto => {
            if io::ivar::looks_like_ivar_tsv(variant_file)
                .map_err(reclassify_generic_as_validation)?
            {
                Ok(VariantInputFormat::Tsv)
            } else {
                Ok(VariantInputFormat::Vcf)
            }
        }
        explicit => Ok(explicit),
    }
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

fn variant_allele_key(variant: &VariantInfo) -> Option<(usize, String, String)> {
    Some((
        *variant.positions.first()?,
        variant.ref_bases.first()?.clone(),
        variant.base_changes.first()?.clone(),
    ))
}

fn has_variant_allele(variants: &[VariantInfo], candidate: &VariantInfo) -> bool {
    let Some(candidate_key) = variant_allele_key(candidate) else {
        return false;
    };
    variants
        .iter()
        .filter_map(variant_allele_key)
        .any(|key| key == candidate_key)
}

fn count_exact_indel_variant_reads(
    state: &mut WorkerState,
    args: &Args,
    contig: &str,
    gene: &Gene,
    variant: &mut VariantInfo,
) -> AppResult<()> {
    let bam = match state.bam.as_mut() {
        Some(b) => b,
        None => {
            return Err(AppError::validation(
                "BAM reader unavailable in worker thread",
            ))
        }
    };
    let bam_header = match state.bam_header.as_ref() {
        Some(h) => h,
        None => {
            return Err(AppError::validation(
                "BAM header unavailable in worker thread",
            ))
        }
    };
    let ref_allele = variant
        .ref_bases
        .first()
        .ok_or_else(|| {
            AppError::validation(format!(
                "Missing REF allele for indel at contig '{}' gene '{}'",
                contig, gene.name
            ))
        })?
        .clone();
    let alt_allele = variant
        .base_changes
        .first()
        .ok_or_else(|| {
            AppError::validation(format!(
                "Missing ALT allele for indel at contig '{}' gene '{}'",
                contig, gene.name
            ))
        })?
        .clone();
    let position = *variant.positions.first().ok_or_else(|| {
        AppError::validation(format!(
            "Missing position for indel at contig '{}' gene '{}'",
            contig, gene.name
        ))
    })?;
    let mut required_components = variant
        .event_components
        .iter()
        .filter_map(|label| variants::parse_component_label(label))
        .collect::<Vec<_>>();
    if required_components.is_empty() {
        required_components =
            variants::decompose_allele(position, &ref_allele, &alt_allele).components;
    }
    let request = read_count::IndelReadCountRequest {
        chrom: contig,
        position,
        ref_allele: &ref_allele,
        alt_allele: &alt_allele,
        required_components: &required_components,
        min_phred_quality: args.min_quality,
        min_mapq: args.min_mapq,
        anchor_depth: args.indel_anchor_depth,
    };
    let summary = read_count::count_indel_reads(bam, bam_header, request).map_err(|e| {
        AppError::validation(format!(
            "Failed counting indel reads for contig '{}' gene '{}' at position {}: {}",
            contig, gene.name, position, e
        ))
    })?;
    // Silent-failure guard: the locus is covered but no read reproduces the
    // exact indel allele/CIGAR. This frequently means the input indel is not
    // left-aligned the same way as the BAM (common in homopolymers/tandem
    // repeats), so the exact-anchor match yields zero support.
    if summary.mnv_total_reads > 0 && summary.mnv_count == 0 {
        log::warn!(
            "Indel at {}:{} {}>{} has {} spanning read(s) but 0 with exact CIGAR support. \
             If it lies in a homopolymer/repeat, left-align the input first \
             (e.g. `bcftools norm -f ref.fa`) so the allele matches the read alignment.",
            contig,
            position,
            ref_allele,
            alt_allele,
            summary.mnv_total_reads
        );
    }
    apply_read_summary(variant, summary);
    Ok(())
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
    let has_indels = variants
        .iter()
        .any(|variant| variant.variant_type == VariantType::Indel);
    if target_positions.is_empty() && !has_indels {
        return Ok((0, 0));
    }

    target_positions.sort_unstable();
    target_positions.dedup();

    let mut cache_hits = 0usize;
    let mut cache_misses = 0usize;
    let cache = if target_positions.is_empty() {
        None
    } else {
        let cache_key = RegionCacheKey {
            contig: contig.to_string(),
            start: gene.start,
            end: gene.end,
            positions: target_positions.clone(),
            min_mapq: args.min_mapq,
        };

        let (cache, hits, misses) = if let Some(cached) = state.region_cache.get(&cache_key) {
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
            let bam_header = match state.bam_header.as_ref() {
                Some(h) => h,
                None => {
                    return Err(AppError::validation(
                        "BAM header unavailable in worker thread",
                    ))
                }
            };
            let built = read_count::build_region_observation_cache(
                bam,
                bam_header,
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
        cache_hits += hits;
        cache_misses += misses;
        Some(cache)
    };

    for variant in variants {
        if variant.variant_type == VariantType::Indel {
            count_exact_indel_variant_reads(state, args, contig, gene, variant)?;
            continue;
        }
        let Some(cache) = cache.as_ref() else {
            continue;
        };
        let summary = read_count::count_reads_from_cache(
            cache,
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

struct PhasedIndelInputs<'a, 'r> {
    contig: &'a str,
    gene: &'a Gene,
    reference: &'a io::Reference<'r>,
    snp_list: &'a [VcfPosition],
    genetic_code: crate::genetic_code::GeneticCode,
}

fn append_supported_phased_indel_haplotypes(
    state: &mut WorkerState,
    args: &Args,
    inputs: PhasedIndelInputs<'_, '_>,
    variants: &mut Vec<VariantInfo>,
) -> AppResult<()> {
    if args.dry_run || state.bam.is_none() {
        return Ok(());
    }

    let candidates = variants::build_phased_indel_haplotype_variants(
        inputs.gene,
        inputs.snp_list,
        inputs.reference,
        inputs.contig,
        inputs.genetic_code,
    );

    for mut candidate in candidates {
        if has_variant_allele(variants, &candidate) {
            continue;
        }
        count_exact_indel_variant_reads(state, args, inputs.contig, inputs.gene, &mut candidate)?;
        let reads = candidate.mnv_reads.unwrap_or(0);
        let depth = candidate.mnv_total_reads.unwrap_or(0);
        let freq = if depth > 0 {
            reads as f64 / depth as f64
        } else {
            0.0
        };
        // Defaults (min_reads = 1, min_freq = 0.0) reproduce the historical
        // "emit if any read supports it" rule; raising either suppresses
        // low-confidence phased haplotypes from dense local windows.
        if reads >= args.phased_indel_min_reads && freq >= args.phased_indel_min_freq {
            variants.push(candidate);
        }
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
