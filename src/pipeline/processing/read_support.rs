//! Per-gene BAM read support: SNP/MNV cache counting, exact indel counting,
//! and supported phased-indel haplotype emission.

use super::{PhasedIndelInputs, RegionCacheKey, WorkerState};
use crate::cli::Args;
use crate::error::{AppError, AppResult};
use crate::io::VcfPosition;
use crate::read_count::{self, ReadCountSummary};
use crate::variants::{self, Gene, VariantInfo, VariantType};
use std::collections::HashSet;

/// Maximum genomic distance between an indel and a SNV for read-based phasing to
/// be attempted. Beyond a single read's reach no record spans both loci, so the
/// query would be wasted; those pairs stay frequency-gated.
const MAX_PHASING_SPAN: usize = 600;

/// BAM-derived phasing for frameshift propagation: the set of `(indel, SNV)`
/// pairs the reads show are in trans. Computed before annotation so it can be
/// threaded into the frameshift gate. Empty without a BAM (or in dry-run), which
/// preserves the pure frequency-gated behaviour.
pub(super) fn compute_frameshift_phasing(
    state: &mut WorkerState,
    args: &Args,
    contig: &str,
    gene: &Gene,
    snp_list: &[VcfPosition],
) -> AppResult<variants::FrameshiftPhasing> {
    if args.dry_run || state.bam.is_none() {
        return Ok(variants::FrameshiftPhasing::default());
    }

    let mut indels: Vec<&VcfPosition> = Vec::new();
    let mut snvs: Vec<(usize, char)> = Vec::new();
    for variant in snp_list
        .iter()
        .filter(|variant| variant.overlaps_interval(gene.start, gene.end))
    {
        if variant.event().class.has_indel_component() {
            indels.push(variant);
        }
        for component in variant.substitution_components() {
            if let Some(alt) = component.alt_allele.chars().next() {
                snvs.push((component.position, alt));
            }
        }
    }
    if indels.is_empty() || snvs.is_empty() {
        return Ok(variants::FrameshiftPhasing::default());
    }
    snvs.sort_unstable();
    snvs.dedup();

    let bam = state
        .bam
        .as_mut()
        .ok_or_else(|| AppError::validation("BAM reader unavailable in worker thread"))?;
    let header = state
        .bam_header
        .as_ref()
        .ok_or_else(|| AppError::validation("BAM header unavailable in worker thread"))?;

    let mut trans_pairs: HashSet<(usize, usize)> = HashSet::new();
    for indel in &indels {
        for &(snv_position, snv_alt) in &snvs {
            let span = indel.position.abs_diff(snv_position);
            if span > MAX_PHASING_SPAN {
                continue;
            }
            if read_count::indel_snv_in_trans(
                bam,
                header,
                contig,
                indel.position,
                &indel.ref_allele,
                &indel.alt_allele,
                snv_position,
                snv_alt,
                args.min_quality,
                args.min_mapq,
            )? {
                trans_pairs.insert((indel.position, snv_position));
            }
        }
    }

    Ok(variants::FrameshiftPhasing::from_trans_pairs(trans_pairs))
}

pub(super) fn apply_read_summary(variant: &mut VariantInfo, summary: ReadCountSummary) {
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

pub(super) fn variant_allele_key(variant: &VariantInfo) -> Option<(usize, String, String)> {
    Some((
        *variant.positions.first()?,
        variant.ref_bases.first()?.clone(),
        variant.base_changes.first()?.clone(),
    ))
}

pub(super) fn has_variant_allele(variants: &[VariantInfo], candidate: &VariantInfo) -> bool {
    let Some(candidate_key) = variant_allele_key(candidate) else {
        return false;
    };
    variants
        .iter()
        .filter_map(variant_allele_key)
        .any(|key| key == candidate_key)
}

pub(super) fn count_exact_indel_variant_reads(
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

pub(super) fn count_gene_variant_reads(
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

pub(super) fn append_supported_phased_indel_haplotypes(
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
