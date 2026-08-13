//! Per-gene BAM read support: SNP/MNV cache counting, exact indel counting,
//! and supported phased-indel haplotype emission.

use super::{PhasedIndelInputs, RegionCacheKey, WorkerState};
use crate::cli::Args;
use crate::error::{AppError, AppResult};
use crate::io::VcfPosition;
use crate::read_count::{self, ReadCountSummary};
use crate::variants::{self, Gene, VariantInfo, VariantType};
use noodles::bam;
use std::collections::HashMap;

/// Maximum genomic distance between an indel and a SNV for read-based phasing to
/// be attempted. Beyond a single read's reach no record spans both loci, so the
/// query would be wasted; those pairs stay frequency-gated.
const MAX_PHASING_SPAN: usize = 600;

/// How many distinct read-observed combinations one window may report. Real
/// mixtures rarely approach this; the bound is there so a noisy pile-up cannot
/// turn one window into thousands of rows.
const MAX_LOCAL_HAPLOTYPES_PER_WINDOW: usize = 64;

/// What the reads said about the upstream indels of a gene: how each one is
/// phased with the downstream SNVs, and what frequency it is actually at.
///
/// The frequency travels with the phasing because both answer the same
/// question, whether this indel has any business shifting the frame of a
/// downstream codon, and both are read from the same BAM pass over the gene.
#[derive(Default)]
pub(super) struct FrameshiftEvidence {
    pub phasing: variants::FrameshiftPhasing,
    /// Keyed by `(position, REF, ALT)`, as the annotation config expects.
    pub observed_indel_freq: HashMap<(usize, String, String), f64>,
}

pub(super) fn compute_frameshift_phasing(
    state: &mut WorkerState,
    args: &Args,
    contig: &str,
    gene: &Gene,
    snp_list: &[VcfPosition],
) -> AppResult<FrameshiftEvidence> {
    if args.dry_run || state.bam.is_none() {
        return Ok(FrameshiftEvidence::default());
    }

    let mut indels: Vec<&VcfPosition> = Vec::new();
    let mut snvs: Vec<(usize, char)> = Vec::new();
    for variant in snp_list
        .iter()
        .filter(|variant| variant.overlaps_interval(gene.start, gene.end))
    {
        if variant.event().class.has_indel_component() {
            // A symbolic ALT (`<DEL>`, `<DUP>`) has no sequence for a read to
            // reproduce, so `observed.allele == alt` can never hold and every
            // spanning molecule would score as informative-and-not-cis, which
            // reads as proof of trans and silently cancels the frameshift the
            // SV should propagate. The reads cannot answer here; leaving the
            // pair out says so.
            if !variant.alt_allele.starts_with('<') {
                indels.push(variant);
            }
        }
        for component in variant.substitution_components() {
            if let Some(alt) = component.alt_allele.chars().next() {
                snvs.push((component.position, alt));
            }
        }
    }
    if indels.is_empty() || snvs.is_empty() {
        return Ok(FrameshiftEvidence::default());
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

    let mut pairs: HashMap<(usize, usize, char), variants::PairLinkage> = HashMap::new();
    for indel in &indels {
        for &(snv_position, snv_alt) in &snvs {
            let span = indel.position.abs_diff(snv_position);
            if span > MAX_PHASING_SPAN {
                continue;
            }
            let linkage = read_count::indel_snv_linkage(
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
                !args.count_mates_separately,
            )?;
            // The cis/trans thresholds belong to the read counting that applies
            // them; the judged answer travels on from here.
            let verdict = if !linkage.is_informative() {
                variants::LinkageVerdict::Unknown
            } else if linkage.is_trans() {
                variants::LinkageVerdict::Trans
            } else {
                variants::LinkageVerdict::Cis
            };
            pairs.insert(
                (indel.position, snv_position, snv_alt.to_ascii_uppercase()),
                variants::PairLinkage {
                    verdict,
                    cis_reads: linkage.cis_reads,
                    informative_reads: linkage.informative_reads,
                },
            );
        }
    }

    // The frequency the reads give this indel, which is what the frameshift gate
    // should weigh. Counted here because the gate runs during annotation, before
    // the per-variant read counting, and this is already the pass that opens the
    // BAM for these same indels.
    let mut observed_indel_freq = HashMap::new();
    for indel in &indels {
        let components =
            variants::decompose_allele(indel.position, &indel.ref_allele, &indel.alt_allele)
                .components;
        let request = read_count::IndelReadCountRequest {
            chrom: contig,
            position: indel.position,
            ref_allele: &indel.ref_allele,
            alt_allele: &indel.alt_allele,
            required_components: &components,
            min_phred_quality: args.min_quality,
            min_mapq: args.min_mapq,
            anchor_depth: args.indel_anchor_depth,
            pair_aware: !args.count_mates_separately,
        };
        let summary = read_count::count_indel_reads(bam, header, request)?;
        // No depth is no measurement. Leaving the entry out lets the gate fall
        // back to what the caller declared instead of reading it as 0.0, which
        // would suppress every propagation from an uncovered indel.
        if summary.mnv_total_reads > 0 {
            observed_indel_freq.insert(
                (
                    indel.position,
                    indel.ref_allele.clone(),
                    indel.alt_allele.clone(),
                ),
                summary.mnv_count as f64 / summary.mnv_total_reads as f64,
            );
        }
    }

    Ok(FrameshiftEvidence {
        phasing: variants::FrameshiftPhasing::from_pairs(pairs),
        observed_indel_freq,
    })
}

/// Carry the input caller's own phase claim onto the multi-position rows, and
/// note where the reads flatly refute it. Run after read counting so both
/// sides of that comparison are available.
pub(super) fn annotate_declared_phase(variants: &mut [VariantInfo], snp_list: &[VcfPosition]) {
    for variant in variants.iter_mut() {
        let Some(mut call) = variants::declared_phase::declared_phase_for_row(
            &variant.positions,
            &variant.base_changes,
            snp_list,
        ) else {
            continue;
        };
        if let (Some(haplotype_reads), Some(informative_reads)) =
            (variant.mnv_reads, variant.mnv_phasing_reads)
        {
            call.contradicted_by_reads = variants::declared_phase::reads_contradict_declared_phase(
                &call,
                haplotype_reads,
                informative_reads,
            );
            if call.contradicted_by_reads {
                log::warn!(
                    "{}:{:?} was declared {} by the input VCF, but {} of {} reads spanning the \
                     site carry the whole haplotype.",
                    variant.chrom,
                    variant.positions,
                    call.verdict.as_str(),
                    haplotype_reads,
                    informative_reads
                );
            }
        }
        variant.annotations.declared_phase = Some(call);
    }
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
    // Reads that span the whole row and carry the rarest constituent SNV: the
    // full-haplotype reads plus the spanning reads that carry that SNV alone.
    // The least-supported SNV is the one with the fewest such reads, and since
    // the full-haplotype reads are common to every position, that is the one
    // with the fewest solo reads.
    variant.mnv_phasing_reads = summary
        .snp_only_informative_counts
        .iter()
        .copied()
        .min()
        .map(|rarest_solo| summary.mnv_count + rarest_solo);
    // How much the substitutions co-occur beyond what their own frequencies
    // predict. The co-occurrence ratio above cannot make that distinction: two
    // alleles on 90% of molecules each meet on 81% of them by arithmetic alone.
    variant.annotations.linkage =
        variants::linkage::codon_linkage(&summary.haplotype_patterns, variant.positions.len());
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
    count_indel_reads_into(bam, bam_header, args, contig, &gene.name, variant)
}

/// Count one indel's exact event support and write it onto the row.
///
/// Split out from the gene path so the intergenic path can use it too: an indel
/// outside every feature used to reach no counter at all, and the row then
/// claimed `Event Reads = 0` at `Event Depth = 0` for an allele every read
/// carried. `context` only labels errors and the warning.
pub(super) fn count_indel_reads_into(
    bam: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    bam_header: &noodles::sam::Header,
    args: &Args,
    contig: &str,
    context: &str,
    variant: &mut VariantInfo,
) -> AppResult<()> {
    let ref_allele = variant
        .ref_bases
        .first()
        .ok_or_else(|| {
            AppError::validation(format!(
                "Missing REF allele for indel at contig '{contig}' {context}"
            ))
        })?
        .clone();
    let alt_allele = variant
        .base_changes
        .first()
        .ok_or_else(|| {
            AppError::validation(format!(
                "Missing ALT allele for indel at contig '{contig}' {context}"
            ))
        })?
        .clone();
    let position = *variant.positions.first().ok_or_else(|| {
        AppError::validation(format!(
            "Missing position for indel at contig '{contig}' {context}"
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
        pair_aware: !args.count_mates_separately,
    };
    let summary = read_count::count_indel_reads(bam, bam_header, request).map_err(|e| {
        AppError::validation(format!(
            "Failed counting indel reads for contig '{contig}' {context} at position {position}: {e}"
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
            !args.count_mates_separately,
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

/// Ask the reads which combinations of a locally-linked group actually travel
/// together. One BAM query per window, whatever the group size: the answer is
/// read off the molecules rather than proposed and checked combination by
/// combination.
fn observe_component_haplotypes(
    state: &mut WorkerState,
    args: &Args,
    contig: &str,
    component: &[&VcfPosition],
) -> AppResult<read_count::LocalHaplotypeObservations> {
    let start = component
        .iter()
        .map(|variant| variant.position)
        .min()
        .unwrap_or(0);
    let end = component
        .iter()
        .map(|variant| variant.position + variant.ref_allele.len().saturating_sub(1))
        .max()
        .unwrap_or(0);
    let component_variants = component
        .iter()
        .map(|variant| {
            let components = variant.event().components;
            read_count::LocalVariant {
                start: variant.position,
                ref_len: variants::observation_ref_len(&variant.ref_allele, &components),
                components,
            }
        })
        .collect::<Vec<_>>();

    let bam = state
        .bam
        .as_mut()
        .ok_or_else(|| AppError::validation("BAM reader unavailable in worker thread"))?;
    let header = state
        .bam_header
        .as_ref()
        .ok_or_else(|| AppError::validation("BAM header unavailable in worker thread"))?;

    let observations = read_count::observe_local_haplotypes(
        bam,
        header,
        read_count::LocalHaplotypeRequest {
            chrom: contig,
            start,
            end,
            variants: &component_variants,
            min_phred_quality: args.min_quality,
            min_mapq: args.min_mapq,
            pair_aware: !args.count_mates_separately,
        },
    )?;

    // A read that stops inside the window witnessed part of a molecule, not a
    // haplotype. Say so when that is where most of the coverage went, since it
    // means the window is wider than the library's fragments.
    if observations.partial_reads > observations.spanning_reads {
        log::debug!(
            "Local haplotype window {}:{}-{} is spanned by {} read(s) but only partly covered by \
             {}; combinations are called from the spanning reads alone.",
            contig,
            start,
            end,
            observations.spanning_reads,
            observations.partial_reads
        );
    }

    Ok(observations)
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

    for component in variants::local_haplotype_components(inputs.gene, inputs.snp_list) {
        let observed = observe_component_haplotypes(state, args, inputs.contig, &component)?;
        let mut haplotypes = observed.haplotypes;
        // Reading the molecules imposes no combinatorial limit, but sequencing
        // error at a called position can still mint a combination of its own.
        // Keep the best-supported ones and say what was set aside rather than
        // quietly truncating.
        haplotypes.sort_by_key(|haplotype| std::cmp::Reverse(haplotype.reads));
        if haplotypes.len() > MAX_LOCAL_HAPLOTYPES_PER_WINDOW {
            let dropped = haplotypes.len() - MAX_LOCAL_HAPLOTYPES_PER_WINDOW;
            let weakest = haplotypes[MAX_LOCAL_HAPLOTYPES_PER_WINDOW].reads;
            log::warn!(
                "Local haplotype window in gene '{}' shows {} distinct combinations on the reads; \
                 reporting the {} best supported and dropping {} at {} read(s) or fewer.",
                inputs.gene.name,
                haplotypes.len(),
                MAX_LOCAL_HAPLOTYPES_PER_WINDOW,
                dropped,
                weakest
            );
            haplotypes.truncate(MAX_LOCAL_HAPLOTYPES_PER_WINDOW);
        }
        // The joint distribution over the window, which is what the linkage
        // statistics read. Discovery already built it to decide which
        // combinations exist; nothing else has to be counted.
        let window_patterns = haplotypes
            .iter()
            .map(|haplotype| (haplotype.carried.clone(), haplotype.reads))
            .collect::<Vec<_>>();

        for haplotype in haplotypes {
            let carried_indices = haplotype
                .carried
                .iter()
                .enumerate()
                .filter_map(|(index, carried)| carried.then_some(index))
                .collect::<Vec<_>>();
            let group = component
                .iter()
                .zip(haplotype.carried.iter())
                .filter_map(|(variant, carried)| carried.then_some(*variant))
                .collect::<Vec<_>>();
            let Some(mut candidate) = variants::phased_haplotype_variant(
                inputs.gene,
                inputs.reference,
                inputs.contig,
                &group,
                inputs.genetic_code,
            ) else {
                continue;
            };
            if has_variant_allele(variants, &candidate) {
                continue;
            }
            count_exact_indel_variant_reads(
                state,
                args,
                inputs.contig,
                inputs.gene,
                &mut candidate,
            )?;
            // The exact count answers "can a molecule reproduce this allele",
            // which is the CIGAR-aware validity check and matters for
            // net-neutral combinations. It cannot answer "how many molecules
            // are this haplotype": it matches over the combination's own span,
            // so a molecule carrying that combination *and more* matches too,
            // and a two-variant row inside a real three-variant species was
            // reported with the species' whole read count at 100% frequency.
            // Which molecules are this combination is what discovery decided,
            // over every variant of the window.
            if candidate.mnv_reads.unwrap_or(0) == 0 {
                continue;
            }
            let depth = candidate.mnv_total_reads.unwrap_or(0);
            candidate.mnv_reads = Some(haplotype.reads);
            candidate.mnv_forward_reads = Some(haplotype.forward_reads);
            candidate.mnv_reverse_reads = Some(haplotype.reverse_reads);
            // Whether this row's own variants travel together, measured the same
            // way as for a codon MNV and over the same molecules its counts come
            // from. Asked about the variants this row carries, not the whole
            // window, since that is what the row claims.
            candidate.annotations.linkage = variants::linkage::codon_linkage(
                &variants::linkage::restricted_to(&window_patterns, &carried_indices),
                carried_indices.len(),
            );
            let reads = haplotype.reads;
            let freq = if depth > 0 {
                reads as f64 / depth as f64
            } else {
                0.0
            };
            if reads >= args.phased_indel_min_reads && freq >= args.phased_indel_min_freq {
                variants.push(candidate);
            }
        }
    }

    Ok(())
}
