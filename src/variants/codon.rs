//! Codon-level variant processing: SNP grouping, MNV detection, and amino
//! acid change calculation.
//!
//! The heavy lifting lives in focused submodules; this file wires them together
//! and hosts the public per-gene / per-transcript entry points.

mod allele_apply;
mod config;
mod gene_path;
mod grouping;
mod indel_effect;
mod nmd;
mod phased;
mod protein;
mod transcript_model;
mod transcript_path;

#[cfg(test)]
mod tests;

pub use config::{
    pair_key, FrameshiftPairKey, FrameshiftPhasing, IndelAnnotationConfig, PairLinkage,
};
pub use gene_path::process_codon;
pub use phased::{local_haplotype_components, phased_haplotype_variant};

use crate::variants::hgvs::CodingPosition;
use crate::variants::types::*;
use std::collections::BTreeMap;

use config::indel_passes_frameshift_gate;
use gene_path::{codon_bounds_for_position, effective_bounds};

// Re-exported so the property tests can check the codon-bounds invariants over
// generated gene shapes rather than one hard-coded feature.
#[cfg(test)]
pub(crate) use gene_path::codon_bounds_for_position as codon_bounds;
use grouping::{merge_snp_into_groups, variant_event_metadata};
use indel_effect::{coding_delta_for_variant, protein_effect_for_indel};
use transcript_model::{
    has_transcript_cds_model, indel_disturbs_codon, indel_is_upstream_of_codon_offset,
    transcript_offsets_for_position, transcript_sequence_for_gene, variant_overlaps_coding_model,
};
use transcript_path::{
    apply_frameshift_labeling, frameshift_ptc_protein_pos, merge_transcript_snp_into_groups,
    process_transcript_codon, transcript_oriented_base, TranscriptSnp,
};

/// The 1-based CDS position of a genomic base, measured the way the substitution
/// paths measure it for `c.`: along the spliced transcript when the gene has a
/// CDS model, otherwise from the phase-adjusted feature start (or end, on the
/// minus strand) plus the feature's protein offset.
fn cds_position(gene: &Gene, position: usize) -> Option<usize> {
    if has_transcript_cds_model(gene) {
        return transcript_offsets_for_position(gene, position)
            .first()
            .map(|offset| offset + 1);
    }
    let (start, end) = effective_bounds(gene);
    if position < start || position > end {
        return None;
    }
    let local = match gene.strand {
        Strand::Plus => position - start,
        Strand::Minus => end - position,
    };
    Some(gene.protein_offset * 3 + local + 1)
}

/// The HGVS coding position of a genomic base: its CDS position, or, in an
/// intron of a spliced transcript, the nearer exon's boundary base and the
/// distance from it (`c.100+2`, `c.101-5`). `None` outside the transcript's
/// coding exons and introns.
fn coding_coordinate(gene: &Gene, position: usize) -> Option<CodingPosition> {
    if let Some(n) = cds_position(gene, position) {
        return Some(CodingPosition::Cds(n));
    }
    if !has_transcript_cds_model(gene) {
        return None;
    }
    // Segments are in transcript order, so each pair brackets one intron with
    // the upstream exon first.
    for pair in gene.cds_segments.windows(2) {
        let (up, down) = (&pair[0], &pair[1]);
        if !gene.intron_separates(up, down) {
            continue;
        }
        let (up_last, down_first, d_up, d_down) = match gene.strand {
            Strand::Plus if position > up.end && position < down.start => {
                (up.end, down.start, position - up.end, down.start - position)
            }
            Strand::Minus if position < up.start && position > down.end => {
                (up.start, down.end, up.start - position, position - down.end)
            }
            _ => continue,
        };
        // HGVS splits an intron down the middle; the middle base of an odd one
        // belongs to the first half.
        return Some(if d_up <= d_down {
            CodingPosition::Intron {
                anchor: cds_position(gene, up_last)?,
                offset: d_up as isize,
            }
        } else {
            CodingPosition::Intron {
                anchor: cds_position(gene, down_first)?,
                offset: -(d_down as isize),
            }
        });
    }
    None
}

/// The coding sequence numbered as `c.` numbers it, for moving an indel 3'-most
/// along it: the spliced transcript, or a feature that starts the protein. A CDS
/// row further into the protein is numbered from the protein's start, which its
/// own sequence does not begin at.
fn numbered_coding_sequence(gene: &Gene, reference: &crate::io::Reference<'_>) -> Option<String> {
    if has_transcript_cds_model(gene) || gene.protein_offset == 0 {
        transcript_model::coding_sequence_for_gene(gene, reference)
    } else {
        None
    }
}

/// HGVS `c.` for an indel in this gene, or `None` when it reaches a UTR or the
/// flank.
pub(super) fn indel_hgvs_c(
    gene: &Gene,
    indel: &crate::io::VcfPosition,
    reference: &crate::io::Reference<'_>,
) -> Option<String> {
    let insertion = crate::variants::hgvs::is_pure_insertion(&indel.ref_allele, &indel.alt_allele);
    let bases = numbered_coding_sequence(gene, reference);
    let cds = bases
        .as_deref()
        .map(|bases| crate::variants::hgvs::CodingSequence {
            bases,
            exon_ends: exon_ends(gene),
        });
    crate::variants::hgvs::coding_indel(
        indel.record_start,
        &indel.ref_allele,
        &indel.alt_allele,
        gene.strand,
        |position| {
            // An insertion after the last base of an internal exon lands at the
            // exon-exon junction of the spliced CDS, which is how the protein
            // annotation reads it too. The base after the anchor is intronic, so
            // its flank is the neighbouring exon's base, one CDS position away.
            if insertion {
                if let Some(anchor) = position.checked_sub(1) {
                    if gene.insertion_at_internal_junction(anchor) {
                        let at = cds_position(gene, anchor)?;
                        return match gene.strand {
                            Strand::Plus => Some(CodingPosition::Cds(at + 1)),
                            Strand::Minus => at.checked_sub(1).map(CodingPosition::Cds),
                        };
                    }
                }
            }
            coding_coordinate(gene, position)
        },
        cds.as_ref(),
    )
}

/// The `c.` number of the last base of every exon an intron follows, in
/// transcript order. A join with no intron (ribosomal slippage) is one reading
/// frame, and the 3' rule may cross it.
fn exon_ends(gene: &Gene) -> Vec<usize> {
    if !has_transcript_cds_model(gene) {
        return Vec::new();
    }
    gene.cds_segments
        .windows(2)
        .filter(|pair| gene.intron_separates(&pair[0], &pair[1]))
        .filter_map(|pair| {
            let last = match gene.strand {
                Strand::Plus => pair[0].end,
                Strand::Minus => pair[0].start,
            };
            cds_position(gene, last)
        })
        .collect()
}

/// The genomic bases this gene object codes, in genomic order: its CDS
/// segments, or its phase-adjusted feature interval.
fn coding_genomic_positions(gene: &Gene) -> Vec<usize> {
    if has_transcript_cds_model(gene) {
        let mut positions: Vec<usize> = gene
            .cds_segments
            .iter()
            .flat_map(|segment| segment.start..=segment.end)
            .collect();
        positions.sort_unstable();
        positions.dedup();
        positions
    } else {
        let (start, end) = effective_bounds(gene);
        (start..=end).collect()
    }
}

/// A gene object's `c.` numbering, for a caller that asks it many times: plain
/// arithmetic for a single feature, which is every gene of a bacterial genome,
/// and a table for a spliced transcript.
pub enum CodingNumbering {
    /// `first` is the `c.` number of the feature's first base in transcript
    /// order, which sits at `start` on the plus strand and at `end` on the minus.
    Linear {
        start: usize,
        end: usize,
        first: usize,
        strand: Strand,
    },
    Spliced(std::collections::BTreeMap<usize, usize>),
}

impl CodingNumbering {
    pub fn of(gene: &Gene) -> Self {
        if has_transcript_cds_model(gene) {
            return Self::Spliced(
                coding_genomic_positions(gene)
                    .into_iter()
                    .filter_map(|position| Some((cds_position(gene, position)?, position)))
                    .collect(),
            );
        }
        let (start, end) = effective_bounds(gene);
        Self::Linear {
            start,
            end,
            first: gene.protein_offset * 3 + 1,
            strand: gene.strand,
        }
    }

    /// The genomic position of coding base `c.n`, if this object holds it.
    pub fn genomic(&self, n: usize) -> Option<usize> {
        match self {
            Self::Linear {
                start,
                end,
                first,
                strand,
            } => {
                let span = end.checked_sub(*start)?;
                let local = n.checked_sub(*first).filter(|&l| l <= span)?;
                Some(match strand {
                    Strand::Plus => start + local,
                    Strand::Minus => end - local,
                })
            }
            Self::Spliced(table) => table.get(&n).copied(),
        }
    }

    /// The highest coding base this object holds, as `c.` number and genomic
    /// position.
    pub fn last(&self) -> Option<(usize, usize)> {
        match self {
            Self::Linear {
                start,
                end,
                first,
                strand,
            } => Some((
                first + end.checked_sub(*start)?,
                match strand {
                    Strand::Plus => *end,
                    Strand::Minus => *start,
                },
            )),
            Self::Spliced(table) => table.last_key_value().map(|(&n, &p)| (n, p)),
        }
    }
}

/// The genomic position of coding base `c.n` of this gene object, or `None`
/// when it does not hold that base (a CDS row further into the protein, say).
pub fn genomic_for_cds(gene: &Gene, n: usize) -> Option<usize> {
    coding_genomic_positions(gene)
        .into_iter()
        .find(|&position| cds_position(gene, position) == Some(n))
}

/// The highest coding position this gene object holds: the protein's last
/// coding base when the object ends the protein.
pub fn last_cds_position(gene: &Gene) -> Option<usize> {
    coding_genomic_positions(gene)
        .into_iter()
        .filter_map(|position| cds_position(gene, position))
        .max()
}

/// The intron of `gene` that holds `position`, as its first and last genomic
/// base, or `None` when the position is in no intron of it.
pub fn intron_span(gene: &Gene, position: usize) -> Option<(usize, usize)> {
    if !has_transcript_cds_model(gene) {
        return None;
    }
    gene.cds_segments.windows(2).find_map(|pair| {
        let (up, down) = (&pair[0], &pair[1]);
        if !gene.intron_separates(up, down) {
            return None;
        }
        let (lo, hi) = match gene.strand {
            Strand::Plus => (up.end + 1, down.start - 1),
            Strand::Minus => (down.end + 1, up.start - 1),
        };
        (lo..=hi).contains(&position).then_some((lo, hi))
    })
}

/// HGVS `c.` for a row with no codon of its own inside `gene`: a change in an
/// intron or at a splice site, numbered from the nearer exon. A deletion or
/// insertion wholly inside one intron is placed 3'-most along the transcript
/// within it, and an insertion repeating the bases before it is a `dup`, as in
/// the coding sequence.
pub fn noncodon_hgvs_c(gene: &Gene, variant: &VariantInfo, sequence: &str) -> Option<String> {
    use crate::variants::hgvs;
    let map = |position: usize| coding_coordinate(gene, position);
    // A row of several substituted bases describes them all, as its `g.` does:
    // spelled as one allele over their span, it decomposes into every one.
    if variant.positions.len() > 1 {
        let (&first, &last) = (
            variant.positions.iter().min()?,
            variant.positions.iter().max()?,
        );
        let mut ref_span = sequence
            .get(first - 1..last)?
            .to_ascii_uppercase()
            .into_bytes();
        let mut alt_span = ref_span.clone();
        for ((&position, r), a) in variant
            .positions
            .iter()
            .zip(&variant.ref_bases)
            .zip(&variant.base_changes)
        {
            let (&[r], &[a]) = (r.as_bytes(), a.as_bytes()) else {
                return None;
            };
            ref_span[position - first] = r.to_ascii_uppercase();
            alt_span[position - first] = a.to_ascii_uppercase();
        }
        let (ref_span, alt_span) = (
            String::from_utf8(ref_span).ok()?,
            String::from_utf8(alt_span).ok()?,
        );
        return hgvs::coding_noncodon(first, &ref_span, &alt_span, gene.strand, map);
    }
    let (pos, ref_allele, alt_allele) = (
        *variant.positions.first()?,
        variant.ref_bases.first()?,
        variant.base_changes.first()?,
    );
    let (start, ref_core, alt_core) = hgvs::indel_core(pos, ref_allele, alt_allele);
    let span = if alt_core.is_empty() && !ref_core.is_empty() {
        let end = start + ref_core.len() - 1;
        intron_span(gene, start).filter(|&(lo, hi)| start >= lo && end <= hi)
    } else if ref_core.is_empty() && !alt_core.is_empty() {
        // Both bases the insertion falls between are intronic.
        intron_span(gene, start).filter(|&(lo, hi)| start > lo && start <= hi)
    } else {
        None
    };
    let Some((lo, hi)) = span else {
        return hgvs::coding_noncodon(pos, ref_allele, alt_allele, gene.strand, map);
    };
    let (pos, ref_allele, alt_allele) =
        hgvs::shift_within(sequence, pos, ref_allele, alt_allele, gene.strand, lo, hi);
    if ref_allele.len() == 1 && alt_allele.len() > 1 {
        // The bases inserted after `pos`, and the copy they would duplicate:
        // the ones before them along the transcript.
        let inserted = &alt_allele[1..];
        let length = inserted.len();
        let copy = match gene.strand {
            Strand::Plus => (pos + 1).checked_sub(length).map(|first| (first, pos)),
            Strand::Minus => Some((pos + 1, pos + length)),
        };
        if let Some((first, last)) =
            copy.filter(|&(first, last)| first >= 1 && last <= sequence.len())
        {
            if sequence[first - 1..last].eq_ignore_ascii_case(inserted) {
                let (a, b) = (map(first)?, map(last)?);
                if a == b {
                    return Some(format!("c.{a}dup"));
                }
                let (low, high) = match gene.strand {
                    Strand::Plus => (a, b),
                    Strand::Minus => (b, a),
                };
                return Some(format!("c.{low}_{high}dup"));
            }
        }
    }
    hgvs::coding_noncodon(pos, &ref_allele, &alt_allele, gene.strand, map)
}

/// Whether the whole of this indel lies upstream of the codon in transcript
/// order, on a gene annotated in genomic coordinates.
///
/// Asked of the bases the record changes, not of the base it starts on. A
/// left-anchored deletion names the base before the gap, so on the minus strand,
/// where the transcript reads from higher coordinates down, a deletion anchored
/// exactly on the codon's last genomic base removes bases that come before that
/// codon and its frame shift has to reach it. Testing the anchor said otherwise,
/// and the codon kept a plain missense label on a frame that had shifted. An
/// insertion places its bases after the anchor, so the gap it opens lies between
/// the anchor and the next base, which is what makes the two strands read their
/// boundary differently.
fn indel_is_upstream_of_codon(
    indel: &crate::io::VcfPosition,
    strand: Strand,
    codon_start: usize,
    codon_end: usize,
) -> bool {
    let components = indel.event().components;
    if components.is_empty() {
        return false;
    }
    components.iter().all(|component| match component.kind {
        crate::variants::AlleleComponentKind::Insertion => match strand {
            Strand::Plus => component.position < codon_start,
            Strand::Minus => component.position >= codon_end,
        },
        _ => {
            let last = component.position + component.ref_allele.len().saturating_sub(1);
            match strand {
                Strand::Plus => last < codon_start,
                Strand::Minus => component.position > codon_end,
            }
        }
    })
}

fn get_mnv_variants_for_transcript(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
    config: &IndelAnnotationConfig,
    phasing: &FrameshiftPhasing,
) -> Vec<VariantInfo> {
    let Some(ref_cds) = transcript_sequence_for_gene(gene, reference) else {
        return Vec::new();
    };

    let mut variants = Vec::new();
    // See merge_snp_into_groups for the multi-allelic alternative interpretation
    // logic. Each codon_start maps to a list of mutually-exclusive groups.
    let mut codon_to_groups: BTreeMap<usize, Vec<Vec<TranscriptSnp>>> = BTreeMap::new();
    let mut indels: Vec<crate::io::VcfPosition> = Vec::new();

    for variant in snp_list
        .iter()
        .filter(|variant| variant_overlaps_coding_model(gene, variant))
    {
        let substitutions = variant.substitution_components();
        if !substitutions.is_empty() {
            for component in substitutions {
                // A base shared by two overlapping CDS rows (a ribosomal-slippage
                // join) is read twice, so it belongs to two codons and must be
                // applied to both. Ordinary models yield exactly one offset here.
                for offset in transcript_offsets_for_position(gene, component.position) {
                    let Some(transcript_alt_base) =
                        transcript_oriented_base(&component.alt_allele, gene.strand)
                    else {
                        continue;
                    };
                    let codon_start = (offset / 3) * 3;
                    let new_snp = TranscriptSnp {
                        snp: Snp {
                            index: component.position,
                            position: component.position,
                            ref_base: component.ref_allele.clone(),
                            base: component.alt_allele.clone(),
                            original_dp: variant.original_dp,
                            original_freq: variant.original_freq,
                            original_info: variant.original_info.clone(),
                        },
                        transcript_offset: offset,
                        transcript_alt_base,
                    };
                    let groups = codon_to_groups.entry(codon_start).or_default();
                    merge_transcript_snp_into_groups(groups, new_snp);
                }
            }
        } else {
            indels.push(variant.clone());
        }
    }

    let placed = coding_placements(gene, reference, &indels);
    let ptc_protein_pos =
        frameshift_ptc_protein_pos(gene, reference, &placed, genetic_code, config);

    let mut codon_starts = codon_to_groups.keys().copied().collect::<Vec<_>>();
    codon_starts.sort_unstable();

    for codon_start in codon_starts {
        let groups = codon_to_groups
            .get(&codon_start)
            .cloned()
            .unwrap_or_default();
        for codon_snps in groups {
            let codon_end = codon_start + 3;
            if codon_end > ref_cds.len() {
                continue;
            }
            let ref_codon = &ref_cds[codon_start..codon_end];
            let mut codon_snps = codon_snps;
            codon_snps.sort_by_key(|s| s.transcript_offset);

            let overlaps_indel = placed
                .iter()
                .any(|indel| indel_disturbs_codon(gene, indel, codon_start, codon_end));

            // Position and ALT together: the read phasing was queried per
            // allele, so a position alone cannot pick the right answer at a
            // multi-allelic site.
            let codon_positions: Vec<(usize, char)> = codon_snps
                .iter()
                .filter_map(|snp| {
                    snp.snp
                        .base
                        .chars()
                        .next()
                        .map(|alt| (snp.snp.position, alt.to_ascii_uppercase()))
                })
                .collect();
            let mut upstream_shift: isize = 0;
            let mut has_symbolic_sv = false;
            let mut frameshift_linkage = Vec::new();
            for (indel, coding) in indels.iter().zip(&placed) {
                if indel_is_upstream_of_codon_offset(gene, coding, codon_start) {
                    // Only let sufficiently frequent upstream indels shift the frame
                    // of downstream codons (see IndelAnnotationConfig).
                    if !indel_passes_frameshift_gate(indel, config) {
                        continue;
                    }
                    // Suppress propagation when the BAM shows this indel is in trans
                    // with the codon's SNV: the codon's molecule does not carry it.
                    // Record what the reads said either way, so the label the row
                    // ends up with can be traced back to the evidence for it.
                    if let Some(linkage) = phasing.linkage_for(indel, &codon_positions) {
                        frameshift_linkage.push(linkage);
                    }
                    if phasing.indel_in_trans_with(indel, &codon_positions) {
                        continue;
                    }
                    if indel.alt_allele.starts_with('<') {
                        has_symbolic_sv = true;
                    } else if let Some(delta) = coding_delta_for_variant(gene, reference, coding) {
                        upstream_shift += delta;
                    } else {
                        has_symbolic_sv = true;
                    }
                }
            }
            // `has_symbolic_sv` means the upstream shift could not be
            // determined, not that there is one. Reading it as a frameshift
            // relabelled every downstream codon of the feature `(fs)` and
            // frameshift_variant/HIGH on the strength of an allele nobody had
            // measured: a plain missense substitution downstream of an <INV>
            // came out HIGH. A codon whose frame cannot be judged is reported
            // the way a codon an indel overlaps already is, as unknown.
            let upstream_shift_unknown = has_symbolic_sv;
            let is_frameshifted = upstream_shift % 3 != 0;

            let mut var_info = process_transcript_codon(
                gene,
                chrom,
                ref_codon,
                codon_start,
                &codon_snps,
                genetic_code,
            );

            if overlaps_indel || upstream_shift_unknown {
                var_info.change_type = ChangeType::IndelOverlap;
                var_info.aa_changes = vec!["Unknown".to_string()];
                var_info.snp_aa_changes =
                    vec!["Unknown".to_string(); var_info.snp_aa_changes.len()];
                var_info.aa_changes_local = vec!["Unknown".to_string()];
                var_info.snp_aa_changes_local =
                    vec!["Unknown".to_string(); var_info.snp_aa_changes_local.len()];
            } else if is_frameshifted {
                apply_frameshift_labeling(&mut var_info, ptc_protein_pos);
            }

            // A nonsense SNV/MNV (still classified Stop gained after any
            // frameshift/overlap relabelling) begins its premature stop at this
            // codon's CDS offset.
            if var_info.change_type == ChangeType::StopGained {
                var_info.annotations.nmd = nmd::snv_mnv_nmd_prediction(gene, codon_start);
            }
            var_info.annotations.frameshift_linkage = frameshift_linkage;

            variants.push(var_info);
        }
    }

    for (indel, coding) in indels.into_iter().zip(placed) {
        let indel_effect::IndelEffect {
            change_type,
            aa_changes,
            aa_changes_local,
            ref_codon,
            alt_codon,
            initiator_lost,
        } = protein_effect_for_indel(gene, reference, &coding, genetic_code);
        let (event_class, event_components) = variant_event_metadata(&indel);
        let nmd = nmd::indel_nmd_prediction(gene, reference, &coding, genetic_code);

        variants.push(VariantInfo {
            chrom: chrom.to_string(),
            gene: gene.name.clone(),
            positions: vec![indel.record_start],
            ref_bases: vec![indel.ref_allele.clone()],
            base_changes: vec![indel.alt_allele.clone()],
            aa_changes,
            snp_aa_changes: vec!["-".to_string()],
            aa_changes_local,
            snp_aa_changes_local: vec!["-".to_string()],
            variant_type: VariantType::Indel,
            change_type,
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
            mnv_phasing_reads: None,
            ref_codon,
            snp_codon: None,
            mnv_codon: alt_codon,
            original_dp: indel.original_dp.map(|value| vec![Some(value)]),
            original_freq: indel.original_freq.map(|value| vec![Some(value)]),
            original_info: indel.original_info.clone(),
            event_class: Some(event_class),
            event_components,
            annotations: crate::variants::VariantAnnotations {
                nmd,
                splice: splice_consequence_for_indel(gene, &indel),
                hgvs_c: indel_hgvs_c(gene, &coding, reference),
                initiator_lost,
                ..Default::default()
            },
        });
    }

    variants
}

/// Where the coding sequence sees each indel (see
/// [`indel_effect::coding_placement`]). The rows keep each record as written,
/// which is what its reads are counted against.
fn coding_placements(
    gene: &Gene,
    reference: &crate::io::Reference,
    indels: &[crate::io::VcfPosition],
) -> Vec<crate::io::VcfPosition> {
    indels
        .iter()
        .map(|indel| {
            indel_effect::coding_placement(gene, reference, indel).unwrap_or_else(|| indel.clone())
        })
        .collect()
}

/// The splice consequence of an indel, taken from the bases it changes.
///
/// The substitution path folds this into the SO term of its coding row. The
/// indel rows left it unset, so the pass that makes each gene answer for itself
/// added a second row for the same event and the same gene, and the two
/// contradicted each other: one said `frameshift_variant` at HIGH with the
/// residues it shifts, the other `splice_region_variant` at LOW with no
/// amino-acid change at all. A substitution at the very same base came back as
/// one row reading `missense_variant&splice_region_variant`.
fn splice_consequence_for_indel(
    gene: &Gene,
    indel: &crate::io::VcfPosition,
) -> Option<crate::variants::SpliceConsequence> {
    indel
        .changed_positions()
        .into_iter()
        .filter_map(|position| {
            crate::variants::splice::splice_consequence_for_position(gene, position)
        })
        .max_by_key(|consequence| consequence.severity())
}

pub fn get_mnv_variants_for_gene(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
) -> Vec<VariantInfo> {
    get_mnv_variants_for_gene_with_config(
        gene,
        snp_list,
        reference,
        chrom,
        genetic_code,
        &IndelAnnotationConfig::default(),
        &FrameshiftPhasing::default(),
    )
}

/// Like [`get_mnv_variants_for_gene`] but with explicit indel-annotation
/// configuration (e.g. the upstream-indel frequency gate for frameshift
/// propagation) and BAM-derived phasing evidence that suppresses propagation
/// from indels shown to be in trans with a codon. The no-config wrapper above
/// preserves historical behaviour.
pub fn get_mnv_variants_for_gene_with_config(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
    config: &IndelAnnotationConfig,
    phasing: &FrameshiftPhasing,
) -> Vec<VariantInfo> {
    if has_transcript_cds_model(gene) {
        return get_mnv_variants_for_transcript(
            gene,
            snp_list,
            reference,
            chrom,
            genetic_code,
            config,
            phasing,
        );
    }

    let mut variants = Vec::new();

    // Each codon_start maps to a list of mutually-exclusive interpretations.
    // After --split-multiallelic, two ALT alleles at the same position
    // produce alternative interpretations of the same codon: each is
    // processed independently to emit one row per multi-allelic combination.
    let mut codon_to_groups: BTreeMap<usize, Vec<Vec<crate::variants::Snp>>> = BTreeMap::new();
    let mut indels: Vec<crate::io::VcfPosition> = Vec::new();

    for variant in snp_list
        .iter()
        .filter(|variant| variant.overlaps_interval(gene.start, gene.end))
    {
        let substitutions = variant.substitution_components();
        if !substitutions.is_empty() {
            for component in substitutions {
                if let Some((codon_start, _)) = codon_bounds_for_position(gene, component.position)
                {
                    let snp = crate::variants::Snp {
                        index: component.position,
                        position: component.position,
                        ref_base: component.ref_allele,
                        base: component.alt_allele,
                        original_dp: variant.original_dp,
                        original_freq: variant.original_freq,
                        original_info: variant.original_info.clone(),
                    };
                    let groups = codon_to_groups.entry(codon_start).or_default();
                    merge_snp_into_groups(groups, snp);
                }
            }
        } else {
            indels.push(variant.clone());
        }
    }

    if codon_to_groups.is_empty() && indels.is_empty() {
        return variants;
    }

    let placed = coding_placements(gene, reference, &indels);
    let ptc_protein_pos =
        frameshift_ptc_protein_pos(gene, reference, &placed, genetic_code, config);

    let mut codon_starts: Vec<usize> = codon_to_groups.keys().copied().collect();
    match gene.strand {
        Strand::Plus => codon_starts.sort_unstable(),
        Strand::Minus => codon_starts.sort_unstable_by(|a, b| b.cmp(a)),
    }

    for codon_start in codon_starts {
        let groups = codon_to_groups
            .get(&codon_start)
            .cloned()
            .unwrap_or_default();
        for codon_snps in groups {
            if codon_snps.is_empty() {
                continue;
            }
            let codon_end = codon_start + 2;
            if codon_end > reference.sequence.len() {
                continue;
            }
            let mut codon_snps = codon_snps;
            codon_snps.sort_by_key(|s| s.position);

            if codon_start == 0 || codon_end > reference.sequence.len() {
                continue;
            }
            let codon_seq = &reference.sequence[(codon_start - 1)..codon_end];

            let overlaps_indel = placed
                .iter()
                .any(|indel| indel.overlaps_interval(codon_start, codon_end));

            let codon_positions: Vec<(usize, char)> = codon_snps
                .iter()
                .filter_map(|snp| {
                    snp.base
                        .chars()
                        .next()
                        .map(|alt| (snp.position, alt.to_ascii_uppercase()))
                })
                .collect();
            let mut upstream_shift: isize = 0;
            let mut has_symbolic_sv = false;
            let mut frameshift_linkage = Vec::new();

            for (indel, coding) in indels.iter().zip(&placed) {
                let is_upstream =
                    indel_is_upstream_of_codon(coding, gene.strand, codon_start, codon_end);

                if is_upstream {
                    // Only sufficiently frequent upstream indels shift the frame of
                    // downstream codons (see IndelAnnotationConfig).
                    if !indel_passes_frameshift_gate(indel, config) {
                        continue;
                    }
                    // Suppress propagation when the BAM shows this indel is in trans
                    // with the codon's SNV: the codon's molecule does not carry it.
                    // Record what the reads said either way, so the label the row
                    // ends up with can be traced back to the evidence for it.
                    if let Some(linkage) = phasing.linkage_for(indel, &codon_positions) {
                        frameshift_linkage.push(linkage);
                    }
                    if phasing.indel_in_trans_with(indel, &codon_positions) {
                        continue;
                    }
                    if indel.alt_allele.starts_with('<') {
                        has_symbolic_sv = true;
                    } else if let Some(delta) = coding_delta_for_variant(gene, reference, coding) {
                        // CDS-clipped coding-length change, consistent with the
                        // transcript model. For indels that span the gene/CDS
                        // boundary this counts only the in-CDS portion, unlike the
                        // raw allele-length difference used previously.
                        upstream_shift += delta;
                    } else {
                        upstream_shift +=
                            (indel.alt_allele.len() as isize) - (indel.ref_allele.len() as isize);
                    }
                }
            }

            // `has_symbolic_sv` means the upstream shift could not be
            // determined, not that there is one. Reading it as a frameshift
            // relabelled every downstream codon of the feature `(fs)` and
            // frameshift_variant/HIGH on the strength of an allele nobody had
            // measured: a plain missense substitution downstream of an <INV>
            // came out HIGH. A codon whose frame cannot be judged is reported
            // the way a codon an indel overlaps already is, as unknown.
            let upstream_shift_unknown = has_symbolic_sv;
            let is_frameshifted = upstream_shift % 3 != 0;

            let (eff_gene_start, eff_gene_end) = effective_bounds(gene);
            let codon_info = CodonInfo {
                codon_list: codon_snps,
                original_codon: codon_seq.to_string(),
                gene_name: gene.name.clone(),
                gene_start: eff_gene_start,
                gene_end: eff_gene_end,
                codon_start,
                codon_end,
                protein_offset: gene.protein_offset,
                starts_protein: gene.starts_protein(),
            };

            let mut var_info = process_codon(codon_info, gene.strand, chrom, genetic_code);

            if overlaps_indel || upstream_shift_unknown {
                var_info.change_type = ChangeType::IndelOverlap;
                var_info.aa_changes = vec!["Unknown".to_string()];
                var_info.snp_aa_changes =
                    vec!["Unknown".to_string(); var_info.snp_aa_changes.len()];
                var_info.aa_changes_local = vec!["Unknown".to_string()];
                var_info.snp_aa_changes_local =
                    vec!["Unknown".to_string(); var_info.snp_aa_changes_local.len()];
            } else if is_frameshifted {
                apply_frameshift_labeling(&mut var_info, ptc_protein_pos);
            }
            var_info.annotations.frameshift_linkage = frameshift_linkage;

            variants.push(var_info);
        }
    }

    for (indel, coding) in indels.into_iter().zip(placed) {
        let indel_effect::IndelEffect {
            change_type,
            aa_changes,
            aa_changes_local,
            ref_codon,
            alt_codon,
            initiator_lost,
        } = protein_effect_for_indel(gene, reference, &coding, genetic_code);
        let (event_class, event_components) = variant_event_metadata(&indel);

        variants.push(VariantInfo {
            chrom: chrom.to_string(),
            gene: gene.name.clone(),
            positions: vec![indel.record_start],
            ref_bases: vec![indel.ref_allele.clone()],
            base_changes: vec![indel.alt_allele.clone()],
            aa_changes,
            snp_aa_changes: vec!["-".to_string()],
            aa_changes_local,
            snp_aa_changes_local: vec!["-".to_string()],
            variant_type: VariantType::Indel,
            change_type,
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
            mnv_phasing_reads: None,
            ref_codon,
            snp_codon: None,
            mnv_codon: alt_codon,
            original_dp: indel.original_dp.map(|value| vec![Some(value)]),
            original_freq: indel.original_freq.map(|value| vec![Some(value)]),
            original_info: indel.original_info.clone(),
            event_class: Some(event_class),
            event_components,
            annotations: crate::variants::VariantAnnotations {
                splice: splice_consequence_for_indel(gene, &indel),
                hgvs_c: indel_hgvs_c(gene, &coding, reference),
                initiator_lost,
                ..Default::default()
            },
        });
    }

    variants
}

/// Build a `VariantInfo` for a variant in the splice region of a gene's intron.
/// It carries the gene name and splice consequence; everything else mirrors an
/// intergenic entry (no codon annotation), so its read support is counted
/// alongside intergenic SNPs and it renders through the normal genic path.
pub fn build_splice_variant(
    chrom: &str,
    vcf_pos: &crate::io::VcfPosition,
    gene_name: &str,
    splice: crate::variants::SpliceConsequence,
) -> VariantInfo {
    let mut variant = build_intergenic_variant(chrom, vcf_pos);
    variant.gene = gene_name.to_string();
    variant.annotations.splice = Some(splice);
    variant
}

/// Build a `VariantInfo` for a variant inside a real intron of a gene, away from
/// its splice sites. It carries the gene name so the variant is reported against
/// the transcript it actually sits in rather than as intergenic, which is what a
/// position inside a gene is not.
pub fn build_intron_variant(
    chrom: &str,
    vcf_pos: &crate::io::VcfPosition,
    gene_name: &str,
) -> VariantInfo {
    build_non_coding_variant(
        chrom,
        vcf_pos,
        gene_name,
        crate::variants::NonCodingReason::Intron,
    )
}

/// Build a `VariantInfo` for a variant inside a gene that carries no codon: an
/// intron, or a transcript that is never translated. It keeps the gene name and
/// the reason, and has no amino-acid annotation because there is none to make.
pub fn build_non_coding_variant(
    chrom: &str,
    vcf_pos: &crate::io::VcfPosition,
    gene_name: &str,
    reason: crate::variants::NonCodingReason,
) -> VariantInfo {
    let mut variant = build_intergenic_variant(chrom, vcf_pos);
    variant.gene = gene_name.to_string();
    variant.annotations.non_coding = Some(reason);
    variant
}

/// Build a `VariantInfo` for a VCF position that falls outside all annotated
/// genes (intergenic).  The entry preserves the original position, alleles and
/// any original metrics but has no codon / amino-acid annotation.
pub fn build_intergenic_variant(chrom: &str, vcf_pos: &crate::io::VcfPosition) -> VariantInfo {
    let event = vcf_pos.event();
    let variant_type = match event.class {
        crate::variants::AlleleEventClass::Snp => VariantType::Snp,
        crate::variants::AlleleEventClass::Mnv => VariantType::Mnv,
        _ => VariantType::Indel,
    };

    // A substitution row carries one entry per base that changes, the same shape
    // the gene path builds. Keeping the record's raw REF and ALT here left the
    // read counter with a multi-base "allele" it could not use: a row outside a
    // CDS got no counts at all, and a threshold then read the absent counts as
    // zero and deleted it. The event already says which bases change.
    let substitutions: Vec<&crate::variants::AlleleComponent> = event
        .components
        .iter()
        .filter(|c| c.kind == crate::variants::AlleleComponentKind::Snp)
        .collect();
    let (positions, ref_bases, base_changes) = if matches!(
        event.class,
        crate::variants::AlleleEventClass::Snp | crate::variants::AlleleEventClass::Mnv
    ) && !substitutions.is_empty()
    {
        (
            substitutions.iter().map(|c| c.position).collect(),
            substitutions.iter().map(|c| c.ref_allele.clone()).collect(),
            substitutions.iter().map(|c| c.alt_allele.clone()).collect(),
        )
    } else {
        (
            vec![vcf_pos.record_start],
            vec![vcf_pos.ref_allele.clone()],
            vec![vcf_pos.alt_allele.clone()],
        )
    };
    let entries = positions.len();

    VariantInfo {
        chrom: chrom.to_string(),
        gene: "intergenic".to_string(),
        positions,
        ref_bases,
        base_changes,
        aa_changes: vec!["-".to_string()],
        snp_aa_changes: vec!["-".to_string()],
        aa_changes_local: vec!["-".to_string()],
        snp_aa_changes_local: vec!["-".to_string()],
        variant_type,
        change_type: ChangeType::Unknown,
        ref_codon: None,
        snp_codon: None,
        mnv_codon: None,
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
        mnv_phasing_reads: None,
        original_dp: vcf_pos.original_dp.map(|dp| vec![Some(dp); entries]),
        original_freq: vcf_pos.original_freq.map(|freq| vec![Some(freq); entries]),
        original_info: vcf_pos.original_info.clone(),
        event_class: Some(event.class.as_str().to_string()),
        event_components: event.component_labels(),
        annotations: crate::variants::VariantAnnotations::default(),
    }
}
