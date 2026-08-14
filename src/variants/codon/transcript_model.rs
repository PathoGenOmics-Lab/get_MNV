//! Spliced-transcript CDS model: offsets, sequence reconstruction, overlaps.

use crate::io::VcfPosition;
use crate::utils::reverse_complement;
use crate::variants::event::AlleleComponentKind;
use crate::variants::{CdsSegment, Gene, Strand};
use std::collections::BTreeSet;

use super::gene_path::effective_bounds;

pub(super) fn has_transcript_cds_model(gene: &Gene) -> bool {
    !gene.cds_segments.is_empty()
}

pub(super) fn transcript_sequence_for_gene(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
) -> Option<String> {
    if !has_transcript_cds_model(gene) {
        return None;
    }

    let mut cds = String::new();
    for segment in &gene.cds_segments {
        if segment.start == 0
            || segment.end == 0
            || segment.start > segment.end
            || segment.end > reference.sequence.len()
        {
            return None;
        }
        let genomic = &reference.sequence[(segment.start - 1)..segment.end];
        match gene.strand {
            Strand::Plus => cds.push_str(genomic),
            Strand::Minus => cds.push_str(&reverse_complement(genomic)),
        }
    }
    Some(cds)
}

/// Every transcript offset a genomic position occupies.
///
/// Normally exactly one. A ribosomal-slippage join, the `join(a..b, b..c)` shape
/// this model supports for ORF1ab and friends, has two CDS rows that share a
/// base, and that base is read twice by the ribosome: it sits at two different
/// offsets in the spliced coding sequence. Taking only the first left the second
/// codon holding the reference base, so a substitution there was annotated with
/// the wrong residue, or produced no row at all when nothing else fell in that
/// codon. The indel path already walked every segment; this is the same walk for
/// a substitution.
pub(super) fn transcript_offsets_for_position(gene: &Gene, position: usize) -> Vec<usize> {
    if !has_transcript_cds_model(gene) {
        return Vec::new();
    }

    let mut offsets = Vec::new();
    let mut offset = 0usize;
    for segment in &gene.cds_segments {
        if position >= segment.start && position <= segment.end {
            offsets.push(match gene.strand {
                Strand::Plus => offset + position.saturating_sub(segment.start),
                Strand::Minus => offset + segment.end.saturating_sub(position),
            });
        }
        offset += segment.end.saturating_sub(segment.start) + 1;
    }
    offsets
}

pub(super) fn transcript_offset_for_position(gene: &Gene, position: usize) -> Option<usize> {
    if !has_transcript_cds_model(gene) {
        return None;
    }

    let mut offset = 0usize;
    for segment in &gene.cds_segments {
        if position >= segment.start && position <= segment.end {
            return Some(match gene.strand {
                Strand::Plus => offset + position.saturating_sub(segment.start),
                Strand::Minus => offset + segment.end.saturating_sub(position),
            });
        }
        offset += segment.end.saturating_sub(segment.start) + 1;
    }
    None
}

/// Whether an insertion anchored at `position` belongs to this CDS `segment`
/// when reconstructing the spliced transcript. True when the anchor is strictly
/// inside the exon, or when it sits exactly on the exon's last base and that base
/// is an internal exon-exon junction: a VCF insertion places its bases *after*
/// the anchor, so an anchor at an internal exon's last base lands the inserted
/// sequence at the junction with the next exon, inside the spliced CDS. An anchor
/// at the CDS terminus is UTR and is excluded (see
/// [`Gene::insertion_at_internal_junction`]).
pub(super) fn insertion_anchor_in_segment(
    gene: &Gene,
    segment: &CdsSegment,
    position: usize,
) -> bool {
    (position >= segment.start && position < segment.end)
        || (position == segment.end && gene.insertion_at_internal_junction(position))
}

pub(super) fn variant_overlaps_coding_model(gene: &Gene, variant: &VcfPosition) -> bool {
    if !has_transcript_cds_model(gene) {
        return variant.overlaps_interval(gene.start, gene.end);
    }

    if gene
        .cds_segments
        .iter()
        .any(|segment| variant.overlaps_interval(segment.start, segment.end))
    {
        return true;
    }
    // An insertion anchored at an internal exon-exon junction is coding even
    // though overlaps_interval (half-open) excludes the exon's last base.
    variant.event().components.iter().any(|component| {
        matches!(component.kind, AlleleComponentKind::Insertion)
            && gene.insertion_at_internal_junction(component.position)
    })
}

pub(super) fn variant_touched_transcript_offsets(gene: &Gene, variant: &VcfPosition) -> Vec<usize> {
    if !has_transcript_cds_model(gene) {
        return Vec::new();
    }

    let mut offsets = BTreeSet::new();
    for component in variant.event().components {
        match component.kind {
            AlleleComponentKind::Snp => {
                offsets.extend(transcript_offsets_for_position(gene, component.position));
            }
            AlleleComponentKind::Insertion => {
                for segment in &gene.cds_segments {
                    if insertion_anchor_in_segment(gene, segment, component.position) {
                        if let Some(offset) =
                            transcript_offset_for_position(gene, component.position)
                        {
                            offsets.insert(offset);
                        }
                    }
                }
            }
            AlleleComponentKind::Deletion
            | AlleleComponentKind::Delins
            | AlleleComponentKind::Symbolic => {
                let component_end =
                    component.position + component.ref_allele.len().saturating_sub(1);
                for segment in &gene.cds_segments {
                    let overlap_start = component.position.max(segment.start);
                    let overlap_end = component_end.min(segment.end);
                    if overlap_start <= overlap_end {
                        for pos in overlap_start..=overlap_end {
                            if let Some(offset) = transcript_offset_for_position(gene, pos) {
                                offsets.insert(offset);
                            }
                        }
                    }
                }
            }
        }
    }

    offsets.into_iter().collect()
}

/// Whether the indel disturbs the codon spanning `[codon_start, codon_end)` in
/// transcript coordinates.
///
/// An insertion places its bases after its anchor, so it falls inside the codon
/// only when the anchor is not the codon's last base: anchored there, the
/// inserted sequence sits between this codon and the next and leaves this
/// codon's three bases intact. That is the boundary rule the documentation
/// states and the genomic path applies through `VcfPosition::overlaps_interval`.
/// Using the anchor's own offset blanked the codon's substitution row to
/// "Indel overlap" and demoted a HIGH start_lost to MODIFIER, so the same
/// reference and the same variants gave two different answers depending only on
/// whether the CDS row carried a Parent attribute.
pub(super) fn indel_disturbs_codon(
    gene: &Gene,
    indel: &VcfPosition,
    codon_start: usize,
    codon_end: usize,
) -> bool {
    let inside = |offset: usize| offset >= codon_start && offset < codon_end;
    for component in indel.event().components {
        match component.kind {
            AlleleComponentKind::Snp => {
                if transcript_offsets_for_position(gene, component.position)
                    .into_iter()
                    .any(inside)
                {
                    return true;
                }
            }
            AlleleComponentKind::Insertion => {
                let anchored = gene
                    .cds_segments
                    .iter()
                    .any(|segment| insertion_anchor_in_segment(gene, segment, component.position));
                if !anchored {
                    continue;
                }
                if let Some(offset) = transcript_offset_for_position(gene, component.position) {
                    // The inserted bases follow the anchor along the genome, so in
                    // transcript order they land after it on the plus strand and
                    // before it on the minus, where the transcript reads from
                    // higher coordinates down. The gap is inside this codon only
                    // when both bases it sits between are, which is the same rule
                    // read from either end.
                    let disturbs = match gene.strand {
                        Strand::Plus => offset >= codon_start && offset + 1 < codon_end,
                        Strand::Minus => offset > codon_start && offset < codon_end,
                    };
                    if disturbs {
                        return true;
                    }
                }
            }
            AlleleComponentKind::Deletion
            | AlleleComponentKind::Delins
            | AlleleComponentKind::Symbolic => {
                let component_end =
                    component.position + component.ref_allele.len().saturating_sub(1);
                for segment in &gene.cds_segments {
                    let overlap_start = component.position.max(segment.start);
                    let overlap_end = component_end.min(segment.end);
                    if overlap_start > overlap_end {
                        continue;
                    }
                    for position in overlap_start..=overlap_end {
                        if transcript_offset_for_position(gene, position).is_some_and(inside) {
                            return true;
                        }
                    }
                }
            }
        }
    }
    false
}

pub(super) fn first_touched_transcript_offset(gene: &Gene, variant: &VcfPosition) -> Option<usize> {
    variant_touched_transcript_offsets(gene, variant)
        .into_iter()
        .min()
}

/// Whether the whole of this indel lies upstream of the codon starting at
/// `codon_start`, in transcript offsets.
///
/// An insertion occupies the boundary between two bases rather than a base, and
/// which two depends on the strand: its bases follow the anchor along the genome,
/// so in transcript order they sit after the anchor on the plus strand and before
/// it on the minus. Testing the anchor's own offset put a minus-strand insertion
/// anchored on a codon's first base inside that codon instead of before it, so a
/// frame that had shifted was reported as a plain codon change.
pub(super) fn indel_is_upstream_of_codon_offset(
    gene: &Gene,
    indel: &VcfPosition,
    codon_start: usize,
) -> bool {
    let Some(first) = first_touched_transcript_offset(gene, indel) else {
        return false;
    };
    let insertion_only = indel
        .event()
        .components
        .iter()
        .all(|component| matches!(component.kind, AlleleComponentKind::Insertion));
    let effective = if insertion_only && gene.strand == Strand::Minus {
        first.saturating_sub(1)
    } else {
        first
    };
    effective < codon_start
}

pub(super) fn coding_sequence_for_gene(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
) -> Option<String> {
    if let Some(transcript_cds) = transcript_sequence_for_gene(gene, reference) {
        return Some(transcript_cds);
    }

    let (eff_start, eff_end) = effective_bounds(gene);
    if eff_start == 0 || eff_end == 0 || eff_start > eff_end || eff_end > reference.sequence.len() {
        return None;
    }
    let genomic = &reference.sequence[(eff_start - 1)..eff_end];
    Some(match gene.strand {
        Strand::Plus => genomic.to_string(),
        Strand::Minus => reverse_complement(genomic),
    })
}
