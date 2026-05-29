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

/// Whether an insertion anchored at `position` falls strictly inside this CDS
/// segment. The upper bound is exclusive (`< segment.end`) on purpose: a VCF
/// insertion places the inserted bases *after* the anchor, so an anchor at the
/// segment's last base sits on the exon|intron boundary, where the inserted
/// sequence is genomically ambiguous between "end of this exon" and "start of
/// the intron" (spliced out). Such boundary insertions are intentionally left
/// out of the spliced-CDS model rather than assigned a possibly-wrong codon;
/// they surface as intergenic instead. Anchors inside an exon are unaffected.
pub(super) fn insertion_anchor_in_segment(segment: &CdsSegment, position: usize) -> bool {
    position >= segment.start && position < segment.end
}

pub(super) fn variant_overlaps_coding_model(gene: &Gene, variant: &VcfPosition) -> bool {
    if !has_transcript_cds_model(gene) {
        return variant.overlaps_interval(gene.start, gene.end);
    }

    gene.cds_segments
        .iter()
        .any(|segment| variant.overlaps_interval(segment.start, segment.end))
}

pub(super) fn variant_touched_transcript_offsets(gene: &Gene, variant: &VcfPosition) -> Vec<usize> {
    if !has_transcript_cds_model(gene) {
        return Vec::new();
    }

    let mut offsets = BTreeSet::new();
    for component in variant.event().components {
        match component.kind {
            AlleleComponentKind::Snp => {
                if let Some(offset) = transcript_offset_for_position(gene, component.position) {
                    offsets.insert(offset);
                }
            }
            AlleleComponentKind::Insertion => {
                for segment in &gene.cds_segments {
                    if insertion_anchor_in_segment(segment, component.position) {
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

pub(super) fn first_touched_transcript_offset(gene: &Gene, variant: &VcfPosition) -> Option<usize> {
    variant_touched_transcript_offsets(gene, variant)
        .into_iter()
        .min()
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
