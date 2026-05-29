//! Apply a variant's allele to a gene's CDS / feature sequence.

use crate::utils::reverse_complement;
use crate::variants::event::{AlleleComponentKind, AlleleEventClass};
use crate::variants::{Gene, Strand};
use std::collections::{BTreeMap, BTreeSet};

use super::gene_path::effective_bounds;
use super::transcript_model::{has_transcript_cds_model, insertion_anchor_in_segment};

pub(super) fn apply_allele_to_transcript(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &crate::io::VcfPosition,
) -> Option<String> {
    if !has_transcript_cds_model(gene) {
        return None;
    }

    let event = variant.event();
    if matches!(event.class, AlleleEventClass::Symbolic) {
        return None;
    }

    let mut touched_transcript = false;
    let mut cds = String::new();

    for segment in &gene.cds_segments {
        if segment.start == 0
            || segment.end == 0
            || segment.start > segment.end
            || segment.end > reference.sequence.len()
        {
            return None;
        }

        let mut snps: BTreeMap<usize, String> = BTreeMap::new();
        let mut insertions: BTreeMap<usize, String> = BTreeMap::new();
        let mut deleted: BTreeSet<usize> = BTreeSet::new();

        for component in &event.components {
            match component.kind {
                AlleleComponentKind::Snp => {
                    if component.position >= segment.start && component.position <= segment.end {
                        match snps.get(&component.position) {
                            Some(existing)
                                if !existing.eq_ignore_ascii_case(&component.alt_allele) =>
                            {
                                return None;
                            }
                            Some(_) => {}
                            None => {
                                snps.insert(component.position, component.alt_allele.clone());
                            }
                        }
                        touched_transcript = true;
                    }
                }
                AlleleComponentKind::Insertion => {
                    if insertion_anchor_in_segment(gene, segment, component.position) {
                        match insertions.get(&component.position) {
                            Some(existing)
                                if !existing.eq_ignore_ascii_case(&component.alt_allele) =>
                            {
                                return None;
                            }
                            Some(_) => {}
                            None => {
                                insertions.insert(component.position, component.alt_allele.clone());
                            }
                        }
                        touched_transcript = true;
                    }
                }
                AlleleComponentKind::Deletion => {
                    let del_end = component.position + component.ref_allele.len().saturating_sub(1);
                    let overlap_start = component.position.max(segment.start);
                    let overlap_end = del_end.min(segment.end);
                    if overlap_start <= overlap_end {
                        for pos in overlap_start..=overlap_end {
                            deleted.insert(pos);
                        }
                        touched_transcript = true;
                    }
                }
                AlleleComponentKind::Delins | AlleleComponentKind::Symbolic => return None,
            }
        }

        if snps.keys().any(|pos| deleted.contains(pos))
            || insertions.keys().any(|pos| deleted.contains(pos))
        {
            return None;
        }

        let mut genomic = String::with_capacity(segment.end - segment.start + 1);
        for pos in segment.start..=segment.end {
            if !deleted.contains(&pos) {
                if let Some(alt_base) = snps.get(&pos) {
                    genomic.push_str(alt_base);
                } else {
                    genomic.push(reference.sequence.as_bytes()[pos - 1] as char);
                }
            }
            if let Some(inserted) = insertions.get(&pos) {
                genomic.push_str(inserted);
            }
        }

        match gene.strand {
            Strand::Plus => cds.push_str(&genomic),
            Strand::Minus => cds.push_str(&reverse_complement(&genomic)),
        }
    }

    touched_transcript.then_some(cds)
}

pub(super) fn apply_allele_to_feature(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &crate::io::VcfPosition,
) -> Option<String> {
    if let Some(transcript_alt) = apply_allele_to_transcript(gene, reference, variant) {
        return Some(transcript_alt);
    }

    let (eff_start, eff_end) = effective_bounds(gene);
    if eff_start == 0 || eff_end == 0 || eff_start > eff_end || eff_end > reference.sequence.len() {
        return None;
    }

    let event = variant.event();
    if matches!(event.class, AlleleEventClass::Symbolic) {
        return None;
    }

    let mut snps: BTreeMap<usize, String> = BTreeMap::new();
    let mut insertions: BTreeMap<usize, String> = BTreeMap::new();
    let mut deleted: BTreeSet<usize> = BTreeSet::new();
    let mut touched_feature = false;

    for component in event.components {
        match component.kind {
            AlleleComponentKind::Snp => {
                if component.position >= eff_start && component.position <= eff_end {
                    match snps.get(&component.position) {
                        Some(existing) if !existing.eq_ignore_ascii_case(&component.alt_allele) => {
                            return None;
                        }
                        Some(_) => {}
                        None => {
                            snps.insert(component.position, component.alt_allele);
                        }
                    }
                    touched_feature = true;
                }
            }
            AlleleComponentKind::Insertion => {
                if component.position >= eff_start && component.position < eff_end {
                    match insertions.get(&component.position) {
                        Some(existing) if !existing.eq_ignore_ascii_case(&component.alt_allele) => {
                            return None;
                        }
                        Some(_) => {}
                        None => {
                            insertions.insert(component.position, component.alt_allele);
                        }
                    }
                    touched_feature = true;
                }
            }
            AlleleComponentKind::Deletion => {
                let del_end = component.position + component.ref_allele.len().saturating_sub(1);
                let overlap_start = component.position.max(eff_start);
                let overlap_end = del_end.min(eff_end);
                if overlap_start <= overlap_end {
                    for pos in overlap_start..=overlap_end {
                        deleted.insert(pos);
                    }
                    touched_feature = true;
                }
            }
            AlleleComponentKind::Delins | AlleleComponentKind::Symbolic => return None,
        }
    }

    if !touched_feature {
        return None;
    }

    if snps.keys().any(|pos| deleted.contains(pos))
        || insertions.keys().any(|pos| deleted.contains(pos))
    {
        return None;
    }

    let mut genomic = String::with_capacity(reference.sequence[(eff_start - 1)..eff_end].len());
    for pos in eff_start..=eff_end {
        if !deleted.contains(&pos) {
            if let Some(alt_base) = snps.get(&pos) {
                genomic.push_str(alt_base);
            } else {
                genomic.push(reference.sequence.as_bytes()[pos - 1] as char);
            }
        }
        if let Some(inserted) = insertions.get(&pos) {
            genomic.push_str(inserted);
        }
    }

    Some(match gene.strand {
        Strand::Plus => genomic,
        Strand::Minus => reverse_complement(&genomic),
    })
}
