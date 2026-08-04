//! Protein-level effect of indels (frameshift, in-frame, stop gained/lost).

use crate::io::VcfPosition;
use crate::variants::{ChangeType, Gene, Strand};

use super::allele_apply::apply_allele_to_feature;
use super::gene_path::effective_bounds;
use super::protein::{describe_protein_change, translate_cds};
use super::transcript_model::{coding_sequence_for_gene, first_touched_transcript_offset};

/// Detect whether an indel destroys the initiator codon, i.e. the reference
/// protein begins with Met and the alternate no longer does.
///
/// Substitutions already report this through
/// [`crate::utils::determine_change_type`], so without this an in-frame indel
/// that removes the translation start was labelled a plain `inframe_deletion` /
/// `inframe_insertion` (`MODERATE`) while the equivalent SNV was `start_lost`
/// (`HIGH`), hiding it from anyone filtering on impact.
///
/// `None` when the reference protein does not start with Met, which is how a
/// partial or phase-shifted CDS model presents itself: there is no initiator to
/// lose. A downstream Met left at position 1 is not reported as a loss.
pub(super) fn indel_start_effect(ref_protein: &str, alt_protein: &str) -> Option<ChangeType> {
    if !ref_protein.starts_with('M') {
        return None;
    }
    (!alt_protein.starts_with('M')).then_some(ChangeType::StartLost)
}

/// Detect whether an in-frame indel creates or removes a stop codon by
/// comparing the number of stop residues (`*`) in the reference and alternate
/// translations. Counting stops (rather than comparing positions) avoids a
/// false "stop gained" for ordinary in-frame deletions, whose terminal stop
/// simply shifts to a lower index in a shorter protein.
pub(super) fn indel_stop_effect(ref_protein: &str, alt_protein: &str) -> Option<ChangeType> {
    let ref_stops = ref_protein.matches('*').count();
    let alt_stops = alt_protein.matches('*').count();
    if alt_stops > ref_stops {
        Some(ChangeType::StopGained)
    } else if ref_stops > 0 && alt_stops < ref_stops {
        Some(ChangeType::StopLost)
    } else {
        None
    }
}

pub(super) fn protein_effect_for_indel(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
    genetic_code: crate::genetic_code::GeneticCode,
) -> (
    ChangeType,
    Vec<String>,
    Vec<String>,
    Option<String>,
    Option<String>,
) {
    // Base classification (frameshift vs in-frame, symbolic handling). When the
    // alternate CDS can be reconstructed we may refine this to Stop gained/lost.
    let base_change_type = indel_change_type_for_variant(gene, reference, variant);

    let Some(ref_cds) = coding_sequence_for_gene(gene, reference) else {
        return (
            base_change_type,
            vec!["Unknown".to_string()],
            vec!["Unknown".to_string()],
            None,
            None,
        );
    };
    let Some(alt_cds) = apply_allele_to_feature(gene, reference, variant) else {
        return (
            base_change_type,
            vec!["Unknown".to_string()],
            vec!["Unknown".to_string()],
            None,
            None,
        );
    };

    let frameshift = (alt_cds.len() as isize - ref_cds.len() as isize) % 3 != 0;
    let ref_protein = translate_cds(&ref_cds, genetic_code);
    let alt_protein = translate_cds(&alt_cds, genetic_code);

    // For in-frame indels, refine the change type when the event destroys the
    // initiator codon or creates/removes a stop (otherwise these high-impact
    // events are hidden under the generic "In-frame Indel" label). Losing the
    // start outranks any stop change: without an initiator the ORF is not
    // translated from here at all. Frameshifts keep the frameshift label because
    // they almost always introduce a downstream stop, so flagging them as "stop
    // gained" would be uninformative, and they are already HIGH impact.
    let change_type = if frameshift {
        base_change_type
    } else {
        indel_start_effect(&ref_protein, &alt_protein)
            .or_else(|| indel_stop_effect(&ref_protein, &alt_protein))
            .unwrap_or(base_change_type)
    };

    let (protein_change, local_change) =
        describe_protein_change(&ref_protein, &alt_protein, gene.protein_offset, frameshift);

    let transcript_codon = first_touched_transcript_offset(gene, variant).and_then(|offset| {
        let codon_start = (offset / 3) * 3;
        let codon_end = codon_start + 3;
        let ref_codon = ref_cds.get(codon_start..codon_end).map(str::to_string)?;
        let alt_codon = alt_cds.get(codon_start..codon_end).map(str::to_string);
        Some((ref_codon, alt_codon))
    });

    if let Some((ref_codon, alt_codon)) = transcript_codon {
        return (
            change_type,
            vec![protein_change],
            vec![local_change],
            Some(ref_codon),
            alt_codon,
        );
    }

    // Legacy path (no transcript CDS model): map the first affected genomic base
    // to its offset in the coding sequence (which is reverse-complemented for the
    // minus strand) and show the codon at that offset in BOTH the reference and
    // alternate CDS. The previous code derived the offset from reference
    // coordinates and reused it to slice the length-changed alt CDS, which on the
    // minus strand selected the wrong reference codon and ran off the end of the
    // shorter alt CDS (yielding mnv_codon = None for every minus-strand deletion).
    let (eff_start, eff_end) = effective_bounds(gene);
    let event = variant.event();
    let touched_lo = event.affected_start.max(eff_start);
    let touched_hi = event.affected_end.min(eff_end);
    let (ref_codon, alt_codon) = if touched_lo <= touched_hi {
        let first_offset = match gene.strand {
            Strand::Plus => touched_lo - eff_start,
            Strand::Minus => eff_end - touched_hi,
        };
        let codon_start = (first_offset / 3) * 3;
        let ref_codon = ref_cds
            .get(codon_start..codon_start + 3)
            .map(str::to_string);
        let alt_codon = alt_cds
            .get(codon_start..codon_start + 3)
            .map(str::to_string);
        (ref_codon, alt_codon)
    } else {
        (None, None)
    };

    (
        change_type,
        vec![protein_change],
        vec![local_change],
        ref_codon,
        alt_codon,
    )
}

pub(super) fn coding_delta_for_variant(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
) -> Option<isize> {
    let ref_cds = coding_sequence_for_gene(gene, reference)?;
    let alt_cds = apply_allele_to_feature(gene, reference, variant)?;
    Some(alt_cds.len() as isize - ref_cds.len() as isize)
}

pub(super) fn indel_change_type_for_variant(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
) -> ChangeType {
    let frameshift = if variant.alt_allele.starts_with('<') {
        true
    } else if let Some(delta) = coding_delta_for_variant(gene, reference, variant) {
        delta % 3 != 0
    } else {
        (variant.ref_allele.len() as isize - variant.alt_allele.len() as isize) % 3 != 0
    };

    if frameshift {
        ChangeType::FrameshiftIndel
    } else {
        ChangeType::InFrameIndel
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn losing_the_initiator_is_start_lost() {
        // The whole start codon deleted, and the start codon disrupted in place.
        assert_eq!(
            indel_start_effect("MKLA*", "KLA*"),
            Some(ChangeType::StartLost)
        );
        assert_eq!(
            indel_start_effect("MKLA*", "RVKLA*"),
            Some(ChangeType::StartLost)
        );
    }

    #[test]
    fn an_intact_initiator_is_not_start_lost() {
        // Insertion and deletion downstream of the start codon.
        assert_eq!(indel_start_effect("MKLA*", "MKLGGA*"), None);
        assert_eq!(indel_start_effect("MKLA*", "MA*"), None);
        // A Met left at position 1 still initiates translation.
        assert_eq!(indel_start_effect("MKLA*", "MLA*"), None);
    }

    #[test]
    fn a_cds_model_without_an_initiator_is_left_alone() {
        // A partial or phase-shifted CDS does not begin at the initiator, so
        // there is no start to lose and the generic in-frame label stands.
        assert_eq!(indel_start_effect("KLA*", "LA*"), None);
        assert_eq!(indel_start_effect("", "KLA*"), None);
    }

    #[test]
    fn start_loss_outranks_a_stop_change() {
        // Deleting the start codon of a protein whose stop also moves must not
        // be reported as a stop change.
        assert_eq!(
            indel_start_effect("MKLA*", "KLA"),
            Some(ChangeType::StartLost)
        );
        assert_eq!(
            indel_stop_effect("MKLA*", "KLA"),
            Some(ChangeType::StopLost)
        );
    }
}
