//! Protein-level effect of indels (frameshift, in-frame, stop gained/lost).

use crate::io::VcfPosition;
use crate::variants::{ChangeType, Gene, Strand};

use super::allele_apply::apply_allele_to_feature;
use super::gene_path::effective_bounds;
use super::protein::{describe_protein_change, translate_cds};
use super::transcript_model::{coding_sequence_for_gene, first_touched_transcript_offset};

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

    // Refinement of the base label (initiator lost, stop gained/lost) lives with
    // the substitution rules in `variants::consequence`, so the two paths cannot
    // drift apart again.
    let change_type = crate::variants::consequence::indel_change_type(
        &ref_protein,
        &alt_protein,
        frameshift,
        base_change_type,
    );

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
