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

    // Whether the frame shifts is a property of the whole feature, so it is
    // measured there rather than on the phase-trimmed sequences translated
    // below: those omit the leading coding bases and undercounted a deletion
    // that reached into them.
    let frameshift = coding_delta_for_variant(gene, reference, variant)
        .map(|delta| delta % 3 != 0)
        .unwrap_or((alt_cds.len() as isize - ref_cds.len() as isize) % 3 != 0);
    let ref_protein = translate_cds(&ref_cds, genetic_code);
    let alt_protein = translate_cds(&alt_cds, genetic_code);
    // What the ribosome actually makes, which stops at the first stop codon.
    // Everything past it is translated by nobody, so two alternates that differ
    // only there produce the same protein.
    let protein_unchanged = super::protein::translated_part(&ref_protein)
        == super::protein::translated_part(&alt_protein);

    // Refinement of the base label (initiator lost, stop gained/lost) lives with
    // the substitution rules in `variants::consequence`, so the two paths cannot
    // drift apart again.
    let mut change_type = crate::variants::consequence::indel_change_type(
        &ref_protein,
        &alt_protein,
        frameshift,
        base_change_type,
    );
    // An in-frame indel that leaves the made protein byte for byte the same
    // changed nothing. An insertion inside the terminal stop codon that keeps a
    // stop there is the case: it was reported as a residue inserted one past the
    // protein's last one, at MODERATE impact, on a row whose own reference and
    // alternate codons were both a stop, while the substitution producing that
    // very same stop codon was correctly called synonymous.
    let protein_change_is_past_the_stop =
        !frameshift && protein_unchanged && change_type == base_change_type;
    if protein_change_is_past_the_stop {
        change_type = ChangeType::Synonymous;
    }

    let (protein_change, local_change) = if protein_change_is_past_the_stop {
        ("Synonymous".to_string(), "Synonymous".to_string())
    } else if touches_phase_skipped_bases(gene, variant) {
        ("Unknown".to_string(), "Unknown".to_string())
    } else {
        describe_protein_change(&ref_protein, &alt_protein, gene.protein_offset, frameshift)
    };

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

/// Whether the variant reaches into the bases the declared phase skips.
///
/// Those bases code, but their codon began in a neighbouring exon, so a single
/// CDS row cannot rebuild it. An indel that reaches into them changes the length
/// of a codon this row cannot see, and translating the phase-trimmed window then
/// describes a protein change that never happened: a three-base deletion came out
/// as a five-residue delins because the window lost two bases rather than three.
/// The length change is still known, so the row keeps its in-frame or frameshift
/// classification and says the residues themselves are unknown.
fn touches_phase_skipped_bases(gene: &Gene, variant: &VcfPosition) -> bool {
    if gene.phase == 0 || !gene.cds_segments.is_empty() {
        return false;
    }
    let phase = gene.phase as usize;
    let (skip_start, skip_end) = match gene.strand {
        crate::variants::Strand::Plus => (gene.start, gene.start + phase - 1),
        crate::variants::Strand::Minus => (gene.end.saturating_sub(phase - 1), gene.end),
    };
    // The interval as it stands. Widening it by one to make room for an
    // insertion's anchor also widened it for a deletion, so a deletion beginning
    // on the first translatable base was reported as touching bases it leaves
    // alone and its residues came out unknown. An insertion is already handled
    // here: its bases land in the gap after the anchor, and `overlaps_interval`
    // counts an anchor only while that gap is still inside the interval.
    variant.overlaps_interval(skip_start, skip_end)
}

/// How many coding bases the variant adds to or removes from the feature.
///
/// Measured over the whole feature, phase included. The phase-trimmed window
/// exists because the codon those leading bases belong to began in a neighbouring
/// exon and cannot be rebuilt from this row, which is a reason not to translate
/// them, not a reason to call them non-coding: the declared phase is precisely
/// the statement that they code. Measuring the delta over the trimmed window
/// counted a deletion of three coding bases as two whenever it reached into that
/// region, so an in-frame deletion was reported as a frameshift and dragged
/// `(fs)` onto every downstream codon of the feature.
pub(super) fn coding_delta_for_variant(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
) -> Option<isize> {
    let mut whole_feature = gene.clone();
    whole_feature.phase = 0;
    let ref_cds = coding_sequence_for_gene(&whole_feature, reference)?;
    let alt_cds = apply_allele_to_feature(&whole_feature, reference, variant)?;
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
