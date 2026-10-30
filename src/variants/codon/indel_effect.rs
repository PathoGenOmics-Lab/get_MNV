//! Protein-level effect of indels (frameshift, in-frame, stop gained/lost).

use crate::io::VcfPosition;
use crate::variants::event::AlleleComponentKind;
use crate::variants::{ChangeType, Gene, Strand};

use super::allele_apply::apply_allele_to_feature;
use super::gene_path::effective_bounds;
use super::protein::{describe_protein_change, translate_cds, translate_from_start};
use super::transcript_model::{coding_sequence_for_gene, first_touched_transcript_offset};

/// The genomic bases this variant's components actually touch.
///
/// `affected_start` and `affected_end` bound the record's own reference span,
/// and for an insertion that span is the padding rather than the junction:
/// `102 TAA>TAAG` and `104 A>AG` describe the same insertion after base 104,
/// and the first begins two bases earlier. The deletion branch of the
/// decomposition sets those bounds to the first deleted base and agrees with
/// its component; the insertion branch does not. Keying the displayed codon on
/// the record's span therefore put a padded insertion in the codon before the
/// one it changes, a codon it leaves alone, on a row whose own amino-acid
/// column named the right one. The components carry the anchor, which is the
/// ruler the transcript path already reads.
fn component_touched_bounds(variant: &VcfPosition) -> Option<(usize, usize)> {
    let changed = variant.changed_positions();
    let lo = *changed.iter().min()?;
    let hi = *changed.iter().max()?;
    Some((lo, hi))
}

/// What an indel does to the protein.
pub(super) struct IndelEffect {
    pub change_type: ChangeType,
    pub aa_changes: Vec<String>,
    pub aa_changes_local: Vec<String>,
    pub ref_codon: Option<String>,
    pub alt_codon: Option<String>,
    /// The indel removes the protein's initiator. Said separately because a
    /// frameshift that does so keeps its frameshift label.
    pub initiator_lost: bool,
}

pub(super) fn protein_effect_for_indel(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
    genetic_code: crate::genetic_code::GeneticCode,
) -> IndelEffect {
    // Base classification (frameshift vs in-frame, symbolic handling). When the
    // alternate CDS can be reconstructed we may refine this to Stop gained/lost.
    let base_change_type = indel_change_type_for_variant(gene, reference, variant);

    let Some(ref_cds) = coding_sequence_for_gene(gene, reference) else {
        return IndelEffect {
            change_type: base_change_type,
            aa_changes: vec!["Unknown".to_string()],
            aa_changes_local: vec!["Unknown".to_string()],
            ref_codon: None,
            alt_codon: None,
            initiator_lost: false,
        };
    };
    let Some(alt_cds) = apply_allele_to_feature(gene, reference, variant) else {
        return IndelEffect {
            change_type: base_change_type,
            aa_changes: vec!["Unknown".to_string()],
            aa_changes_local: vec!["Unknown".to_string()],
            ref_codon: None,
            alt_codon: None,
            initiator_lost: false,
        };
    };

    // Whether the frame shifts is a property of the whole feature, so it is
    // measured there rather than on the phase-trimmed sequences translated
    // below: those omit the leading coding bases and undercounted a deletion
    // that reached into them.
    let frameshift = coding_delta_for_variant(gene, reference, variant)
        .map(|delta| delta % 3 != 0)
        .unwrap_or((alt_cds.len() as isize - ref_cds.len() as isize) % 3 != 0);
    // A feature that begins at the protein's first codon is read from its
    // initiator, so a GTG start is Met and removing it is seen as such. A CDS row
    // further into the protein, or one that starts mid-codon, has no initiator.
    let starts_protein = gene.starts_protein();
    let translate = |seq: &str| {
        if starts_protein {
            translate_from_start(seq, genetic_code)
        } else {
            translate_cds(seq, genetic_code)
        }
    };
    let ref_protein = translate(&ref_cds);
    // The alternate is read from its initiator only while it still has one. A
    // deletion that takes the start codon away puts another codon first, and
    // reading that one as an initiator too made a GTG or TTG from inside the
    // gene Met again: fgd1 losing its GTG start, the 32 bases before it and 134
    // more codons came out as the in-frame `Ala2_Val136del`. Nothing is
    // predicted to start there, so the alternate makes no protein from this
    // start, which the rules below read as a lost initiator and write `Met1?`.
    let alt_protein = if starts_protein
        && ref_protein.starts_with('M')
        && removes_initiation_codon(gene, variant, &ref_cds, &alt_cds)
    {
        String::new()
    } else {
        translate(&alt_cds)
    };
    let initiator_lost = starts_protein
        && crate::variants::consequence::initiator_lost_between(&ref_protein, &alt_protein);
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
    // very same stop codon was correctly called synonymous. A stop past the
    // first one counts for nothing either: an insertion that copies the terminal
    // stop codon (`TAA` to `TAATAA`) adds a stop no ribosome reaches, and
    // counting the stops of the whole translation called it a stop gained at
    // HIGH impact, which a catalogue's `LoF` entry then matched.
    let protein_change_is_past_the_stop = !frameshift && protein_unchanged;
    if protein_change_is_past_the_stop {
        change_type = ChangeType::Synonymous;
    }

    let (protein_change, local_change) = if initiator_lost {
        // HGVS: with the initiator gone, where translation starts, if it does,
        // cannot be predicted, so no residues are named.
        ("Met1?".to_string(), "Met1?".to_string())
    } else if protein_change_is_past_the_stop {
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
        return IndelEffect {
            change_type,
            aa_changes: vec![protein_change],
            aa_changes_local: vec![local_change],
            ref_codon: Some(ref_codon),
            alt_codon,
            initiator_lost,
        };
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
    let (touched_start, touched_end) =
        component_touched_bounds(variant).unwrap_or((event.affected_start, event.affected_end));
    let touched_lo = touched_start.max(eff_start);
    let touched_hi = touched_end.min(eff_end);
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

    IndelEffect {
        change_type,
        aa_changes: vec![protein_change],
        aa_changes_local: vec![local_change],
        ref_codon,
        alt_codon,
        initiator_lost,
    }
}

/// The genomic positions of the protein's first codon, in the order the
/// ribosome reads them. Empty for a feature shorter than a codon.
fn first_codon_positions(gene: &Gene) -> Vec<usize> {
    let segments = if gene.cds_segments.is_empty() {
        vec![effective_bounds(gene)]
    } else {
        gene.cds_segments
            .iter()
            .map(|segment| (segment.start, segment.end))
            .collect()
    };
    let mut positions = Vec::with_capacity(3);
    for (start, end) in segments {
        let wanted = 3 - positions.len();
        match gene.strand {
            Strand::Plus => positions.extend((start..=end).take(wanted)),
            Strand::Minus => positions.extend((start..=end).rev().take(wanted)),
        }
        if positions.len() == 3 {
            return positions;
        }
    }
    Vec::new()
}

/// Whether the variant takes the protein's initiation codon away.
///
/// It does when it deletes one of the codon's bases or inserts between them and
/// the codon now first is not the same three bases. Deleting `GTGGCA` from the
/// front of `GTGGCAGTG` is deleting the `GCAGTG` behind the start, which leaves
/// the start where it was and is written that way (`Ala2_Val3del`). A deletion
/// that also removes the base before the start codon takes it away whatever
/// codon comes up in its place: the sequence the ribosome finds there is no
/// longer the one it started from. A substitution leaves the codon where it was,
/// and the initiator rules read whatever it makes of it.
fn removes_initiation_codon(
    gene: &Gene,
    variant: &VcfPosition,
    ref_cds: &str,
    alt_cds: &str,
) -> bool {
    let codon = first_codon_positions(gene);
    let Some(&first) = codon.first() else {
        return false;
    };
    let before = match gene.strand {
        Strand::Plus => first.checked_sub(1),
        Strand::Minus => first.checked_add(1),
    };
    let mut touched = false;
    let mut reaches_before = false;
    for component in variant.event().components {
        match component.kind {
            AlleleComponentKind::Deletion => {
                let last = component.position + component.ref_allele.len().saturating_sub(1);
                let deleted = component.position..=last;
                if codon.iter().any(|position| deleted.contains(position)) {
                    touched = true;
                    reaches_before |= before.is_some_and(|base| deleted.contains(&base));
                }
            }
            // Inserted bases follow their anchor along the genome, so they land
            // inside the codon when the anchor and the base after it are both
            // the codon's.
            AlleleComponentKind::Insertion => {
                touched |= codon.contains(&component.position)
                    && codon.contains(&(component.position + 1));
            }
            _ => {}
        }
    }
    let same_codon_first = match (ref_cds.get(..3), alt_cds.get(..3)) {
        (Some(reference), Some(alternate)) => reference.eq_ignore_ascii_case(alternate),
        _ => false,
    };
    touched && (reaches_before || !same_codon_first)
}

/// Where the coding sequence sees a deletion its record writes across the
/// protein's start: at its 3'-most place along the transcript within the first
/// exon, as HGVS places it. `None` for any other record, which is seen where it
/// is written.
///
/// A deletion in a repeat can be written anywhere along the repeat, and the
/// feature clips a deletion at its edges. Left-aligned, as `bcftools norm` and
/// `--left-align` write it, deleting the second ATG of a gene that starts
/// `ATGATG` right after a G begins on that G, so only two of its three bases
/// were coding: a frameshift that lost the start, where the same deletion
/// written inside the gene was `Met2del`. Placed 3'-most, a deletion lies
/// inside the coding sequence whenever any spelling of it does, and every
/// spelling of it is read the same. One that still reaches the base before the
/// start codon there really removes the start.
pub(super) fn coding_placement(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
) -> Option<VcfPosition> {
    if !gene.starts_protein() {
        return None;
    }
    let event = variant.event();
    let [component] = event.components.as_slice() else {
        return None;
    };
    if component.kind != AlleleComponentKind::Deletion {
        return None;
    }
    let first = *first_codon_positions(gene).first()?;
    let before = match gene.strand {
        Strand::Plus => first.checked_sub(1)?,
        Strand::Minus => first.checked_add(1)?,
    };
    let deleted = component.position..component.position + component.ref_allele.len();
    if !(deleted.contains(&before) && deleted.contains(&first)) {
        return None;
    }
    let (exon_start, exon_end) = match gene.cds_segments.first() {
        Some(segment) => (segment.start, segment.end),
        None => effective_bounds(gene),
    };
    let (pos, ref_allele, alt_allele) = crate::variants::hgvs::shift_within(
        reference.sequence,
        variant.record_start,
        &variant.ref_allele,
        &variant.alt_allele,
        gene.strand,
        exon_start,
        exon_end,
    );
    Some(VcfPosition {
        record_start: pos,
        ref_allele,
        alt_allele,
        ..variant.clone()
    })
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
    // A symbolic allele names a structural event without spelling out its
    // sequence, and get_MNV does not read SVTYPE, END or SVLEN, so how many
    // coding bases it adds or removes is not known here. Answering "the frame
    // shifts" turned that absence into a positive claim: an <INV>, which
    // preserves length exactly and can never shift a frame, was reported as
    // frameshift_variant at HIGH, and a <DEL> was too whatever its SVLEN said.
    if variant.alt_allele.starts_with('<') {
        return ChangeType::Unknown;
    }
    let frameshift = if let Some(delta) = coding_delta_for_variant(gene, reference, variant) {
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
