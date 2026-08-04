//! The one place where a protein-level change becomes a [`ChangeType`].
//!
//! Substitutions and indels reach this question from different directions: a
//! substitution knows the residue before and after at a known protein position,
//! while an indel only knows the reference and alternate translations. They were
//! previously classified by two separate implementations, and the initiator rule
//! existed on the substitution side only, so an in-frame indel that deleted the
//! start codon was reported as a plain `inframe_deletion` at `MODERATE` impact
//! while the equivalent SNV was `start_lost` at `HIGH`.
//!
//! Both rule sets live here so that a rule stated once cannot drift, and
//! `both_paths_agree_on_the_initiator_rule` holds them to it.

use crate::variants::ChangeType;

/// Whether the initiator methionine is gone.
///
/// The single definition of "the start codon is lost", shared by both entry
/// points. `reference` and `alternate` are the residues at protein position 1.
/// A stop at position 1 is deliberately excluded: it is reported as a stop
/// gained, which is the more specific call, and both paths must agree on that.
fn initiator_lost(reference: char, alternate: char) -> bool {
    reference == 'M' && alternate != 'M' && alternate != '*'
}

/// First residue of a translated codon, `X` when there is none.
///
/// Callers hold the translation as a string; this is the one place that turns it
/// into the residue the classifier works on, so no caller has to remember that
/// an empty translation means "ambiguous" rather than "no change".
pub fn first_residue(translation: &str) -> char {
    translation.chars().next().unwrap_or('X')
}

/// Consequence of a single-residue substitution at 1-based `protein_position`.
///
/// Residues are one-letter codes, with `*` for a stop and `X` for anything
/// ambiguous or untranslatable.
pub fn substitution_change_type(
    reference: char,
    alternate: char,
    protein_position: usize,
) -> ChangeType {
    if reference == 'X' || alternate == 'X' {
        return ChangeType::Unknown;
    }
    if protein_position == 1 && initiator_lost(reference, alternate) {
        return ChangeType::StartLost;
    }
    if reference == alternate {
        ChangeType::Synonymous
    } else if alternate == '*' {
        ChangeType::StopGained
    } else if reference == '*' {
        ChangeType::StopLost
    } else {
        ChangeType::NonSynonymous
    }
}

/// Whether an indel destroys the initiator codon.
///
/// `None` when the reference protein does not start with Met, which is how a
/// partial or phase-shifted CDS model presents itself: there is no initiator to
/// lose. A Met left at position 1 still initiates translation.
fn indel_start_effect(reference_protein: &str, alternate_protein: &str) -> Option<ChangeType> {
    let reference = reference_protein.chars().next()?;
    // An alternate that translates to nothing has lost its initiator too.
    let alternate = alternate_protein.chars().next().unwrap_or('-');
    initiator_lost(reference, alternate).then_some(ChangeType::StartLost)
}

/// Whether an in-frame indel creates or removes a stop codon, by comparing the
/// number of stop residues in the two translations. Counting stops (rather than
/// comparing positions) avoids a false "stop gained" for an ordinary in-frame
/// deletion, whose terminal stop simply shifts to a lower index.
fn indel_stop_effect(reference_protein: &str, alternate_protein: &str) -> Option<ChangeType> {
    let reference_stops = reference_protein.matches('*').count();
    let alternate_stops = alternate_protein.matches('*').count();
    if alternate_stops > reference_stops {
        Some(ChangeType::StopGained)
    } else if reference_stops > 0 && alternate_stops < reference_stops {
        Some(ChangeType::StopLost)
    } else {
        None
    }
}

/// Consequence of an indel, refining `base` (frameshift or in-frame) from the
/// reference and alternate translations.
///
/// Losing the start outranks any stop change: without an initiator the ORF is
/// not translated from here at all. A frameshift keeps `base`, because it almost
/// always introduces a downstream stop, so flagging it as "stop gained" would be
/// uninformative, and it is already `HIGH` impact.
pub(crate) fn indel_change_type(
    reference_protein: &str,
    alternate_protein: &str,
    frameshift: bool,
    base: ChangeType,
) -> ChangeType {
    if frameshift {
        return base;
    }
    indel_start_effect(reference_protein, alternate_protein)
        .or_else(|| indel_stop_effect(reference_protein, alternate_protein))
        .unwrap_or(base)
}

#[cfg(test)]
mod tests {
    use super::*;

    const RESIDUES: [char; 21] = [
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
        'W', 'Y', '*',
    ];

    /// The rule that broke: it existed on the substitution side only, so an
    /// in-frame indel removing the start codon came out as a generic in-frame
    /// indel. Both entry points must now answer identically for every residue.
    #[test]
    fn both_paths_agree_on_the_initiator_rule() {
        for &alternate in RESIDUES.iter() {
            let by_substitution = substitution_change_type('M', alternate, 1);
            // The same event seen by the indel path: an in-frame change that
            // leaves `alternate` where the initiator Met was.
            let by_indel = indel_change_type(
                &format!("M{}", "KLA*"),
                &format!("{alternate}{}", "KLA*"),
                false,
                ChangeType::InFrameIndel,
            );
            let substitution_says_start_lost = by_substitution == ChangeType::StartLost;
            let indel_says_start_lost = by_indel == ChangeType::StartLost;
            assert_eq!(
                substitution_says_start_lost, indel_says_start_lost,
                "M -> {alternate} at position 1: substitution said {by_substitution:?}, \
                 indel said {by_indel:?}"
            );
        }
    }

    #[test]
    fn a_stop_at_position_one_is_a_stop_gained_on_both_paths() {
        // Documented precedence: Met1 -> stop is the more specific stop gained,
        // not a start lost, and both paths must make the same call.
        assert_eq!(
            substitution_change_type('M', '*', 1),
            ChangeType::StopGained
        );
        assert_eq!(
            indel_change_type("MKLA*", "*KLA*", false, ChangeType::InFrameIndel),
            ChangeType::StopGained
        );
    }

    #[test]
    fn the_initiator_rule_only_applies_at_position_one() {
        assert_eq!(substitution_change_type('M', 'T', 1), ChangeType::StartLost);
        assert_eq!(
            substitution_change_type('M', 'T', 50),
            ChangeType::NonSynonymous
        );
        assert_eq!(
            substitution_change_type('M', 'M', 1),
            ChangeType::Synonymous
        );
    }

    #[test]
    fn substitutions_classify_the_ordinary_cases() {
        assert_eq!(
            substitution_change_type('A', 'T', 15),
            ChangeType::NonSynonymous
        );
        assert_eq!(
            substitution_change_type('G', 'G', 92),
            ChangeType::Synonymous
        );
        assert_eq!(
            substitution_change_type('E', '*', 112),
            ChangeType::StopGained
        );
        assert_eq!(substitution_change_type('*', 'Q', 50), ChangeType::StopLost);
        assert_eq!(substitution_change_type('X', 'T', 15), ChangeType::Unknown);
        assert_eq!(substitution_change_type('A', 'X', 15), ChangeType::Unknown);
    }

    #[test]
    fn indels_keep_the_frameshift_label() {
        // A frameshift is already HIGH impact and almost always introduces a
        // downstream stop, so it is not relabelled.
        assert_eq!(
            indel_change_type("MKLA*", "KLA*", true, ChangeType::FrameshiftIndel),
            ChangeType::FrameshiftIndel
        );
    }

    #[test]
    fn indels_detect_stop_changes_and_leave_the_rest_alone() {
        assert_eq!(
            indel_change_type("MKLA*", "MKL*A*", false, ChangeType::InFrameIndel),
            ChangeType::StopGained
        );
        assert_eq!(
            indel_change_type("MKLA*", "MKLA", false, ChangeType::InFrameIndel),
            ChangeType::StopLost
        );
        assert_eq!(
            indel_change_type("MKLA*", "MKLGGA*", false, ChangeType::InFrameIndel),
            ChangeType::InFrameIndel
        );
    }

    #[test]
    fn a_cds_model_without_an_initiator_is_left_alone() {
        // A partial or phase-shifted CDS does not begin at the initiator, so
        // there is no start to lose and the generic label stands.
        assert_eq!(
            indel_change_type("KLA*", "LA*", false, ChangeType::InFrameIndel),
            ChangeType::InFrameIndel
        );
    }
}
