//! Translation and protein-change description helpers.

use crate::utils::aa_three_letter;

pub(super) fn translate_cds(seq: &str, genetic_code: crate::genetic_code::GeneticCode) -> String {
    let mut protein = String::with_capacity(seq.len() / 3);
    for codon in seq.as_bytes().as_chunks::<3>().0 {
        let codon = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        protein.push(genetic_code.translate(&codon));
    }
    protein
}

/// [`translate_cds`] for a sequence that begins at the protein's first codon,
/// read as the ribosome reads it: an initiation codon there is methionine. A
/// gene that starts on GTG otherwise began with valine, and an indel could not
/// be seen to remove its start.
pub(super) fn translate_from_start(
    seq: &str,
    genetic_code: crate::genetic_code::GeneticCode,
) -> String {
    let mut protein = translate_cds(seq, genetic_code);
    if let Some(first) = seq.as_bytes().get(0..3) {
        let codon = [
            first[0].to_ascii_uppercase(),
            first[1].to_ascii_uppercase(),
            first[2].to_ascii_uppercase(),
        ];
        if genetic_code.is_start(&codon) {
            protein.replace_range(0..1, "M");
        }
    }
    protein
}

/// The part of a translation that is actually made: everything up to and
/// including the first stop.
///
/// `translate_cds` renders every codon of the sequence it is given, stop codons
/// included, which is what the stop-position searches want. Comparing two such
/// translations base for base compares residues downstream of the stop, which no
/// ribosome ever reaches: an insertion inside a gene's terminal stop codon that
/// leaves a stop in place came out as an inserted residue one past the protein's
/// last one, on a row whose own reference and alternate codons were both TAA.
pub(super) fn translated_part(protein: &str) -> &str {
    match protein.find('*') {
        Some(stop) => &protein[..=stop],
        None => protein,
    }
}

pub(super) fn aa_segment_three_letter(segment: &[char]) -> String {
    segment
        .iter()
        .map(|aa| aa_three_letter(*aa))
        .collect::<Vec<_>>()
        .join("")
}

pub(super) fn describe_protein_change(
    ref_protein: &str,
    alt_protein: &str,
    protein_offset: usize,
    frameshift: bool,
) -> (String, String) {
    if ref_protein == alt_protein {
        if frameshift {
            let ref_chars: Vec<char> = ref_protein.chars().collect();
            let idx = ref_chars.len().saturating_sub(1);
            let local_pos = idx + 1;
            let protein_pos = protein_offset + local_pos;
            let ref_aa = ref_chars.get(idx).copied().unwrap_or('X');
            return (
                format!("{}{}fs", aa_three_letter(ref_aa), protein_pos),
                format!("{}{}fs", aa_three_letter(ref_aa), local_pos),
            );
        }
        return ("Synonymous".to_string(), "Synonymous".to_string());
    }

    let ref_chars: Vec<char> = ref_protein.chars().collect();
    let alt_chars: Vec<char> = alt_protein.chars().collect();
    let mut prefix = 0usize;
    while prefix < ref_chars.len()
        && prefix < alt_chars.len()
        && ref_chars[prefix] == alt_chars[prefix]
    {
        prefix += 1;
    }

    let local_pos = prefix + 1;
    let protein_pos = protein_offset + local_pos;
    let ref_aa = ref_chars.get(prefix).copied().unwrap_or('X');
    let alt_aa = alt_chars.get(prefix).copied().unwrap_or('X');

    if frameshift {
        let local = format!(
            "{}{}{}fs",
            aa_three_letter(ref_aa),
            local_pos,
            aa_three_letter(alt_aa)
        );
        let protein = format!(
            "{}{}{}fs",
            aa_three_letter(ref_aa),
            protein_pos,
            aa_three_letter(alt_aa)
        );
        return (protein, local);
    }

    let mut suffix = 0usize;
    while suffix + prefix < ref_chars.len()
        && suffix + prefix < alt_chars.len()
        && ref_chars[ref_chars.len() - 1 - suffix] == alt_chars[alt_chars.len() - 1 - suffix]
    {
        suffix += 1;
    }

    let ref_mid_end = ref_chars.len().saturating_sub(suffix);
    let alt_mid_end = alt_chars.len().saturating_sub(suffix);
    let ref_mid = &ref_chars[prefix..ref_mid_end];
    let alt_mid = &alt_chars[prefix..alt_mid_end];

    let build = |start_pos: usize| -> String {
        if ref_mid.len() == 1 && alt_mid.len() == 1 {
            return format!(
                "{}{}{}",
                aa_three_letter(ref_mid[0]),
                start_pos,
                aa_three_letter(alt_mid[0])
            );
        }

        let end_pos = start_pos + ref_mid.len().saturating_sub(1);
        let ref_start = ref_mid.first().copied().unwrap_or(ref_aa);
        let ref_end = ref_mid.last().copied().unwrap_or(ref_aa);
        let ref_range = if start_pos == end_pos {
            format!("{}{}", aa_three_letter(ref_start), start_pos)
        } else {
            format!(
                "{}{}_{}{}",
                aa_three_letter(ref_start),
                start_pos,
                aa_three_letter(ref_end),
                end_pos
            )
        };

        if ref_mid.is_empty() && repeats_residues_before(&ref_chars, prefix, alt_mid) {
            // HGVS writes an insertion that repeats the residues just before it
            // as a duplication of them, `Phe433dup` or `His68_Leu70dup`, as the
            // `c.` and `g.` descriptors of the same row already do. `prefix` is
            // the longest common prefix, so the insertion already sits at its
            // 3'-most place and the copy it repeats is the 3'-most copy.
            let last = start_pos - 1;
            let first = start_pos - alt_mid.len();
            let first_aa = aa_three_letter(alt_mid[0]);
            if first == last {
                format!("{first_aa}{last}dup")
            } else {
                let last_aa = aa_three_letter(alt_mid[alt_mid.len() - 1]);
                format!("{first_aa}{first}_{last_aa}{last}dup")
            }
        } else if ref_mid.is_empty() {
            // HGVS protein insertions name the two flanking (unchanged) residues
            // bracketing the insertion point, e.g. `Lys2_Phe3insGly`, matching
            // the residue-aware form already used for `del` / `delins`. The
            // flanks are the last residue of the common prefix (`prefix - 1`) and
            // the first of the common suffix (`prefix`). A `wrapping_sub` on a
            // zero prefix yields `usize::MAX`, so `get` returns `None` and we fall
            // back to bare positions for the degenerate N-/C-terminal insertions.
            let anchor = start_pos.saturating_sub(1);
            let inserted = aa_segment_three_letter(alt_mid);
            let left = ref_chars.get(prefix.wrapping_sub(1)).copied();
            let right = ref_chars.get(prefix).copied();
            match (left, right) {
                (Some(left), Some(right)) => format!(
                    "{}{anchor}_{}{start_pos}ins{inserted}",
                    aa_three_letter(left),
                    aa_three_letter(right)
                ),
                _ => format!("{anchor}_{start_pos}ins{inserted}"),
            }
        } else if alt_mid.is_empty() {
            format!("{ref_range}del")
        } else {
            format!("{ref_range}delins{}", aa_segment_three_letter(alt_mid))
        }
    };

    (build(protein_pos), build(local_pos))
}

/// Whether residues inserted after the first `prefix` residues repeat the ones
/// just before them, which HGVS writes as a duplication of those residues. A
/// residue that could not be read (`X`, a codon with an `N`) repeats nothing.
fn repeats_residues_before(ref_chars: &[char], prefix: usize, inserted: &[char]) -> bool {
    !inserted.is_empty()
        && !inserted.contains(&'X')
        && prefix >= inserted.len()
        && ref_chars[prefix - inserted.len()..prefix] == *inserted
}

pub(super) fn complement_base(base: char) -> char {
    match base.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        'N' => 'N',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::describe_protein_change;

    #[test]
    fn inframe_insertion_names_both_flanking_residues() {
        // Ref MKF, alt inserts Gly between Lys2 and Phe3 (-> MKGF). Proper HGVS
        // protein insertion names both flanking residues, like `del`/`delins`.
        let (protein, local) = describe_protein_change("MKF", "MKGF", 0, false);
        assert_eq!(protein, "Lys2_Phe3insGly");
        assert_eq!(local, "Lys2_Phe3insGly");
    }

    #[test]
    fn inframe_insertion_applies_protein_offset_to_positions_only() {
        // The exon protein offset shifts the numbering but not the residue names.
        let (protein, local) = describe_protein_change("MKF", "MKGF", 100, false);
        assert_eq!(protein, "Lys102_Phe103insGly");
        assert_eq!(local, "Lys2_Phe3insGly");
    }

    #[test]
    fn multi_residue_insertion_names_flanks_once() {
        // Insert Gly-Ser between Lys2 and Phe3.
        let (protein, _) = describe_protein_change("MKF", "MKGSF", 0, false);
        assert_eq!(protein, "Lys2_Phe3insGlySer");
    }

    #[test]
    fn c_terminal_insertion_falls_back_to_bare_positions() {
        // Insertion after the last residue has no right-hand flank; keep the
        // bare-position form rather than fabricating a residue.
        let (protein, _) = describe_protein_change("MKF", "MKFG", 0, false);
        assert_eq!(protein, "3_4insGly");
    }

    #[test]
    fn insertion_repeating_the_residues_before_it_is_a_duplication() {
        // rpoB 761102 A>ATTC: a second Phe after Phe433, which HGVS and the WHO
        // catalogue write p.Phe433dup rather than p.Phe433_Met434insPhe.
        let (protein, local) = describe_protein_change("MKFM", "MKFFM", 430, false);
        assert_eq!(protein, "Phe433dup");
        assert_eq!(local, "Phe3dup");
        // A run of residues: His68_Leu70dup.
        let (run, _) = describe_protein_change("MHQLA", "MHQLHQLA", 66, false);
        assert_eq!(run, "His68_Leu70dup");
        // Beside the stop, and with no residue after it at all, the copy is
        // still a duplication: no flank needs naming.
        let (before_stop, _) = describe_protein_change("MKF*", "MKFF*", 0, false);
        assert_eq!(before_stop, "Phe3dup");
        let (c_terminal, _) = describe_protein_change("MKF", "MKFF", 0, false);
        assert_eq!(c_terminal, "Phe3dup");
    }

    #[test]
    fn duplication_in_a_repeat_names_its_three_prime_copy() {
        // Phe3 and Phe4 are both Phe: the copy is the 3'-most, Phe4.
        let (single, _) = describe_protein_change("MKFFM", "MKFFFM", 0, false);
        assert_eq!(single, "Phe4dup");
        // His-Gln twice: the copy is the second pair, His4_Gln5.
        let (pair, _) = describe_protein_change("MHQHQA", "MHQHQHQA", 0, false);
        assert_eq!(pair, "His4_Gln5dup");
        // Phe-Met put in front of Phe3-Met4 makes the same protein as Phe-Met
        // put after them, so it too is their duplication.
        let (in_front, _) = describe_protein_change("MKFMA", "MKFMFMA", 0, false);
        assert_eq!(in_front, "Phe3_Met4dup");
    }

    #[test]
    fn insertion_of_residues_not_just_before_it_stays_an_insertion() {
        // Lys is in the protein, but not just before the insertion point.
        let (protein, _) = describe_protein_change("MKFM", "MKFKM", 0, false);
        assert_eq!(protein, "Phe3_Met4insLys");
        // Phe-Phe after a single Phe repeats more than the residues before it.
        let (longer, _) = describe_protein_change("MKFM", "MKFFFM", 0, false);
        assert_eq!(longer, "Phe3_Met4insPhePhe");
        // Two residues that could not be read are not known to be the same.
        let (unread, _) = describe_protein_change("MKXM", "MKXXM", 0, false);
        assert_eq!(unread, "X3_Met4insX");
    }

    #[test]
    fn inframe_deletion_and_delins_keep_residue_aware_hgvs() {
        // Regression guard: the residue-aware del/delins forms are unchanged.
        let (del, _) = describe_protein_change("MKLF", "MKF", 0, false);
        assert_eq!(del, "Leu3del");
        let (delins, _) = describe_protein_change("MKLF", "MKWWF", 0, false);
        assert_eq!(delins, "Leu3delinsTrpTrp");
    }
}
