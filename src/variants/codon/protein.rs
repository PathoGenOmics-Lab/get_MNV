//! Translation and protein-change description helpers.

use crate::utils::aa_three_letter;

pub(super) fn translate_cds(seq: &str, genetic_code: crate::genetic_code::GeneticCode) -> String {
    let mut protein = String::with_capacity(seq.len() / 3);
    for codon in seq.as_bytes().chunks_exact(3) {
        let codon = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        protein.push(genetic_code.translate(&codon));
    }
    protein
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

        if ref_mid.is_empty() {
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
    fn inframe_deletion_and_delins_keep_residue_aware_hgvs() {
        // Regression guard: the residue-aware del/delins forms are unchanged.
        let (del, _) = describe_protein_change("MKLF", "MKF", 0, false);
        assert_eq!(del, "Leu3del");
        let (delins, _) = describe_protein_change("MKLF", "MKWWF", 0, false);
        assert_eq!(delins, "Leu3delinsTrpTrp");
    }
}
