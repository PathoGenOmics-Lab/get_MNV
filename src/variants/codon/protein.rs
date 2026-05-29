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
            let anchor = start_pos.saturating_sub(1);
            format!(
                "{anchor}_{start_pos}ins{}",
                aa_segment_three_letter(alt_mid)
            )
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
