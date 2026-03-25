use protein_translate::translate;

pub fn determine_change_type(aa_change: &str) -> String {
    if aa_change.is_empty() {
        return "Unknown".to_string();
    }

    let original = aa_change.chars().next().unwrap_or('X');
    let mutated = aa_change.chars().last().unwrap_or('X');

    if original == mutated {
        "Synonymous".to_string()
    } else if mutated == '*' {
        "Stop gained".to_string()
    } else if original == '*' {
        "Stop lost".to_string()
    } else {
        "Non-synonymous".to_string()
    }
}

fn aa_three_letter(code: char) -> &'static str {
    match code {
        'A' => "Ala",
        'C' => "Cys",
        'D' => "Asp",
        'E' => "Glu",
        'F' => "Phe",
        'G' => "Gly",
        'H' => "His",
        'I' => "Ile",
        'K' => "Lys",
        'L' => "Leu",
        'M' => "Met",
        'N' => "Asn",
        'P' => "Pro",
        'Q' => "Gln",
        'R' => "Arg",
        'S' => "Ser",
        'T' => "Thr",
        'V' => "Val",
        'W' => "Trp",
        'Y' => "Tyr",
        '*' => "*",
        _ => "X",
    }
}

pub fn iupac_aa(aa_change_1_letter: &str) -> String {
    if aa_change_1_letter.len() < 3 {
        return "Invalid_format".to_string();
    }

    let first_char = aa_change_1_letter.chars().next().unwrap_or('X');
    let last_char = aa_change_1_letter.chars().last().unwrap_or('X');

    let first = aa_three_letter(first_char);
    let last = aa_three_letter(last_char);

    let position = &aa_change_1_letter[1..aa_change_1_letter.len() - 1];

    format!("{}{}{}", first, position, last)
}

pub fn reverse_complement(seq: &str) -> String {
    let rc = bio::alphabets::dna::revcomp(seq.as_bytes());
    String::from_utf8(rc)
        .expect("Critical error: reverse-complement sequence contains invalid UTF-8.")
}

pub fn process_translate(seq: &[u8]) -> String {
    translate(seq)
}

// ==========================================
// Unit tests
// Run with: cargo test
// ==========================================
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determine_change_type() {
        // Cover all major amino-acid change categories.
        assert_eq!(determine_change_type("A15T"), "Non-synonymous");
        assert_eq!(determine_change_type("G92G"), "Synonymous");
        assert_eq!(determine_change_type("E112*"), "Stop gained");
        assert_eq!(determine_change_type("*50Q"), "Stop lost");
        assert_eq!(determine_change_type(""), "Unknown");
    }

    #[test]
    fn test_iupac_aa_expansion() {
        // Validate conversion from one-letter to three-letter amino-acid notation.
        assert_eq!(iupac_aa("A15T"), "Ala15Thr");
        assert_eq!(iupac_aa("G92*"), "Gly92*");
        assert_eq!(iupac_aa("E112Q"), "Glu112Gln");
        assert_eq!(iupac_aa(""), "Invalid_format");
    }

    #[test]
    fn test_reverse_complement() {
        // Validate DNA reverse-complement behavior.
        assert_eq!(reverse_complement("ATGC"), "GCAT");
        assert_eq!(reverse_complement("CCGGAATT"), "AATTCCGG");
    }

    #[test]
    fn test_process_translate() {
        // ATG = Methionine (M), TAA = stop codon (*)
        assert_eq!(process_translate(b"ATG"), "M");
        assert_eq!(process_translate(b"TAA"), "*");
    }
}
