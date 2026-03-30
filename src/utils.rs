//! Utility functions: amino acid change classification, IUPAC three-letter
//! conversion, reverse complement, and codon translation.

/// Translate a DNA codon (3 bytes) to a single-character amino acid.
/// Uses the standard genetic code (NCBI table 1). Returns 'X' for
/// unknown/ambiguous codons.
fn translate_codon(codon: &[u8; 3]) -> char {
    match codon {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}

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

    let position = if aa_change_1_letter.len() > 2 {
        &aa_change_1_letter[1..aa_change_1_letter.len() - 1]
    } else {
        "?"
    };

    format!("{first}{position}{last}")
}

pub fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            // IUPAC ambiguity codes
            b'R' | b'r' => b'Y',
            b'Y' | b'y' => b'R',
            b'S' | b's' => b'S',
            b'W' | b'w' => b'W',
            b'K' | b'k' => b'M',
            b'M' | b'm' => b'K',
            b'B' | b'b' => b'V',
            b'V' | b'v' => b'B',
            b'D' | b'd' => b'H',
            b'H' | b'h' => b'D',
            other => other,
        })
        .map(|b| b as char)
        .collect()
}

/// Translate a DNA sequence to a protein string. Only the first codon
/// (3 bases) is translated — matching the original single-codon usage
/// throughout the codebase.
pub fn process_translate(seq: &[u8]) -> String {
    if seq.len() < 3 {
        return "X".to_string();
    }
    let codon: [u8; 3] = [
        seq[0].to_ascii_uppercase(),
        seq[1].to_ascii_uppercase(),
        seq[2].to_ascii_uppercase(),
    ];
    translate_codon(&codon).to_string()
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
