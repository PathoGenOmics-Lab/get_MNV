//! Utility functions: amino acid change classification, IUPAC three-letter
//! conversion, reverse complement, and codon translation.

/// Translate a DNA codon (3 bytes) to a single-character amino acid.
/// Uses the standard/bacterial genetic code (NCBI tables 1/11 — identical
/// for internal sense codons; alternative start codons are not relevant
/// here since we only translate individual codons, not initiation sites).
/// Returns 'X' for unknown/ambiguous codons.
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

    // If either amino acid is unknown/ambiguous, we cannot classify the change.
    if original == 'X' || mutated == 'X' {
        return "Unknown".to_string();
    }

    // Start lost: the initiator Met (protein position 1) changed to another
    // residue. A Met1 -> stop is left classified as "Stop gained".
    if original == 'M' && mutated != 'M' && mutated != '*' {
        let position = aa_change
            .get(1..aa_change.len().saturating_sub(1))
            .and_then(|s| s.parse::<usize>().ok());
        if position == Some(1) {
            return "Start lost".to_string();
        }
    }

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

pub(crate) fn aa_three_letter(code: char) -> &'static str {
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

/// Grantham (1974) physicochemical properties of the 20 standard amino acids:
/// `(composition, polarity, molecular volume)`.
fn grantham_properties(aa: char) -> Option<(f64, f64, f64)> {
    let props = match aa.to_ascii_uppercase() {
        'S' => (1.42, 9.2, 32.0),
        'R' => (0.65, 10.5, 124.0),
        'L' => (0.00, 4.9, 111.0),
        'P' => (0.39, 8.0, 32.5),
        'T' => (0.71, 8.6, 61.0),
        'A' => (0.00, 8.1, 31.0),
        'V' => (0.00, 5.9, 84.0),
        'G' => (0.74, 9.0, 3.0),
        'I' => (0.00, 5.2, 111.0),
        'F' => (0.00, 5.2, 132.0),
        'Y' => (0.20, 6.2, 136.0),
        'C' => (2.75, 5.5, 55.0),
        'H' => (0.58, 10.4, 96.0),
        'Q' => (0.89, 10.5, 85.0),
        'N' => (1.33, 11.6, 56.0),
        'K' => (0.33, 11.3, 119.0),
        'D' => (1.38, 13.0, 54.0),
        'E' => (0.92, 12.3, 83.0),
        'M' => (0.00, 5.7, 105.0),
        'W' => (0.13, 5.4, 170.0),
        _ => return None,
    };
    Some(props)
}

/// Grantham distance between two amino acids (single-letter codes), rounded to
/// the nearest integer. Returns `None` unless both are standard amino acids.
///
/// Uses Grantham's original formula `D = ρ · sqrt(α·Δc² + β·Δp² + γ·Δv²)` with
/// `α = 1.833`, `β = 0.1018`, `γ = 0.000399`, `ρ = 50.723`, reproducing the
/// published matrix (Ser↔Arg ≈ 110, Leu↔Ile ≈ 5, Cys↔Trp ≈ 215).
pub fn grantham_distance(a: char, b: char) -> Option<u16> {
    let (ca, pa, va) = grantham_properties(a)?;
    let (cb, pb, vb) = grantham_properties(b)?;
    const ALPHA: f64 = 1.833;
    const BETA: f64 = 0.1018;
    const GAMMA: f64 = 0.000399;
    const RHO: f64 = 50.723;
    let d = RHO
        * (ALPHA * (ca - cb).powi(2) + BETA * (pa - pb).powi(2) + GAMMA * (va - vb).powi(2)).sqrt();
    Some(d.round() as u16)
}

/// Grantham's conservation category for a distance: conservative (0–50),
/// moderately conservative (51–100), moderately radical (101–150) or
/// radical (>150).
pub fn grantham_category(distance: u16) -> &'static str {
    match distance {
        0..=50 => "conservative",
        51..=100 => "moderately conservative",
        101..=150 => "moderately radical",
        _ => "radical",
    }
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

/// Translate a DNA sequence to a protein string using the standard genetic
/// code (NCBI table 1/11). Only the first codon (3 bases) is translated.
///
/// **Prefer [`GeneticCode::translate_seq`](crate::genetic_code::GeneticCode::translate_seq)**
/// which supports all NCBI translation tables via `--translation-table`.
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
    fn test_grantham_distance_and_category() {
        // Published Grantham (1974) values, within ±1.
        assert!((grantham_distance('S', 'R').unwrap() as i32 - 110).abs() <= 1);
        assert!((grantham_distance('L', 'I').unwrap() as i32 - 5).abs() <= 1);
        assert!((grantham_distance('C', 'W').unwrap() as i32 - 215).abs() <= 1);
        assert_eq!(grantham_distance('D', 'K'), grantham_distance('K', 'D'));
        assert_eq!(grantham_distance('*', 'A'), None);
        assert_eq!(grantham_category(5), "conservative");
        assert_eq!(grantham_category(110), "moderately radical");
        assert_eq!(grantham_category(215), "radical");
    }

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
    fn test_determine_change_type_start_lost() {
        // Initiator Met (protein position 1) changed to another residue.
        assert_eq!(determine_change_type("M1T"), "Start lost");
        assert_eq!(determine_change_type("M1L"), "Start lost");
        // Internal Met (not the start codon) is an ordinary substitution.
        assert_eq!(determine_change_type("M50T"), "Non-synonymous");
        // Synonymous start, and start->stop, keep their usual classification.
        assert_eq!(determine_change_type("M1M"), "Synonymous");
        assert_eq!(determine_change_type("M1*"), "Stop gained");
    }

    #[test]
    fn test_determine_change_type_ambiguous_x() {
        // When either amino acid is X (ambiguous codon), classify as Unknown
        // rather than accidentally matching first==last.
        assert_eq!(determine_change_type("X15X"), "Unknown");
        assert_eq!(determine_change_type("A15X"), "Unknown");
        assert_eq!(determine_change_type("X15A"), "Unknown");
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
