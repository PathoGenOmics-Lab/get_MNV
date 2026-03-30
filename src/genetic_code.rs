//! NCBI genetic code translation tables.
//!
//! Supports tables commonly used in genomics:
//! - 1: Standard (default for eukaryotes)
//! - 2: Vertebrate Mitochondrial
//! - 3: Yeast Mitochondrial
//! - 4: Mold, Protozoan, Coelenterate Mitochondrial; Mycoplasma/Spiroplasma
//! - 5: Invertebrate Mitochondrial
//! - 6: Ciliate, Dasycladacean, Hexamita Nuclear
//! - 11: Bacterial, Archaeal, Plant Plastid
//! - 12: Alternative Yeast Nuclear
//! - 25: Candidate Division SR1 and Gracilibacteria
//!
//! Reference: <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>

/// A genetic code translation table identified by NCBI table number.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeneticCode {
    table: u8,
}

/// All supported NCBI table numbers.
pub const SUPPORTED_TABLES: &[u8] = &[1, 2, 3, 4, 5, 6, 11, 12, 25];

impl GeneticCode {
    /// Create a genetic code from an NCBI table number.
    /// Returns `None` if the table is not supported.
    ///
    /// # Examples
    /// ```
    /// use get_mnv::genetic_code::GeneticCode;
    /// assert!(GeneticCode::new(11).is_some());
    /// assert!(GeneticCode::new(99).is_none());
    /// ```
    pub fn new(table: u8) -> Option<Self> {
        if SUPPORTED_TABLES.contains(&table) {
            Some(Self { table })
        } else {
            None
        }
    }

    /// Return the NCBI table number.
    pub fn table_number(self) -> u8 {
        self.table
    }

    /// Translate a DNA codon (3 uppercase bytes) to a single-character amino
    /// acid using this genetic code.  Returns `'X'` for ambiguous/unknown
    /// codons.
    ///
    /// # Examples
    /// ```
    /// use get_mnv::genetic_code::GeneticCode;
    /// let standard = GeneticCode::new(1).unwrap();
    /// assert_eq!(standard.translate(b"ATG"), 'M');
    /// assert_eq!(standard.translate(b"TAA"), '*');
    ///
    /// let vert_mito = GeneticCode::new(2).unwrap();
    /// assert_eq!(vert_mito.translate(b"TGA"), 'W'); // stop in table 1
    /// assert_eq!(vert_mito.translate(b"AGA"), '*'); // Arg in table 1
    /// ```
    pub fn translate(self, codon: &[u8; 3]) -> char {
        // Apply table-specific overrides first, then fall through to standard.
        match self.table {
            2 => match codon {
                // Vertebrate Mitochondrial
                b"TGA" => return 'W',
                b"AGA" | b"AGG" => return '*',
                _ => {}
            },
            3 => match codon {
                // Yeast Mitochondrial
                b"TGA" => return 'W',
                b"CTT" | b"CTC" | b"CTA" | b"CTG" => return 'T',
                _ => {}
            },
            4 => {
                // Mold, Protozoan, Coelenterate Mitochondrial; Mycoplasma
                if codon == b"TGA" {
                    return 'W';
                }
            }
            5 => match codon {
                // Invertebrate Mitochondrial
                b"TGA" => return 'W',
                b"AGA" | b"AGG" => return 'S',
                _ => {}
            },
            6 => match codon {
                // Ciliate Nuclear
                b"TAA" | b"TAG" => return 'Q',
                _ => {}
            },
            12 => {
                // Alternative Yeast Nuclear
                if codon == b"CTG" {
                    return 'S';
                }
            }
            25 => {
                // Candidate Division SR1 and Gracilibacteria
                if codon == b"TGA" {
                    return 'G';
                }
            }
            // Tables 1 and 11 have identical sense-codon translations.
            1 | 11 => {}
            _ => {}
        }
        // Standard genetic code (shared by tables 1, 11, and as fallback).
        translate_standard(codon)
    }

    /// Translate a DNA sequence (first 3 bases) to a single-character amino
    /// acid string.  Returns `"X"` if the sequence is shorter than 3 bases.
    pub fn translate_seq(self, seq: &[u8]) -> String {
        if seq.len() < 3 {
            return "X".to_string();
        }
        let codon: [u8; 3] = [
            seq[0].to_ascii_uppercase(),
            seq[1].to_ascii_uppercase(),
            seq[2].to_ascii_uppercase(),
        ];
        self.translate(&codon).to_string()
    }
}

impl Default for GeneticCode {
    /// Default is table 11 (Bacterial/Archaeal/Plant Plastid) since the
    /// primary use case is prokaryotic genomes.
    fn default() -> Self {
        Self { table: 11 }
    }
}

impl std::fmt::Display for GeneticCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self.table {
            1 => "Standard",
            2 => "Vertebrate Mitochondrial",
            3 => "Yeast Mitochondrial",
            4 => "Mold/Protozoan/Coelenterate Mito; Mycoplasma",
            5 => "Invertebrate Mitochondrial",
            6 => "Ciliate Nuclear",
            11 => "Bacterial/Archaeal/Plant Plastid",
            12 => "Alternative Yeast Nuclear",
            25 => "SR1/Gracilibacteria",
            _ => "Unknown",
        };
        write!(f, "{} (NCBI table {})", name, self.table)
    }
}

/// Standard genetic code translation (NCBI table 1).
fn translate_standard(codon: &[u8; 3]) -> char {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_supported_tables_construct() {
        for &table in SUPPORTED_TABLES {
            assert!(GeneticCode::new(table).is_some(), "table {table}");
        }
    }

    #[test]
    fn test_unsupported_table() {
        assert!(GeneticCode::new(99).is_none());
        assert!(GeneticCode::new(0).is_none());
    }

    #[test]
    fn test_default_is_table_11() {
        assert_eq!(GeneticCode::default().table_number(), 11);
    }

    #[test]
    fn test_standard_codons_shared() {
        // Tables 1 and 11 should produce identical results for all 64 codons.
        let t1 = GeneticCode::new(1).unwrap();
        let t11 = GeneticCode::new(11).unwrap();
        let bases = [b'A', b'T', b'C', b'G'];
        for &a in &bases {
            for &b in &bases {
                for &c in &bases {
                    let codon = [a, b, c];
                    assert_eq!(
                        t1.translate(&codon),
                        t11.translate(&codon),
                        "Mismatch for codon {}{}{}",
                        a as char,
                        b as char,
                        c as char
                    );
                }
            }
        }
    }

    #[test]
    fn test_table2_vertebrate_mito() {
        let gc = GeneticCode::new(2).unwrap();
        assert_eq!(gc.translate(b"TGA"), 'W');
        assert_eq!(gc.translate(b"AGA"), '*');
        assert_eq!(gc.translate(b"AGG"), '*');
        // Standard codons still work
        assert_eq!(gc.translate(b"ATG"), 'M');
        assert_eq!(gc.translate(b"TAA"), '*');
    }

    #[test]
    fn test_table3_yeast_mito() {
        let gc = GeneticCode::new(3).unwrap();
        assert_eq!(gc.translate(b"TGA"), 'W');
        assert_eq!(gc.translate(b"CTT"), 'T');
        assert_eq!(gc.translate(b"CTC"), 'T');
        assert_eq!(gc.translate(b"CTA"), 'T');
        assert_eq!(gc.translate(b"CTG"), 'T');
    }

    #[test]
    fn test_table4_mold_mito() {
        let gc = GeneticCode::new(4).unwrap();
        assert_eq!(gc.translate(b"TGA"), 'W');
        // AGA still Arg (not changed)
        assert_eq!(gc.translate(b"AGA"), 'R');
    }

    #[test]
    fn test_table5_invertebrate_mito() {
        let gc = GeneticCode::new(5).unwrap();
        assert_eq!(gc.translate(b"TGA"), 'W');
        assert_eq!(gc.translate(b"AGA"), 'S');
        assert_eq!(gc.translate(b"AGG"), 'S');
    }

    #[test]
    fn test_table6_ciliate() {
        let gc = GeneticCode::new(6).unwrap();
        assert_eq!(gc.translate(b"TAA"), 'Q');
        assert_eq!(gc.translate(b"TAG"), 'Q');
        // TGA is still stop
        assert_eq!(gc.translate(b"TGA"), '*');
    }

    #[test]
    fn test_table12_alt_yeast() {
        let gc = GeneticCode::new(12).unwrap();
        assert_eq!(gc.translate(b"CTG"), 'S');
        // Other CTN still Leu
        assert_eq!(gc.translate(b"CTT"), 'L');
    }

    #[test]
    fn test_table25_sr1() {
        let gc = GeneticCode::new(25).unwrap();
        assert_eq!(gc.translate(b"TGA"), 'G');
    }

    #[test]
    fn test_translate_seq() {
        let gc = GeneticCode::default();
        assert_eq!(gc.translate_seq(b"ATG"), "M");
        assert_eq!(gc.translate_seq(b"atg"), "M"); // handles lowercase
        assert_eq!(gc.translate_seq(b"AT"), "X"); // too short
    }

    #[test]
    fn test_display() {
        let gc = GeneticCode::new(11).unwrap();
        assert_eq!(gc.to_string(), "Bacterial/Archaeal/Plant Plastid (NCBI table 11)");
    }
}
