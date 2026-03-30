//! Core domain types: variant enums, gene/SNP/codon structs, and display traits.

use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum VariantType {
    Snp,
    Mnv,
    SnpMnv,
    Indel,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum ChangeType {
    Unknown,
    Synonymous,
    NonSynonymous,
    StopGained,
    StopLost,
    IndelOverlap,
    FrameshiftSynonymous,
    FrameshiftNonSynonymous,
    FrameshiftStopGained,
    FrameshiftStopLost,
    FrameshiftUnknown,
    FrameshiftIndel,
    InFrameIndel,
}

impl ChangeType {
    /// Parse a human-readable label into a `ChangeType`.
    ///
    /// # Examples
    /// ```
    /// use get_mnv::variants::ChangeType;
    /// assert_eq!(ChangeType::from_label("Synonymous"), ChangeType::Synonymous);
    /// assert_eq!(ChangeType::from_label("unknown_label"), ChangeType::Unknown);
    /// ```
    pub fn from_label(label: &str) -> Self {
        match label {
            "Synonymous" => ChangeType::Synonymous,
            "Non-synonymous" => ChangeType::NonSynonymous,
            "Stop gained" => ChangeType::StopGained,
            "Stop lost" => ChangeType::StopLost,
            "Unknown" => ChangeType::Unknown,
            "Indel overlap" => ChangeType::IndelOverlap,
            "Synonymous (frameshift)" => ChangeType::FrameshiftSynonymous,
            "Non-synonymous (frameshift)" => ChangeType::FrameshiftNonSynonymous,
            "Stop gained (frameshift)" => ChangeType::FrameshiftStopGained,
            "Stop lost (frameshift)" => ChangeType::FrameshiftStopLost,
            "Unknown (frameshift)" => ChangeType::FrameshiftUnknown,
            "Frameshift Indel" => ChangeType::FrameshiftIndel,
            "In-frame Indel" => ChangeType::InFrameIndel,
            _ => ChangeType::Unknown,
        }
    }

    /// Convert to the frameshift variant of this change type.
    ///
    /// # Examples
    /// ```
    /// use get_mnv::variants::ChangeType;
    /// assert_eq!(
    ///     ChangeType::Synonymous.with_frameshift(),
    ///     ChangeType::FrameshiftSynonymous
    /// );
    /// ```
    pub fn with_frameshift(self) -> Self {
        match self {
            ChangeType::Synonymous => ChangeType::FrameshiftSynonymous,
            ChangeType::NonSynonymous => ChangeType::FrameshiftNonSynonymous,
            ChangeType::StopGained => ChangeType::FrameshiftStopGained,
            ChangeType::StopLost => ChangeType::FrameshiftStopLost,
            ChangeType::Unknown => ChangeType::FrameshiftUnknown,
            ChangeType::IndelOverlap => ChangeType::IndelOverlap,
            ChangeType::FrameshiftSynonymous => ChangeType::FrameshiftSynonymous,
            ChangeType::FrameshiftNonSynonymous => ChangeType::FrameshiftNonSynonymous,
            ChangeType::FrameshiftStopGained => ChangeType::FrameshiftStopGained,
            ChangeType::FrameshiftStopLost => ChangeType::FrameshiftStopLost,
            ChangeType::FrameshiftUnknown => ChangeType::FrameshiftUnknown,
            ChangeType::FrameshiftIndel => ChangeType::FrameshiftIndel,
            ChangeType::InFrameIndel => ChangeType::InFrameIndel,
        }
    }
}

impl VariantType {
    /// Return the short display label.
    ///
    /// # Examples
    /// ```
    /// use get_mnv::variants::VariantType;
    /// assert_eq!(VariantType::Snp.as_str(), "SNP");
    /// assert_eq!(VariantType::Mnv.as_str(), "MNV");
    /// ```
    pub fn as_str(self) -> &'static str {
        match self {
            VariantType::Snp => "SNP",
            VariantType::Mnv => "MNV",
            VariantType::SnpMnv => "SNP/MNV",
            VariantType::Indel => "INDEL",
        }
    }
}

impl Display for VariantType {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.write_str(self.as_str())
    }
}

impl Display for ChangeType {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        let value = match self {
            ChangeType::Unknown => "Unknown",
            ChangeType::Synonymous => "Synonymous",
            ChangeType::NonSynonymous => "Non-synonymous",
            ChangeType::StopGained => "Stop gained",
            ChangeType::StopLost => "Stop lost",
            ChangeType::IndelOverlap => "Indel overlap",
            ChangeType::FrameshiftSynonymous => "Synonymous (frameshift)",
            ChangeType::FrameshiftNonSynonymous => "Non-synonymous (frameshift)",
            ChangeType::FrameshiftStopGained => "Stop gained (frameshift)",
            ChangeType::FrameshiftStopLost => "Stop lost (frameshift)",
            ChangeType::FrameshiftUnknown => "Unknown (frameshift)",
            ChangeType::FrameshiftIndel => "Frameshift Indel",
            ChangeType::InFrameIndel => "In-frame Indel",
        };
        f.write_str(value)
    }
}

impl FromStr for Strand {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Plus),
            "-" => Ok(Strand::Minus),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Gene {
    pub name: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Snp {
    pub index: usize,
    pub position: usize,
    pub ref_base: String,
    pub base: String,
    pub original_dp: Option<usize>,
    pub original_freq: Option<f64>,
    pub original_info: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CodonInfo {
    pub codon_list: Vec<Snp>,
    pub original_codon: String,
    pub gene_name: String,
    pub gene_start: usize,
    pub gene_end: usize,
    pub codon_start: usize,
    pub codon_end: usize,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct VariantInfo {
    pub chrom: String,
    pub gene: String,
    pub positions: Vec<usize>,
    pub ref_bases: Vec<String>,
    pub base_changes: Vec<String>,
    pub aa_changes: Vec<String>,
    pub snp_aa_changes: Vec<String>,
    pub variant_type: VariantType,
    pub change_type: ChangeType,
    pub snp_reads: Option<Vec<usize>>,
    pub snp_forward_reads: Option<Vec<usize>>,
    pub snp_reverse_reads: Option<Vec<usize>>,
    pub mnv_reads: Option<usize>,
    pub mnv_forward_reads: Option<usize>,
    pub mnv_reverse_reads: Option<usize>,
    pub mnv_total_reads: Option<usize>,
    pub total_reads: Option<Vec<usize>>,
    pub total_forward_reads: Option<Vec<usize>>,
    pub total_reverse_reads: Option<Vec<usize>>,
    pub mnv_total_forward_reads: Option<usize>,
    pub mnv_total_reverse_reads: Option<usize>,
    pub ref_codon: Option<String>,
    pub snp_codon: Option<String>,
    pub mnv_codon: Option<String>,
    pub original_dp: Option<Vec<usize>>,
    pub original_freq: Option<Vec<f64>>,
    pub original_info: Option<String>,
}


#[cfg(test)]
mod tests {
    use super::*;

    // ---- VariantType ----

    #[test]
    fn test_variant_type_display() {
        assert_eq!(VariantType::Snp.to_string(), "SNP");
        assert_eq!(VariantType::Mnv.to_string(), "MNV");
        assert_eq!(VariantType::SnpMnv.to_string(), "SNP/MNV");
        assert_eq!(VariantType::Indel.to_string(), "INDEL");
    }

    #[test]
    fn test_variant_type_as_str() {
        assert_eq!(VariantType::Snp.as_str(), "SNP");
    }

    // ---- ChangeType ----

    #[test]
    fn test_change_type_display_roundtrip() {
        let types = [
            ChangeType::Synonymous,
            ChangeType::NonSynonymous,
            ChangeType::StopGained,
            ChangeType::StopLost,
            ChangeType::Unknown,
            ChangeType::IndelOverlap,
            ChangeType::FrameshiftSynonymous,
            ChangeType::FrameshiftNonSynonymous,
            ChangeType::FrameshiftStopGained,
            ChangeType::FrameshiftStopLost,
            ChangeType::FrameshiftUnknown,
            ChangeType::FrameshiftIndel,
            ChangeType::InFrameIndel,
        ];
        for ct in &types {
            let label = ct.to_string();
            assert_eq!(ChangeType::from_label(&label), *ct, "roundtrip failed for {label}");
        }
    }

    #[test]
    fn test_change_type_from_unknown_label() {
        assert_eq!(ChangeType::from_label("garbage"), ChangeType::Unknown);
    }

    #[test]
    fn test_with_frameshift() {
        assert_eq!(ChangeType::Synonymous.with_frameshift(), ChangeType::FrameshiftSynonymous);
        assert_eq!(ChangeType::NonSynonymous.with_frameshift(), ChangeType::FrameshiftNonSynonymous);
        assert_eq!(ChangeType::StopGained.with_frameshift(), ChangeType::FrameshiftStopGained);
        assert_eq!(ChangeType::StopLost.with_frameshift(), ChangeType::FrameshiftStopLost);
        assert_eq!(ChangeType::Unknown.with_frameshift(), ChangeType::FrameshiftUnknown);
        // Frameshift of frameshift returns same
        assert_eq!(ChangeType::FrameshiftSynonymous.with_frameshift(), ChangeType::FrameshiftSynonymous);
    }

    // ---- Strand ----

    #[test]
    fn test_strand_from_str() {
        assert_eq!("+".parse::<Strand>().unwrap(), Strand::Plus);
        assert_eq!("-".parse::<Strand>().unwrap(), Strand::Minus);
        assert!("x".parse::<Strand>().is_err());
        assert!("".parse::<Strand>().is_err());
    }

    // ---- ChangeType ordering ----

    #[test]
    fn test_change_type_ord() {
        // PartialOrd is derived, verify it doesn't panic
        assert!(ChangeType::Synonymous < ChangeType::NonSynonymous);
    }

    // ---- VariantType ordering ----

    #[test]
    fn test_variant_type_ord() {
        assert!(VariantType::Snp < VariantType::Mnv);
    }
}
