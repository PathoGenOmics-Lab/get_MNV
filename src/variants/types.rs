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
    StartLost,
    StopGained,
    StopLost,
    IndelOverlap,
    FrameshiftSynonymous,
    FrameshiftNonSynonymous,
    FrameshiftStopGained,
    FrameshiftStopLost,
    FrameshiftDownstreamOfStop,
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
            "Start lost" => ChangeType::StartLost,
            "Stop gained" => ChangeType::StopGained,
            "Stop lost" => ChangeType::StopLost,
            "Unknown" => ChangeType::Unknown,
            "Indel overlap" => ChangeType::IndelOverlap,
            "Synonymous (frameshift)" => ChangeType::FrameshiftSynonymous,
            "Non-synonymous (frameshift)" => ChangeType::FrameshiftNonSynonymous,
            "Stop gained (frameshift)" => ChangeType::FrameshiftStopGained,
            "Stop lost (frameshift)" => ChangeType::FrameshiftStopLost,
            "Downstream of premature stop" => ChangeType::FrameshiftDownstreamOfStop,
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
            // The initiator codon (protein position 1) has nothing upstream, so a
            // start-lost is never reached through downstream frameshift propagation.
            ChangeType::StartLost => ChangeType::StartLost,
            ChangeType::StopGained => ChangeType::FrameshiftStopGained,
            ChangeType::StopLost => ChangeType::FrameshiftStopLost,
            ChangeType::Unknown => ChangeType::FrameshiftUnknown,
            ChangeType::IndelOverlap => ChangeType::IndelOverlap,
            ChangeType::FrameshiftSynonymous => ChangeType::FrameshiftSynonymous,
            ChangeType::FrameshiftNonSynonymous => ChangeType::FrameshiftNonSynonymous,
            ChangeType::FrameshiftStopGained => ChangeType::FrameshiftStopGained,
            ChangeType::FrameshiftStopLost => ChangeType::FrameshiftStopLost,
            ChangeType::FrameshiftDownstreamOfStop => ChangeType::FrameshiftDownstreamOfStop,
            ChangeType::FrameshiftUnknown => ChangeType::FrameshiftUnknown,
            ChangeType::FrameshiftIndel => ChangeType::FrameshiftIndel,
            ChangeType::InFrameIndel => ChangeType::InFrameIndel,
        }
    }
}

/// How the amino-acid consequence of a *combined* MNV compares with the
/// consequence of its individual constituent SNVs. Only meaningful for
/// MNV / SNP-MNV records; single SNVs are always `NotApplicable`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum ConsequenceShift {
    #[default]
    NotApplicable,
    /// Same severity as the most severe individual SNV.
    Concordant,
    /// More severe than any single SNV alone (e.g. two individually-synonymous
    /// SNVs producing a non-synonymous residue) — what per-SNV annotation misses.
    Gained,
    /// Less severe than the most severe single SNV (e.g. a nonsense SNV rescued
    /// to a missense by its neighbour).
    Masked,
}

impl ConsequenceShift {
    pub fn as_str(self) -> &'static str {
        match self {
            ConsequenceShift::NotApplicable => "-",
            ConsequenceShift::Concordant => "Concordant",
            ConsequenceShift::Gained => "MNV-gained",
            ConsequenceShift::Masked => "MNV-masked",
        }
    }
}

impl Display for ConsequenceShift {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.write_str(self.as_str())
    }
}

/// Nonsense-mediated decay (NMD) prediction for a premature termination codon,
/// following the 50-nucleotide rule: a PTC located more than 50 nt upstream of
/// the last exon-exon junction of a multi-exon transcript triggers NMD; a PTC in
/// the last exon or within 50 nt of the last junction escapes it. Single-exon
/// transcripts have no junction and are never annotated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NmdPrediction {
    Triggering,
    Escaping,
}

impl NmdPrediction {
    pub fn as_str(self) -> &'static str {
        match self {
            NmdPrediction::Triggering => "NMD-triggering",
            NmdPrediction::Escaping => "NMD-escaping",
        }
    }
}

impl Display for NmdPrediction {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.write_str(self.as_str())
    }
}

/// Splice consequence of a variant near an internal exon-exon junction, using
/// the standard Sequence Ontology terms. The two intronic bases at each end of
/// an intron are the essential donor / acceptor sites (HIGH impact); the exon's
/// first/last 3 bases and the intron's 3rd-8th bases are the wider splice region
/// (LOW impact). Only defined for the multi-exon transcript CDS model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SpliceConsequence {
    Donor,
    Acceptor,
    Region,
}

impl SpliceConsequence {
    /// Sequence Ontology term.
    pub fn as_str(self) -> &'static str {
        match self {
            SpliceConsequence::Donor => "splice_donor_variant",
            SpliceConsequence::Acceptor => "splice_acceptor_variant",
            SpliceConsequence::Region => "splice_region_variant",
        }
    }

    /// Predicted impact: the essential donor / acceptor sites are HIGH, the
    /// wider splice region is LOW.
    pub fn impact(self) -> &'static str {
        match self {
            SpliceConsequence::Donor | SpliceConsequence::Acceptor => "HIGH",
            SpliceConsequence::Region => "LOW",
        }
    }

    /// Rank used to keep the most severe match when a position is near more than
    /// one junction (donor / acceptor outrank the wider region).
    pub fn severity(self) -> u8 {
        match self {
            SpliceConsequence::Region => 1,
            SpliceConsequence::Donor | SpliceConsequence::Acceptor => 2,
        }
    }
}

/// Extra, biologist-facing annotations layered on top of the core variant call
/// (Grantham distance, MNV-vs-SNV consequence shift). Grouped in one `Default`
/// struct so it can be added to `VariantInfo` without touching every
/// constructor.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct VariantAnnotations {
    /// Grantham distance of a genuine missense change (None otherwise).
    #[serde(default)]
    pub grantham: Option<u16>,
    /// How the combined MNV consequence compares with its individual SNVs.
    #[serde(default)]
    pub consequence_shift: ConsequenceShift,
    /// COSMIC-style doublet-base-substitution class for an MNV of two adjacent
    /// single-base substitutions, e.g. `CC>TT` (reverse-complement collapsed).
    /// `None` for single SNVs, indels, and non-adjacent or >2-SNV MNVs.
    #[serde(default)]
    pub dbs_class: Option<String>,
    /// 50-nt-rule nonsense-mediated decay prediction for a premature stop.
    /// `None` unless the variant introduces a premature stop in a multi-exon
    /// transcript.
    #[serde(default)]
    pub nmd: Option<NmdPrediction>,
    /// HGVS coding (`c.`) descriptor for a coding substitution (SNV or MNV),
    /// e.g. `c.30A>G` or `c.[28G>A;30T>C]`. `None` for indels and non-coding
    /// variants; the genomic (`g.`) descriptor is derived at output time.
    #[serde(default)]
    pub hgvs_c: Option<String>,
    /// Splice consequence near an internal exon-exon junction (multi-exon
    /// transcripts only). Folded into the SO term / impact rather than reported
    /// as a separate column: an exonic coding change is combined (e.g.
    /// `missense_variant&splice_region_variant`), an intronic site stands alone.
    #[serde(default)]
    pub splice: Option<SpliceConsequence>,
    /// Set when the variant lies inside a real intron of the named gene but
    /// outside its splice sites, which makes it an `intron_variant`. Splice-site
    /// positions are intronic too, but they carry the more specific `splice`
    /// annotation instead, so the two are mutually exclusive by construction.
    #[serde(default)]
    pub intronic: bool,
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
            ChangeType::StartLost => "Start lost",
            ChangeType::StopGained => "Stop gained",
            ChangeType::StopLost => "Stop lost",
            ChangeType::IndelOverlap => "Indel overlap",
            ChangeType::FrameshiftSynonymous => "Synonymous (frameshift)",
            ChangeType::FrameshiftNonSynonymous => "Non-synonymous (frameshift)",
            ChangeType::FrameshiftStopGained => "Stop gained (frameshift)",
            ChangeType::FrameshiftStopLost => "Stop lost (frameshift)",
            ChangeType::FrameshiftDownstreamOfStop => "Downstream of premature stop",
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
pub struct CdsSegment {
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone)]
pub struct Gene {
    pub name: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    /// GFF phase (column 8) for CDS features: 0, 1 or 2.
    /// Indicates how many bases must be removed from the start of the feature
    /// (or from the end, for minus strand features) to reach the first base of
    /// the first complete codon. Defaults to 0 for features without phase
    /// information (e.g. TSV gene files, gene/exon features).
    pub phase: u8,
    /// Amino-acid position of the first complete codon of this feature within
    /// the full protein coded by the parent transcript (0-based). For
    /// multi-exon eukaryotic CDS this is the cumulative count of complete
    /// codons contributed by all prior exons of the same transcript. For
    /// single-feature inputs (TSV gene files, non-CDS features, prokaryotic
    /// annotations) this is always 0 and the historical per-feature numbering
    /// is preserved.
    pub protein_offset: usize,
    /// Parent transcript identifier for GFF/GTF CDS-derived records. `None`
    /// for TSV annotations and non-transcript features.
    pub transcript_id: Option<String>,
    /// CDS segments belonging to the parent transcript in transcript order.
    /// Empty for legacy per-feature annotations. When present, codon grouping
    /// and indel protein effects are computed on the spliced CDS sequence.
    pub cds_segments: Vec<CdsSegment>,
}

impl Gene {
    /// Whether a real intron separates two consecutive CDS segments given in
    /// transcript order.
    ///
    /// Consecutive segments do not always imply splicing. A CDS annotated as a
    /// ribosomal-slippage join, such as SARS-CoV-2 ORF1ab
    /// (`join(266..13468,13468..21555)`), is one continuous reading frame whose
    /// segments abut or overlap: the ribosome re-reads the shared base rather
    /// than splicing anything out. Only a genuine gap between the segments is an
    /// intron, and only an intron carries donor/acceptor sites and an exon-exon
    /// junction.
    pub(crate) fn intron_separates(&self, earlier: &CdsSegment, later: &CdsSegment) -> bool {
        match self.strand {
            Strand::Plus => later.start > earlier.end + 1,
            Strand::Minus => earlier.start > later.end + 1,
        }
    }

    /// Whether an insertion anchored at genomic `position` lands at an internal
    /// exon-exon junction of the spliced CDS — i.e. immediately after the last
    /// base of an internal coding exon. The inserted bases then fall between two
    /// exons in the spliced transcript, so the insertion is coding, even though
    /// the exon's terminal base is excluded from the half-open exon interval used
    /// elsewhere. An anchor at the CDS terminus (after the last exon on the plus
    /// strand, before the first on the minus strand) is UTR, not an internal
    /// junction. `cds_segments` is in transcript order, so the terminal exon is
    /// the last entry on the plus strand and the first on the minus strand.
    pub(crate) fn insertion_at_internal_junction(&self, position: usize) -> bool {
        let segment_count = self.cds_segments.len();
        self.cds_segments.iter().enumerate().any(|(idx, segment)| {
            position == segment.end
                && match self.strand {
                    Strand::Plus => idx + 1 < segment_count,
                    Strand::Minus => idx > 0,
                }
        })
    }
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
    /// Amino-acid offset of the current feature within the full protein
    /// (0-based). Non-zero only for CDS exons of a multi-exon transcript; the
    /// final `aa_pos` reported in the output is computed as
    /// `protein_offset + local_aa_pos`.
    #[serde(default)]
    pub protein_offset: usize,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct VariantInfo {
    pub chrom: String,
    pub gene: String,
    pub positions: Vec<usize>,
    pub ref_bases: Vec<String>,
    pub base_changes: Vec<String>,
    /// Amino acid change in **full-protein** numbering. For multi-exon
    /// eukaryotic CDS this is the position relative to the protein N-terminus
    /// (compatible with VEP/ANNOVAR/SnpEff/UniProt). For prokaryotes,
    /// single-exon features and TSV gene files this is identical to
    /// `aa_changes_local` because the protein offset is 0.
    pub aa_changes: Vec<String>,
    pub snp_aa_changes: Vec<String>,
    /// Amino acid change in **per-feature** (exon-local) numbering — what
    /// get_MNV reported up to and including version 1.1.1, where each CDS
    /// row of the GFF was treated as an independent feature and the codon
    /// counter was reset to 1 at every exon. Kept alongside `aa_changes` so
    /// existing pipelines/scripts that rely on the old numbering can still
    /// consume both columns. Identical to `aa_changes` whenever the feature
    /// has no transcript context (`protein_offset == 0`).
    #[serde(default)]
    pub aa_changes_local: Vec<String>,
    #[serde(default)]
    pub snp_aa_changes_local: Vec<String>,
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
    /// Canonical allele event class derived from REF/ALT, e.g. `mnv`,
    /// `insertion`, `deletion`, `complex_indel`. Kept separate from
    /// `variant_type` so length-changing events that also contain SNV/MNV
    /// components can be represented without inventing false downstream MNVs.
    #[serde(default)]
    pub event_class: Option<String>,
    /// Human-readable decomposition of the original allele into event
    /// components such as `SNV:10:A>G`, `INS:10:+T` or `DEL:11-12:TG`.
    #[serde(default)]
    pub event_components: Vec<String>,
    /// Biologist-facing annotations (Grantham distance, MNV-vs-SNV consequence
    /// shift). Defaulted so older serialized records still deserialize.
    #[serde(default)]
    pub annotations: VariantAnnotations,
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
            assert_eq!(
                ChangeType::from_label(&label),
                *ct,
                "roundtrip failed for {label}"
            );
        }
    }

    #[test]
    fn test_change_type_from_unknown_label() {
        assert_eq!(ChangeType::from_label("garbage"), ChangeType::Unknown);
    }

    #[test]
    fn test_with_frameshift() {
        assert_eq!(
            ChangeType::Synonymous.with_frameshift(),
            ChangeType::FrameshiftSynonymous
        );
        assert_eq!(
            ChangeType::NonSynonymous.with_frameshift(),
            ChangeType::FrameshiftNonSynonymous
        );
        assert_eq!(
            ChangeType::StopGained.with_frameshift(),
            ChangeType::FrameshiftStopGained
        );
        assert_eq!(
            ChangeType::StopLost.with_frameshift(),
            ChangeType::FrameshiftStopLost
        );
        assert_eq!(
            ChangeType::Unknown.with_frameshift(),
            ChangeType::FrameshiftUnknown
        );
        // Frameshift of frameshift returns same
        assert_eq!(
            ChangeType::FrameshiftSynonymous.with_frameshift(),
            ChangeType::FrameshiftSynonymous
        );
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
