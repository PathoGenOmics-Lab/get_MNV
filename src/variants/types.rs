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
    /// SNVs producing a non-synonymous residue), which is what per-SNV annotation misses.
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
    /// Why a variant that lies inside a gene has no codon annotation. Splice-site
    /// positions are intronic too, but they carry the more specific `splice`
    /// annotation instead, so the two are mutually exclusive by construction.
    #[serde(default)]
    pub non_coding: Option<NonCodingReason>,
    /// What the reads said about each upstream indel sharing molecules with
    /// this codon. Empty when no BAM was given or no upstream indel was in
    /// range. Without it a codon that is not labelled frameshifted looks the
    /// same whether the reads proved the indel is on other molecules or nobody
    /// ever asked.
    #[serde(default)]
    pub frameshift_linkage: Vec<FrameshiftLinkage>,
    /// What the input VCF's own phasing said about this row's alleles, when it
    /// said anything. A claim by the caller, kept separate from the read
    /// evidence rather than merged into it.
    #[serde(default)]
    pub declared_phase: Option<DeclaredPhaseCall>,
    /// Read-level linkage disequilibrium between this codon's substitutions.
    /// `None` when no molecule observed every position, or when one of the
    /// alleles is carried by all of them or none, which leaves nothing to
    /// correlate. See [`crate::variants::linkage`].
    #[serde(default)]
    pub linkage: Option<crate::variants::linkage::CodonLinkage>,
}

/// The input caller's phase claim for one row, and whether the reads refute it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct DeclaredPhaseCall {
    pub verdict: LinkageVerdict,
    /// The `PS` phase set the claim came from.
    pub phase_set: Option<usize>,
    /// Set when the reads leave the claim no room: a declared cis that no read
    /// carries whole, or a declared trans that every informative read carries.
    #[serde(default)]
    pub contradicted_by_reads: bool,
}

impl std::fmt::Display for DeclaredPhaseCall {
    /// `cis:12345`, or just `cis` when the caller phased without a phase set,
    /// with `|contradicted-by-reads` appended when the BAM refutes it.
    ///
    /// The separator is `|` rather than `;` because `;` ends an INFO field: a
    /// semicolon here has to be percent-encoded, leaving the reader with
    /// `cis:42%3Bcontradicted-by-reads` in the one field whose purpose is to
    /// catch a human's eye. `|` passes through intact, as it does for `COMP`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.phase_set {
            Some(phase_set) => write!(f, "{}:{}", self.verdict.as_str(), phase_set)?,
            None => write!(f, "{}", self.verdict.as_str())?,
        }
        if self.contradicted_by_reads {
            write!(f, "|contradicted-by-reads")?;
        }
        Ok(())
    }
}

/// What the reads say about two variants sharing molecules.
///
/// The thresholds behind these live with the read counting that applies them
/// (`read_count::IndelSnvLinkage`); this is the judged answer, not a rule to be
/// re-derived. Ordered by how much the answer constrains interpretation, which
/// is the order used to pick the pair worth reporting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum LinkageVerdict {
    /// Reads settle it: the two are on different molecules.
    Trans,
    /// Reads settle it: the two travel together.
    Cis,
    /// Too few reads span both loci to say either way.
    Unknown,
}

impl LinkageVerdict {
    pub fn as_str(self) -> &'static str {
        match self {
            LinkageVerdict::Trans => "trans",
            LinkageVerdict::Cis => "cis",
            LinkageVerdict::Unknown => "unknown",
        }
    }
}

/// Read-level evidence about one upstream indel and this codon.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FrameshiftLinkage {
    pub indel_position: usize,
    /// Reads carrying this codon's substitution that also span the indel locus.
    /// Only these can answer the question.
    pub informative_reads: usize,
    /// Of those, how many also carry the indel.
    pub cis_reads: usize,
    pub verdict: LinkageVerdict,
}

impl std::fmt::Display for FrameshiftLinkage {
    /// `trans:1234:0/18`: the verdict, the indel's position, and the reads
    /// behind it.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}/{}",
            self.verdict.as_str(),
            self.indel_position,
            self.cis_reads,
            self.informative_reads
        )
    }
}

/// Why a variant inside an annotated gene carries no amino-acid change. Both
/// cases are `MODIFIER`: the variant is in the transcript but not in a codon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NonCodingReason {
    /// Inside a real intron, away from the splice sites.
    Intron,
    /// Inside a transcript that is not translated at all (a non-coding RNA gene
    /// or a pseudogene).
    NonCodingTranscript,
    /// Inside coding sequence, but no whole codon could be built around it: the
    /// annotated feature's length is not a multiple of three, or its end runs
    /// past the contig, so the position falls in a trailing partial codon.
    ///
    /// Such a variant used to be dropped from every output with no warning,
    /// because the gene path could not annotate it and the intergenic fallback
    /// considered it covered by the gene. It is real and it is inside a gene;
    /// only its amino-acid effect is unknown.
    IncompleteCodon,
}

impl NonCodingReason {
    /// Sequence Ontology term.
    pub fn as_str(self) -> &'static str {
        match self {
            NonCodingReason::Intron => "intron_variant",
            NonCodingReason::NonCodingTranscript => "non_coding_transcript_exon_variant",
            NonCodingReason::IncompleteCodon => "coding_sequence_variant",
        }
    }
}

/// Whether an annotated feature is translated.
///
/// The GFF path selects CDS features and is therefore coding by construction.
/// The TSV format carries this in its optional sixth column: without it a
/// non-coding RNA gene is translated as though it were a protein, which invents
/// amino-acid changes for something that is never translated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Biotype {
    #[default]
    ProteinCoding,
    NonCoding,
}

impl Biotype {
    pub fn is_coding(self) -> bool {
        matches!(self, Biotype::ProteinCoding)
    }

    /// Parse a biotype from the TSV annotation. The vocabulary is closed on
    /// purpose: an unrecognised value is an error rather than a silent guess,
    /// because guessing wrong either invents a protein or hides a real one.
    pub fn parse(value: &str) -> Result<Self, String> {
        const CODING: [&str; 4] = ["protein_coding", "coding", "CDS", "mRNA"];
        const NON_CODING: [&str; 13] = [
            "ncRNA",
            "rRNA",
            "tRNA",
            "tmRNA",
            "miRNA",
            "snRNA",
            "snoRNA",
            "misc_RNA",
            "antisense_RNA",
            "SRP_RNA",
            "RNase_P_RNA",
            "non_coding",
            "pseudogene",
        ];
        let trimmed = value.trim();
        if CODING.iter().any(|c| c.eq_ignore_ascii_case(trimmed)) {
            return Ok(Biotype::ProteinCoding);
        }
        if NON_CODING.iter().any(|c| c.eq_ignore_ascii_case(trimmed)) {
            return Ok(Biotype::NonCoding);
        }
        Err(format!(
            "unknown biotype '{trimmed}'. Coding: {}. Non-coding: {}",
            CODING.join(", "),
            NON_CODING.join(", ")
        ))
    }
}

#[cfg(test)]
mod biotype_tests {
    use super::Biotype;

    #[test]
    fn coding_and_non_coding_vocabularies_are_recognised() {
        assert_eq!(Biotype::parse("protein_coding"), Ok(Biotype::ProteinCoding));
        assert_eq!(Biotype::parse("CDS"), Ok(Biotype::ProteinCoding));
        assert_eq!(Biotype::parse("ncRNA"), Ok(Biotype::NonCoding));
        assert_eq!(Biotype::parse("tRNA"), Ok(Biotype::NonCoding));
        assert_eq!(Biotype::parse("pseudogene"), Ok(Biotype::NonCoding));
        // Case and surrounding whitespace are not meaningful.
        assert_eq!(Biotype::parse("  rRNA  "), Ok(Biotype::NonCoding));
        assert_eq!(Biotype::parse("Protein_Coding"), Ok(Biotype::ProteinCoding));
    }

    #[test]
    fn an_unknown_biotype_is_an_error_rather_than_a_guess() {
        // A typo must not silently decide whether a feature gets translated: it
        // would either invent a protein or hide a real one.
        let err =
            Biotype::parse("protein-coding").expect_err("hyphen is not the accepted spelling");
        assert!(err.contains("unknown biotype"), "{err}");
        assert!(
            err.contains("protein_coding"),
            "the message lists what is accepted: {err}"
        );
        assert!(Biotype::parse("").is_err());
    }

    #[test]
    fn the_default_is_coding_so_existing_annotations_are_unchanged() {
        // Four-column TSV files predate this column and must keep translating.
        assert_eq!(Biotype::default(), Biotype::ProteinCoding);
        assert!(Biotype::default().is_coding());
        assert!(!Biotype::NonCoding.is_coding());
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
    /// Whether this feature is translated. Always coding on the GFF path, which
    /// selects CDS features; set from the TSV annotation's optional sixth
    /// column otherwise.
    pub biotype: Biotype,
}

/// A real intron of a spliced transcript: the genomic gap between two CDS
/// segments, together with where that exon-exon junction falls in the spliced
/// CDS and which exon bases flank it.
///
/// Derived once by [`Gene::introns`] so that splice-site classification,
/// intronic-position tests and the NMD junction all read the same geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Intron {
    /// Lowest genomic base of the intron.
    pub start: usize,
    /// Highest genomic base of the intron.
    pub end: usize,
    /// Spliced-CDS offset at which the exon after this intron begins, i.e. the
    /// exon-exon junction in transcript coordinates.
    pub junction_offset: usize,
    /// Terminal base of the exon before the intron, in transcript order: the
    /// donor boundary.
    pub exon_before_end: usize,
    /// First base of the exon after the intron, in transcript order: the
    /// acceptor boundary.
    pub exon_after_start: usize,
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

    /// The transcript's real introns, in transcript order.
    ///
    /// Splice sites, intronic positions and the NMD junction are all questions
    /// about the same geometry, and each used to re-derive it from
    /// `cds_segments` on its own. Deriving it once means the three answers
    /// cannot disagree about where an intron is, or about whether a pair of
    /// segments has one at all.
    ///
    /// Yields nothing for an unspliced transcript: a single CDS segment, or
    /// segments that abut or overlap in a ribosomal-slippage join.
    pub(crate) fn introns(&self) -> impl Iterator<Item = Intron> + '_ {
        self.cds_segments
            .windows(2)
            .scan(0usize, |offset_after, pair| {
                // `offset_after` accumulates the spliced length up to and
                // including `earlier`, which is where the next exon begins in
                // transcript coordinates.
                let (earlier, later) = (&pair[0], &pair[1]);
                *offset_after += earlier.end.saturating_sub(earlier.start) + 1;
                Some((*offset_after, earlier, later))
            })
            .filter_map(|(junction_offset, earlier, later)| {
                if !self.intron_separates(earlier, later) {
                    return None;
                }
                // `cds_segments` is in transcript order, so on the minus strand
                // the intron lies below the earlier exon rather than above it.
                let (start, end) = match self.strand {
                    Strand::Plus => (earlier.end + 1, later.start - 1),
                    Strand::Minus => (later.end + 1, earlier.start - 1),
                };
                Some(Intron {
                    start,
                    end,
                    junction_offset,
                    exon_before_end: match self.strand {
                        Strand::Plus => earlier.end,
                        Strand::Minus => earlier.start,
                    },
                    exon_after_start: match self.strand {
                        Strand::Plus => later.start,
                        Strand::Minus => later.end,
                    },
                })
            })
    }

    /// Whether an insertion anchored at genomic `position` lands at an internal
    /// exon-exon junction of the spliced CDS, i.e. immediately after the last
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
    /// Reads entitled to testify about linkage: those that observe every
    /// position of the row and carry the least-supported constituent SNV,
    /// whether or not they also carry the rest of the haplotype. This is the
    /// denominator of `MNV Phasing Support`; `None` (not zero) when no read
    /// spans the whole codon, which is a lack of evidence rather than evidence
    /// of absence.
    #[serde(default)]
    pub mnv_phasing_reads: Option<usize>,
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

#[cfg(test)]
mod intron_tests {
    use crate::test_support::{joined_gene, spliced_gene, transcript_gene};
    use crate::variants::Strand;

    #[test]
    fn a_spliced_transcript_yields_its_introns_in_transcript_order() {
        // exon1 801-900 (100 nt), intron 901-1000, exon2 1001-1200.
        let gene = spliced_gene("t", Strand::Plus, &[(801, 900), (1001, 1200)]);
        let introns: Vec<_> = gene.introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!((introns[0].start, introns[0].end), (901, 1000));
        // The exon after the intron begins at spliced-CDS offset 100.
        assert_eq!(introns[0].junction_offset, 100);
        assert_eq!(introns[0].exon_before_end, 900);
        assert_eq!(introns[0].exon_after_start, 1001);
    }

    #[test]
    fn minus_strand_introns_keep_transcript_orientation() {
        // Transcript order is descending genomic: exon 1001-1200 then 801-900,
        // so the donor sits at the low end of the intron and the acceptor high.
        let gene = spliced_gene("t", Strand::Minus, &[(1001, 1200), (801, 900)]);
        let introns: Vec<_> = gene.introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!((introns[0].start, introns[0].end), (901, 1000));
        assert_eq!(introns[0].junction_offset, 200);
        assert_eq!(introns[0].exon_before_end, 1001);
        assert_eq!(introns[0].exon_after_start, 900);
    }

    #[test]
    fn an_unspliced_transcript_has_no_introns() {
        let single = transcript_gene("t", Strand::Plus, &[(801, 900)]);
        assert_eq!(single.introns().count(), 0);
        let abutting = joined_gene("t", Strand::Plus, &[(801, 900), (901, 1200)]);
        assert_eq!(abutting.introns().count(), 0);
        // SARS-CoV-2 ORF1ab: join(266..13468, 13468..21555).
        let slippage = joined_gene("orf1ab", Strand::Plus, &[(266, 13468), (13468, 21555)]);
        assert_eq!(slippage.introns().count(), 0);
    }

    #[test]
    fn junction_offsets_accumulate_across_several_exons() {
        // 30 nt, intron, 60 nt, intron, 90 nt: junctions at offsets 30 and 90.
        let gene = spliced_gene("t", Strand::Plus, &[(1, 30), (101, 160), (301, 390)]);
        let offsets: Vec<_> = gene.introns().map(|i| i.junction_offset).collect();
        assert_eq!(offsets, vec![30, 90]);
    }

    #[test]
    fn a_slippage_join_after_a_real_intron_does_not_move_the_junction() {
        // exon 1-120, intron, exon 201-260, then a joined segment 261-320. The
        // only junction is the one after the real intron, at offset 120.
        let gene = transcript_gene("t", Strand::Plus, &[(1, 120), (201, 260), (261, 320)]);
        let introns: Vec<_> = gene.introns().collect();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0].junction_offset, 120);
        assert_eq!((introns[0].start, introns[0].end), (121, 200));
    }
}
