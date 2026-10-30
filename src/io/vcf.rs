//! VCF loading for plain and BGZF-compressed VCF files and for BCF, metrics
//! extraction, allele normalisation, and original INFO field preservation.
//!
//! BCF is decoded to VCF text as it is read, so every variant file reaches the
//! same line parser and a record means the same thing whatever it is stored as.
//!
//! Plain `.vcf` files use the fast text parser in `vcf_fast.rs`.

use super::validation::validate_vcf_allele;
use crate::error::{AppError, AppResult};
use crate::variants::{
    decompose_allele, substitution_components, AlleleComponent, AlleleComponentKind, AlleleEvent,
};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Read, Write};

/// INFO keys get_MNV writes itself. An input's values for them are dropped, so
/// a run over get_MNV's own output does not carry the old ones forward.
pub(crate) const GET_MNV_INFO_TAGS: &[&str] = &[
    "GENE", "AA", "CT", "TYPE", "ODP", "OFREQ", "SR", "SRF", "SRR", "MR", "MRF", "MRR", "DP",
    "FREQ", "SBP", "MSBP", "EC", "COMP", "ER", "ERF", "ERR", "EDP", "EFREQ", "SO", "IMPACT", "GD",
    "MNVSHIFT", "DBS", "MNVPS", "MNVPR", "FSPH", "DPHASE", "LD", "LDP", "NMD", "HGVSG", "HGVSC",
    "SBS96", "ADJSUB", "CAT", "ID83",
];

#[derive(Debug, Clone)]
pub struct VcfPosition {
    /// The record's POS: where the record begins, which is not always where it
    /// acts. A record may be padded with reference bases on its left, and a
    /// left-anchored indel names the base before the change. Anything asking
    /// where this variant *acts* wants [`VcfPosition::changed_positions`]; this
    /// field is for writing the record back out, sorting, and messages.
    pub record_start: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub original_dp: Option<usize>,
    pub original_freq: Option<f64>,
    pub original_info: Option<String>,
    /// Phase the input caller declared for this allele, from a `|`-separated
    /// `GT` and the `PS` phase set. `None` when the record is unphased, which
    /// is the usual case for haploid pathogen callers. get_MNV never phases on
    /// its own, so this is the caller's claim, not an observation.
    pub declared_phase: Option<DeclaredPhase>,
}

/// What the input VCF says about which haplotype carries this allele.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeclaredPhase {
    /// The `PS` value: only alleles sharing it were phased with each other.
    /// `None` when the genotype is phased but no phase set was given, in which
    /// case nothing here can be compared across records.
    pub phase_set: Option<usize>,
    /// Which haplotype slots of a `|`-separated `GT` carry this record's ALT.
    /// `1|0` gives `[0]`, `0|1` gives `[1]`, `1|1` gives `[0, 1]`.
    pub alt_haplotypes: Vec<u8>,
    /// Whether any slot is a no-call (`.`). Such a slot might hold this ALT, so
    /// two half-called records cannot be shown to be on different haplotypes
    /// even when their called slots are disjoint.
    pub has_unknown_slot: bool,
}

impl DeclaredPhase {
    /// Whether the caller placed both alleles on at least one shared
    /// haplotype. Comparable only within a phase set: two records phased in
    /// different sets say nothing about each other, and neither does an
    /// unphased one.
    pub fn shares_haplotype_with(&self, other: &DeclaredPhase) -> Option<bool> {
        if self.phase_set != other.phase_set {
            return None;
        }
        // Phased with no phase set names no group, so two such records make no
        // claim about each other.
        self.phase_set?;
        if self
            .alt_haplotypes
            .iter()
            .any(|slot| other.alt_haplotypes.contains(slot))
        {
            return Some(true);
        }
        // Disjoint called slots settle it as trans only when both genotypes
        // were fully called. A `.` could be holding this ALT.
        if self.has_unknown_slot || other.has_unknown_slot {
            return None;
        }
        Some(false)
    }
}

/// Read `GT` and `PS` for one sample into a [`DeclaredPhase`].
///
/// Only a `|`-separated genotype counts: `/` means the caller did not phase,
/// and treating it as phase would invent linkage the input never claimed. The
/// ALT of interest is allele index `alt_index + 1`, since index 0 is REF.
pub fn parse_declared_phase(
    genotype: Option<&str>,
    phase_set: Option<&str>,
    alt_index: usize,
) -> Option<DeclaredPhase> {
    let genotype = genotype?;
    if !genotype.contains('|') {
        return None;
    }
    let wanted = (alt_index + 1).to_string();
    let alt_haplotypes = genotype
        .split('|')
        .enumerate()
        .filter_map(|(slot, allele)| (allele == wanted).then_some(slot as u8))
        .collect::<Vec<_>>();
    if alt_haplotypes.is_empty() {
        return None;
    }
    Some(DeclaredPhase {
        phase_set: phase_set.and_then(|value| value.parse::<usize>().ok()),
        alt_haplotypes,
        has_unknown_slot: genotype.split('|').any(|allele| allele == "."),
    })
}

/// Whether the selected sample's genotype carries ALT number `alt_index`.
///
/// A VCF record lists every ALT seen at that site across the whole cohort, so a
/// multi-sample file names alleles a given sample does not have. Annotating all
/// of them for every sample made each sample carry every variant: the cohort
/// matrix of `--sample all` showed a variant for a sample whose genotype reads
/// `0/0`. A genotype that is absent, empty or partly a no-call says nothing
/// about carriage and keeps the allele, since unknown is not absence, and so
/// does one this parser cannot read.
pub fn sample_carries_alt(genotype: Option<&str>, alt_index: usize) -> bool {
    let Some(genotype) = genotype else {
        return true;
    };
    let genotype = genotype.trim();
    if genotype.is_empty() || genotype == "." {
        return true;
    }
    let wanted = (alt_index + 1).to_string();
    let mut fully_called = false;
    for allele in genotype.split(['/', '|']) {
        if allele == wanted {
            return true;
        }
        if allele == "." || allele.parse::<usize>().is_err() {
            return true;
        }
        fully_called = true;
    }
    !fully_called
}

impl VcfPosition {
    pub fn event(&self) -> AlleleEvent {
        decompose_allele(self.record_start, &self.ref_allele, &self.alt_allele)
    }

    /// Every genomic base this record's own components change.
    ///
    /// A record may start on a base it does not touch: `40 GCTGCTG>GCTGCTA`
    /// changes base 46 and nothing else, and a left-anchored deletion
    /// `90 AGT>A` removes 91 and 92 while starting at 90. Anything that asks
    /// where a record acts, rather than where it begins, asks this.
    pub fn changed_positions(&self) -> Vec<usize> {
        use crate::variants::AlleleComponentKind;

        let mut positions = Vec::new();
        for component in self.event().components {
            match component.kind {
                AlleleComponentKind::Snp | AlleleComponentKind::Insertion => {
                    positions.push(component.position);
                }
                AlleleComponentKind::Deletion
                | AlleleComponentKind::Delins
                | AlleleComponentKind::Symbolic => {
                    let end = component.position + component.ref_allele.len().saturating_sub(1);
                    positions.extend(component.position..=end);
                }
            }
        }
        if positions.is_empty() {
            positions.push(self.record_start);
        }
        positions.sort_unstable();
        positions.dedup();
        positions
    }

    pub fn substitution_components(&self) -> Vec<AlleleComponent> {
        substitution_components(self.record_start, &self.ref_allele, &self.alt_allele)
    }

    pub fn is_single_nucleotide_substitution(&self) -> bool {
        let components = self.substitution_components();
        components.len() == 1
            && self.ref_allele.chars().count() == 1
            && self.alt_allele.chars().count() == 1
    }

    pub fn overlaps_interval(&self, start: usize, end: usize) -> bool {
        if start == 0 || end == 0 || start > end {
            return false;
        }
        let event = self.event();
        event
            .components
            .iter()
            .any(|component| match component.kind {
                AlleleComponentKind::Snp => {
                    component.position >= start && component.position <= end
                }
                AlleleComponentKind::Insertion => {
                    component.position >= start && component.position < end
                }
                AlleleComponentKind::Deletion
                | AlleleComponentKind::Delins
                | AlleleComponentKind::Symbolic => {
                    let component_end = component
                        .position
                        .saturating_add(component.ref_allele.len().saturating_sub(1));
                    component.position <= end && component_end >= start
                }
            })
    }
}

pub(crate) fn parse_optional_freq_token(raw_token: &str) -> Option<f64> {
    let token = raw_token.trim();
    if token.is_empty() || token == "." {
        return None;
    }
    let has_percent = token.ends_with('%');
    let numeric = token.trim_end_matches('%').parse::<f64>().ok()?;
    if has_percent || (numeric > 1.0 && numeric <= 100.0) {
        Some(numeric / 100.0)
    } else {
        Some(numeric)
    }
}

/// The counts of an `AD` or `AO` field, each in its own place. A missing value
/// stays a gap instead of closing up: dropped before indexing, `AD=.,5,3` gave
/// the first ALT the second ALT's 3 reads.
fn allele_counts(raw: &str) -> Vec<Option<i64>> {
    raw.split(',')
        .map(|value| value.trim().parse::<i64>().ok().filter(|count| *count >= 0))
        .collect()
}

/// ALT number `alt_index`'s share of the reads in an `AD` field, which lists
/// REF first. Unknown when this ALT's count or REF's is missing. Another ALT's
/// missing count adds nothing, since that is what `bcftools merge` writes for
/// an allele this sample was not called with.
pub(crate) fn freq_from_ad(raw: &str, alt_index: usize) -> Option<f64> {
    let counts = allele_counts(raw);
    let alt = (*counts.get(alt_index + 1)?)?;
    let total = counts[0]? + counts[1..].iter().flatten().sum::<i64>();
    (total > 0).then(|| alt as f64 / total as f64)
}

/// ALT number `alt_index`'s share from FreeBayes' `AO`, one count per ALT, over
/// the depth when it is known and over `AO + RO` otherwise.
pub(crate) fn freq_from_ao_ro(
    ao_raw: &str,
    ro_raw: Option<&str>,
    original_dp: Option<usize>,
    alt_index: usize,
) -> Option<f64> {
    let counts = allele_counts(ao_raw);
    let alt = (*counts.get(alt_index)?)?;
    let total = match original_dp {
        Some(dp) => dp as i64,
        None => {
            let ro = ro_raw
                .and_then(|raw| raw.trim().parse::<i64>().ok())
                .unwrap_or(0);
            counts.iter().flatten().sum::<i64>() + ro
        }
    };
    (total > 0).then(|| alt as f64 / total as f64)
}

pub(crate) fn parse_optional_freq_index(raw: &str, alt_index: usize) -> Option<f64> {
    let values = raw.split(',').map(str::trim).collect::<Vec<_>>();
    if values.is_empty() {
        return None;
    }
    let token = if alt_index < values.len() {
        values[alt_index]
    } else if values.len() == 1 {
        values[0]
    } else {
        return None;
    };
    parse_optional_freq_token(token)
}

pub(crate) fn parse_optional_depth(raw: &str) -> Option<usize> {
    let raw_value = raw.trim();
    let first = raw_value.split(',').next()?;
    if first.is_empty() || first == "." {
        return None;
    }
    if let Ok(value) = first.parse::<usize>() {
        return Some(value);
    }
    first
        .parse::<f64>()
        .ok()
        .filter(|value| *value >= 0.0)
        .map(|value| value.round() as usize)
}

/// The longest structural event [`spell_out_symbolic`] writes out. Past this a
/// record stays symbolic: its alleles would run to tens of kilobases in every
/// output row, and an event that size is a matter for an SV annotator.
pub(crate) const MAX_SPELLED_SYMBOLIC_SPAN: usize = 10_000;

/// A symbolic deletion or duplication written out as the explicit alleles the
/// reference gives it, so it is annotated like any other indel: frameshift or
/// in frame, its protein change, HGVS, and the frame shift it passes to the
/// codons after it. Without this the tool knew it was a deletion but not of how
/// many bases, and reported the row and everything downstream of it as unknown.
///
/// `<DEL>` (and its `<DEL:...>` subtypes) removes POS+1 to END. `<DUP>` and
/// `<DUP:TANDEM>` are read as a tandem copy of POS+1 to END, written as an
/// insertion of that segment right after END. The extent comes from `END`, or
/// from `SVLEN` (this ALT's entry, sign ignored) when `END` is absent.
///
/// `None` leaves the record symbolic: another symbolic type, no extent, an
/// extent past the contig or over [`MAX_SPELLED_SYMBOLIC_SPAN`], no reference
/// for the contig, or a padding base the reference does not have.
pub(crate) fn spell_out_symbolic(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    info: &str,
    alt_index: usize,
    sequence: Option<&str>,
) -> Option<(usize, String, String)> {
    let sequence = sequence?;
    let kind = alt_allele.strip_prefix('<')?.strip_suffix('>')?;
    let deletion = kind == "DEL" || kind.starts_with("DEL:");
    let duplication = kind == "DUP" || kind == "DUP:TANDEM";
    if !(deletion || duplication) || ref_allele.len() != 1 {
        return None;
    }
    let end = match get_info_value(info, "END").and_then(|value| value.parse::<usize>().ok()) {
        Some(end) => end,
        None => {
            let svlen = get_info_value(info, "SVLEN")?
                .split(',')
                .nth(alt_index)?
                .trim()
                .parse::<i64>()
                .ok()?;
            pos.checked_add(usize::try_from(svlen.unsigned_abs()).ok()?)?
        }
    };
    if end <= pos || end > sequence.len() || end - pos > MAX_SPELLED_SYMBOLIC_SPAN {
        return None;
    }
    let padding = sequence.get(pos - 1..pos)?;
    if !padding.eq_ignore_ascii_case(ref_allele) {
        return None;
    }
    if deletion {
        Some((
            pos,
            sequence[pos - 1..end].to_ascii_uppercase(),
            padding.to_ascii_uppercase(),
        ))
    } else {
        let anchor = sequence[end - 1..end].to_ascii_uppercase();
        let copy = sequence[pos..end].to_ascii_uppercase();
        Some((end, anchor.clone(), format!("{anchor}{copy}")))
    }
}

/// How the readers rewrite each record's alleles before anything annotates
/// them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AlleleNormalization {
    /// Kept as written.
    #[default]
    AsWritten,
    /// `--normalize-alleles`: the context REF and ALT share is trimmed wherever
    /// that leaves the event where it was ([`normalize_ref_alt`]).
    Trim,
    /// `--left-align`: trimmed and moved to the leftmost place that describes
    /// the same change, as `bcftools norm -f` writes it. A record left-alignment
    /// leaves as written is still trimmed, so this includes [`Self::Trim`].
    LeftAlign,
}

impl AlleleNormalization {
    /// The record as this setting writes it, and whether that moved its event
    /// rather than only respelling it, which is what [`log_left_aligned`]
    /// counts. `sequence` is the record's contig, which left-alignment walks;
    /// without it a record is only trimmed.
    pub(crate) fn apply(
        self,
        pos: usize,
        ref_allele: String,
        alt_allele: String,
        sequence: Option<&str>,
    ) -> (usize, String, String, bool) {
        let trimmed = |ref_allele: &str, alt_allele: &str| {
            let (pos, ref_allele, alt_allele) = normalize_ref_alt(pos, ref_allele, alt_allele);
            (pos, ref_allele, alt_allele, false)
        };
        match self {
            Self::AsWritten => (pos, ref_allele, alt_allele, false),
            Self::Trim => trimmed(&ref_allele, &alt_allele),
            Self::LeftAlign => {
                let aligned = sequence.and_then(|seq| {
                    super::left_align::left_align(pos, &ref_allele, &alt_allele, seq)
                });
                match aligned {
                    Some((new_pos, new_ref, new_alt)) => {
                        let moved =
                            crate::variants::decompose_allele(pos, &ref_allele, &alt_allele)
                                .components
                                != crate::variants::decompose_allele(new_pos, &new_ref, &new_alt)
                                    .components;
                        (new_pos, new_ref, new_alt, moved)
                    }
                    None => trimmed(&ref_allele, &alt_allele),
                }
            }
        }
    }
}

pub(crate) fn log_left_aligned(moved: usize) {
    if moved > 0 {
        log::info!(
            "Left-aligned {moved} indel(s) against the reference: each now sits at the \
             leftmost position that describes the same change, where aligners place it"
        );
    }
}

pub(crate) fn normalize_ref_alt(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> (usize, String, String) {
    let is_symbolic = alt_allele.starts_with('<') && alt_allele.ends_with('>');
    if is_symbolic {
        return (pos, ref_allele.to_string(), alt_allele.to_string());
    }

    let untouched = || (pos, ref_allele.to_string(), alt_allele.to_string());
    let before = crate::variants::decompose_allele(pos, ref_allele, alt_allele);

    // Trimming shared context must not move the event. Inside a repeat it can:
    // `29 CT>CTGCT` trims canonically to `29 C>CTGC`, the same insertion
    // anchored one base earlier. get_MNV matches indel support against the exact
    // CIGAR the aligner used, so the relocated allele found none of its own
    // reads and a row with 20 of 20 supporting reads was reported at 0, under a
    // warning that blamed the input for a move this flag had made. Placing an
    // indel canonically would need reference-aware left-alignment of both sides,
    // which this flag does not promise and cannot do from the alleles alone.
    //
    // So: take the most trimming that leaves the event where it was. The
    // canonical suffix-then-prefix trim first, since that is the form other
    // tools expect; the prefix alone when the suffix step is what moved it; and
    // the allele untouched when neither is safe.
    for trim_suffix in [true, false] {
        if let Some(result) = trimmed(pos, ref_allele, alt_allele, trim_suffix) {
            let after = crate::variants::decompose_allele(result.0, &result.1, &result.2);
            if after.components == before.components {
                return result;
            }
        }
    }
    untouched()
}

/// One trimming attempt: the shared suffix (optionally) then the shared prefix,
/// always leaving at least one base in each allele. `None` when there is nothing
/// left to describe.
fn trimmed(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    trim_suffix: bool,
) -> Option<(usize, String, String)> {
    let ref_chars: Vec<char> = ref_allele.chars().collect();
    let alt_chars: Vec<char> = alt_allele.chars().collect();
    let mut start = 0usize;
    let mut ref_end = ref_chars.len();
    let mut alt_end = alt_chars.len();

    if trim_suffix {
        while ref_end - start > 1
            && alt_end - start > 1
            && ref_chars[ref_end - 1] == alt_chars[alt_end - 1]
        {
            ref_end -= 1;
            alt_end -= 1;
        }
    }
    while ref_end - start > 1 && alt_end - start > 1 && ref_chars[start] == alt_chars[start] {
        start += 1;
    }

    let norm_ref: String = ref_chars[start..ref_end].iter().collect();
    let norm_alt: String = alt_chars[start..alt_end].iter().collect();
    if norm_ref.is_empty() || norm_alt.is_empty() {
        return None;
    }
    Some((pos + start, norm_ref, norm_alt))
}

/// Parse a VCF line into fields. Returns None if the line is a header or empty.
/// Whether a line is the `#CHROM` line that names the columns, in whichever
/// case it is written. Every reader asks this one question, so a lowercase
/// `#chrom` cannot be a header to one reader and nothing to another.
pub(crate) fn is_column_header(line: &str) -> bool {
    line.get(..6)
        .is_some_and(|start| start.eq_ignore_ascii_case("#CHROM"))
}

fn parse_vcf_line(line: &str) -> Option<Vec<&str>> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    Some(line.split('\t').collect())
}

/// Extract INFO field value by key from a raw INFO string.
fn get_info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    for field in info.split(';') {
        if let Some((k, v)) = field.split_once('=') {
            if k == key {
                return Some(v);
            }
        } else if field == key {
            return Some("");
        }
    }
    None
}

/// Extract FORMAT field value for a given sample index.
fn get_format_value<'a>(format_keys: &[&str], sample_field: &'a str, key: &str) -> Option<&'a str> {
    let idx = format_keys.iter().position(|k| *k == key)?;
    let values: Vec<&str> = sample_field.split(':').collect();
    values
        .get(idx)
        .copied()
        .filter(|v| !v.is_empty() && *v != ".")
}

fn parse_original_metrics_from_fields(
    info: &str,
    format_keys: &[&str],
    sample_field: Option<&str>,
    alt_index: usize,
) -> (Option<usize>, Option<f64>) {
    let mut original_dp: Option<usize> = None;
    let mut original_freq: Option<f64> = None;

    // Try FORMAT fields first (sample-specific)
    if let Some(sample) = sample_field {
        // DP
        if let Some(dp_str) = get_format_value(format_keys, sample, "DP") {
            original_dp = parse_optional_depth(dp_str);
        }
        // FREQ / AF
        for tag in &["FREQ", "AF"] {
            if original_freq.is_none() {
                if let Some(val) = get_format_value(format_keys, sample, tag) {
                    original_freq = parse_optional_freq_index(val, alt_index);
                }
            }
        }
        // AD → derive freq
        if original_freq.is_none() {
            if let Some(ad_str) = get_format_value(format_keys, sample, "AD") {
                original_freq = freq_from_ad(ad_str, alt_index);
            }
        }
        // AO/RO → derive freq
        if original_freq.is_none() {
            if let Some(ao_str) = get_format_value(format_keys, sample, "AO") {
                let ro = get_format_value(format_keys, sample, "RO");
                original_freq = freq_from_ao_ro(ao_str, ro, original_dp, alt_index);
            }
        }
    }

    // Fall back to INFO fields
    if original_dp.is_none() {
        if let Some(dp_str) = get_info_value(info, "DP") {
            original_dp = parse_optional_depth(dp_str);
        }
    }
    for tag in &["AF", "FREQ"] {
        if original_freq.is_none() {
            if let Some(val) = get_info_value(info, tag) {
                original_freq = parse_optional_freq_index(val, alt_index);
            }
        }
    }
    if original_freq.is_none() {
        if let Some(ad_str) = get_info_value(info, "AD") {
            original_freq = freq_from_ad(ad_str, alt_index);
        }
    }
    if original_freq.is_none() {
        if let Some(ao_str) = get_info_value(info, "AO") {
            let ro = get_info_value(info, "RO");
            original_freq = freq_from_ao_ro(ao_str, ro, original_dp, alt_index);
        }
    }

    (original_dp, original_freq)
}

/// Read a `key=...` attribute value out of a `##INFO` header line, up to the
/// next `,` or `>`.
fn header_attr(line: &str, key: &str) -> Option<String> {
    let start = line.find(key)? + key.len();
    let rest = &line[start..];
    let end = rest.find([',', '>']).unwrap_or(rest.len());
    Some(rest[..end].to_string())
}

/// Parse `Number=` from `##INFO` header lines into a map of field ID → its
/// per-allele cardinality code (`A`, `R` or `G`). Fixed-Number fields (1, 0,
/// ., a constant) are omitted because they need no per-allele subsetting.
pub(crate) fn per_allele_info_numbers(header_lines: &[String]) -> HashMap<String, char> {
    let mut map = HashMap::new();
    for line in header_lines {
        if let (Some(id), Some(number)) = (header_attr(line, "ID="), header_attr(line, "Number=")) {
            match number.as_str() {
                "A" => drop(map.insert(id, 'A')),
                "R" => drop(map.insert(id, 'R')),
                "G" => drop(map.insert(id, 'G')),
                _ => {}
            }
        }
    }
    map
}

/// Subset a per-allele INFO field value to a single ALT so it is valid on a
/// single-ALT output record. `Number=A` keeps the value at `alt_idx`;
/// `Number=R` keeps `[REF, ALT@alt_idx]`. `Number=G` and arrays too short for
/// `alt_idx` are dropped (cannot be safely attributed to one allele).
fn subset_per_allele_field(key: &str, value: &str, kind: char, alt_idx: usize) -> Option<String> {
    let elems: Vec<&str> = value.split(',').collect();
    match kind {
        'A' => elems.get(alt_idx).map(|v| format!("{key}={v}")),
        'R' => match (elems.first(), elems.get(alt_idx + 1)) {
            (Some(r), Some(a)) => Some(format!("{key}={r},{a}")),
            _ => None,
        },
        _ => None,
    }
}

/// Filter get_mnv's own INFO tags out of an original record's INFO string and
/// subset any per-allele field to `alt_idx`, so the result is valid when copied
/// onto a single-ALT output record (`--keep-original-info` + multiallelic
/// input previously copied whole `Number=A/R` arrays, producing invalid VCF).
pub(crate) fn filter_and_subset_original_info<F: Fn(&str) -> bool>(
    info: &str,
    is_own_tag: F,
    per_allele: &HashMap<String, char>,
    alt_idx: usize,
) -> Option<String> {
    if info.is_empty() || info == "." {
        return None;
    }
    let parts: Vec<String> = info
        .split(';')
        .filter_map(|field| {
            let (key, value) = match field.split_once('=') {
                Some((k, v)) => (k, Some(v)),
                None => (field, None),
            };
            if is_own_tag(key) {
                return None;
            }
            match (per_allele.get(key), value) {
                (Some(&kind), Some(v)) => subset_per_allele_field(key, v, kind, alt_idx),
                _ => Some(field.to_string()),
            }
        })
        .collect();
    if parts.is_empty() {
        None
    } else {
        Some(parts.join(";"))
    }
}

/// An input `ANN` kept for one allele of a record, pointed at the allele its
/// output record will carry.
///
/// snpEff writes each entry against the record's ALT as the record spells it.
/// get_MNV writes a padded SNV as the one base that changes, and a normalized
/// or spelled-out allele as such, so an entry left as it was named an allele no
/// output record has and was never shown. Entries for the record's other ALTs
/// go: each split record carries its own.
pub(crate) fn retarget_input_ann(
    info: Option<String>,
    record_alt: &str,
    position: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<String> {
    let info = info?;
    if !info.split(';').any(|field| field.starts_with("ANN=")) {
        return Some(info);
    }
    let event = crate::variants::decompose_allele(position, ref_allele, alt_allele);
    let output_alt = match event.components.as_slice() {
        [one] if one.kind == crate::variants::AlleleComponentKind::Snp => one.alt_allele.clone(),
        _ => alt_allele.to_string(),
    };
    let fields: Vec<String> = info
        .split(';')
        .filter_map(|field| {
            let Some(value) = field.strip_prefix("ANN=") else {
                return Some(field.to_string());
            };
            let entries: Vec<String> = value
                .split(',')
                .filter_map(|entry| {
                    let (allele, rest) = entry.split_once('|')?;
                    allele
                        .eq_ignore_ascii_case(record_alt)
                        .then(|| format!("{output_alt}|{rest}"))
                })
                .collect();
            (!entries.is_empty()).then(|| format!("ANN={}", entries.join(",")))
        })
        .collect();
    (!fields.is_empty()).then(|| fields.join(";"))
}

/// Whether a path names a BCF file, the binary form of VCF.
pub fn is_bcf_path(path: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".bcf")
}

/// Opens a variant file as VCF text, whatever it is stored as: plain, BGZF or
/// BCF. BCF is decoded on a background thread and handed over through a pipe,
/// so a large file is never held in memory as text and a caller that stops at
/// the header stops the decoder too.
pub(crate) fn open_vcf_text(vcf_file: &str) -> AppResult<Box<dyn BufRead>> {
    if is_bcf_path(vcf_file) {
        return Ok(Box::new(DecodedBcf::spawn(vcf_file)?));
    }
    let mut file = std::fs::File::open(vcf_file)?;
    if vcf_file.ends_with(".gz") {
        ensure_bgzf_eof(&mut file, vcf_file)?;
        Ok(Box::new(BufReader::new(noodles::bgzf::io::Reader::new(
            file,
        ))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// The empty block every BGZF writer ends a file with.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Checks that a BGZF file, a `.vcf.gz` or a BCF, ends with its end-of-file
/// block. The reader stops cleanly at the end of any block, so a file cut
/// between two blocks, by a copy or a download that did not finish, read as a
/// VCF with fewer records and the run succeeded on part of the input. A file
/// that does not start like BGZF is left to the reader, which says what it is.
fn ensure_bgzf_eof(file: &mut std::fs::File, path: &str) -> AppResult<()> {
    use std::io::{Seek, SeekFrom};

    // A pipe cannot be checked without reading it whole, so it is read as it
    // comes.
    if !file.metadata()?.is_file() {
        return Ok(());
    }
    let mut start = [0u8; 4];
    let is_bgzf = file.read_exact(&mut start).is_ok() && start == BGZF_EOF[..4];
    let mut end = [0u8; 28];
    let complete = !is_bgzf
        || (file.seek(SeekFrom::End(-28)).is_ok()
            && file.read_exact(&mut end).is_ok()
            && end == BGZF_EOF);
    file.seek(SeekFrom::Start(0))?;
    if complete {
        Ok(())
    } else {
        Err(AppError::validation(format!(
            "'{path}' does not end with the BGZF end-of-file block, so it was most likely cut \
             short, by a copy or a download that did not finish or by a writer that was \
             stopped. Reading it would annotate only the calls before the cut. Check the \
             file, or write it again with bgzip or bcftools."
        )))
    }
}

/// VCF text decoded from a BCF file by a background thread.
///
/// The end of the text is only the end of the file if the decoder finished
/// cleanly, so reaching it joins the thread and turns a decoding error into a
/// read error. Without that, a damaged BCF would read as a VCF that simply has
/// fewer records.
struct DecodedBcf {
    text: BufReader<std::io::PipeReader>,
    decoder: Option<std::thread::JoinHandle<std::io::Result<()>>>,
}

impl DecodedBcf {
    fn spawn(path: &str) -> AppResult<Self> {
        // Open here rather than in the thread, so a missing file fails with
        // its own message before anything is read.
        let mut file =
            std::fs::File::open(path).map_err(|e| format!("Cannot open BCF file '{path}': {e}"))?;
        ensure_bgzf_eof(&mut file, path)?;
        let (text, sink) = std::io::pipe()?;
        let decoder = std::thread::spawn(move || decode_bcf_as_vcf_text(file, sink));
        Ok(Self {
            text: BufReader::with_capacity(64 * 1024, text),
            decoder: Some(decoder),
        })
    }

    fn finish(&mut self) -> std::io::Result<()> {
        match self.decoder.take() {
            Some(handle) => handle.join().unwrap_or_else(|_| {
                Err(std::io::Error::other(
                    "the BCF decoder stopped unexpectedly",
                ))
            }),
            None => Ok(()),
        }
    }
}

impl Read for DecodedBcf {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let read = self.text.read(buf)?;
        if read == 0 && !buf.is_empty() {
            self.finish()?;
        }
        Ok(read)
    }
}

impl BufRead for DecodedBcf {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        if self.text.fill_buf()?.is_empty() {
            self.finish()?;
        }
        self.text.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        self.text.consume(amount);
    }
}

fn decode_bcf_as_vcf_text(file: std::fs::File, sink: std::io::PipeWriter) -> std::io::Result<()> {
    use noodles::vcf::variant::io::Write as _;

    let mut reader = noodles::bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|e| std::io::Error::new(e.kind(), format!("invalid BCF header: {e}")))?;
    let mut sink = std::io::BufWriter::new(sink);
    let mut writer = noodles::vcf::io::Writer::new(Vec::new());
    writer.write_header(&header)?;
    sink.write_all(writer.get_ref())?;
    writer.get_mut().clear();
    for (index, result) in reader.records().enumerate() {
        let record = result.map_err(|e| {
            std::io::Error::new(e.kind(), format!("invalid BCF record {}: {e}", index + 1))
        })?;
        writer.write_variant_record(&header, &record)?;
        // The VCF writer percent-encodes a `%` inside a string value as `%25`,
        // as VCF 4.3 asks. A text VCF from VarScan writes `FREQ=20%` with a bare
        // `%`, and that is what the parsers read: left encoded, every such
        // value came out unreadable from a BCF and fine from the same VCF. The
        // other encodings stand for characters a VCF line cannot hold bare.
        let line = writer.get_mut();
        sink.write_all(&decode_percent_sign(line))?;
        line.clear();
    }
    sink.flush()
}

/// `%25` back to `%`, and nothing else.
fn decode_percent_sign(line: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    if !line.windows(3).any(|w| w == b"%25") {
        return std::borrow::Cow::Borrowed(line);
    }
    let mut out = Vec::with_capacity(line.len());
    let mut i = 0;
    while i < line.len() {
        if line[i..].starts_with(b"%25") {
            out.push(b'%');
            i += 3;
        } else {
            out.push(line[i]);
            i += 1;
        }
    }
    std::borrow::Cow::Owned(out)
}

pub fn list_vcf_samples(vcf_file: &str) -> AppResult<Vec<String>> {
    let reader = open_vcf_text(vcf_file)?;

    // Parse header to find #CHROM line
    for line_result in reader.lines() {
        let line = line_result?;
        if is_column_header(&line) {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() > 9 {
                return Ok(fields[9..].iter().map(|s| s.to_string()).collect());
            }
            return Ok(Vec::new());
        }
        if !line.starts_with('#') {
            break;
        }
    }
    Ok(Vec::new())
}

fn resolve_sample_index(samples: &[String], sample_name: Option<&str>) -> AppResult<Option<usize>> {
    if let Some(name) = sample_name {
        if samples.is_empty() {
            return Err(format!("Requested sample '{name}' but VCF has no sample columns").into());
        }
        if let Some(index) = samples.iter().position(|value| value == name) {
            return Ok(Some(index));
        }
        return Err(format!(
            "Sample '{}' not found in VCF header. Available samples: {}",
            name,
            samples.join(", ")
        )
        .into());
    }

    if samples.is_empty() {
        return Ok(None);
    }
    if samples.len() > 1 {
        log::warn!(
            "Multi-sample VCF detected ({} samples). Using first sample '{}' for original FORMAT metrics. Use --sample to select another sample.",
            samples.len(),
            samples[0]
        );
    }
    Ok(Some(0))
}

/// Whether an input `##INFO` header gives way to the one get_MNV writes: one of
/// its own keys, or `ANN`, whose input entries are kept after get_MNV's own.
pub(crate) fn header_written_by_get_mnv(id: &str) -> bool {
    GET_MNV_INFO_TAGS.contains(&id) || id == "ANN"
}

pub fn extract_original_info_headers(vcf_file: &str) -> AppResult<Vec<String>> {
    let mut lines = Vec::new();

    let reader = open_vcf_text(vcf_file)?;

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("##INFO=") {
            // Extract ID from the header line
            if let Some(id_start) = line.find("ID=") {
                let rest = &line[id_start + 3..];
                if let Some(id_end) = rest.find([',', '>']) {
                    let id = &rest[..id_end];
                    if !header_written_by_get_mnv(id) {
                        lines.push(line.clone());
                    }
                }
            }
        }
        if !line.starts_with('#') {
            break;
        }
    }
    Ok(lines)
}

/// Loads a variant file with the reader its container takes: the plain-text
/// one for `.vcf`, the BGZF-aware one for `.vcf.gz` and BCF. Everything that
/// loads calls goes through here, so nothing reads a file another way than a
/// run does.
pub fn load_vcf(
    vcf_file: &str,
    sample_name: Option<&str>,
    split_multiallelic: bool,
    normalization: AlleleNormalization,
    keep_original_info: bool,
    references: Option<&crate::io::ReferenceMap>,
) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    if super::vcf_fast::use_fast_parser(vcf_file) {
        super::vcf_fast::load_vcf_text(
            vcf_file,
            sample_name,
            split_multiallelic,
            normalization,
            keep_original_info,
            references,
        )
    } else {
        load_vcf_positions_by_contig(
            vcf_file,
            sample_name,
            split_multiallelic,
            normalization,
            keep_original_info,
            references,
        )
    }
}

pub fn load_vcf_positions_by_contig(
    vcf_file: &str,
    sample_name: Option<&str>,
    split_multiallelic: bool,
    normalization: AlleleNormalization,
    keep_original_info: bool,
    references: Option<&crate::io::ReferenceMap>,
) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    log::info!("Loading VCF: {vcf_file}");

    let samples = list_vcf_samples(vcf_file)?;
    let sample_index = resolve_sample_index(&samples, sample_name)?;
    let own_tags: HashSet<&str> = GET_MNV_INFO_TAGS.iter().copied().collect();
    let per_allele_info = if keep_original_info {
        per_allele_info_numbers(&extract_original_info_headers(vcf_file)?)
    } else {
        HashMap::new()
    };

    let reader = open_vcf_text(vcf_file)?;

    let mut positions_by_contig: HashMap<String, Vec<VcfPosition>> = HashMap::new();
    let mut split_count = 0usize;
    let mut not_carried = 0usize;
    let mut overlapped = 0usize;
    let mut spelled_out = 0usize;
    let mut left_aligned = 0usize;
    let mut record_idx = 0usize;
    let mut header_seen = false;

    for line_result in reader.lines() {
        let line = line_result?;
        if is_column_header(&line) {
            header_seen = true;
        }
        let fields = match parse_vcf_line(&line) {
            Some(f) => f,
            None => continue,
        };
        // Skip any data line appearing before the #CHROM header, matching the
        // plain-text fast parser. A well-formed VCF has no records before it.
        if !header_seen {
            continue;
        }
        record_idx += 1;

        if fields.len() < 8 {
            return Err(format!("VCF record {} has fewer than 8 fields", record_idx).into());
        }

        let chrom = fields[0];
        let pos: usize = fields[1].parse().map_err(|_| {
            format!(
                "Invalid VCF position at record {}: {}",
                record_idx, fields[1]
            )
        })?;
        if pos == 0 {
            return Err(format!(
                "Invalid VCF position 0 at record {record_idx} (VCF coordinates are 1-based)"
            )
            .into());
        }
        let ref_allele = fields[3];
        let alt_field = fields[4];
        let info = fields[7];

        if ref_allele.is_empty() {
            return Err(format!(
                "Invalid VCF allele at record {} (pos {}): empty REF",
                record_idx, pos
            )
            .into());
        }
        validate_vcf_allele(ref_allele, record_idx, chrom, pos, "REF")?;

        let alts: Vec<&str> = alt_field.split(',').collect();
        if alts.len() > 1 && !split_multiallelic {
            return Err(format!(
                "Multiallelic VCF record {} at {}:{} is not supported. Split multiallelic sites first (e.g. bcftools norm -m -).",
                record_idx, chrom, pos
            ).into());
        }

        // Get FORMAT and sample fields
        let format_keys: Vec<&str> = if fields.len() > 8 {
            fields[8].split(':').collect()
        } else {
            Vec::new()
        };
        let sample_field = sample_index.and_then(|idx| fields.get(9 + idx).copied());
        let genotype = sample_field.and_then(|sample| get_format_value(&format_keys, sample, "GT"));
        let phase_set =
            sample_field.and_then(|sample| get_format_value(&format_keys, sample, "PS"));

        // Register the contig as soon as a record for it is read, before any
        // allele is kept or skipped. An input whose every allele is skipped, a
        // cohort sample that carries none of them, then still reports the contig
        // it covered, which is what separates "this run has nothing to annotate"
        // from "this file holds no records at all". Only the second is an error.
        positions_by_contig.entry(chrom.to_string()).or_default();

        let mut pushed_here = 0usize;
        for (alt_idx, alt_allele) in alts.iter().enumerate() {
            if alt_allele.is_empty() || *alt_allele == "." {
                return Err(format!(
                    "Invalid VCF allele at record {} (pos {}): REF='{}' ALT='{}'",
                    record_idx, pos, ref_allele, alt_allele
                )
                .into());
            }
            if !sample_carries_alt(genotype, alt_idx) {
                not_carried += 1;
                continue;
            }
            if *alt_allele == "*" {
                overlapped += 1;
                continue;
            }
            let sequence = references
                .and_then(|refs| refs.get(chrom))
                .map(String::as_str);
            let (start, explicit_ref, explicit_alt) =
                match spell_out_symbolic(pos, ref_allele, alt_allele, info, alt_idx, sequence) {
                    Some(spelled) => {
                        spelled_out += 1;
                        spelled
                    }
                    None => (pos, ref_allele.to_string(), alt_allele.to_string()),
                };
            let (normalized_pos, normalized_ref, normalized_alt, moved) =
                normalization.apply(start, explicit_ref, explicit_alt, sequence);
            left_aligned += usize::from(moved);
            validate_vcf_allele(&normalized_ref, record_idx, chrom, normalized_pos, "REF")?;
            validate_vcf_allele(&normalized_alt, record_idx, chrom, normalized_pos, "ALT")?;

            let (original_dp, original_freq) =
                parse_original_metrics_from_fields(info, &format_keys, sample_field, alt_idx);
            let original_info = if keep_original_info {
                retarget_input_ann(
                    filter_and_subset_original_info(
                        info,
                        |k| own_tags.contains(k),
                        &per_allele_info,
                        alt_idx,
                    ),
                    alt_allele,
                    normalized_pos,
                    &normalized_ref,
                    &normalized_alt,
                )
            } else {
                None
            };

            positions_by_contig
                .entry(chrom.to_string())
                .or_default()
                .push(VcfPosition {
                    record_start: normalized_pos,
                    ref_allele: normalized_ref,
                    alt_allele: normalized_alt,
                    original_dp,
                    original_freq,
                    original_info,
                    declared_phase: parse_declared_phase(genotype, phase_set, alt_idx),
                });
            pushed_here += 1;
        }
        split_count += pushed_here.saturating_sub(1);
    }

    if !header_seen {
        return Err("No #CHROM header line found in VCF".into());
    }

    for values in positions_by_contig.values_mut() {
        values.sort_by_key(|v| v.record_start);
    }

    if split_multiallelic && split_count > 0 {
        log::info!("Split {split_count} additional ALT alleles from multiallelic VCF records");
    }
    if not_carried > 0 {
        log::info!(
            "Skipped {not_carried} ALT alleles the selected sample's genotype does not carry"
        );
    }
    log_overlapped(overlapped);
    log_spelled_out(spelled_out);
    log_left_aligned(left_aligned);

    Ok(positions_by_contig)
}

/// A `*` ALT is not a change at its position. It says this allele is missing
/// there because a deletion that starts upstream covers it, and that deletion
/// is annotated at its own record. GATK writes one beside a real ALT whenever
/// one sample's deletion overlaps another sample's SNV; with splitting the
/// default, such a record stopped a whole cohort at the first sample that
/// carried the `*`.
pub(crate) fn log_overlapped(overlapped: usize) {
    if overlapped > 0 {
        log::info!(
            "Skipped {overlapped} '*' allele(s): each stands for a deletion that starts \
             upstream, which is annotated at its own record"
        );
    }
}

pub(crate) fn log_spelled_out(spelled_out: usize) {
    if spelled_out > 0 {
        log::info!(
            "Wrote out {spelled_out} symbolic deletion/duplication allele(s) from the reference, \
             so they are annotated as explicit indels"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::{spell_out_symbolic, VcfPosition, MAX_SPELLED_SYMBOLIC_SPAN};

    fn variant(position: usize, ref_allele: &str, alt_allele: &str) -> VcfPosition {
        VcfPosition {
            record_start: position,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
            declared_phase: None,
        }
    }

    #[test]
    fn test_sample_carries_alt_reads_the_genotype() {
        use super::sample_carries_alt;
        // Sin genotipo no se sabe nada, asi que se conserva el alelo.
        assert!(sample_carries_alt(None, 0));
        assert!(sample_carries_alt(Some("./."), 0));
        assert!(sample_carries_alt(Some("."), 0));
        assert!(sample_carries_alt(Some(""), 0));
        // Una llamada completa sin ALT dice que la muestra no lo lleva.
        assert!(!sample_carries_alt(Some("0/0"), 0));
        assert!(!sample_carries_alt(Some("0|0"), 0));
        assert!(!sample_carries_alt(Some("0"), 0));
        // Y una que si lo lleva.
        assert!(sample_carries_alt(Some("0/1"), 0));
        assert!(sample_carries_alt(Some("1|1"), 0));
        assert!(sample_carries_alt(Some("1"), 0));
        // Multialelico: 1/2 lleva los dos, 1/1 solo el primero.
        assert!(sample_carries_alt(Some("1/2"), 0));
        assert!(sample_carries_alt(Some("1/2"), 1));
        assert!(sample_carries_alt(Some("1/1"), 0));
        assert!(!sample_carries_alt(Some("1/1"), 1));
        // Un genotipo a medias no descarta nada.
        assert!(sample_carries_alt(Some("./1"), 1));
        assert!(sample_carries_alt(Some("0/."), 0));
        assert!(sample_carries_alt(Some("raro"), 0));
    }

    #[test]
    fn test_list_vcf_samples_accepts_lowercase_chrom_header() {
        use std::io::Write;
        let path =
            std::env::temp_dir().join(format!("get_mnv_lc_samples_{}.vcf", std::process::id()));
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "#chrom\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A"
        )
        .unwrap();
        drop(f);
        let samples = super::list_vcf_samples(path.to_str().unwrap()).unwrap();
        assert_eq!(samples, vec!["SAMPLE_A".to_string()]);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_load_vcf_accepts_lowercase_chrom_and_skips_preheader_data() {
        use std::io::Write;
        let path = std::env::temp_dir().join(format!("get_mnv_lc_load_{}.vcf", std::process::id()));
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        // A stray data-shaped line before the header must be ignored, and a
        // lowercase #chrom header must be accepted (parity with the fast parser).
        writeln!(f, "chr1\t1\t.\tA\tG\t.\tPASS\t.").unwrap();
        writeln!(f, "#chrom\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(f, "chr1\t5\t.\tA\tT\t.\tPASS\t.").unwrap();
        drop(f);
        let by_contig = super::load_vcf_positions_by_contig(
            path.to_str().unwrap(),
            None,
            false,
            super::AlleleNormalization::AsWritten,
            false,
            None,
        )
        .unwrap();
        assert_eq!(by_contig.get("chr1").map(std::vec::Vec::len), Some(1));
        assert_eq!(by_contig["chr1"][0].record_start, 5);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_per_allele_info_numbers_parses_number() {
        let headers = vec![
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"alt count\">".to_string(),
            "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"allele depth\">".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">".to_string(),
        ];
        let map = super::per_allele_info_numbers(&headers);
        assert_eq!(map.get("AC"), Some(&'A'));
        assert_eq!(map.get("AD"), Some(&'R'));
        assert_eq!(map.get("DP"), None); // fixed-Number field omitted
    }

    #[test]
    fn test_filter_and_subset_original_info_multiallelic() {
        let headers = vec![
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">".to_string(),
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">".to_string(),
            "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"x\">".to_string(),
        ];
        let map = super::per_allele_info_numbers(&headers);
        let info = "AC=5,7;AF=0.5,0.7;AD=10,5,7;DP=22;FOO=bar";
        // Split allele index 1: per-allele fields subset to that ALT.
        let out = super::filter_and_subset_original_info(info, |k| k == "GENE", &map, 1).unwrap();
        assert_eq!(out, "AC=7;AF=0.7;AD=10,7;DP=22;FOO=bar");
    }

    #[test]
    fn test_filter_and_subset_original_info_biallelic_unchanged() {
        let headers = vec![
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">".to_string(),
            "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"x\">".to_string(),
        ];
        let map = super::per_allele_info_numbers(&headers);
        let out =
            super::filter_and_subset_original_info("AC=5;AD=10,5", |_| false, &map, 0).unwrap();
        assert_eq!(out, "AC=5;AD=10,5");
    }

    #[test]
    fn insertion_overlaps_only_between_interval_bases() {
        let insertion = variant(10, "A", "AT");
        assert!(insertion.overlaps_interval(9, 11));
        assert!(insertion.overlaps_interval(10, 11));
        assert!(!insertion.overlaps_interval(11, 12));
        assert!(!insertion.overlaps_interval(10, 10));
    }

    #[test]
    fn anchored_deletion_overlaps_deleted_reference_span() {
        let deletion = variant(10, "AT", "A");
        assert!(deletion.overlaps_interval(11, 11));
        assert!(deletion.overlaps_interval(10, 11));
        assert!(!deletion.overlaps_interval(10, 10));
    }

    #[test]
    fn symbolic_deletions_and_duplications_are_spelled_out_from_the_reference() {
        //               1234567890
        let sequence = "acgtacgtac";
        // END removes 4-6; the reference's lowercase is written in uppercase.
        assert_eq!(
            spell_out_symbolic(3, "G", "<DEL>", "SVTYPE=DEL;END=6", 0, Some(sequence)),
            Some((3, "GTAC".to_string(), "G".to_string()))
        );
        // SVLEN is read per ALT, its sign ignored, when END is missing.
        assert_eq!(
            spell_out_symbolic(3, "G", "<DEL:ME:ALU>", "SVLEN=-5,-2", 1, Some(sequence)),
            Some((3, "GTA".to_string(), "G".to_string()))
        );
        // A tandem copy of 4-6 goes in after 6.
        assert_eq!(
            spell_out_symbolic(3, "G", "<DUP:TANDEM>", "END=6", 0, Some(sequence)),
            Some((6, "C".to_string(), "CTAC".to_string()))
        );
    }

    #[test]
    fn symbolic_alleles_that_cannot_be_spelled_out_stay_symbolic() {
        let sequence = "acgtacgtac";
        let spell = |alt: &str, ref_base: &str, info: &str, seq: Option<&str>| {
            spell_out_symbolic(3, ref_base, alt, info, 0, seq)
        };
        assert_eq!(spell("<INV>", "G", "END=6", Some(sequence)), None);
        assert_eq!(spell("<DEL>", "G", "SVTYPE=DEL", Some(sequence)), None);
        assert_eq!(spell("<DEL>", "G", "END=11", Some(sequence)), None);
        assert_eq!(spell("<DEL>", "G", "END=3", Some(sequence)), None);
        assert_eq!(spell("<DEL>", "T", "END=6", Some(sequence)), None);
        assert_eq!(spell("<DEL>", "G", "END=6", None), None);
        assert_eq!(spell("A", "G", "END=6", Some(sequence)), None);

        let long = "acgt".repeat(5_000);
        let over = format!("END={}", 3 + MAX_SPELLED_SYMBOLIC_SPAN + 1);
        assert_eq!(spell("<DEL>", "G", &over, Some(&long)), None);
        let at_limit = format!("END={}", 3 + MAX_SPELLED_SYMBOLIC_SPAN);
        assert!(spell("<DEL>", "G", &at_limit, Some(&long)).is_some());
    }

    /// Trimming must not move an indel. Suffix trimming on a length-changing
    /// allele re-anchors it: `29 CT>CTGCT` became `29 C>CTGC`, the same
    /// insertion at a different placement inside a repeat, which no longer
    /// matched the CIGAR the aligner used. A row with 20 of 20 supporting reads
    /// was then reported at 0, with a warning blaming the input file.
    #[test]
    fn normalising_an_indel_keeps_its_anchor() {
        let (pos, ref_allele, alt_allele) = super::normalize_ref_alt(29, "CT", "CTGCT");
        assert_eq!(
            (pos, ref_allele.as_str(), alt_allele.as_str()),
            (30, "T", "TGCT"),
            "the insertion must stay anchored after position 30"
        );

        // A deletion written with trailing shared context, same rule.
        let (pos, ref_allele, alt_allele) = super::normalize_ref_alt(29, "CTGCT", "CT");
        assert_eq!(
            (pos, ref_allele.as_str(), alt_allele.as_str()),
            (30, "TGCT", "T")
        );

        // Found by `normalising_preserves_the_decomposed_event`: here the
        // canonical suffix-then-prefix trim IS safe, and refusing it would leave
        // a longer allele than it needs to be. The rule is not "never trim an
        // indel", it is "keep the trim that leaves the event where it was".
        let (pos, ref_allele, alt_allele) = super::normalize_ref_alt(5, "CG", "CTG");
        assert_eq!(
            (pos, ref_allele.as_str(), alt_allele.as_str()),
            (5, "C", "CT"),
            "an unambiguous insertion should reach its minimal anchored form"
        );
    }

    /// A substitution has no anchor to lose, so it is trimmed from both ends as
    /// the canonical form requires.
    #[test]
    fn normalising_a_substitution_trims_both_ends() {
        let (pos, ref_allele, alt_allele) = super::normalize_ref_alt(10, "CTA", "CGA");
        assert_eq!(
            (pos, ref_allele.as_str(), alt_allele.as_str()),
            (11, "T", "G")
        );
    }
}
