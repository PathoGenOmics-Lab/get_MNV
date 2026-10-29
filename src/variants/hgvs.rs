//! HGVS genomic (`g.`) and coding (`c.`) sequence-variant descriptors.
//!
//! Genomic descriptors are built from the variant's genomic REF/ALT and are
//! computed at output time. Coding descriptors need the gene strand and the CDS
//! offset, so they are built during annotation (in the codon paths) and carried
//! on [`crate::variants::VariantAnnotations::hgvs_c`].
//!
//! Multi-nucleotide substitutions (MNVs) use the HGVS allele-bracket form
//! `g.[POS1REF>ALT;POS2REF>ALT]`, which represents co-occurring substitutions
//! without needing the unchanged bases between them. Coding descriptors cover
//! substitutions and indels, with inserted bases on the coding strand. Following
//! HGVS, a deletion or insertion inside a repeat is placed at its 3'-most
//! equivalent position (in the genomic sequence for `g.`, along the coding
//! sequence for `c.`, never beyond the CDS), and an insertion that repeats the
//! bases just before it is written as a duplication (`dup`). A base in an intron
//! of a spliced transcript is numbered from the nearest exon (`c.100+2`,
//! `c.101-5`); UTR and flanking positions (`c.-15`, `c.*3`) are not built, and
//! a change reaching them gets no `c.` descriptor. Neither form carries a
//! reference-sequence accession prefix; the Chromosome and Gene columns supply
//! that context.

use crate::utils::reverse_complement;
use crate::variants::{Strand, VariantInfo, VariantType};

/// Common prefix and suffix lengths shared by a REF/ALT pair, with the two
/// windows kept disjoint so REF/ALT split cleanly into prefix, core, and suffix.
/// Bases compare without regard to case, as they do in a VCF: a lowercase REF
/// beside an uppercase ALT shares its anchor all the same.
fn trim_common(ref_bases: &[u8], alt_bases: &[u8]) -> (usize, usize) {
    let max = ref_bases.len().min(alt_bases.len());
    let mut prefix = 0;
    while prefix < max && ref_bases[prefix].eq_ignore_ascii_case(&alt_bases[prefix]) {
        prefix += 1;
    }
    let mut suffix = 0;
    while suffix < max - prefix
        && ref_bases[ref_bases.len() - 1 - suffix]
            .eq_ignore_ascii_case(&alt_bases[alt_bases.len() - 1 - suffix])
    {
        suffix += 1;
    }
    (prefix, suffix)
}

/// HGVS genomic descriptor for an indel from its VCF-style anchored REF/ALT at
/// genomic `pos` (1-based). Emits `del` / `ins` / `delins` from the trimmed core.
fn indel_genomic(pos: usize, ref_allele: &str, alt_allele: &str) -> String {
    let (prefix, suffix) = trim_common(ref_allele.as_bytes(), alt_allele.as_bytes());
    let ref_core = &ref_allele[prefix..ref_allele.len() - suffix];
    let alt_core = alt_allele[prefix..alt_allele.len() - suffix].to_ascii_uppercase();
    let start = pos + prefix; // genomic position of the first changed REF base
    if ref_core.is_empty() {
        // Pure insertion, placed between start-1 and start.
        format!("g.{}_{}ins{}", start.saturating_sub(1), start, alt_core)
    } else {
        let end = start + ref_core.len() - 1;
        let range = if ref_core.len() == 1 {
            format!("{start}")
        } else {
            format!("{start}_{end}")
        };
        if alt_core.is_empty() {
            format!("g.{range}del")
        } else {
            format!("g.{range}delins{alt_core}")
        }
    }
}

fn substitution_genomic(variant: &VariantInfo) -> Option<String> {
    // `>` describes exactly one nucleotide in HGVS. An annotated codon row
    // already carries one entry per changed base, but a fallback row (intergenic,
    // intron, splice site, a position inside a gene that could not be given a
    // codon) carries the record's whole REF and ALT, and formatting those with
    // `>` produced `g.28GC>GA`: not a valid descriptor, naming a base that did
    // not change, at a coordinate the change did not happen at. Decompose such
    // an entry into the bases that actually differ, which is what the row's own
    // Event Components column says.
    let mut items: Vec<(usize, String, String)> = Vec::new();
    for ((&pos, r), a) in variant
        .positions
        .iter()
        .zip(variant.ref_bases.iter())
        .zip(variant.base_changes.iter())
    {
        if r.chars().count() == 1 && a.chars().count() == 1 {
            items.push((pos, r.clone(), a.clone()));
            continue;
        }
        for component in crate::variants::decompose_allele(pos, r, a).components {
            if component.kind == crate::variants::AlleleComponentKind::Snp {
                items.push((
                    component.position,
                    component.ref_allele,
                    component.alt_allele,
                ));
            }
        }
    }
    if items.is_empty() {
        return None;
    }
    items.sort_by_key(|(pos, _, _)| *pos);
    if items.len() == 1 {
        let (pos, r, a) = &items[0];
        Some(format!("g.{pos}{r}>{a}"))
    } else {
        let inner = items
            .iter()
            .map(|(pos, r, a)| format!("{pos}{r}>{a}"))
            .collect::<Vec<_>>()
            .join(";");
        Some(format!("g.[{inner}]"))
    }
}

/// HGVS genomic descriptor for a variant, or `None` when it cannot be built
/// (no positions, or REF/ALT vectors whose lengths disagree with the positions).
pub fn genomic(variant: &VariantInfo) -> Option<String> {
    if variant.positions.is_empty()
        || variant.positions.len() != variant.ref_bases.len()
        || variant.positions.len() != variant.base_changes.len()
    {
        return None;
    }
    match variant.variant_type {
        VariantType::Indel => {
            let pos = *variant.positions.first()?;
            let alt_allele = variant.base_changes.first()?;
            if is_symbolic(alt_allele) {
                return None;
            }
            Some(indel_genomic(pos, variant.ref_bases.first()?, alt_allele))
        }
        _ => substitution_genomic(variant),
    }
}

/// An ALT that names no sequence (`<INV>`, a breakend, the `*` of an
/// overlapping deletion): there are no bases to write a descriptor from. A
/// symbolic deletion or duplication is spelled out from the reference before it
/// gets here, and that one has its bases.
fn is_symbolic(allele: &str) -> bool {
    allele.starts_with('<') || allele.contains('[') || allele.contains(']') || allele == "*"
}

/// Build the HGVS coding-substitution descriptor from per-SNV
/// `(cds_position_1based, coding_strand_ref, coding_strand_alt)` entries. A
/// single SNV yields `c.POSREF>ALT`; multiple yield the allele bracket
/// `c.[...;...]`, ordered by ascending coding position. `None` for an empty slice.
pub(crate) fn coding_substitution(entries: &[(usize, char, char)]) -> Option<String> {
    if entries.is_empty() {
        return None;
    }
    let mut sorted = entries.to_vec();
    sorted.sort_by_key(|&(pos, _, _)| pos);
    if sorted.len() == 1 {
        let (pos, r, a) = sorted[0];
        Some(format!("c.{pos}{r}>{a}"))
    } else {
        let inner = sorted
            .iter()
            .map(|&(pos, r, a)| format!("{pos}{r}>{a}"))
            .collect::<Vec<_>>()
            .join(";");
        Some(format!("c.[{inner}]"))
    }
}

/// Whether REF/ALT, once their shared bases are trimmed, only add bases.
pub(crate) fn is_pure_insertion(ref_allele: &str, alt_allele: &str) -> bool {
    let (prefix, suffix) = trim_common(ref_allele.as_bytes(), alt_allele.as_bytes());
    prefix + suffix == ref_allele.len() && alt_allele.len() > ref_allele.len()
}

/// A position in HGVS coding numbering: a coding base, or an intronic base
/// counted from the nearest exon (`c.100+2` after an exon's last base, `c.101-5`
/// before an exon's first).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CodingPosition {
    Cds(usize),
    Intron { anchor: usize, offset: isize },
}

impl CodingPosition {
    /// Transcript order: `c.100`, `c.100+1`, ..., `c.101-1`, `c.101`.
    fn order(self) -> (usize, isize) {
        match self {
            Self::Cds(n) => (n, 0),
            Self::Intron { anchor, offset } => (anchor, offset),
        }
    }

    /// Whether `next` is the position right after this one along the transcript.
    fn precedes(self, next: Self) -> bool {
        match (self, next) {
            (Self::Cds(a), Self::Cds(b)) => b == a + 1,
            (Self::Cds(a), Self::Intron { anchor, offset }) => anchor == a && offset == 1,
            (Self::Intron { anchor, offset }, Self::Cds(b)) => anchor == b && offset == -1,
            (
                Self::Intron {
                    anchor: a,
                    offset: x,
                },
                Self::Intron {
                    anchor: b,
                    offset: y,
                },
            ) => (a == b && y == x + 1) || (b > a && x > 0 && y < 0),
        }
    }
}

impl std::fmt::Display for CodingPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::Cds(n) => write!(f, "{n}"),
            Self::Intron { anchor, offset } if offset > 0 => write!(f, "{anchor}+{offset}"),
            Self::Intron { anchor, offset } => write!(f, "{anchor}{offset}"),
        }
    }
}

/// A deletion of 1-based `start..=end` moved to its 3'-most equivalent place in
/// `seq`, no further than base `limit`: while the base after it equals its first
/// base, deleting one base further on removes the same sequence.
pub(crate) fn shift_deletion(
    seq: &[u8],
    mut start: usize,
    mut end: usize,
    limit: usize,
) -> (usize, usize) {
    while end < limit.min(seq.len()) && seq[start - 1].eq_ignore_ascii_case(&seq[end]) {
        start += 1;
        end += 1;
    }
    (start, end)
}

/// An insertion of `inserted` after 1-based base `after`, moved to its 3'-most
/// equivalent place in `seq`, no further than after base `limit`, its bases
/// rotated as it moves.
pub(crate) fn shift_insertion(
    seq: &[u8],
    mut after: usize,
    inserted: &str,
    limit: usize,
) -> (usize, String) {
    let mut bases = inserted.to_ascii_uppercase().into_bytes();
    while !bases.is_empty()
        && after < limit.min(seq.len())
        && seq[after].eq_ignore_ascii_case(&bases[0])
    {
        bases.rotate_left(1);
        after += 1;
    }
    (after, String::from_utf8(bases).unwrap_or_default())
}

/// Whether an insertion after base `after` repeats the bases just before it,
/// which HGVS writes as a duplication of those bases.
fn duplicates(seq: &[u8], after: usize, inserted: &str) -> bool {
    let length = inserted.len();
    length > 0
        && after >= length
        && seq[after - length..after].eq_ignore_ascii_case(inserted.as_bytes())
}

fn range(start: CodingPosition, end: CodingPosition) -> String {
    if start == end {
        format!("{start}")
    } else {
        format!("{start}_{end}")
    }
}

/// The coding sequence an indel is moved along, and where each of its exons
/// ends. HGVS moves an indel 3'-most within its exon only: the rule does not
/// cross an exon-exon junction.
pub(crate) struct CodingSequence<'a> {
    pub bases: &'a str,
    /// The `c.` number of the last base of every exon followed by an intron,
    /// ascending. The end of the coding sequence closes the last one.
    pub exon_ends: Vec<usize>,
}

impl CodingSequence<'_> {
    /// The last base of the exon holding coding base `n`.
    fn exon_end(&self, n: usize) -> usize {
        self.exon_ends
            .iter()
            .copied()
            .find(|&end| end >= n)
            .unwrap_or(self.bases.len())
            .min(self.bases.len())
    }

    /// A coding position as `c.` writes it: past the stop codon it is `*N`.
    fn number(&self, n: usize) -> String {
        match n.checked_sub(self.bases.len()) {
            Some(past) if past > 0 => format!("*{past}"),
            _ => n.to_string(),
        }
    }
}

/// HGVS coding (`c.`) descriptor for an indel from its VCF-style anchored
/// REF/ALT at genomic `pos`, given the gene strand, a map from a genomic position
/// to its coding position, and the coding sequence when it is known.
///
/// `None` when the alleles describe no change, or when a deleted base, or either
/// base an insertion falls between, has no coding position: the change then
/// reaches a UTR or the flank. With the coding sequence, a deletion or insertion
/// that lies wholly in coding bases is placed 3'-most along it, within its exon,
/// and an insertion that repeats the bases before it becomes a `dup`. A symbolic
/// ALT has no bases to describe and gets none.
pub(crate) fn coding_indel(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    strand: Strand,
    coding_position: impl Fn(usize) -> Option<CodingPosition>,
    cds: Option<&CodingSequence<'_>>,
) -> Option<String> {
    if is_symbolic(alt_allele) {
        return None;
    }
    let (prefix, suffix) = trim_common(ref_allele.as_bytes(), alt_allele.as_bytes());
    let ref_core = &ref_allele[prefix..ref_allele.len() - suffix];
    let alt_core = &alt_allele[prefix..alt_allele.len() - suffix];
    let start = pos + prefix;
    let inserted = match strand {
        Strand::Plus => alt_core.to_ascii_uppercase(),
        Strand::Minus => reverse_complement(&alt_core.to_ascii_uppercase()),
    };

    if ref_core.is_empty() {
        if alt_core.is_empty() {
            return None;
        }
        let (a, b) = (
            coding_position(start.checked_sub(1)?)?,
            coding_position(start)?,
        );
        let (before, after) = if a.order() <= b.order() {
            (a, b)
        } else {
            (b, a)
        };
        if !before.precedes(after) {
            return None;
        }
        if let (CodingPosition::Cds(after_base), Some(cds)) = (before, cds) {
            let limit = cds.exon_end(after_base);
            let (after_base, shifted) =
                shift_insertion(cds.bases.as_bytes(), after_base, &inserted, limit);
            if duplicates(cds.bases.as_bytes(), after_base, &shifted) {
                let first = CodingPosition::Cds(after_base + 1 - shifted.len());
                return Some(format!(
                    "c.{}dup",
                    range(first, CodingPosition::Cds(after_base))
                ));
            }
            return Some(format!(
                "c.{after_base}_{}ins{shifted}",
                cds.number(after_base + 1)
            ));
        }
        return Some(format!("c.{before}_{after}ins{inserted}"));
    }

    let mut mapped = Vec::with_capacity(ref_core.len());
    for genomic in start..start + ref_core.len() {
        mapped.push(coding_position(genomic)?);
    }
    mapped.sort_by_key(|p| p.order());
    // Consecutive genomic bases are consecutive along the transcript; a gap means
    // the map gave one base two places, which a slippage join can do.
    if !mapped.windows(2).all(|pair| pair[0].precedes(pair[1])) {
        return None;
    }
    let (mut low, mut high) = (*mapped.first()?, *mapped.last()?);
    if alt_core.is_empty() {
        if let (CodingPosition::Cds(a), CodingPosition::Cds(b), Some(cds)) = (low, high, cds) {
            let (a, b) = shift_deletion(cds.bases.as_bytes(), a, b, cds.exon_end(b));
            low = CodingPosition::Cds(a);
            high = CodingPosition::Cds(b);
        }
        Some(format!("c.{}del", range(low, high)))
    } else {
        Some(format!("c.{}delins{inserted}", range(low, high)))
    }
}

/// The changed core of an anchored REF/ALT: where it starts and what it removes
/// and adds once the shared bases are trimmed.
pub(crate) fn indel_core(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> (usize, String, String) {
    let (prefix, suffix) = trim_common(ref_allele.as_bytes(), alt_allele.as_bytes());
    (
        pos + prefix,
        ref_allele[prefix..ref_allele.len() - suffix].to_ascii_uppercase(),
        alt_allele[prefix..alt_allele.len() - suffix].to_ascii_uppercase(),
    )
}

/// A pure deletion or insertion moved 3'-most along the transcript, rightward
/// on the plus strand and leftward on the minus, without leaving the genomic
/// span `lo..=hi`, and returned anchored as a VCF record would be. Anything else
/// comes back as it was given.
pub(crate) fn shift_within(
    seq: &str,
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    strand: Strand,
    lo: usize,
    hi: usize,
) -> (usize, String, String) {
    let unchanged = || (pos, ref_allele.to_string(), alt_allele.to_string());
    let seq = seq.as_bytes();
    let (start, ref_core, alt_core) = indel_core(pos, ref_allele, alt_allele);
    let base = |p: usize| seq[p - 1] as char;
    if alt_core.is_empty() && !ref_core.is_empty() {
        let (mut s, mut e) = (start, start + ref_core.len() - 1);
        if s < 2 || e > seq.len() {
            return unchanged();
        }
        match strand {
            Strand::Plus => {
                while e < hi && e < seq.len() && seq[s - 1].eq_ignore_ascii_case(&seq[e]) {
                    s += 1;
                    e += 1;
                }
            }
            Strand::Minus => {
                while s > lo && s >= 3 && seq[s - 2].eq_ignore_ascii_case(&seq[e - 1]) {
                    s -= 1;
                    e -= 1;
                }
            }
        }
        let anchored: String = (s - 1..=e).map(base).collect();
        return (
            s - 1,
            anchored.to_ascii_uppercase(),
            base(s - 1).to_ascii_uppercase().to_string(),
        );
    }
    if ref_core.is_empty() && !alt_core.is_empty() {
        let mut after = start - 1;
        let mut bases = alt_core.into_bytes();
        if after == 0 || after > seq.len() {
            return unchanged();
        }
        match strand {
            Strand::Plus => {
                while after < hi && after < seq.len() && seq[after].eq_ignore_ascii_case(&bases[0])
                {
                    bases.rotate_left(1);
                    after += 1;
                }
            }
            Strand::Minus => {
                while after > lo && seq[after - 1].eq_ignore_ascii_case(&bases[bases.len() - 1]) {
                    bases.rotate_right(1);
                    after -= 1;
                }
            }
        }
        let anchor = base(after).to_ascii_uppercase().to_string();
        let inserted = String::from_utf8(bases).unwrap_or_default();
        return (after, anchor.clone(), format!("{anchor}{inserted}"));
    }
    unchanged()
}

/// HGVS coding descriptor for a change at a position with no codon of its own,
/// such as a substitution or indel in an intron, from its anchored REF/ALT.
/// Substitutions are numbered like any other (`c.100+5G>A`), several in a
/// bracket; a change mixing an indel with substitutions gets none.
pub(crate) fn coding_noncodon(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    strand: Strand,
    coding_position: impl Fn(usize) -> Option<CodingPosition>,
) -> Option<String> {
    let event = crate::variants::decompose_allele(pos, ref_allele, alt_allele);
    let substitutions: Vec<_> = event
        .components
        .iter()
        .filter(|c| c.kind == crate::variants::AlleleComponentKind::Snp)
        .collect();
    if substitutions.len() != event.components.len() {
        if event.components.len() == 1 {
            return coding_indel(pos, ref_allele, alt_allele, strand, coding_position, None);
        }
        return None;
    }
    let oriented = |base: &str| match strand {
        Strand::Plus => base.to_ascii_uppercase(),
        Strand::Minus => reverse_complement(&base.to_ascii_uppercase()),
    };
    let mut entries: Vec<(CodingPosition, String, String)> = Vec::new();
    for component in substitutions {
        entries.push((
            coding_position(component.position)?,
            oriented(&component.ref_allele),
            oriented(&component.alt_allele),
        ));
    }
    entries.sort_by_key(|(p, _, _)| p.order());
    let items: Vec<String> = entries
        .iter()
        .map(|(p, r, a)| format!("{p}{r}>{a}"))
        .collect();
    match items.as_slice() {
        [] => None,
        [one] => Some(format!("c.{one}")),
        many => Some(format!("c.[{}]", many.join(";"))),
    }
}

/// [`genomic`] with indels placed as HGVS asks: 3'-most along `sequence`, and an
/// insertion that repeats the bases before it written as a `dup`.
pub fn genomic_normalized(variant: &VariantInfo, sequence: &str) -> Option<String> {
    if variant.variant_type != VariantType::Indel {
        return genomic(variant);
    }
    let pos = *variant.positions.first()?;
    let ref_allele = variant.ref_bases.first()?;
    let alt_allele = variant.base_changes.first()?;
    if alt_allele.starts_with('<') {
        return genomic(variant);
    }
    let (prefix, suffix) = trim_common(ref_allele.as_bytes(), alt_allele.as_bytes());
    let ref_core = &ref_allele[prefix..ref_allele.len() - suffix];
    let alt_core = &alt_allele[prefix..alt_allele.len() - suffix];
    let start = pos + prefix;
    let seq = sequence.as_bytes();
    if ref_core.is_empty() && !alt_core.is_empty() && start >= 1 && start - 1 <= seq.len() {
        let (after, shifted) = shift_insertion(seq, start - 1, alt_core, seq.len());
        if duplicates(seq, after, &shifted) {
            let first = after + 1 - shifted.len();
            return Some(if first == after {
                format!("g.{after}dup")
            } else {
                format!("g.{first}_{after}dup")
            });
        }
        return Some(format!("g.{after}_{}ins{shifted}", after + 1));
    }
    if alt_core.is_empty() && !ref_core.is_empty() && start + ref_core.len() - 1 <= seq.len() {
        let (first, last) = shift_deletion(seq, start, start + ref_core.len() - 1, seq.len());
        return Some(if first == last {
            format!("g.{first}del")
        } else {
            format!("g.{first}_{last}del")
        });
    }
    genomic(variant)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::{ChangeType, VariantAnnotations};

    fn variant(
        variant_type: VariantType,
        positions: Vec<usize>,
        ref_bases: &[&str],
        base_changes: &[&str],
    ) -> VariantInfo {
        VariantInfo {
            chrom: "chr1".to_string(),
            gene: "geneA".to_string(),
            positions,
            ref_bases: ref_bases.iter().map(|s| s.to_string()).collect(),
            base_changes: base_changes.iter().map(|s| s.to_string()).collect(),
            aa_changes: vec!["-".to_string()],
            snp_aa_changes: vec!["-".to_string()],
            aa_changes_local: vec!["-".to_string()],
            snp_aa_changes_local: vec!["-".to_string()],
            variant_type,
            change_type: ChangeType::Unknown,
            snp_reads: None,
            snp_forward_reads: None,
            snp_reverse_reads: None,
            mnv_reads: None,
            mnv_forward_reads: None,
            mnv_reverse_reads: None,
            mnv_total_reads: None,
            total_reads: None,
            total_forward_reads: None,
            total_reverse_reads: None,
            mnv_total_forward_reads: None,
            mnv_total_reverse_reads: None,
            mnv_phasing_reads: None,
            ref_codon: None,
            snp_codon: None,
            mnv_codon: None,
            original_dp: None,
            original_freq: None,
            original_info: None,
            event_class: None,
            event_components: Vec::new(),
            annotations: VariantAnnotations::default(),
        }
    }

    #[test]
    fn genomic_snv() {
        let v = variant(VariantType::Snp, vec![100], &["A"], &["G"]);
        assert_eq!(genomic(&v).as_deref(), Some("g.100A>G"));
    }

    #[test]
    fn genomic_mnv_uses_allele_bracket_sorted_by_position() {
        let v = variant(VariantType::SnpMnv, vec![30, 28], &["T", "G"], &["A", "T"]);
        assert_eq!(genomic(&v).as_deref(), Some("g.[28G>T;30T>A]"));
    }

    #[test]
    fn genomic_indels() {
        // Deletion: REF=ATG ALT=A removes TG at 101-102.
        let del = variant(VariantType::Indel, vec![100], &["ATG"], &["A"]);
        assert_eq!(genomic(&del).as_deref(), Some("g.101_102del"));
        // Single-base deletion.
        let del1 = variant(VariantType::Indel, vec![100], &["AT"], &["A"]);
        assert_eq!(genomic(&del1).as_deref(), Some("g.101del"));
        // Insertion: REF=A ALT=ATG inserts TG between 100 and 101.
        let ins = variant(VariantType::Indel, vec![100], &["A"], &["ATG"]);
        assert_eq!(genomic(&ins).as_deref(), Some("g.100_101insTG"));
        // Delins: REF=ATG ALT=AC replaces TG (101-102) with C.
        let delins = variant(VariantType::Indel, vec![100], &["ATG"], &["AC"]);
        assert_eq!(genomic(&delins).as_deref(), Some("g.101_102delinsC"));
    }

    #[test]
    fn coding_substitution_single_and_bracket() {
        assert_eq!(
            coding_substitution(&[(30, 'A', 'G')]).as_deref(),
            Some("c.30A>G")
        );
        // Bracket, ordered by ascending coding position (input out of order).
        assert_eq!(
            coding_substitution(&[(90, 'T', 'C'), (28, 'G', 'A')]).as_deref(),
            Some("c.[28G>A;90T>C]")
        );
        assert_eq!(coding_substitution(&[]), None);
    }

    /// A plus-strand CDS running from genomic 101: CDS position = genomic - 100.
    fn plus_cds(position: usize) -> Option<CodingPosition> {
        (101..=400)
            .contains(&position)
            .then(|| CodingPosition::Cds(position - 100))
    }

    /// A minus-strand CDS whose first coding base is genomic 400.
    fn minus_cds(position: usize) -> Option<CodingPosition> {
        (101..=400)
            .contains(&position)
            .then(|| CodingPosition::Cds(401 - position))
    }

    /// Two exons, 101-200 (c.1-100) and 301-400 (c.101-200), on the plus strand,
    /// with the intron between them numbered from the nearer exon.
    fn spliced(position: usize) -> Option<CodingPosition> {
        match position {
            101..=200 => Some(CodingPosition::Cds(position - 100)),
            301..=400 => Some(CodingPosition::Cds(position - 200)),
            201..=250 => Some(CodingPosition::Intron {
                anchor: 100,
                offset: (position - 200) as isize,
            }),
            251..=300 => Some(CodingPosition::Intron {
                anchor: 101,
                offset: position as isize - 301,
            }),
            _ => None,
        }
    }

    #[test]
    fn coding_indels_on_the_plus_strand() {
        // REF=ATG ALT=A at 110 deletes genomic 111-112, CDS 11-12.
        assert_eq!(
            coding_indel(110, "ATG", "A", Strand::Plus, plus_cds, None).as_deref(),
            Some("c.11_12del")
        );
        assert_eq!(
            coding_indel(110, "AT", "A", Strand::Plus, plus_cds, None).as_deref(),
            Some("c.11del")
        );
        // Insertion after genomic 110 sits between CDS 10 and 11.
        assert_eq!(
            coding_indel(110, "A", "ATG", Strand::Plus, plus_cds, None).as_deref(),
            Some("c.10_11insTG")
        );
        assert_eq!(
            coding_indel(110, "ATG", "AC", Strand::Plus, plus_cds, None).as_deref(),
            Some("c.11_12delinsC")
        );
    }

    #[test]
    fn coding_indels_on_the_minus_strand_count_from_the_other_end() {
        // Genomic 111-112 are CDS 290-289; the range is written low to high.
        assert_eq!(
            coding_indel(110, "ATG", "A", Strand::Minus, minus_cds, None).as_deref(),
            Some("c.289_290del")
        );
        // Inserted TG reads CA on the coding strand, between CDS 290 and 291.
        assert_eq!(
            coding_indel(110, "A", "ATG", Strand::Minus, minus_cds, None).as_deref(),
            Some("c.290_291insCA")
        );
    }

    #[test]
    fn an_indel_reaching_outside_the_cds_gets_no_coding_descriptor() {
        // Deletes genomic 400-401: the second base is past the CDS end.
        assert_eq!(
            coding_indel(399, "ATG", "A", Strand::Plus, plus_cds, None),
            None
        );
        // Inserted after the last coding base, so between CDS 300 and nothing.
        assert_eq!(
            coding_indel(400, "A", "AT", Strand::Plus, plus_cds, None),
            None
        );
        assert_eq!(
            coding_indel(110, "A", "A", Strand::Plus, plus_cds, None),
            None
        );
    }

    #[test]
    fn an_indel_in_an_intron_is_numbered_from_the_nearer_exon() {
        // Deletes 200-202: the exon's last base and the first two intronic ones.
        assert_eq!(
            coding_indel(199, "AGGG", "A", Strand::Plus, spliced, None).as_deref(),
            Some("c.100_100+2del")
        );
        // Deep in the intron's second half: counted back from the next exon.
        assert_eq!(
            coding_indel(294, "ACG", "A", Strand::Plus, spliced, None).as_deref(),
            Some("c.101-6_101-5del")
        );
        // Right after the exon's last base: between c.100 and c.100+1.
        assert_eq!(
            coding_indel(200, "A", "AT", Strand::Plus, spliced, None).as_deref(),
            Some("c.100_100+1insT")
        );
        // An intronic substitution.
        assert_eq!(
            coding_noncodon(205, "G", "A", Strand::Plus, spliced).as_deref(),
            Some("c.100+5G>A")
        );
    }

    #[test]
    fn indels_in_a_repeat_move_to_the_three_prime_end_and_copies_are_duplications() {
        // c.1-12 = ATG CAG CAG CAG: deleting the first CAG (c.4-6) is, by the
        // 3' rule, deleting the last one (c.10-12).
        let cds = CodingSequence {
            bases: "ATGCAGCAGCAG",
            exon_ends: Vec::new(),
        };
        let map = |position: usize| {
            (101..=112)
                .contains(&position)
                .then(|| CodingPosition::Cds(position - 100))
        };
        assert_eq!(
            coding_indel(103, "GCAG", "G", Strand::Plus, map, Some(&cds)).as_deref(),
            Some("c.10_12del")
        );
        // Inserting CAG after c.3 copies the CAG that follows: a duplication of
        // the last copy, c.10_12dup.
        assert_eq!(
            coding_indel(103, "G", "GCAG", Strand::Plus, map, Some(&cds)).as_deref(),
            Some("c.10_12dup")
        );
        // A single base inserted next to its twin is a one-base dup.
        assert_eq!(
            coding_indel(101, "A", "AT", Strand::Plus, map, Some(&cds)).as_deref(),
            Some("c.2dup")
        );
        // Not a repeat: an ordinary insertion.
        assert_eq!(
            coding_indel(101, "A", "AC", Strand::Plus, map, Some(&cds)).as_deref(),
            Some("c.1_2insC")
        );
    }

    /// The 3' rule stops at an exon-exon junction (HGVS's own example: c.3921del,
    /// not c.3922del), and a position past the stop codon is written `*N`.
    #[test]
    fn indels_move_three_prime_within_their_exon_only() {
        // Exon 1 is c.1-6 ATGGAT, exon 2 c.7-12 TCATAA: the T that ends exon 1
        // is repeated by the T that starts exon 2.
        let cds = CodingSequence {
            bases: "ATGGATTCATAA",
            exon_ends: vec![6],
        };
        let map = |position: usize| {
            (101..=112)
                .contains(&position)
                .then(|| CodingPosition::Cds(position - 100))
        };
        assert_eq!(
            coding_indel(105, "AT", "A", Strand::Plus, map, Some(&cds)).as_deref(),
            Some("c.6del")
        );
        assert_eq!(
            coding_indel(106, "T", "TT", Strand::Plus, map, Some(&cds)).as_deref(),
            Some("c.6dup")
        );
        // Without the junction the same deletion moves to c.7.
        let joined = CodingSequence {
            bases: "ATGGATTCATAA",
            exon_ends: Vec::new(),
        };
        assert_eq!(
            coding_indel(105, "AT", "A", Strand::Plus, map, Some(&joined)).as_deref(),
            Some("c.7del")
        );
        // Inserting AAG before the last two bases of TAA moves past the stop:
        // after c.12 comes c.*1, never a c.13.
        assert_eq!(
            coding_indel(110, "T", "TAAG", Strand::Plus, map, Some(&joined)).as_deref(),
            Some("c.12_*1insGAA")
        );
    }

    #[test]
    fn a_symbolic_allele_and_a_lowercase_one_get_what_they_can() {
        let map = |position: usize| Some(CodingPosition::Cds(position));
        assert_eq!(
            coding_indel(40, "T", "<INV>", Strand::Plus, map, None),
            None
        );
        assert_eq!(
            coding_indel(40, "T", "<INV>", Strand::Minus, map, None),
            None
        );
        let inversion = variant(VariantType::Indel, vec![40], &["T"], &["<INV>"]);
        assert_eq!(genomic(&inversion), None);
        // A lowercase REF shares its anchor with an uppercase ALT.
        let lower = variant(VariantType::Indel, vec![13], &["g"], &["GAA"]);
        assert_eq!(genomic(&lower).as_deref(), Some("g.13_14insAA"));
        assert_eq!(
            genomic_normalized(&lower, "TTTTTTTTTTTTGAAAAC").as_deref(),
            Some("g.16_17dup")
        );
    }

    #[test]
    fn genomic_indels_are_normalized_against_the_sequence() {
        //               1234567890
        let sequence = "ACAGCAGCAT";
        // Deleting CAG at 2-4 leaves ACAGCAT, and so does deleting GCA at 7-9:
        // the 3' rule runs through the partial copy at the end of the repeat.
        let del = variant(VariantType::Indel, vec![1], &["ACAG"], &["A"]);
        assert_eq!(
            genomic_normalized(&del, sequence).as_deref(),
            Some("g.7_9del")
        );
        let dup = variant(VariantType::Indel, vec![4], &["G"], &["GCAG"]);
        assert_eq!(
            genomic_normalized(&dup, sequence).as_deref(),
            Some("g.7_9dup")
        );
        let ins = variant(VariantType::Indel, vec![9], &["A"], &["AGG"]);
        assert_eq!(
            genomic_normalized(&ins, sequence).as_deref(),
            Some("g.9_10insGG")
        );
        let snv = variant(VariantType::Snp, vec![2], &["C"], &["T"]);
        assert_eq!(
            genomic_normalized(&snv, sequence).as_deref(),
            Some("g.2C>T")
        );
    }

    /// `>` is a single-nucleotide form in HGVS. A fallback row (intergenic,
    /// intron, splice site, or a position inside a gene with no codon) carries
    /// the record's whole REF and ALT, and formatting those with `>` produced
    /// `g.1250AA>AT`: not valid, naming an unchanged base at a coordinate the
    /// change did not happen at, and contradicting the row's own components.
    #[test]
    fn a_multi_base_substitution_never_uses_the_single_base_form() {
        let left_padded = variant(VariantType::Snp, vec![1250], &["AA"], &["AT"]);
        assert_eq!(genomic(&left_padded).as_deref(), Some("g.1251A>T"));

        // Two bases changed in one record: the allele bracket, one entry each.
        let both = variant(VariantType::Mnv, vec![100], &["AC"], &["GT"]);
        assert_eq!(genomic(&both).as_deref(), Some("g.[100A>G;101C>T]"));
    }
}
