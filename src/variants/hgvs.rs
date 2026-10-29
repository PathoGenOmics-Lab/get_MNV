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
//! substitutions and, when every base they touch lies in the CDS, indels as
//! `del`, `ins` and `delins`, with the inserted bases on the coding strand. An
//! indel that reaches into an intron or a UTR gets no `c.` descriptor rather
//! than an offset form (`c.123+5del`) this module does not build. Neither form
//! is 3'-shifted, an insertion is never rewritten as `dup`, and neither carries
//! a reference-sequence accession prefix; the Chromosome and Gene columns supply
//! that context.

use crate::utils::reverse_complement;
use crate::variants::{Strand, VariantInfo, VariantType};

/// Common prefix and suffix lengths shared by a REF/ALT pair, with the two
/// windows kept disjoint so REF/ALT split cleanly into prefix, core, and suffix.
fn trim_common(ref_bases: &[u8], alt_bases: &[u8]) -> (usize, usize) {
    let max = ref_bases.len().min(alt_bases.len());
    let mut prefix = 0;
    while prefix < max && ref_bases[prefix] == alt_bases[prefix] {
        prefix += 1;
    }
    let mut suffix = 0;
    while suffix < max - prefix
        && ref_bases[ref_bases.len() - 1 - suffix] == alt_bases[alt_bases.len() - 1 - suffix]
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
    let alt_core = &alt_allele[prefix..alt_allele.len() - suffix];
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
            Some(indel_genomic(
                pos,
                variant.ref_bases.first()?,
                variant.base_changes.first()?,
            ))
        }
        _ => substitution_genomic(variant),
    }
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

/// HGVS coding (`c.`) descriptor for an indel from its VCF-style anchored
/// REF/ALT at genomic `pos`, given the gene strand and a map from a genomic
/// position to its 1-based CDS position.
///
/// `None` when the alleles describe no change, or when a deleted base, or
/// either base an insertion falls between, has no CDS position: the change then
/// reaches outside the coding sequence. A deletion must also map onto
/// consecutive CDS positions, which is what tells a deletion that crosses an
/// intron, and so removes a splice junction, from one inside an exon.
pub(crate) fn coding_indel(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    strand: Strand,
    cds_position: impl Fn(usize) -> Option<usize>,
) -> Option<String> {
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
        let before = cds_position(start.checked_sub(1)?)?;
        let after = cds_position(start)?;
        if before.abs_diff(after) != 1 {
            return None;
        }
        return Some(format!(
            "c.{}_{}ins{inserted}",
            before.min(after),
            before.max(after)
        ));
    }

    let mut mapped = Vec::with_capacity(ref_core.len());
    for genomic in start..start + ref_core.len() {
        mapped.push(cds_position(genomic)?);
    }
    let low = *mapped.iter().min()?;
    let high = *mapped.iter().max()?;
    if high - low + 1 != mapped.len() {
        return None;
    }
    let range = if low == high {
        format!("{low}")
    } else {
        format!("{low}_{high}")
    };
    if alt_core.is_empty() {
        Some(format!("c.{range}del"))
    } else {
        Some(format!("c.{range}delins{inserted}"))
    }
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
    fn plus_cds(position: usize) -> Option<usize> {
        (101..=400).contains(&position).then(|| position - 100)
    }

    /// A minus-strand CDS whose first coding base is genomic 400.
    fn minus_cds(position: usize) -> Option<usize> {
        (101..=400).contains(&position).then(|| 401 - position)
    }

    #[test]
    fn coding_indels_on_the_plus_strand() {
        // REF=ATG ALT=A at 110 deletes genomic 111-112, CDS 11-12.
        assert_eq!(
            coding_indel(110, "ATG", "A", Strand::Plus, plus_cds).as_deref(),
            Some("c.11_12del")
        );
        assert_eq!(
            coding_indel(110, "AT", "A", Strand::Plus, plus_cds).as_deref(),
            Some("c.11del")
        );
        // Insertion after genomic 110 sits between CDS 10 and 11.
        assert_eq!(
            coding_indel(110, "A", "ATG", Strand::Plus, plus_cds).as_deref(),
            Some("c.10_11insTG")
        );
        assert_eq!(
            coding_indel(110, "ATG", "AC", Strand::Plus, plus_cds).as_deref(),
            Some("c.11_12delinsC")
        );
    }

    #[test]
    fn coding_indels_on_the_minus_strand_count_from_the_other_end() {
        // Genomic 111-112 are CDS 290-289; the range is written low to high.
        assert_eq!(
            coding_indel(110, "ATG", "A", Strand::Minus, minus_cds).as_deref(),
            Some("c.289_290del")
        );
        // Inserted TG reads CA on the coding strand, between CDS 290 and 291.
        assert_eq!(
            coding_indel(110, "A", "ATG", Strand::Minus, minus_cds).as_deref(),
            Some("c.290_291insCA")
        );
    }

    #[test]
    fn an_indel_reaching_outside_the_cds_gets_no_coding_descriptor() {
        // Deletes genomic 400-401: the second base is past the CDS end.
        assert_eq!(coding_indel(399, "ATG", "A", Strand::Plus, plus_cds), None);
        // Inserted after the last coding base, so between CDS 300 and nothing.
        assert_eq!(coding_indel(400, "A", "AT", Strand::Plus, plus_cds), None);
        // Two exons, 101-200 and 301-400: a deletion over 199-302 maps to CDS
        // positions that are not consecutive, because it removes the intron
        // between them as well.
        let spliced = |position: usize| match position {
            101..=200 => Some(position - 100),
            301..=400 => Some(position - 200),
            _ => None,
        };
        assert_eq!(
            coding_indel(198, "AGG", "A", Strand::Plus, spliced).as_deref(),
            Some("c.99_100del")
        );
        assert_eq!(coding_indel(199, "AGGG", "A", Strand::Plus, spliced), None);
        // After the exon's last base the next genomic base is intronic, so the
        // map alone cannot place the insertion; the codon layer knows the
        // junction rule and supplies that flank itself.
        assert_eq!(coding_indel(200, "A", "AT", Strand::Plus, spliced), None);
        assert_eq!(coding_indel(110, "A", "A", Strand::Plus, plus_cds), None);
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
