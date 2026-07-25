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
//! substitutions only; the protein-level HGVS (`p.`, in the AA columns) already
//! conveys indel consequences. Neither form is 3'-shifted, and neither carries a
//! reference-sequence accession prefix — the Chromosome and Gene columns supply
//! that context.

use crate::variants::{VariantInfo, VariantType};

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
    let mut items: Vec<(usize, &str, &str)> = variant
        .positions
        .iter()
        .zip(variant.ref_bases.iter())
        .zip(variant.base_changes.iter())
        .map(|((&pos, r), a)| (pos, r.as_str(), a.as_str()))
        .collect();
    if items.is_empty() {
        return None;
    }
    items.sort_by_key(|&(pos, _, _)| pos);
    if items.len() == 1 {
        let (pos, r, a) = items[0];
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
}
