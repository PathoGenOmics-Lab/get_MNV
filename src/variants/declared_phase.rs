//! What the input VCF's own phasing says about a multi-position row.
//!
//! get_MNV does not phase. When the caller did, that claim is worth carrying
//! through: long-read pipelines (WhatsHap, HiPhase, Clair3) publish `GT` with
//! `|` and a `PS` phase set, and until now get_MNV ignored it, so a codon MNV
//! the caller had already resolved arrived unannotated. It is a claim, not an
//! observation, which is why it is reported beside the read evidence rather
//! than folded into it.

use crate::io::{DeclaredPhase, VcfPosition};
use crate::variants::event::AlleleComponentKind;

use super::types::{DeclaredPhaseCall, LinkageVerdict};

/// The phase the input declared for the alleles of one row.
///
/// Only records inside the same phase set can be compared, so a row whose
/// positions were phased separately, or not at all, gets nothing. A single
/// pair placed on disjoint haplotypes settles the row as trans: the row claims
/// one molecule carries every one of its alleles, and that pair says otherwise.
pub fn declared_phase_for_row(
    positions: &[usize],
    alternates: &[String],
    snp_list: &[VcfPosition],
) -> Option<DeclaredPhaseCall> {
    if positions.len() < 2 || positions.len() != alternates.len() {
        return None;
    }

    let phases = positions
        .iter()
        .zip(alternates.iter())
        .map(|(position, alternate)| declared_phase_at(*position, alternate, snp_list))
        .collect::<Vec<_>>();

    // Every position must carry a claim, and all of them must sit in one phase
    // set. A row is a statement about all of its alleles together: publishing
    // `cis` because two of three positions were phased together would speak for
    // an allele the caller never placed.
    let mut phase_set = None;
    for phase in &phases {
        let Some(phase) = phase else { return None };
        match (phase_set, phase.phase_set) {
            (_, None) => return None,
            (None, Some(set)) => phase_set = Some(set),
            (Some(current), Some(set)) if current != set => return None,
            (Some(_), Some(_)) => {}
        }
    }

    let mut comparisons = 0usize;
    let mut all_share = true;
    for (index, left) in phases.iter().enumerate() {
        let Some(left) = left else { continue };
        for right in phases.iter().skip(index + 1).flatten() {
            let shares = left.shares_haplotype_with(right)?;
            comparisons += 1;
            all_share &= shares;
        }
    }

    if comparisons == 0 {
        return None;
    }
    Some(DeclaredPhaseCall {
        verdict: if all_share {
            LinkageVerdict::Cis
        } else {
            LinkageVerdict::Trans
        },
        phase_set,
        // Set by the caller once the read evidence is in.
        contradicted_by_reads: false,
    })
}

/// The declared phase of the input record that carries *this* allele at this
/// position.
///
/// Matching on position alone is not enough. A multiallelic site splits into one
/// record per ALT, and a genotype like `2|0` places the second ALT on a
/// haplotype and the first on none at all; a position-only lookup hands the
/// second ALT's claim to the first, inventing a claim the caller never made.
/// The record is matched by an allele component that agrees on both the
/// coordinate and the base, which also covers an MNP decomposed across several
/// row positions.
///
/// The component must also be a substitution. An insertion is anchored at the
/// base it follows and its component carries the *inserted* sequence, so the
/// record `28 G>GT` yields the component `(28, "T")`, byte for byte the one the
/// SNV `28 G>T` yields. Without the kind test the insertion's `GT`/`PS` was
/// published as the SNV's declared phase, so a row claimed a haplotype the input
/// never placed it on, and which of the two records came first in the file
/// decided the answer.
fn declared_phase_at<'a>(
    position: usize,
    alternate: &str,
    snp_list: &'a [VcfPosition],
) -> Option<&'a DeclaredPhase> {
    snp_list
        .iter()
        .find(|variant| {
            variant.declared_phase.is_some()
                && variant.event().components.iter().any(|component| {
                    matches!(component.kind, AlleleComponentKind::Snp)
                        && component.position == position
                        && component.alt_allele.eq_ignore_ascii_case(alternate)
                })
        })
        .and_then(|variant| variant.declared_phase.as_ref())
}

/// Whether the reads flatly contradict what the caller declared.
///
/// Only the two extremes count, because they are facts rather than judgements:
/// not one read carrying both alleles refutes a declared cis, and every
/// informative read carrying both refutes a declared trans. Anything between
/// is a matter of degree, and the phasing-support column already reports it.
pub fn reads_contradict_declared_phase(
    call: &DeclaredPhaseCall,
    haplotype_reads: usize,
    informative_reads: usize,
) -> bool {
    if informative_reads == 0 {
        return false;
    }
    match call.verdict {
        LinkageVerdict::Cis => haplotype_reads == 0,
        LinkageVerdict::Trans => haplotype_reads == informative_reads,
        LinkageVerdict::Unknown => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alts() -> Vec<String> {
        vec!["C".to_string(), "G".to_string()]
    }

    fn phased(position: usize, alt: &str, genotype: &str, phase_set: Option<&str>) -> VcfPosition {
        VcfPosition {
            position,
            ref_allele: "A".to_string(),
            alt_allele: alt.to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
            declared_phase: crate::io::parse_declared_phase(Some(genotype), phase_set, 0),
        }
    }

    #[test]
    fn same_haplotype_in_one_phase_set_is_cis() {
        let snps = vec![
            phased(10, "C", "1|0", Some("7")),
            phased(12, "G", "1|0", Some("7")),
        ];
        let call = declared_phase_for_row(&[10, 12], &alts(), &snps).expect("comparable");
        assert_eq!(call.verdict, LinkageVerdict::Cis);
        assert_eq!(call.phase_set, Some(7));
    }

    #[test]
    fn opposite_haplotypes_in_one_phase_set_are_trans() {
        let snps = vec![
            phased(10, "C", "1|0", Some("7")),
            phased(12, "G", "0|1", Some("7")),
        ];
        let call = declared_phase_for_row(&[10, 12], &alts(), &snps).expect("comparable");
        assert_eq!(call.verdict, LinkageVerdict::Trans);
    }

    #[test]
    fn a_row_whose_other_position_was_never_phased_claims_nothing() {
        // Only position 10 carries a claim. Publishing `cis` off the back of it
        // would speak for the allele at 12, which the caller never placed.
        let mut unphased = phased(12, "G", "1|0", Some("7"));
        unphased.declared_phase = None;
        let snps = vec![phased(10, "C", "1|0", Some("7")), unphased];
        assert!(declared_phase_for_row(&[10, 12], &alts(), &snps).is_none());
    }

    #[test]
    fn different_phase_sets_say_nothing_about_each_other() {
        let snps = vec![
            phased(10, "C", "1|0", Some("7")),
            phased(12, "G", "1|0", Some("9")),
        ];
        assert!(declared_phase_for_row(&[10, 12], &alts(), &snps).is_none());
    }

    #[test]
    fn an_unphased_genotype_declares_nothing() {
        // `/` is the caller saying it did not resolve the phase. Reading it as
        // phase would invent linkage the input never claimed.
        let snps = vec![
            phased(10, "C", "1/0", Some("7")),
            phased(12, "G", "1/0", Some("7")),
        ];
        assert!(declared_phase_for_row(&[10, 12], &alts(), &snps).is_none());
    }

    #[test]
    fn phased_without_a_phase_set_is_not_comparable_across_records() {
        let snps = vec![phased(10, "C", "1|0", None), phased(12, "G", "0|1", None)];
        assert!(declared_phase_for_row(&[10, 12], &alts(), &snps).is_none());
    }

    #[test]
    fn a_split_multiallelic_site_does_not_lend_its_claim_to_the_other_allele() {
        // The site is G>T,C and the genotype is `2|0`: the second ALT (the C)
        // sits on haplotype 0 and the first (the T) on no haplotype at all.
        // Splitting produces two records at position 10, and matching on the
        // coordinate alone would hand the C's claim to the T, inventing a
        // linkage the caller never declared.
        let mut carries_c = phased(10, "C", "1|0", Some("7"));
        carries_c.ref_allele = "G".to_string();
        let mut carries_t = phased(10, "T", "1|0", Some("7"));
        carries_t.ref_allele = "G".to_string();
        carries_t.declared_phase = None; // allele 1 appears nowhere in `2|0`
        let partner = phased(12, "G", "1|0", Some("7"));

        let snps = vec![carries_t, carries_c, partner];
        assert_eq!(
            declared_phase_for_row(&[10, 12], &["C".to_string(), "G".to_string()], &snps)
                .map(|call| call.verdict),
            Some(LinkageVerdict::Cis)
        );
        assert!(
            declared_phase_for_row(&[10, 12], &["T".to_string(), "G".to_string()], &snps).is_none(),
            "the T carries no declared phase, so its row must claim none"
        );
    }

    #[test]
    fn only_the_extremes_count_as_a_contradiction() {
        let cis = DeclaredPhaseCall {
            verdict: LinkageVerdict::Cis,
            phase_set: Some(7),
            contradicted_by_reads: false,
        };
        assert!(reads_contradict_declared_phase(&cis, 0, 20));
        assert!(!reads_contradict_declared_phase(&cis, 1, 20));
        // No informative read is not a contradiction, it is silence.
        assert!(!reads_contradict_declared_phase(&cis, 0, 0));

        let trans = DeclaredPhaseCall {
            verdict: LinkageVerdict::Trans,
            phase_set: Some(7),
            contradicted_by_reads: false,
        };
        assert!(reads_contradict_declared_phase(&trans, 20, 20));
        assert!(!reads_contradict_declared_phase(&trans, 19, 20));
    }
}
