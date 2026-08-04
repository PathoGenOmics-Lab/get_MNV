//! What the input VCF's own phasing says about a multi-position row.
//!
//! get_MNV does not phase. When the caller did, that claim is worth carrying
//! through: long-read pipelines (WhatsHap, HiPhase, Clair3) publish `GT` with
//! `|` and a `PS` phase set, and until now get_MNV ignored it, so a codon MNV
//! the caller had already resolved arrived unannotated. It is a claim, not an
//! observation, which is why it is reported beside the read evidence rather
//! than folded into it.

use crate::io::{DeclaredPhase, VcfPosition};

use super::types::{DeclaredPhaseCall, LinkageVerdict};

/// The phase the input declared for the alleles of one row.
///
/// Only records inside the same phase set can be compared, so a row whose
/// positions were phased separately, or not at all, gets nothing. A single
/// pair placed on disjoint haplotypes settles the row as trans: the row claims
/// one molecule carries every one of its alleles, and that pair says otherwise.
pub fn declared_phase_for_row(
    positions: &[usize],
    snp_list: &[VcfPosition],
) -> Option<DeclaredPhaseCall> {
    if positions.len() < 2 {
        return None;
    }

    let phases = positions
        .iter()
        .map(|position| declared_phase_at(*position, snp_list))
        .collect::<Vec<_>>();

    let mut phase_set = None;
    let mut comparisons = 0usize;
    let mut all_share = true;
    for (index, left) in phases.iter().enumerate() {
        let Some(left) = left else { continue };
        for right in phases.iter().skip(index + 1).flatten() {
            let Some(shares) = left.shares_haplotype_with(right) else {
                continue;
            };
            comparisons += 1;
            phase_set = left.phase_set;
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

/// The declared phase of whichever input record covers this position. A record
/// may be an MNP that was decomposed into several row positions, so the match
/// is by overlap rather than by an exact start.
fn declared_phase_at(position: usize, snp_list: &[VcfPosition]) -> Option<&DeclaredPhase> {
    snp_list
        .iter()
        .find(|variant| {
            variant.declared_phase.is_some() && variant.overlaps_interval(position, position)
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
        let call = declared_phase_for_row(&[10, 12], &snps).expect("comparable");
        assert_eq!(call.verdict, LinkageVerdict::Cis);
        assert_eq!(call.phase_set, Some(7));
    }

    #[test]
    fn opposite_haplotypes_in_one_phase_set_are_trans() {
        let snps = vec![
            phased(10, "C", "1|0", Some("7")),
            phased(12, "G", "0|1", Some("7")),
        ];
        let call = declared_phase_for_row(&[10, 12], &snps).expect("comparable");
        assert_eq!(call.verdict, LinkageVerdict::Trans);
    }

    #[test]
    fn different_phase_sets_say_nothing_about_each_other() {
        let snps = vec![
            phased(10, "C", "1|0", Some("7")),
            phased(12, "G", "1|0", Some("9")),
        ];
        assert!(declared_phase_for_row(&[10, 12], &snps).is_none());
    }

    #[test]
    fn an_unphased_genotype_declares_nothing() {
        // `/` is the caller saying it did not resolve the phase. Reading it as
        // phase would invent linkage the input never claimed.
        let snps = vec![
            phased(10, "C", "1/0", Some("7")),
            phased(12, "G", "1/0", Some("7")),
        ];
        assert!(declared_phase_for_row(&[10, 12], &snps).is_none());
    }

    #[test]
    fn phased_without_a_phase_set_is_not_comparable_across_records() {
        let snps = vec![phased(10, "C", "1|0", None), phased(12, "G", "0|1", None)];
        assert!(declared_phase_for_row(&[10, 12], &snps).is_none());
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
