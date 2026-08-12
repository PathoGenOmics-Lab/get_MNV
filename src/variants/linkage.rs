//! Linkage disequilibrium between the substitutions of one codon, measured on
//! the molecules that carry them.
//!
//! A bare co-occurrence ratio cannot tell a haplotype from an accident of
//! frequency. Two substitutions each present on 95% of molecules will be found
//! together on about 90% of them whether or not they have anything to do with
//! each other, and a ratio reports that as near-perfect linkage. The question
//! worth asking is not "how often do they co-occur" but "how much more often
//! than their own frequencies predict", which is exactly what linkage
//! disequilibrium measures.
//!
//! `D'` normalises that excess by the most it could have been given the two
//! frequencies, so it lands in `[-1, 1]` and is comparable between sites:
//!
//! - `+1`: the alleles travel together as far as their frequencies allow. One
//!   haplotype, and the codon-level MNV is real.
//! - `~0`: they co-occur exactly as often as chance predicts. Two variants that
//!   happen to share a codon, which is not the same thing as an MNV.
//! - `-1`: they exclude each other. Both are present, never on one molecule.
//!   In a haploid population that is two competing lineages, the signature of a
//!   mixed infection or of distinct sub-populations rather than of one variant.
//!
//! The Fisher exact p-value says whether the departure from independence is
//! more than the depth could produce by chance, so a `D'` of 1.0 from four
//! molecules is not read as a `D'` of 1.0 from four hundred.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

use crate::output::stats::fisher_exact_two_tailed;

/// Read-level linkage between the substitutions of one codon.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CodonLinkage {
    /// Normalised linkage disequilibrium, in `[-1, 1]`.
    pub d_prime: f64,
    /// Squared correlation of the two alleles across molecules, in `[0, 1]`.
    pub r_squared: f64,
    /// Two-tailed Fisher exact p-value for the 2x2 table of molecule counts.
    pub p_value: f64,
    /// Molecules the table was built from: those observing every position of
    /// the codon.
    pub molecules: usize,
}

/// The 2x2 table of molecules for one pair of positions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PairTable {
    /// Molecules carrying both alternates.
    pub both: usize,
    /// Molecules carrying the first alternate only.
    pub first_only: usize,
    /// Molecules carrying the second alternate only.
    pub second_only: usize,
    /// Molecules carrying neither.
    pub neither: usize,
}

impl PairTable {
    pub fn total(&self) -> usize {
        self.both + self.first_only + self.second_only + self.neither
    }
}

/// Linkage for one pair of positions.
///
/// `None` when the question is not answerable rather than when the answer is
/// zero: no molecule observed both positions, or one of the two alleles is
/// carried by every molecule or by none, in which case there is no variation to
/// correlate and `D'` has no value (its denominator is zero). Reporting `0.0`
/// there would claim independence was measured and found.
pub fn pair_linkage(table: PairTable) -> Option<CodonLinkage> {
    let total = table.total();
    if total == 0 {
        return None;
    }
    let n = total as f64;
    let p_first = (table.both + table.first_only) as f64 / n;
    let p_second = (table.both + table.second_only) as f64 / n;

    // A locus every molecule shares, or none does, carries no information about
    // linkage with anything.
    let variance_first = p_first * (1.0 - p_first);
    let variance_second = p_second * (1.0 - p_second);
    if variance_first <= 0.0 || variance_second <= 0.0 {
        return None;
    }

    let observed_both = table.both as f64 / n;
    let disequilibrium = observed_both - p_first * p_second;

    // The largest excess (or deficit) these two frequencies could produce.
    let max_disequilibrium = if disequilibrium >= 0.0 {
        (p_first * (1.0 - p_second)).min((1.0 - p_first) * p_second)
    } else {
        (p_first * p_second).min((1.0 - p_first) * (1.0 - p_second))
    };
    if max_disequilibrium <= 0.0 {
        return None;
    }

    Some(CodonLinkage {
        d_prime: (disequilibrium / max_disequilibrium).clamp(-1.0, 1.0),
        r_squared: ((disequilibrium * disequilibrium) / (variance_first * variance_second))
            .clamp(0.0, 1.0),
        p_value: fisher_exact_two_tailed(
            table.both,
            table.first_only,
            table.second_only,
            table.neither,
        ),
        molecules: total,
    })
}

/// Linkage for a whole codon, from the observed haplotype patterns.
///
/// Each entry pairs a pattern (which alternates a molecule carries, in the
/// row's position order) with the number of molecules showing it; only
/// molecules observing every position are included.
///
/// A codon-level MNV asserts that one molecule carries *all* of its
/// substitutions, so the claim is only as good as its weakest pair: the pair
/// with the lowest `D'` is the one reported. With two positions there is one
/// pair and this is simply that pair.
pub fn codon_linkage(patterns: &[(Vec<bool>, usize)], positions: usize) -> Option<CodonLinkage> {
    if positions < 2 {
        return None;
    }
    let mut weakest: Option<CodonLinkage> = None;
    for first in 0..positions {
        for second in (first + 1)..positions {
            let mut table = PairTable::default();
            for (pattern, count) in patterns {
                if pattern.len() != positions {
                    continue;
                }
                let slot = match (pattern[first], pattern[second]) {
                    (true, true) => &mut table.both,
                    (true, false) => &mut table.first_only,
                    (false, true) => &mut table.second_only,
                    (false, false) => &mut table.neither,
                };
                *slot += count;
            }
            let Some(linkage) = pair_linkage(table) else {
                continue;
            };
            if weakest.is_none_or(|current| linkage.d_prime < current.d_prime) {
                weakest = Some(linkage);
            }
        }
    }
    weakest
}

/// The joint distribution restricted to some of the window's variants.
///
/// A local haplotype row is usually a subset of the variants its window holds,
/// and the linkage question for that row is about *its* variants: whether they
/// travel together. Projecting the observed combinations onto those positions
/// and merging the ones that become identical asks exactly that, over the
/// molecules that observed the whole window.
pub fn restricted_to(
    patterns: &[(Vec<bool>, usize)],
    indices: &[usize],
) -> Vec<(Vec<bool>, usize)> {
    let mut merged: BTreeMap<Vec<bool>, usize> = BTreeMap::new();
    for (pattern, count) in patterns {
        if indices.iter().any(|index| *index >= pattern.len()) {
            continue;
        }
        let projected = indices.iter().map(|index| pattern[*index]).collect();
        *merged.entry(projected).or_default() += count;
    }
    merged.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(value: f64, expected: f64) {
        assert!(
            (value - expected).abs() < 1e-9,
            "expected {expected}, got {value}"
        );
    }

    #[test]
    fn perfect_co_occurrence_is_complete_linkage() {
        let linkage = pair_linkage(PairTable {
            both: 20,
            first_only: 0,
            second_only: 0,
            neither: 20,
        })
        .expect("both alleles vary");
        approx(linkage.d_prime, 1.0);
        approx(linkage.r_squared, 1.0);
        assert!(linkage.p_value < 1e-6);
    }

    #[test]
    fn alleles_that_never_share_a_molecule_are_complete_repulsion() {
        // Both present on half the molecules, never together. In a haploid
        // population this is two competing lineages, not one variant.
        let linkage = pair_linkage(PairTable {
            both: 0,
            first_only: 20,
            second_only: 20,
            neither: 0,
        })
        .expect("both alleles vary");
        approx(linkage.d_prime, -1.0);
        assert!(linkage.p_value < 1e-6);
    }

    #[test]
    fn two_common_alleles_that_merely_coincide_are_not_linked() {
        // Each allele is on 90% of molecules and they co-occur on 81% of them,
        // exactly what independence predicts. A co-occurrence ratio calls this
        // 0.9 and reads as strong linkage; the excess over chance is zero.
        let linkage = pair_linkage(PairTable {
            both: 81,
            first_only: 9,
            second_only: 9,
            neither: 1,
        })
        .expect("both alleles vary");
        approx(linkage.d_prime, 0.0);
        approx(linkage.r_squared, 0.0);
        approx(linkage.p_value, 1.0);
    }

    #[test]
    fn a_rare_allele_riding_inside_a_common_one_is_completely_linked() {
        // Three molecules carry the rare allele and all three also carry the
        // common one. D' is 1: given that the rare allele is rare, it could not
        // be more associated than this. r squared is small, because the two
        // frequencies are far apart, which is why both are reported.
        let linkage = pair_linkage(PairTable {
            both: 3,
            first_only: 0,
            second_only: 60,
            neither: 37,
        })
        .expect("both alleles vary");
        approx(linkage.d_prime, 1.0);
        assert!(linkage.r_squared < 0.1, "got {}", linkage.r_squared);
    }

    #[test]
    fn an_asymmetric_table_pins_the_arithmetic() {
        // 100 molecules: 30 carry both, 10 the first only, 20 the second only,
        // 40 neither. So p1 = 0.4, p2 = 0.5 and they co-occur on 0.30 where
        // independence predicts 0.20. The excess is 0.10 and the most it could
        // have been is min(0.4*0.5, 0.6*0.5) = 0.20, so D' is exactly 0.5, and
        // r squared is 0.01 / (0.4*0.6*0.5*0.5) = 1/6. Neither number survives
        // a changed operator, which the symmetric cases above do.
        let linkage = pair_linkage(PairTable {
            both: 30,
            first_only: 10,
            second_only: 20,
            neither: 40,
        })
        .expect("both alleles vary");
        approx(linkage.d_prime, 0.5);
        approx(linkage.r_squared, 1.0 / 6.0);
        assert_eq!(linkage.molecules, 100);
    }

    #[test]
    fn an_asymmetric_deficit_pins_the_other_branch() {
        // The mirror image: they co-occur on 0.10 where 0.20 is predicted. The
        // deficit is normalised by min(p1*p2, (1-p1)*(1-p2)) = 0.20 instead,
        // which is the branch the positive case never reaches.
        let linkage = pair_linkage(PairTable {
            both: 10,
            first_only: 30,
            second_only: 40,
            neither: 20,
        })
        .expect("both alleles vary");
        approx(linkage.d_prime, -0.5);
        approx(linkage.r_squared, 1.0 / 6.0);
    }

    #[test]
    fn the_other_half_of_the_maximum_decides_when_it_is_the_smaller() {
        // p1 = 0.6, p2 = 0.5, so the two candidates for the largest possible
        // excess are 0.6*0.5 = 0.30 and 0.4*0.5 = 0.20 and the second wins.
        // The tables above are all shaped so the first wins, which leaves the
        // second expression untested and any change to it invisible.
        // Observed co-occurrence 0.35 against 0.30 predicted, so D' = 0.25.
        let linkage = pair_linkage(PairTable {
            both: 35,
            first_only: 25,
            second_only: 15,
            neither: 25,
        })
        .expect("both alleles vary");
        approx(linkage.d_prime, 0.25);
    }

    #[test]
    fn the_other_half_of_the_deficit_maximum_decides_too() {
        // p1 = 0.6, p2 = 0.7: the candidates are 0.42 and 0.4*0.3 = 0.12, and
        // the second wins. Observed 0.35 against 0.42 predicted, so the deficit
        // is 0.07 and D' = -0.07/0.12.
        let linkage = pair_linkage(PairTable {
            both: 35,
            first_only: 25,
            second_only: 35,
            neither: 5,
        })
        .expect("both alleles vary");
        approx(linkage.d_prime, -0.07 / 0.12);
    }

    #[test]
    fn a_two_position_codon_is_its_only_pair() {
        // The ordinary case, and the one the module exists for. Nothing here
        // reached codon_linkage with two positions, so a guard that turned
        // exactly that case away went unnoticed.
        let patterns = vec![
            (vec![true, true], 30),
            (vec![true, false], 10),
            (vec![false, true], 20),
            (vec![false, false], 40),
        ];
        let linkage = codon_linkage(&patterns, 2).expect("both alleles vary");
        approx(linkage.d_prime, 0.5);
        assert_eq!(linkage.molecules, 100);
    }

    #[test]
    fn a_locus_carried_by_every_molecule_has_no_linkage_to_report() {
        // Nothing varies at the second locus, so there is no correlation to
        // measure. That is unknown, not independent.
        assert!(pair_linkage(PairTable {
            both: 20,
            first_only: 0,
            second_only: 20,
            neither: 0,
        })
        .is_none());
    }

    #[test]
    fn no_molecule_spanning_the_pair_yields_nothing() {
        assert!(pair_linkage(PairTable::default()).is_none());
    }

    #[test]
    fn the_weakest_pair_decides_a_three_position_codon() {
        // Positions 0 and 1 travel together perfectly; position 2 is
        // independent of both. The codon claims one molecule carries all three,
        // and that claim is only as good as the weakest pair.
        let patterns = vec![
            (vec![true, true, true], 25),
            (vec![true, true, false], 25),
            (vec![false, false, true], 25),
            (vec![false, false, false], 25),
        ];
        let linkage = codon_linkage(&patterns, 3).expect("some pair varies");
        approx(linkage.d_prime, 0.0);
    }

    #[test]
    fn restricting_merges_the_combinations_that_become_identical() {
        // Three variants observed; asking about the first two only. The two
        // combinations that differ solely in the third become one.
        let patterns = vec![
            (vec![true, true, true], 12),
            (vec![true, true, false], 8),
            (vec![false, false, true], 5),
        ];
        assert_eq!(
            restricted_to(&patterns, &[0, 1]),
            vec![(vec![false, false], 5), (vec![true, true], 20)]
        );
    }

    #[test]
    fn a_haplotype_inside_another_is_still_perfectly_linked() {
        // 12 molecules carry all three variants and 8 carry only the first two.
        // Asking about the first two: all 20 carry both, so they are as linked
        // as they could be, even though only 8 molecules are the two-variant
        // species. The read count and the linkage answer different questions
        // and both are right.
        let patterns = vec![(vec![true, true, true], 12), (vec![true, true, false], 8)];
        assert!(codon_linkage(&restricted_to(&patterns, &[0, 1]), 2).is_none());
    }

    #[test]
    fn a_single_position_row_has_no_pair() {
        assert!(codon_linkage(&[(vec![true], 10)], 1).is_none());
    }
}
