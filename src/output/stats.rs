//! Strand-bias statistics: Fisher exact two-tailed test for SNP and MNV
//! read count distributions.

use crate::variants::VariantInfo;

pub(crate) fn fisher_exact_two_tailed(a: usize, b: usize, c: usize, d: usize) -> f64 {
    let row1 = a + b;
    let row2 = c + d;
    let col1 = a + c;
    let total = row1 + row2;
    if total == 0 {
        return 1.0;
    }

    // Tabulate ln(n!) once. Computing it per term made the whole test quadratic
    // in the molecule count, which is fine for a handful of reads and is not
    // fine at the depths intra-host sequencing reaches.
    let mut ln_factorials = Vec::with_capacity(total + 1);
    let mut running = 0.0f64;
    ln_factorials.push(running);
    for n in 1..=total {
        running += (n as f64).ln();
        ln_factorials.push(running);
    }
    let choose = |n: usize, k: usize| -> f64 {
        if k > n {
            f64::NEG_INFINITY
        } else {
            ln_factorials[n] - ln_factorials[k] - ln_factorials[n - k]
        }
    };
    let probability = |candidate: usize| -> f64 {
        (choose(row1, candidate) + choose(row2, col1.saturating_sub(candidate))
            - choose(total, col1))
        .exp()
    };

    /// Slack for the "as extreme as" comparison, in units of the observed
    /// probability, so it stays meaningful however small that probability is.
    const RELATIVE_TOLERANCE: f64 = 1e-7;

    let observed = probability(a);
    let min_a = col1.saturating_sub(row2);
    let max_a = row1.min(col1);
    let mut p_value = 0.0f64;
    for candidate in min_a..=max_a {
        let p = probability(candidate);
        // "At least as extreme as the observed table" has to be judged
        // relatively. An absolute slack swallows every candidate once the
        // observed probability itself drops below it, so any p-value under
        // roughly 1e-12 was floored rather than computed, and two records could
        // come out with the stronger strand bias carrying the larger p-value.
        if p <= observed * (1.0 + RELATIVE_TOLERANCE) {
            p_value += p;
        }
    }
    p_value.min(1.0)
}

pub(crate) fn snp_strand_bias_p_value(variant: &VariantInfo, index: usize) -> Option<f64> {
    let snp_forward = *variant.snp_forward_reads.as_ref()?.get(index)?;
    let snp_reverse = *variant.snp_reverse_reads.as_ref()?.get(index)?;
    let total_forward = *variant.total_forward_reads.as_ref()?.get(index)?;
    let total_reverse = *variant.total_reverse_reads.as_ref()?.get(index)?;
    let ref_forward = total_forward.saturating_sub(snp_forward);
    let ref_reverse = total_reverse.saturating_sub(snp_reverse);
    Some(fisher_exact_two_tailed(
        snp_forward,
        snp_reverse,
        ref_forward,
        ref_reverse,
    ))
}

pub(crate) fn mnv_strand_bias_p_value(variant: &VariantInfo) -> Option<f64> {
    let mnv_forward = variant.mnv_forward_reads?;
    let mnv_reverse = variant.mnv_reverse_reads?;
    let total_forward = variant.mnv_total_forward_reads?;
    let total_reverse = variant.mnv_total_reverse_reads?;
    let ref_forward = total_forward.saturating_sub(mnv_forward);
    let ref_reverse = total_reverse.saturating_sub(mnv_reverse);
    Some(fisher_exact_two_tailed(
        mnv_forward,
        mnv_reverse,
        ref_forward,
        ref_reverse,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fisher_exact_all_zeros() {
        let p = fisher_exact_two_tailed(0, 0, 0, 0);
        assert!((p - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_fisher_exact_symmetric() {
        let p1 = fisher_exact_two_tailed(5, 5, 5, 5);
        let p2 = fisher_exact_two_tailed(5, 5, 5, 5);
        assert!((p1 - p2).abs() < 1e-10);
    }

    #[test]
    fn test_fisher_exact_extreme_bias() {
        let p = fisher_exact_two_tailed(10, 0, 0, 10);
        assert!(p < 0.01);
    }

    #[test]
    fn test_fisher_exact_no_bias() {
        let p = fisher_exact_two_tailed(5, 5, 5, 5);
        assert!(p > 0.5);
    }

    #[test]
    fn test_fisher_exact_p_value_bounded() {
        for a in 0..=5 {
            for b in 0..=5 {
                for c in 0..=5 {
                    for d in 0..=5 {
                        let p = fisher_exact_two_tailed(a, b, c, d);
                        assert!(
                            (0.0..=1.0 + 1e-10).contains(&p),
                            "p={p} for ({a},{b},{c},{d})"
                        );
                    }
                }
            }
        }
    }

    /// Against scipy's `fisher_exact`, including tables extreme enough that the
    /// p-value falls far below any absolute slack.
    ///
    /// The "at least as extreme" test used an absolute tolerance, so once the
    /// observed probability dropped under it every candidate table qualified and
    /// the p-value was floored near 1e-12 instead of computed. A 400/0/0/400
    /// split came out ~1e-12 where the answer is ~1e-239.
    #[test]
    fn fisher_matches_an_independent_implementation_on_extreme_tables() {
        // (a, b, c, d, scipy's two-tailed p)
        let cases: &[(usize, usize, usize, usize, f64)] = &[
            (12, 1, 10, 10, 2.159082e-02),
            (30, 0, 0, 30, 1.691123e-17),
            (60, 0, 0, 60, 2.070074e-35),
            (100, 0, 0, 100, 2.208761e-59),
            (200, 1, 1, 200, 1.967060e-115),
            (5, 5, 5, 5, 1.0),
            (50, 2, 3, 45, 3.031762e-22),
            // The pair that made --min-strand-bias-p contradict itself: the
            // stronger bias was floored to 1.57e-13 and passed a 1e-13
            // threshold that the weaker one, at 1.58e-14, was tagged by.
            (25, 0, 0, 25, 1.582146e-14),
            (40, 0, 0, 40, 1.860340e-23),
            (400, 0, 0, 400, 1.063590e-239),
        ];
        for &(a, b, c, d, expected) in cases {
            let got = fisher_exact_two_tailed(a, b, c, d);
            let relative = (got - expected).abs() / expected;
            assert!(
                relative < 1e-4,
                "fisher({a},{b},{c},{d}) = {got:e}, expected {expected:e} (relative error {relative:e})"
            );
        }
    }

    /// The ordering the strand-bias filter depends on: a more extreme table must
    /// never report a larger p-value than a less extreme one.
    #[test]
    fn a_stronger_bias_never_reports_a_weaker_p_value() {
        let weak = fisher_exact_two_tailed(25, 0, 0, 25);
        let strong = fisher_exact_two_tailed(40, 0, 0, 40);
        // Both sit below a 1e-13 strand-bias threshold, so both must be tagged;
        // the floor let the stronger one through.
        assert!(
            weak < 1e-13 && strong < 1e-13,
            "weak={weak:e} strong={strong:e}"
        );
        assert!(
            strong < weak,
            "stronger bias gave p={strong:e} against p={weak:e} for the weaker one"
        );
    }
}
