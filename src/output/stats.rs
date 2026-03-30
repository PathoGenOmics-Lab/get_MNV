//! Strand-bias statistics: Fisher exact two-tailed test for SNP and MNV
//! read count distributions.

use crate::variants::VariantInfo;

fn ln_factorial(value: usize) -> f64 {
    if value <= 1 {
        0.0
    } else {
        (2..=value).map(|n| (n as f64).ln()).sum()
    }
}

fn ln_choose(n: usize, k: usize) -> f64 {
    if k > n {
        f64::NEG_INFINITY
    } else {
        ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
    }
}

fn hypergeometric_probability(
    a: usize,
    row1: usize,
    row2: usize,
    col1: usize,
    total: usize,
) -> f64 {
    (ln_choose(row1, a) + ln_choose(row2, col1.saturating_sub(a)) - ln_choose(total, col1)).exp()
}

pub(crate) fn fisher_exact_two_tailed(a: usize, b: usize, c: usize, d: usize) -> f64 {
    let row1 = a + b;
    let row2 = c + d;
    let col1 = a + c;
    let total = row1 + row2;
    if total == 0 {
        return 1.0;
    }

    let observed = hypergeometric_probability(a, row1, row2, col1, total);
    let min_a = col1.saturating_sub(row2);
    let max_a = row1.min(col1);
    let mut p_value = 0.0f64;
    for candidate in min_a..=max_a {
        let p = hypergeometric_probability(candidate, row1, row2, col1, total);
        if p <= observed + 1e-12 {
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
}
