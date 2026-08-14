//! SNP-to-codon grouping and per-variant event metadata.

use crate::variants::Snp;
use std::collections::BTreeSet;

/// Merge a new SNP into a set of alternative codon interpretations.
///
/// `groups` is a list of mutually-exclusive interpretations of the SNVs at
/// a single codon. Each interpretation is a Vec<Snp> of compatible SNVs.
///
/// For an incoming `snp`:
///   - For each existing interpretation that has no SNP at `snp.position`,
///     add `snp` to that interpretation (always compatible).
///   - For each existing interpretation that already has `snp.position` with
///     the SAME alt base, keep it unchanged (true duplicate).
///   - For each existing interpretation that already has `snp.position` with
///     a DIFFERENT alt base (multi-allelic alternative), KEEP the original
///     and ADD a new interpretation that replaces that position with `snp`.
///
/// The result is the Cartesian product over multi-allelic positions, with
/// interpretations dedup'd by the sorted (pos, alt) key.
pub(super) fn merge_snp_into_groups(groups: &mut Vec<Vec<Snp>>, snp: Snp) {
    if groups.is_empty() {
        groups.push(vec![snp]);
        return;
    }
    let mut new_groups: Vec<Vec<Snp>> = Vec::with_capacity(groups.len() * 2);
    let mut seen: BTreeSet<Vec<(usize, String)>> = BTreeSet::new();

    let push_if_new =
        |g: Vec<Snp>, new_groups: &mut Vec<Vec<Snp>>, seen: &mut BTreeSet<Vec<(usize, String)>>| {
            let mut key: Vec<(usize, String)> =
                g.iter().map(|s| (s.position, s.base.clone())).collect();
            key.sort();
            if seen.insert(key) {
                new_groups.push(g);
            }
        };

    for group in groups.iter() {
        if let Some(idx) = group.iter().position(|s| s.position == snp.position) {
            if group[idx].base == snp.base {
                push_if_new(group.clone(), &mut new_groups, &mut seen);
            } else {
                // Multi-allelic alternative: keep original and add a clone
                // with the alt replaced. Both are valid interpretations.
                push_if_new(group.clone(), &mut new_groups, &mut seen);
                let mut alt_group = group.clone();
                alt_group[idx] = snp.clone();
                push_if_new(alt_group, &mut new_groups, &mut seen);
            }
        } else {
            let mut new_group = group.clone();
            new_group.push(snp.clone());
            push_if_new(new_group, &mut new_groups, &mut seen);
        }
    }
    *groups = new_groups;
}
/// Representative `original_info` for the SNPs that form a codon group.
///
/// When the SNVs come from VCF records with *different* INFO, two independent
/// INFO field-sets cannot be combined into one valid INFO column: joining them
/// (the previous behaviour) produced malformed entries such as `AF=0.5|DP=12`
/// and duplicate keys once re-split on `;`/`=`. Use the first record's INFO as
/// the representative for the combined row; when every SNV shares the same INFO
/// (the common case) this is exactly that shared string.
pub(super) fn merge_original_info(snps: &[Snp]) -> Option<String> {
    snps.iter()
        .filter_map(|s| s.original_info.as_deref())
        .next()
        .map(str::to_string)
}

/// The per-SNV values of a codon group, each entry keeping its own presence.
///
/// These used to discard the whole vector the moment one entry was missing, so
/// a substitution that declared its depth and frequency lost both because the
/// other substitution in its codon declared neither. Whether a measurement the
/// input supplied survived depended on an unrelated neighbour landing in the
/// same codon, and the loss reached the per-SNV rows of the constituent that
/// did declare them. The vector is absent only when no entry has a value.
pub(super) fn collect_all_usize(
    values: impl Iterator<Item = Option<usize>>,
) -> Option<Vec<Option<usize>>> {
    let out: Vec<Option<usize>> = values.collect();
    out.iter().any(Option::is_some).then_some(out)
}

pub(super) fn collect_all_f64(
    values: impl Iterator<Item = Option<f64>>,
) -> Option<Vec<Option<f64>>> {
    let out: Vec<Option<f64>> = values.collect();
    out.iter().any(Option::is_some).then_some(out)
}

pub(super) fn event_metadata(
    position: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> (String, Vec<String>) {
    let event = crate::variants::decompose_allele(position, ref_allele, alt_allele);
    (event.class.as_str().to_string(), event.component_labels())
}

pub(super) fn variant_event_metadata(variant: &crate::io::VcfPosition) -> (String, Vec<String>) {
    event_metadata(
        variant.record_start,
        &variant.ref_allele,
        &variant.alt_allele,
    )
}

#[cfg(test)]
mod tests {
    use super::merge_original_info;
    use crate::variants::Snp;

    fn snp(original_info: Option<&str>) -> Snp {
        Snp {
            index: 1,
            position: 1,
            ref_base: "A".to_string(),
            base: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: original_info.map(str::to_string),
        }
    }

    #[test]
    fn test_merge_original_info_shared_value() {
        let snps = vec![snp(Some("DP=10;AF=0.5")), snp(Some("DP=10;AF=0.5"))];
        assert_eq!(merge_original_info(&snps).as_deref(), Some("DP=10;AF=0.5"));
    }

    #[test]
    fn test_merge_original_info_divergent_keeps_first_without_pipe() {
        // Regression: differing INFO from two source records must not be joined
        // with '|' (which corrupts the INFO column once re-split on ';'/'='); the
        // first record's INFO is used as the representative for the MNV row.
        let snps = vec![snp(Some("DP=10;AF=0.5")), snp(Some("DP=12;AF=0.6"))];
        let merged = merge_original_info(&snps).expect("some info");
        assert_eq!(merged, "DP=10;AF=0.5");
        assert!(!merged.contains('|'));
    }

    #[test]
    fn test_merge_original_info_none_when_absent() {
        assert_eq!(merge_original_info(&[snp(None), snp(None)]), None);
    }
}
