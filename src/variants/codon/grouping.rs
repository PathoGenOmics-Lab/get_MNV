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
/// Merge `original_info` from all SNPs in a codon group.
/// When all SNPs share the same info string, return that single string.
/// When they differ, concatenate unique info strings with `|` as separator
/// so no original INFO data is lost for MNV records.
pub(super) fn merge_original_info(snps: &[Snp]) -> Option<String> {
    let infos: Vec<&str> = snps
        .iter()
        .filter_map(|s| s.original_info.as_deref())
        .collect();
    if infos.is_empty() {
        return None;
    }
    // Deduplicate while preserving order
    let mut seen = std::collections::HashSet::new();
    let unique: Vec<&str> = infos.into_iter().filter(|s| seen.insert(*s)).collect();
    Some(unique.join("|"))
}

pub(super) fn collect_all_usize(values: impl Iterator<Item = Option<usize>>) -> Option<Vec<usize>> {
    let mut out = Vec::new();
    for value in values {
        out.push(value?);
    }
    Some(out)
}

pub(super) fn collect_all_f64(values: impl Iterator<Item = Option<f64>>) -> Option<Vec<f64>> {
    let mut out = Vec::new();
    for value in values {
        out.push(value?);
    }
    Some(out)
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
    event_metadata(variant.position, &variant.ref_allele, &variant.alt_allele)
}
