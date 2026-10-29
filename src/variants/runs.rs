//! Runs of adjacent single-base substitutions, wherever they fall.
//!
//! A codon row only sees the substitutions inside its own codon, so two changed
//! bases on either side of a codon boundary came out as two unrelated SNVs, with
//! no doublet class, even when the caller wrote them as a single record. In
//! COSMIC's terms adjacent substitutions are one event: two make a doublet
//! (DBS78), three or more a multi-base substitution (MBS), and neither belongs
//! in the single-base spectrum. This pass looks across all of a contig's rows at
//! once, which is the only place a run that straddles two rows can be seen.
//!
//! Runs are read from positions, as SigProfiler reads them: two adjacent
//! substitutions in one sample are one event. Read support for the run as a
//! whole is not counted; a codon row's own MNV columns still say what the reads
//! showed inside that codon.

use super::{parse_component_label, AlleleComponentKind, VariantInfo, VariantType};
use std::collections::{BTreeMap, HashMap};

/// A maximal stretch of consecutive substituted bases, at least two long.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubstitutionRun {
    pub start: usize,
    pub end: usize,
    pub ref_bases: String,
    pub alt_bases: String,
}

impl SubstitutionRun {
    pub fn length(&self) -> usize {
        self.end - self.start + 1
    }

    /// HGVS genomic descriptor of the whole run: adjacent substitutions on one
    /// allele are a deletion-insertion in HGVS, not a list of substitutions.
    pub fn descriptor(&self) -> String {
        format!("g.{}_{}delins{}", self.start, self.end, self.alt_bases)
    }

    /// COSMIC DBS78 class, for a run of exactly two bases.
    pub fn dbs_class(&self) -> Option<String> {
        if self.length() != 2 {
            return None;
        }
        crate::utils::dbs_class(&self.ref_bases, &self.alt_bases)
    }
}

/// What the spectrum calls a base that belongs to a run: `DBS` inside a doublet,
/// `MBS` inside a longer run.
pub fn run_marker(length: usize) -> &'static str {
    if length == 2 {
        "DBS"
    } else {
        "MBS"
    }
}

fn substitutions(variant: &VariantInfo) -> impl Iterator<Item = (usize, char, char)> + '_ {
    variant
        .event_components
        .iter()
        .filter_map(|label| parse_component_label(label))
        .filter(|component| component.kind == AlleleComponentKind::Snp)
        .filter_map(|component| {
            let mut r = component.ref_allele.chars();
            let mut a = component.alt_allele.chars();
            match (r.next(), r.next(), a.next(), a.next()) {
                (Some(r), None, Some(a), None) => Some((
                    component.position,
                    r.to_ascii_uppercase(),
                    a.to_ascii_uppercase(),
                )),
                _ => None,
            }
        })
}

/// Every run of two or more adjacent substitutions among `variants`.
///
/// Indel rows are left out: a complex indel's substitution parts are part of an
/// indel event, and a substitution that also stands on its own has its own row.
/// A position with two different alternates (a split multiallelic site) is not
/// one base of a run, since the two cannot both be on the molecule.
pub fn find_runs(variants: &[VariantInfo]) -> Vec<SubstitutionRun> {
    let mut by_position: BTreeMap<usize, Option<(char, char)>> = BTreeMap::new();
    for variant in variants
        .iter()
        .filter(|v| v.variant_type != VariantType::Indel)
    {
        for (position, r, a) in substitutions(variant) {
            by_position
                .entry(position)
                .and_modify(|seen| {
                    if seen.is_some_and(|(_, alt)| alt != a) {
                        *seen = None;
                    }
                })
                .or_insert(Some((r, a)));
        }
    }

    let mut runs = Vec::new();
    let mut current: Option<SubstitutionRun> = None;
    for (&position, base) in &by_position {
        let Some((r, a)) = *base else {
            runs.extend(current.take());
            continue;
        };
        match current.as_mut() {
            Some(run) if run.end + 1 == position => {
                run.end = position;
                run.ref_bases.push(r);
                run.alt_bases.push(a);
            }
            _ => {
                runs.extend(current.take());
                current = Some(SubstitutionRun {
                    start: position,
                    end: position,
                    ref_bases: r.to_string(),
                    alt_bases: a.to_string(),
                });
            }
        }
    }
    runs.extend(current);
    runs.retain(|run| run.length() >= 2);
    runs
}

/// Annotate every row with the runs its substitutions belong to, and give the
/// doublet class to every row of a doublet. Returns, by position, the length of
/// the run each base belongs to, for the single-base spectrum to leave out.
pub fn annotate_adjacent_substitutions(variants: &mut [VariantInfo]) -> HashMap<usize, usize> {
    let runs = find_runs(variants);
    let mut run_at: HashMap<usize, usize> = HashMap::new();
    for (index, run) in runs.iter().enumerate() {
        for position in run.start..=run.end {
            run_at.insert(position, index);
        }
    }

    for variant in variants.iter_mut() {
        if variant.variant_type == VariantType::Indel {
            continue;
        }
        let mut touched: Vec<usize> = Vec::new();
        for (position, _, _) in substitutions(variant) {
            if let Some(&index) = run_at.get(&position) {
                if !touched.contains(&index) {
                    touched.push(index);
                }
            }
        }
        touched.sort_by_key(|&index| runs[index].start);
        variant.annotations.adjacent_substitutions = touched
            .iter()
            .map(|&index| runs[index].descriptor())
            .collect();
        // The doublet class comes from the run, not from the codon: two
        // substituted bases of one codon next to a third in the neighbouring
        // codon are a run of three, which COSMIC does not call a doublet, and a
        // doublet split by a codon boundary is one.
        let classes: Vec<String> = touched
            .iter()
            .filter_map(|&index| runs[index].dbs_class())
            .collect();
        let in_a_run = !touched.is_empty();
        if in_a_run {
            variant.annotations.dbs_class = (!classes.is_empty()).then(|| classes.join(", "));
        } else if variant.annotations.dbs_class.is_some() {
            // A codon's own doublet that is not a run here: one of its two bases
            // also carries another alternate, so the two are not one molecule.
            variant.annotations.dbs_class = None;
        }
    }

    run_at
        .into_iter()
        .map(|(position, index)| (position, runs[index].length()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::{ChangeType, VariantAnnotations};

    fn row(variant_type: VariantType, components: &[&str]) -> VariantInfo {
        VariantInfo {
            chrom: "chr1".to_string(),
            gene: "geneA".to_string(),
            positions: vec![1],
            ref_bases: vec!["A".to_string()],
            base_changes: vec!["C".to_string()],
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
            event_components: components.iter().map(|c| c.to_string()).collect(),
            annotations: VariantAnnotations::default(),
        }
    }

    #[test]
    fn a_doublet_across_a_codon_boundary_is_one_event_on_both_rows() {
        let mut rows = vec![
            row(VariantType::Snp, &["SNV:103:G>A"]),
            row(VariantType::Snp, &["SNV:104:G>T"]),
            row(VariantType::Snp, &["SNV:200:C>T"]),
        ];
        let members = annotate_adjacent_substitutions(&mut rows);
        for r in &rows[..2] {
            assert_eq!(
                r.annotations.adjacent_substitutions,
                vec!["g.103_104delinsAT"]
            );
            // GG>AT is CC>AT read from the other strand.
            assert_eq!(r.annotations.dbs_class.as_deref(), Some("CC>AT"));
        }
        assert!(rows[2].annotations.adjacent_substitutions.is_empty());
        assert_eq!(members.get(&103), Some(&2));
        assert_eq!(members.get(&200), None);
    }

    #[test]
    fn three_in_a_row_are_a_multi_base_substitution_not_a_doublet() {
        // A row's own doublet (28, 29) next to a substituted base on another
        // row (30): three adjacent bases, one multi-base substitution.
        let mut rows = vec![
            row(VariantType::SnpMnv, &["SNV:28:G>T", "SNV:29:C>A"]),
            row(VariantType::Snp, &["SNV:30:T>G"]),
        ];
        rows[0].annotations.dbs_class = Some("GC>TA".to_string());
        let members = annotate_adjacent_substitutions(&mut rows);
        assert_eq!(
            rows[0].annotations.adjacent_substitutions,
            vec!["g.28_30delinsTAG"]
        );
        assert_eq!(rows[0].annotations.dbs_class, None);
        assert_eq!(rows[1].annotations.dbs_class, None);
        assert_eq!(run_marker(members[&29]), "MBS");
    }

    #[test]
    fn a_split_multiallelic_base_breaks_the_run() {
        let mut rows = vec![
            row(VariantType::Snp, &["SNV:10:A>C"]),
            row(VariantType::Snp, &["SNV:10:A>G"]),
            row(VariantType::Snp, &["SNV:11:T>A"]),
        ];
        annotate_adjacent_substitutions(&mut rows);
        assert!(rows
            .iter()
            .all(|r| r.annotations.adjacent_substitutions.is_empty()));
    }

    #[test]
    fn indel_rows_neither_join_nor_receive_a_run() {
        let mut rows = vec![
            row(VariantType::Indel, &["DEL:9:AT", "SNV:12:C>T"]),
            row(VariantType::Snp, &["SNV:13:G>A"]),
        ];
        annotate_adjacent_substitutions(&mut rows);
        assert!(rows
            .iter()
            .all(|r| r.annotations.adjacent_substitutions.is_empty()));
    }
}
