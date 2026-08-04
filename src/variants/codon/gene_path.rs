//! Genomic-coordinate codon construction, processing, and codon-bounds math.

use crate::utils::{grantham_distance, iupac_aa, reverse_complement};
use crate::variants::consequence::{first_residue, substitution_change_type};

/// Coarse severity of an amino-acid change: stop (3) > missense (2) >
/// synonymous (1) > ambiguous (0).
fn aa_severity(orig: char, mutated: char) -> u8 {
    if orig == 'X' || mutated == 'X' {
        0
    } else if orig == mutated {
        // Synonymous, including a retained stop (`*` -> `*`). Must precede the
        // stop branch so a stop-retained change is not scored as stop severity.
        1
    } else if orig == '*' || mutated == '*' {
        3
    } else {
        2
    }
}

/// Build the biologist annotations (Grantham distance for a genuine missense,
/// and the MNV-vs-SNV consequence shift) from the translated amino acids.
/// Shared by the genomic and transcript codon paths.
pub(super) fn compute_annotations(
    orig_aa: &str,
    mut_aa: &str,
    single_aas: &[char],
    change_type: ChangeType,
) -> VariantAnnotations {
    let orig_c = orig_aa.chars().next().unwrap_or('X');
    let mut_c = mut_aa.chars().next().unwrap_or('X');
    VariantAnnotations {
        grantham: if change_type == ChangeType::NonSynonymous {
            grantham_distance(orig_c, mut_c)
        } else {
            None
        },
        consequence_shift: compute_consequence_shift(orig_c, mut_c, single_aas),
        // Set by the caller, which has the raw per-SNV alleles and positions.
        dbs_class: None,
        // Set by the caller for premature stops in a multi-exon transcript.
        nmd: None,
        // Set by the caller, which has the strand and CDS offsets.
        hgvs_c: None,
        // Set by the transcript caller for exonic near-junction substitutions.
        splice: None,
        // This path only builds coding annotations, which by definition have a codon.
        non_coding: None,
    }
}

/// COSMIC-style doublet-base-substitution class for a codon MNV, when it is
/// exactly two genomically-adjacent single-base substitutions. Bases are taken
/// on the reference (genomic) strand; the DBS canonicalisation collapses the two
/// strands, so no strand orientation is applied here. `None` otherwise.
pub(super) fn dbs_class_for_codon(snvs: &[&Snp]) -> Option<String> {
    if snvs.len() != 2 {
        return None;
    }
    let mut pair: Vec<(usize, char, char)> = snvs
        .iter()
        .filter_map(|snv| {
            Some((
                snv.position,
                single_base(&snv.ref_base)?,
                single_base(&snv.base)?,
            ))
        })
        .collect();
    if pair.len() != 2 {
        return None;
    }
    pair.sort_by_key(|(pos, _, _)| *pos);
    let (p0, r0, a0) = pair[0];
    let (p1, r1, a1) = pair[1];
    if p1 != p0 + 1 {
        return None;
    }
    crate::utils::dbs_class(&format!("{r0}{r1}"), &format!("{a0}{a1}"))
}

/// HGVS coding (`c.`) descriptor for a genomic-coordinate codon's substitutions.
/// The CDS position is measured from the (phase-adjusted) feature start on the
/// plus strand and from the feature end on the minus strand, plus the feature's
/// protein offset; bases are taken in the coding-strand orientation.
pub(super) fn coding_substitution_for_codon(
    codon_info: &CodonInfo,
    strand: Strand,
) -> Option<String> {
    let entries: Vec<(usize, char, char)> = codon_info
        .codon_list
        .iter()
        .filter_map(|snp| {
            let local_nt = match strand {
                Strand::Plus => snp.position.checked_sub(codon_info.gene_start)?,
                Strand::Minus => codon_info.gene_end.checked_sub(snp.position)?,
            };
            let cds_pos = codon_info.protein_offset * 3 + local_nt + 1;
            let cref = super::transcript_path::transcript_oriented_base(&snp.ref_base, strand)?;
            let calt = super::transcript_path::transcript_oriented_base(&snp.base, strand)?;
            Some((cds_pos, cref, calt))
        })
        .collect();
    crate::variants::hgvs::coding_substitution(&entries)
}

/// The single uppercase base of a one-character allele, or `None` for empty or
/// multi-base alleles (which cannot form a doublet substitution).
fn single_base(allele: &str) -> Option<char> {
    let mut chars = allele.chars();
    let base = chars.next()?;
    if chars.next().is_some() {
        return None;
    }
    Some(base.to_ascii_uppercase())
}

/// Compare the combined MNV consequence against its individual SNVs.
fn compute_consequence_shift(
    orig: char,
    combined_mut: char,
    single_muts: &[char],
) -> ConsequenceShift {
    if single_muts.len() <= 1 || orig == 'X' {
        return ConsequenceShift::NotApplicable;
    }
    let combined = aa_severity(orig, combined_mut);
    let max_individual = single_muts
        .iter()
        .map(|&m| aa_severity(orig, m))
        .max()
        .unwrap_or(0);
    match combined.cmp(&max_individual) {
        std::cmp::Ordering::Greater => ConsequenceShift::Gained,
        std::cmp::Ordering::Less => ConsequenceShift::Masked,
        std::cmp::Ordering::Equal => ConsequenceShift::Concordant,
    }
}
use crate::variants::{
    ChangeType, CodonInfo, ConsequenceShift, Gene, Snp, Strand, VariantAnnotations, VariantInfo,
    VariantType,
};

use super::grouping::{collect_all_f64, collect_all_usize, merge_original_info};

pub(super) fn construct_codon(codon_info: &CodonInfo, target_snps: &[&Snp]) -> String {
    let mut codon = String::with_capacity(3);
    for i in 0..3 {
        let current_pos = codon_info.codon_start + i;
        if let Some(snp) = target_snps.iter().find(|&&s| s.position == current_pos) {
            // Uppercase ALT alleles to match the reference (ensures the codon
            // lookup table always matches).
            for c in snp.base.chars() {
                codon.push(c.to_ascii_uppercase());
            }
        } else {
            codon.push(codon_info.original_codon.chars().nth(i).unwrap_or('N'));
        }
    }
    codon
}

pub fn process_codon(
    codon_info: CodonInfo,
    strand: Strand,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
) -> VariantInfo {
    // process_codon is only ever called with at least one SNP in the codon
    // (the BTreeMap groups in `get_mnv_variants_for_gene` only ever contain
    // non-empty Vec<Snp>). We still guard here so that an accidental
    // direct caller does not panic via `codon_list[0]` or produce a row
    // claiming a SNP at "position 0".
    if codon_info.codon_list.is_empty() {
        log::error!(
            "process_codon called with empty codon_list for gene '{}' at codon {}-{}; \
             this is a logic bug, please file an issue.",
            codon_info.gene_name,
            codon_info.codon_start,
            codon_info.codon_end
        );
        return VariantInfo {
            chrom: chrom.to_string(),
            gene: codon_info.gene_name,
            positions: Vec::new(),
            ref_bases: Vec::new(),
            base_changes: Vec::new(),
            aa_changes: vec!["-".to_string()],
            snp_aa_changes: vec!["-".to_string()],
            aa_changes_local: vec!["-".to_string()],
            snp_aa_changes_local: vec!["-".to_string()],
            variant_type: VariantType::Snp,
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
            ref_codon: Some(codon_info.original_codon),
            snp_codon: None,
            mnv_codon: None,
            original_dp: None,
            original_freq: None,
            original_info: None,
            event_class: None,
            event_components: Vec::new(),
            annotations: crate::variants::VariantAnnotations::default(),
        };
    }
    let ref_codon = codon_info.original_codon.clone();

    let mnv_snps: Vec<&Snp> = codon_info.codon_list.iter().collect();
    let mnv_codon = construct_codon(&codon_info, &mnv_snps);

    let snp_codons: Vec<String> = codon_info
        .codon_list
        .iter()
        .map(|snp| construct_codon(&codon_info, &[snp]))
        .collect();
    let snp_codon = snp_codons.join(", ");

    let (orig_aa, mut_aa) = match strand {
        Strand::Minus => (
            genetic_code.translate_seq(reverse_complement(&ref_codon).as_bytes()),
            genetic_code.translate_seq(reverse_complement(&mnv_codon).as_bytes()),
        ),
        Strand::Plus => (
            genetic_code.translate_seq(ref_codon.as_bytes()),
            genetic_code.translate_seq(mnv_codon.as_bytes()),
        ),
    };

    let local_aa_pos = match strand {
        Strand::Plus => {
            codon_info.codon_list[0]
                .position
                .saturating_sub(codon_info.gene_start)
                / 3
                + 1
        }
        Strand::Minus => {
            codon_info
                .gene_end
                .saturating_sub(codon_info.codon_list[0].position)
                / 3
                + 1
        }
    };
    let aa_pos = codon_info.protein_offset + local_aa_pos;

    let combined_change = format!("{orig_aa}{aa_pos}{mut_aa}");
    let combined_change_local = format!("{orig_aa}{local_aa_pos}{mut_aa}");
    let combined_aa = iupac_aa(&combined_change);
    let combined_aa_local = iupac_aa(&combined_change_local);
    let change_type =
        substitution_change_type(first_residue(&orig_aa), first_residue(&mut_aa), aa_pos);

    // Single-residue result of each SNV on its own (reused for the per-SNP
    // columns and the MNV-vs-SNV consequence comparison).
    let single_aas: Vec<char> = codon_info
        .codon_list
        .iter()
        .map(|snp| {
            let single_codon = construct_codon(&codon_info, &[snp]);
            let single = match strand {
                Strand::Minus => reverse_complement(&single_codon),
                Strand::Plus => single_codon,
            };
            genetic_code
                .translate_seq(single.as_bytes())
                .chars()
                .next()
                .unwrap_or('X')
        })
        .collect();
    let snp_changes: Vec<String> = single_aas
        .iter()
        .map(|aa| iupac_aa(&format!("{orig_aa}{aa_pos}{aa}")))
        .collect();
    let snp_changes_local: Vec<String> = single_aas
        .iter()
        .map(|aa| iupac_aa(&format!("{orig_aa}{local_aa_pos}{aa}")))
        .collect();

    let mut annotations = compute_annotations(&orig_aa, &mut_aa, &single_aas, change_type);
    annotations.dbs_class = dbs_class_for_codon(&mnv_snps);
    annotations.hgvs_c = coding_substitution_for_codon(&codon_info, strand);

    VariantInfo {
        chrom: chrom.to_string(),
        gene: codon_info.gene_name,
        positions: codon_info.codon_list.iter().map(|s| s.position).collect(),
        ref_bases: codon_info
            .codon_list
            .iter()
            .map(|s| s.ref_base.clone())
            .collect(),
        base_changes: codon_info
            .codon_list
            .iter()
            .map(|s| s.base.clone())
            .collect(),
        aa_changes: vec![combined_aa],
        snp_aa_changes: snp_changes,
        aa_changes_local: vec![combined_aa_local],
        snp_aa_changes_local: snp_changes_local,
        variant_type: if codon_info.codon_list.len() == 1 {
            VariantType::Snp
        } else {
            VariantType::SnpMnv
        },
        change_type,
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
        ref_codon: Some(ref_codon),
        snp_codon: Some(snp_codon),
        mnv_codon: Some(mnv_codon),
        original_dp: collect_all_usize(codon_info.codon_list.iter().map(|s| s.original_dp)),
        original_freq: collect_all_f64(codon_info.codon_list.iter().map(|s| s.original_freq)),
        original_info: merge_original_info(&codon_info.codon_list),
        event_class: Some(if codon_info.codon_list.len() == 1 {
            "snp".to_string()
        } else {
            "mnv".to_string()
        }),
        event_components: codon_info
            .codon_list
            .iter()
            .map(|s| format!("SNV:{}:{}>{}", s.position, s.ref_base, s.base))
            .collect(),
        annotations,
    }
}

/// Return the effective (phase-adjusted) start/end of the gene for codon math.
///
/// GFF column 8 (phase) tells us how many bases must be skipped from the start
/// of the CDS to reach the first complete codon. For features on the plus
/// strand the skip is applied at `gene.start`; for the minus strand, where the
/// biological "start" is `gene.end`, the skip is applied at `gene.end`.
/// Features without phase information (TSV gene files, gene/exon/UTR rows)
/// carry `phase = 0` and behave exactly as before.
pub(super) fn effective_bounds(gene: &Gene) -> (usize, usize) {
    let phase = gene.phase as usize;
    match gene.strand {
        Strand::Plus => (gene.start.saturating_add(phase), gene.end),
        Strand::Minus => (gene.start, gene.end.saturating_sub(phase)),
    }
}

// Process one gene at a time to keep memory usage low.
pub(super) fn codon_bounds_for_position(gene: &Gene, position: usize) -> Option<(usize, usize)> {
    let (eff_start, eff_end) = effective_bounds(gene);
    // The variant fell inside the GFF feature interval but outside the
    // phase-adjusted region (the first `phase` bases of a plus-strand CDS or
    // the last `phase` bases of a minus-strand CDS belong to a codon that
    // physically spans into the *previous* exon — we cannot reconstruct that
    // codon from this single CDS row, so the variant has to be dropped). Warn
    // explicitly so the user knows: silently dropping was the trap that
    // hid the codon-grouping bug behind issue #12 for so long.
    if position < eff_start || position > eff_end {
        if gene.phase > 0 && position >= gene.start && position <= gene.end {
            log::warn!(
                "Variant at {}:{} falls in the phase-skipped region of CDS '{}' \
                 (phase={}, exon {}-{}); the codon spans into a neighbouring exon \
                 and cannot be reconstructed from a single GFF row. Variant skipped.",
                gene.name,
                position,
                gene.name,
                gene.phase,
                gene.start,
                gene.end
            );
        }
        return None;
    }
    let incomplete_codon_log = |reason: &str| {
        log::debug!(
            "SNP at position {} in gene '{}' falls in incomplete codon ({}; gene length {} not multiple of 3, phase {})",
            position, gene.name, reason, eff_end - eff_start + 1, gene.phase
        );
    };
    match gene.strand {
        Strand::Plus => {
            let offset = position - eff_start;
            let codon_start = eff_start + (offset / 3) * 3;
            let codon_end = codon_start + 2;
            if codon_end <= eff_end {
                Some((codon_start, codon_end))
            } else {
                incomplete_codon_log("plus-strand codon end past CDS");
                None
            }
        }
        Strand::Minus => {
            let offset = eff_end - position;
            let codon_end = eff_end.saturating_sub((offset / 3) * 3);
            if codon_end < eff_start + 2 {
                incomplete_codon_log("minus-strand codon would underflow CDS start");
                return None;
            }
            let codon_start = codon_end - 2;
            if codon_start >= eff_start {
                Some((codon_start, codon_end))
            } else {
                incomplete_codon_log("minus-strand codon start before effective start");
                None
            }
        }
    }
}

#[cfg(test)]
mod annotation_tests {
    use super::{aa_severity, compute_consequence_shift, dbs_class_for_codon};
    use crate::variants::{ConsequenceShift, Snp};

    fn snp(position: usize, ref_base: &str, base: &str) -> Snp {
        Snp {
            index: position,
            position,
            ref_base: ref_base.to_string(),
            base: base.to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }
    }

    #[test]
    fn test_dbs_class_for_adjacent_pair_is_order_independent() {
        let a = snp(10, "C", "T");
        let b = snp(11, "C", "T");
        assert_eq!(dbs_class_for_codon(&[&a, &b]).as_deref(), Some("CC>TT"));
        assert_eq!(dbs_class_for_codon(&[&b, &a]).as_deref(), Some("CC>TT"));
    }

    #[test]
    fn test_dbs_class_none_for_non_adjacent_or_single() {
        let a = snp(10, "C", "T");
        let far = snp(12, "C", "T");
        // Same codon, positions 1 and 3: not a genomic doublet.
        assert_eq!(dbs_class_for_codon(&[&a, &far]), None);
        // A single SNV is not a doublet.
        assert_eq!(dbs_class_for_codon(&[&a]), None);
        // Three SNVs are a triplet, not a DBS.
        let b = snp(11, "C", "T");
        assert_eq!(dbs_class_for_codon(&[&a, &b, &far]), None);
    }

    #[test]
    fn test_aa_severity_retained_stop_is_synonymous() {
        // A retained stop (`*` -> `*`) is synonymous (1), not stop severity (3).
        assert_eq!(aa_severity('*', '*'), 1);
        assert_eq!(aa_severity('A', 'A'), 1);
        assert_eq!(aa_severity('A', '*'), 3); // stop gained
        assert_eq!(aa_severity('*', 'A'), 3); // stop lost
        assert_eq!(aa_severity('M', 'A'), 2); // missense
        assert_eq!(aa_severity('X', 'A'), 0);
    }

    #[test]
    fn test_consequence_shift_masked_when_mnv_retains_stop() {
        // Reference stop codon MNV: one SNV alone loses the stop (*->W), the
        // combined MNV retains it (*->*). The MNV masks the stop-loss call.
        assert_eq!(
            compute_consequence_shift('*', '*', &['W', '*']),
            ConsequenceShift::Masked
        );
    }

    #[test]
    fn test_consequence_shift_gained_two_synonymous() {
        // Two individually-synonymous SNVs producing a non-synonymous residue.
        assert_eq!(
            compute_consequence_shift('L', 'F', &['L', 'L']),
            ConsequenceShift::Gained
        );
    }

    #[test]
    fn test_consequence_shift_single_snv_not_applicable() {
        assert_eq!(
            compute_consequence_shift('M', 'A', &['A']),
            ConsequenceShift::NotApplicable
        );
    }
}
