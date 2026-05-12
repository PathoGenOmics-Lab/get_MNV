//! Codon-level variant processing: SNP grouping, MNV detection, and amino
//! acid change calculation.

use super::types::*;
use crate::utils::{determine_change_type, iupac_aa, reverse_complement};
use std::collections::BTreeMap;

/// Merge `original_info` from all SNPs in a codon group.
/// When all SNPs share the same info string, return that single string.
/// When they differ, concatenate unique info strings with `|` as separator
/// so no original INFO data is lost for MNV records.
fn merge_original_info(snps: &[Snp]) -> Option<String> {
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

fn collect_all_usize(values: impl Iterator<Item = Option<usize>>) -> Option<Vec<usize>> {
    let mut out = Vec::new();
    for value in values {
        out.push(value?);
    }
    Some(out)
}

fn collect_all_f64(values: impl Iterator<Item = Option<f64>>) -> Option<Vec<f64>> {
    let mut out = Vec::new();
    for value in values {
        out.push(value?);
    }
    Some(out)
}

fn construct_codon(codon_info: &CodonInfo, target_snps: &[&Snp]) -> String {
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
    let change_type = ChangeType::from_label(&determine_change_type(&combined_change));

    let snp_changes: Vec<String> = codon_info
        .codon_list
        .iter()
        .map(|snp| {
            let single_codon = construct_codon(&codon_info, &[snp]);
            let single = match strand {
                Strand::Minus => reverse_complement(&single_codon),
                Strand::Plus => single_codon,
            };
            let single_aa = genetic_code.translate_seq(single.as_bytes());
            iupac_aa(&format!("{orig_aa}{aa_pos}{single_aa}"))
        })
        .collect();
    let snp_changes_local: Vec<String> = codon_info
        .codon_list
        .iter()
        .map(|snp| {
            let single_codon = construct_codon(&codon_info, &[snp]);
            let single = match strand {
                Strand::Minus => reverse_complement(&single_codon),
                Strand::Plus => single_codon,
            };
            let single_aa = genetic_code.translate_seq(single.as_bytes());
            iupac_aa(&format!("{orig_aa}{local_aa_pos}{single_aa}"))
        })
        .collect();

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
fn effective_bounds(gene: &Gene) -> (usize, usize) {
    let phase = gene.phase as usize;
    match gene.strand {
        Strand::Plus => (gene.start.saturating_add(phase), gene.end),
        Strand::Minus => (gene.start, gene.end.saturating_sub(phase)),
    }
}

// Process one gene at a time to keep memory usage low.
fn codon_bounds_for_position(gene: &Gene, position: usize) -> Option<(usize, usize)> {
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

pub fn get_mnv_variants_for_gene(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
) -> Vec<VariantInfo> {
    let mut variants = Vec::new();

    let mut codon_to_snps: BTreeMap<usize, Vec<crate::variants::Snp>> = BTreeMap::new();
    let mut indels: Vec<crate::io::VcfPosition> = Vec::new();

    for variant in snp_list
        .iter()
        .filter(|snp| snp.position >= gene.start && snp.position <= gene.end)
    {
        let is_snp = variant.ref_allele.len() == 1
            && variant.alt_allele.len() == 1
            && !variant.alt_allele.starts_with('<');
        if is_snp {
            if let Some((codon_start, _)) = codon_bounds_for_position(gene, variant.position) {
                let group = codon_to_snps.entry(codon_start).or_default();
                // After --split-multiallelic, two ALT alleles at the same
                // position produce separate VCF records.  Only keep the first
                // ALT per position within each codon group — otherwise
                // construct_codon would overwrite the same base twice, and the
                // MNV annotation would be incorrect.
                if !group.iter().any(|s| s.position == variant.position) {
                    group.push(crate::variants::Snp {
                        index: variant.position,
                        position: variant.position,
                        ref_base: variant.ref_allele.clone(),
                        base: variant.alt_allele.clone(),
                        original_dp: variant.original_dp,
                        original_freq: variant.original_freq,
                        original_info: variant.original_info.clone(),
                    });
                } else {
                    log::debug!(
                        "Skipping duplicate ALT at position {} in gene '{}' (multi-allelic split)",
                        variant.position,
                        gene.name
                    );
                }
            }
        } else {
            indels.push(variant.clone());
        }
    }

    if codon_to_snps.is_empty() && indels.is_empty() {
        return variants;
    }

    let mut codon_starts: Vec<usize> = codon_to_snps.keys().copied().collect();
    match gene.strand {
        Strand::Plus => codon_starts.sort_unstable(),
        Strand::Minus => codon_starts.sort_unstable_by(|a, b| b.cmp(a)),
    }

    for codon_start in codon_starts {
        let codon_snps = codon_to_snps.get(&codon_start).cloned().unwrap_or_default();
        if codon_snps.is_empty() {
            continue;
        }
        let codon_end = codon_start + 2;
        if codon_end > reference.sequence.len() {
            continue;
        }
        let mut codon_snps = codon_snps;
        codon_snps.sort_by_key(|s| s.position);

        if codon_start == 0 || codon_end > reference.sequence.len() {
            continue;
        }
        let codon_seq = &reference.sequence[(codon_start - 1)..codon_end];

        let overlaps_indel = indels.iter().any(|indel| {
            let indel_end = indel.position + indel.ref_allele.len() - 1;
            indel.position <= codon_end && indel_end >= codon_start
        });

        let mut upstream_shift: isize = 0;
        let mut has_symbolic_sv = false;

        for indel in &indels {
            let is_upstream = match gene.strand {
                Strand::Plus => indel.position < codon_start,
                Strand::Minus => indel.position > codon_end,
            };

            if is_upstream {
                if indel.alt_allele.starts_with('<') {
                    has_symbolic_sv = true;
                } else {
                    upstream_shift +=
                        (indel.alt_allele.len() as isize) - (indel.ref_allele.len() as isize);
                }
            }
        }

        let is_frameshifted = has_symbolic_sv || upstream_shift % 3 != 0;

        let (eff_gene_start, eff_gene_end) = effective_bounds(gene);
        let codon_info = CodonInfo {
            codon_list: codon_snps,
            original_codon: codon_seq.to_string(),
            gene_name: gene.name.clone(),
            gene_start: eff_gene_start,
            gene_end: eff_gene_end,
            codon_start,
            codon_end,
            protein_offset: gene.protein_offset,
        };

        let mut var_info = process_codon(codon_info, gene.strand, chrom, genetic_code);

        if overlaps_indel {
            var_info.change_type = ChangeType::IndelOverlap;
            var_info.aa_changes = vec!["Unknown".to_string()];
            var_info.snp_aa_changes = vec!["Unknown".to_string(); var_info.snp_aa_changes.len()];
            var_info.aa_changes_local = vec!["Unknown".to_string()];
            var_info.snp_aa_changes_local =
                vec!["Unknown".to_string(); var_info.snp_aa_changes_local.len()];
        } else if is_frameshifted {
            var_info.change_type = var_info.change_type.with_frameshift();
            var_info.aa_changes = var_info
                .aa_changes
                .into_iter()
                .map(|s| format!("{s} (fs)"))
                .collect();
            var_info.snp_aa_changes = var_info
                .snp_aa_changes
                .into_iter()
                .map(|s| format!("{s} (fs)"))
                .collect();
            var_info.aa_changes_local = var_info
                .aa_changes_local
                .into_iter()
                .map(|s| format!("{s} (fs)"))
                .collect();
            var_info.snp_aa_changes_local = var_info
                .snp_aa_changes_local
                .into_iter()
                .map(|s| format!("{s} (fs)"))
                .collect();
        }

        variants.push(var_info);
    }

    for indel in indels {
        let is_frameshift = if indel.alt_allele.starts_with('<') {
            true
        } else {
            (indel.ref_allele.len() as isize - indel.alt_allele.len() as isize) % 3 != 0
        };

        let change_type = if is_frameshift {
            ChangeType::FrameshiftIndel
        } else {
            ChangeType::InFrameIndel
        };

        variants.push(VariantInfo {
            chrom: chrom.to_string(),
            gene: gene.name.clone(),
            positions: vec![indel.position],
            ref_bases: vec![indel.ref_allele.clone()],
            base_changes: vec![indel.alt_allele.clone()],
            aa_changes: vec!["-".to_string()],
            snp_aa_changes: vec!["-".to_string()],
            aa_changes_local: vec!["-".to_string()],
            snp_aa_changes_local: vec!["-".to_string()],
            variant_type: VariantType::Indel,
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
            ref_codon: None,
            snp_codon: None,
            mnv_codon: None,
            original_dp: indel.original_dp.map(|value| vec![value]),
            original_freq: indel.original_freq.map(|value| vec![value]),
            original_info: indel.original_info.clone(),
        });
    }

    variants
}

/// Build a `VariantInfo` for a VCF position that falls outside all annotated
/// genes (intergenic).  The entry preserves the original position, alleles and
/// any original metrics but has no codon / amino-acid annotation.
pub fn build_intergenic_variant(chrom: &str, vcf_pos: &crate::io::VcfPosition) -> VariantInfo {
    let is_snp = vcf_pos.ref_allele.len() == 1
        && vcf_pos.alt_allele.len() == 1
        && !vcf_pos.alt_allele.starts_with('<');

    let variant_type = if is_snp {
        VariantType::Snp
    } else {
        VariantType::Indel
    };

    VariantInfo {
        chrom: chrom.to_string(),
        gene: "intergenic".to_string(),
        positions: vec![vcf_pos.position],
        ref_bases: vec![vcf_pos.ref_allele.clone()],
        base_changes: vec![vcf_pos.alt_allele.clone()],
        aa_changes: vec!["-".to_string()],
        snp_aa_changes: vec!["-".to_string()],
        aa_changes_local: vec!["-".to_string()],
        snp_aa_changes_local: vec!["-".to_string()],
        variant_type,
        change_type: ChangeType::Unknown,
        ref_codon: None,
        snp_codon: None,
        mnv_codon: None,
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
        original_dp: vcf_pos.original_dp.map(|dp| vec![dp]),
        original_freq: vcf_pos.original_freq.map(|freq| vec![freq]),
        original_info: vcf_pos.original_info.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::{codon_bounds_for_position, process_codon};
    use crate::utils::reverse_complement;
    use crate::variants::{ChangeType, CodonInfo, Gene, Snp, Strand, VariantType};

    fn next_u64(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        *seed
    }

    fn random_base(seed: &mut u64) -> char {
        match next_u64(seed) % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        }
    }

    #[test]
    fn test_property_reverse_complement_involution_for_random_sequences() {
        let mut seed = 123456789u64;
        for _ in 0..200 {
            let mut seq = String::new();
            for _ in 0..30 {
                seq.push(random_base(&mut seed));
            }
            let twice = reverse_complement(&reverse_complement(&seq));
            assert_eq!(seq, twice);
        }
    }

    #[test]
    fn test_property_codon_bounds_plus_strand_cover_position() {
        let gene = Gene {
            name: "gene_plus".to_string(),
            start: 100,
            end: 399,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
        };
        for pos in gene.start..=gene.end {
            let bounds = codon_bounds_for_position(&gene, pos).expect("expected codon bounds");
            assert_eq!(bounds.1 - bounds.0 + 1, 3);
            assert!(bounds.0 <= pos && pos <= bounds.1);
            assert!(bounds.0 >= gene.start && bounds.1 <= gene.end);
        }
    }

    #[test]
    fn test_property_codon_bounds_minus_strand_cover_position() {
        let gene = Gene {
            name: "gene_minus".to_string(),
            start: 100,
            end: 399,
            strand: Strand::Minus,
            phase: 0,
            protein_offset: 0,
        };
        for pos in gene.start..=gene.end {
            let bounds = codon_bounds_for_position(&gene, pos).expect("expected codon bounds");
            assert_eq!(bounds.1 - bounds.0 + 1, 3);
            assert!(bounds.0 <= pos && pos <= bounds.1);
            assert!(bounds.0 >= gene.start && bounds.1 <= gene.end);
        }
    }

    #[test]
    fn test_change_type_from_label_roundtrip() {
        let labels = [
            "Synonymous",
            "Non-synonymous",
            "Stop gained",
            "Stop lost",
            "Unknown",
            "Indel overlap",
            "Synonymous (frameshift)",
            "Non-synonymous (frameshift)",
            "Stop gained (frameshift)",
            "Stop lost (frameshift)",
            "Unknown (frameshift)",
            "Frameshift Indel",
            "In-frame Indel",
        ];
        for label in labels {
            let ct = ChangeType::from_label(label);
            assert_eq!(ct.to_string(), label);
        }
    }

    #[test]
    fn test_change_type_from_label_unknown_fallback() {
        assert_eq!(ChangeType::from_label("garbage"), ChangeType::Unknown);
        assert_eq!(ChangeType::from_label(""), ChangeType::Unknown);
    }

    #[test]
    fn test_with_frameshift_idempotent() {
        let fs = ChangeType::FrameshiftSynonymous;
        assert_eq!(fs.with_frameshift(), ChangeType::FrameshiftSynonymous);
        let base = ChangeType::StopGained;
        assert_eq!(base.with_frameshift(), ChangeType::FrameshiftStopGained);
    }

    #[test]
    fn test_variant_type_display() {
        assert_eq!(VariantType::Snp.to_string(), "SNP");
        assert_eq!(VariantType::Mnv.to_string(), "MNV");
        assert_eq!(VariantType::SnpMnv.to_string(), "SNP/MNV");
        assert_eq!(VariantType::Indel.to_string(), "INDEL");
    }

    #[test]
    fn test_strand_from_str() {
        assert_eq!("+".parse::<Strand>(), Ok(Strand::Plus));
        assert_eq!("-".parse::<Strand>(), Ok(Strand::Minus));
        assert!("?".parse::<Strand>().is_err());
    }

    #[test]
    fn test_codon_bounds_phase_skipped_position_returns_none() {
        // Plus strand, phase=2: positions 100 and 101 are in the
        // phase-skipped region and must NOT be reported as belonging to a
        // codon (the codon they sit in spans into the previous exon).
        let gene = Gene {
            name: "cds_plus_phase2".to_string(),
            start: 100,
            end: 120,
            strand: Strand::Plus,
            phase: 2,
            protein_offset: 0,
        };
        assert!(codon_bounds_for_position(&gene, 100).is_none());
        assert!(codon_bounds_for_position(&gene, 101).is_none());
        // 102 is the first base of the first complete codon.
        assert_eq!(codon_bounds_for_position(&gene, 102), Some((102, 104)));

        // Minus strand symmetric: phase=2 means the LAST 2 bases of the
        // exon are skipped (they belong to a codon ending in the next exon
        // in transcript order). 120 and 119 → None, 118 → first complete
        // codon.
        let minus = Gene {
            name: "cds_minus_phase2".to_string(),
            start: 100,
            end: 120,
            strand: Strand::Minus,
            phase: 2,
            protein_offset: 0,
        };
        assert!(codon_bounds_for_position(&minus, 120).is_none());
        assert!(codon_bounds_for_position(&minus, 119).is_none());
        assert_eq!(codon_bounds_for_position(&minus, 118), Some((116, 118)));
    }

    #[test]
    fn test_codon_bounds_plus_strand_with_phase_1() {
        // Phase=1 means skip 1 base from gene.start: first codon is [start+1, start+3]
        let gene = Gene {
            name: "cds_plus_phase1".to_string(),
            start: 100,
            end: 120,
            strand: Strand::Plus,
            phase: 1,
            protein_offset: 0,
        };
        // Position 100 is in the skipped region, no codon
        assert!(codon_bounds_for_position(&gene, 100).is_none());
        // 101,102,103 → first codon
        assert_eq!(codon_bounds_for_position(&gene, 101), Some((101, 103)));
        assert_eq!(codon_bounds_for_position(&gene, 103), Some((101, 103)));
        // 104 → next codon
        assert_eq!(codon_bounds_for_position(&gene, 104), Some((104, 106)));
    }

    #[test]
    fn test_codon_bounds_minus_strand_phase_1_gnaq_regression() {
        // Regression for issue #12: GNAQ CDS chr9:77794463-77794592, minus strand, phase=1.
        // Without phase fix, 77794516 and 77794517 were grouped together and 77794518 left alone.
        // With phase=1 applied, 77794517 and 77794518 must share a codon, and 77794516 must be on its own.
        let gene = Gene {
            name: "GNAQ_cds".to_string(),
            start: 77_794_463,
            end: 77_794_592,
            strand: Strand::Minus,
            phase: 1,
            protein_offset: 0,
        };
        let b516 = codon_bounds_for_position(&gene, 77_794_516).expect("bounds for 516");
        let b517 = codon_bounds_for_position(&gene, 77_794_517).expect("bounds for 517");
        let b518 = codon_bounds_for_position(&gene, 77_794_518).expect("bounds for 518");
        assert_eq!(b517, b518, "517 and 518 must share a codon under phase=1");
        assert_ne!(
            b516, b517,
            "516 must NOT share a codon with 517/518 under phase=1"
        );
    }

    #[test]
    fn test_codon_bounds_outside_gene_returns_none() {
        let gene = Gene {
            name: "g".to_string(),
            start: 100,
            end: 199,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
        };
        assert!(codon_bounds_for_position(&gene, 50).is_none());
        assert!(codon_bounds_for_position(&gene, 300).is_none());
    }

    #[test]
    fn test_get_mnv_variants_for_gene_mixed_snps_and_indels() {
        use super::get_mnv_variants_for_gene;
        use crate::io::{Reference, VcfPosition};

        let gene = Gene {
            name: "testGene".to_string(),
            start: 1,
            end: 9,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
        };
        let reference = Reference {
            sequence: "ATGATGATG",
        };

        let snps = vec![
            VcfPosition {
                position: 2,
                ref_allele: "T".to_string(),
                alt_allele: "C".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            VcfPosition {
                position: 5,
                ref_allele: "TG".to_string(),
                alt_allele: "T".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ];

        let variants = get_mnv_variants_for_gene(
            &gene,
            &snps,
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );
        assert!(
            variants.len() >= 2,
            "expected SNP + indel, got {}",
            variants.len()
        );
        let has_snp = variants.iter().any(|v| v.variant_type == VariantType::Snp);
        let has_indel = variants
            .iter()
            .any(|v| v.variant_type == VariantType::Indel);
        assert!(has_snp, "missing SNP variant");
        assert!(has_indel, "missing indel variant");
    }

    #[test]
    fn test_build_intergenic_variant_snp() {
        use super::build_intergenic_variant;
        use crate::io::VcfPosition;

        let pos = VcfPosition {
            position: 42,
            ref_allele: "A".to_string(),
            alt_allele: "G".to_string(),
            original_dp: Some(30),
            original_freq: Some(0.8),
            original_info: None,
        };
        let v = build_intergenic_variant("chrX", &pos);
        assert_eq!(v.gene, "intergenic");
        assert_eq!(v.variant_type, VariantType::Snp);
        assert_eq!(v.positions, vec![42]);
        assert_eq!(v.original_dp, Some(vec![30]));
    }

    #[test]
    fn test_build_intergenic_variant_indel() {
        use super::build_intergenic_variant;
        use crate::io::VcfPosition;

        let pos = VcfPosition {
            position: 10,
            ref_allele: "AT".to_string(),
            alt_allele: "A".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        };
        let v = build_intergenic_variant("chr1", &pos);
        assert_eq!(v.variant_type, VariantType::Indel);
    }

    #[test]
    fn test_process_codon_mnv_two_snps_in_codon() {
        let codon_info = CodonInfo {
            codon_list: vec![
                Snp {
                    index: 100,
                    position: 100,
                    ref_base: "A".to_string(),
                    base: "T".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
                Snp {
                    index: 101,
                    position: 101,
                    ref_base: "T".to_string(),
                    base: "C".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
            ],
            original_codon: "ATG".to_string(),
            gene_name: "gene1".to_string(),
            gene_start: 100,
            gene_end: 399,
            codon_start: 100,
            codon_end: 102,
            protein_offset: 0,
        };
        let result = process_codon(
            codon_info,
            Strand::Plus,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );
        assert_eq!(result.variant_type, VariantType::SnpMnv);
        assert_eq!(result.positions.len(), 2);
        assert!(result.mnv_codon.is_some());
    }

    #[test]
    fn test_process_codon_emits_expected_variant_type_for_single_snp() {
        let codon_info = CodonInfo {
            codon_list: vec![Snp {
                index: 101,
                position: 101,
                ref_base: "T".to_string(),
                base: "C".to_string(),
                original_dp: Some(20),
                original_freq: Some(0.5),
                original_info: None,
            }],
            original_codon: "ATG".to_string(),
            gene_name: "gene1".to_string(),
            gene_start: 100,
            gene_end: 399,
            codon_start: 100,
            codon_end: 102,
            protein_offset: 0,
        };

        let result = process_codon(
            codon_info,
            Strand::Plus,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );
        assert_eq!(result.variant_type, VariantType::Snp);
        assert!(matches!(
            result.change_type,
            ChangeType::Synonymous
                | ChangeType::NonSynonymous
                | ChangeType::StopGained
                | ChangeType::StopLost
                | ChangeType::Unknown
        ));
    }

    #[test]
    fn test_construct_codon_uppercase_alt() {
        // ALT alleles should be uppercased to match the codon lookup table.
        let codon_info = CodonInfo {
            codon_list: vec![Snp {
                index: 101,
                position: 101,
                ref_base: "A".to_string(),
                base: "t".to_string(), // lowercase ALT
                original_dp: None,
                original_freq: None,
                original_info: None,
            }],
            original_codon: "ATG".to_string(),
            gene_name: "test".to_string(),
            gene_start: 100,
            gene_end: 120,
            codon_start: 100,
            codon_end: 102,
            protein_offset: 0,
        };
        let result = process_codon(
            codon_info,
            Strand::Plus,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );
        // Position 101 is the 2nd base (index 1). ATG → ATG with pos 101 T→t
        // The codon becomes "ATt" → uppercased to "ATT" → Ile (I), not X.
        assert_ne!(
            result.aa_changes[0], "X",
            "Lowercase ALT should not produce unknown amino acid"
        );
    }

    #[test]
    fn test_duplicate_position_dedup() {
        // After --split-multiallelic, two ALTs at the same position should
        // not both appear in the same codon group.
        use crate::io::{Reference, VcfPosition};
        use crate::variants::get_mnv_variants_for_gene;

        let gene = Gene {
            name: "geneA".to_string(),
            start: 100,
            end: 111,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
        };
        // Two VCF records at position 101 with different ALTs
        let snps = vec![
            VcfPosition {
                position: 101,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            VcfPosition {
                position: 101,
                ref_allele: "A".to_string(),
                alt_allele: "G".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ];
        let seq = "N".repeat(99) + "ATGATGATGATG";
        let reference = Reference { sequence: &seq };
        let variants = get_mnv_variants_for_gene(
            &gene,
            &snps,
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );
        // Should produce exactly 1 variant (first ALT wins, duplicate skipped)
        assert_eq!(
            variants.len(),
            1,
            "Duplicate position should be deduplicated"
        );
        assert_eq!(variants[0].positions.len(), 1);
        assert_eq!(variants[0].base_changes[0], "T"); // first ALT kept
    }
}
