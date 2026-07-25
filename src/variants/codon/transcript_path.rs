//! Spliced-transcript codon path: transcript SNPs, codon processing, frameshift/PTC.

use crate::utils::{determine_change_type, iupac_aa};
use crate::variants::{ChangeType, Gene, Snp, Strand, VariantInfo, VariantType};
use std::collections::BTreeSet;

use super::allele_apply::apply_allele_to_feature;
use super::config::{indel_passes_frameshift_gate, IndelAnnotationConfig};
use super::grouping::{collect_all_f64, collect_all_usize, merge_original_info};
use super::indel_effect::coding_delta_for_variant;
use super::protein::{complement_base, translate_cds};
use super::transcript_model::coding_sequence_for_gene;

#[derive(Debug, Clone)]
pub(super) struct TranscriptSnp {
    pub(super) snp: Snp,
    pub(super) transcript_offset: usize,
    pub(super) transcript_alt_base: char,
}

pub(super) fn transcript_oriented_base(base: &str, strand: Strand) -> Option<char> {
    let base = base.chars().next()?;
    Some(match strand {
        Strand::Plus => base.to_ascii_uppercase(),
        Strand::Minus => complement_base(base),
    })
}

pub(super) fn construct_transcript_codon(
    ref_codon: &str,
    target_snps: &[&TranscriptSnp],
) -> String {
    let mut codon = ref_codon.chars().collect::<Vec<_>>();
    for snp in target_snps {
        let idx = snp.transcript_offset % 3;
        if let Some(slot) = codon.get_mut(idx) {
            *slot = snp.transcript_alt_base.to_ascii_uppercase();
        }
    }
    codon.into_iter().collect()
}

pub(super) fn process_transcript_codon(
    gene: &Gene,
    chrom: &str,
    ref_codon: &str,
    codon_start_offset: usize,
    codon_snps: &[TranscriptSnp],
    genetic_code: crate::genetic_code::GeneticCode,
) -> VariantInfo {
    let mnv_snps = codon_snps.iter().collect::<Vec<_>>();
    let mnv_codon = construct_transcript_codon(ref_codon, &mnv_snps);
    let snp_codons = codon_snps
        .iter()
        .map(|snp| construct_transcript_codon(ref_codon, &[snp]))
        .collect::<Vec<_>>();
    let snp_codon = snp_codons.join(", ");

    let orig_aa = genetic_code.translate_seq(ref_codon.as_bytes());
    let mut_aa = genetic_code.translate_seq(mnv_codon.as_bytes());
    let aa_pos = (codon_start_offset / 3) + 1;
    // Classify from the one-letter form (e.g. "L2L", "M1T"); `determine_change_type`
    // inspects the first/last character, so it must NOT receive the three-letter
    // IUPAC string (which would turn "Leu2Leu" into a spurious Non-synonymous and
    // hide synonymous / start-lost calls). The IUPAC form is kept only for display.
    let combined_change = format!("{orig_aa}{aa_pos}{mut_aa}");
    let combined_aa = iupac_aa(&combined_change);
    let change_type = ChangeType::from_label(&determine_change_type(&combined_change));

    let single_aas: Vec<char> = codon_snps
        .iter()
        .map(|snp| {
            let single_codon = construct_transcript_codon(ref_codon, &[snp]);
            genetic_code
                .translate_seq(single_codon.as_bytes())
                .chars()
                .next()
                .unwrap_or('X')
        })
        .collect();
    let snp_changes = single_aas
        .iter()
        .map(|aa| iupac_aa(&format!("{orig_aa}{aa_pos}{aa}")))
        .collect::<Vec<_>>();
    let mut annotations =
        super::gene_path::compute_annotations(&orig_aa, &mut_aa, &single_aas, change_type);
    let snp_refs: Vec<&Snp> = codon_snps.iter().map(|ts| &ts.snp).collect();
    annotations.dbs_class = super::gene_path::dbs_class_for_codon(&snp_refs);
    let hgvs_entries: Vec<(usize, char, char)> = codon_snps
        .iter()
        .filter_map(|ts| {
            let cref = transcript_oriented_base(&ts.snp.ref_base, gene.strand)?;
            Some((ts.transcript_offset + 1, cref, ts.transcript_alt_base))
        })
        .collect();
    annotations.hgvs_c = crate::variants::hgvs::coding_substitution(&hgvs_entries);

    let raw_snps = codon_snps
        .iter()
        .map(|snp| snp.snp.clone())
        .collect::<Vec<_>>();

    VariantInfo {
        chrom: chrom.to_string(),
        gene: gene.name.clone(),
        positions: codon_snps.iter().map(|s| s.snp.position).collect(),
        ref_bases: codon_snps.iter().map(|s| s.snp.ref_base.clone()).collect(),
        base_changes: codon_snps.iter().map(|s| s.snp.base.clone()).collect(),
        aa_changes: vec![combined_aa.clone()],
        snp_aa_changes: snp_changes.clone(),
        aa_changes_local: vec![combined_aa],
        snp_aa_changes_local: snp_changes,
        variant_type: if codon_snps.len() == 1 {
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
        ref_codon: Some(ref_codon.to_string()),
        snp_codon: Some(snp_codon),
        mnv_codon: Some(mnv_codon),
        original_dp: collect_all_usize(raw_snps.iter().map(|s| s.original_dp)),
        original_freq: collect_all_f64(raw_snps.iter().map(|s| s.original_freq)),
        original_info: merge_original_info(&raw_snps),
        event_class: Some(if codon_snps.len() == 1 {
            "snp".to_string()
        } else {
            "mnv".to_string()
        }),
        event_components: codon_snps
            .iter()
            .map(|s| format!("SNV:{}:{}>{}", s.snp.position, s.snp.ref_base, s.snp.base))
            .collect(),
        annotations,
    }
}

/// Protein position (1-based) of the first premature stop introduced by a single
/// upstream frameshift indel, when it truncates the protein earlier than the
/// natural stop. Returns `None` for the multi-indel case or when the frameshift
/// does not introduce an earlier stop (frameshift propagation is then unchanged).
pub(super) fn frameshift_ptc_protein_pos(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    indels: &[crate::io::VcfPosition],
    genetic_code: crate::genetic_code::GeneticCode,
    config: &IndelAnnotationConfig,
) -> Option<usize> {
    let fs_indels: Vec<&crate::io::VcfPosition> = indels
        .iter()
        .filter(|indel| !indel.alt_allele.starts_with('<'))
        .filter(|indel| indel_passes_frameshift_gate(indel, config))
        .filter(|indel| {
            coding_delta_for_variant(gene, reference, indel)
                .map(|delta| delta % 3 != 0)
                .unwrap_or(false)
        })
        .collect();
    if fs_indels.len() != 1 {
        return None;
    }
    let ref_cds = coding_sequence_for_gene(gene, reference)?;
    let alt_cds = apply_allele_to_feature(gene, reference, fs_indels[0])?;
    let alt_protein = translate_cds(&alt_cds, genetic_code);
    let ref_protein = translate_cds(&ref_cds, genetic_code);
    let alt_stop = alt_protein.find('*')?;
    let ref_stop = ref_protein.find('*').unwrap_or(ref_protein.len());
    (alt_stop < ref_stop).then_some(alt_stop + 1)
}

/// Annotate a frameshifted downstream codon. If it lies past the premature stop
/// (`ptc_protein_pos`) it is reported as untranslated ("downstream of premature
/// stop"); otherwise the previous "(fs)" amino-acid annotation is kept.
pub(super) fn apply_frameshift_labeling(
    var_info: &mut VariantInfo,
    ptc_protein_pos: Option<usize>,
) {
    let codon_pos = var_info.aa_changes.first().and_then(|aa| {
        aa.chars()
            .filter(char::is_ascii_digit)
            .collect::<String>()
            .parse::<usize>()
            .ok()
    });
    let past_ptc = matches!((ptc_protein_pos, codon_pos), (Some(ptc), Some(pos)) if pos > ptc);
    if past_ptc {
        const LABEL: &str = "downstream of premature stop";
        var_info.change_type = ChangeType::FrameshiftDownstreamOfStop;
        var_info.aa_changes = vec![LABEL.to_string()];
        var_info.snp_aa_changes = vec![LABEL.to_string(); var_info.snp_aa_changes.len()];
        var_info.aa_changes_local = vec![LABEL.to_string()];
        var_info.snp_aa_changes_local =
            vec![LABEL.to_string(); var_info.snp_aa_changes_local.len()];
        return;
    }
    var_info.change_type = var_info.change_type.with_frameshift();
    var_info.aa_changes = var_info
        .aa_changes
        .iter()
        .map(|s| format!("{s} (fs)"))
        .collect();
    var_info.snp_aa_changes = var_info
        .snp_aa_changes
        .iter()
        .map(|s| format!("{s} (fs)"))
        .collect();
    var_info.aa_changes_local = var_info
        .aa_changes_local
        .iter()
        .map(|s| format!("{s} (fs)"))
        .collect();
    var_info.snp_aa_changes_local = var_info
        .snp_aa_changes_local
        .iter()
        .map(|s| format!("{s} (fs)"))
        .collect();
}
/// Transcript-coordinates variant of merge_snp_into_groups.
pub(super) fn merge_transcript_snp_into_groups(
    groups: &mut Vec<Vec<TranscriptSnp>>,
    snp: TranscriptSnp,
) {
    if groups.is_empty() {
        groups.push(vec![snp]);
        return;
    }
    let mut new_groups: Vec<Vec<TranscriptSnp>> = Vec::with_capacity(groups.len() * 2);
    let mut seen: BTreeSet<Vec<(usize, String)>> = BTreeSet::new();

    let push_if_new = |g: Vec<TranscriptSnp>,
                       new_groups: &mut Vec<Vec<TranscriptSnp>>,
                       seen: &mut BTreeSet<Vec<(usize, String)>>| {
        let mut key: Vec<(usize, String)> = g
            .iter()
            .map(|s| (s.transcript_offset, s.transcript_alt_base.to_string()))
            .collect();
        key.sort();
        if seen.insert(key) {
            new_groups.push(g);
        }
    };

    for group in groups.iter() {
        if let Some(idx) = group
            .iter()
            .position(|s| s.transcript_offset == snp.transcript_offset)
        {
            if group[idx].transcript_alt_base == snp.transcript_alt_base {
                push_if_new(group.clone(), &mut new_groups, &mut seen);
            } else {
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
