//! Codon-level variant processing: SNP grouping, MNV detection, and amino
//! acid change calculation.
//!
//! The heavy lifting lives in focused submodules; this file wires them together
//! and hosts the public per-gene / per-transcript entry points.

mod allele_apply;
mod config;
mod gene_path;
mod grouping;
mod indel_effect;
mod nmd;
mod phased;
mod protein;
mod transcript_model;
mod transcript_path;

#[cfg(test)]
mod tests;

pub use config::{FrameshiftPhasing, IndelAnnotationConfig, PairLinkage};
pub use gene_path::process_codon;
pub use phased::{local_haplotype_components, phased_haplotype_variant};

use crate::variants::types::*;
use std::collections::BTreeMap;

use config::indel_passes_frameshift_gate;
use gene_path::{codon_bounds_for_position, effective_bounds};

// Re-exported so the property tests can check the codon-bounds invariants over
// generated gene shapes rather than one hard-coded feature.
#[cfg(test)]
pub(crate) use gene_path::codon_bounds_for_position as codon_bounds;
use grouping::{merge_snp_into_groups, variant_event_metadata};
use indel_effect::{coding_delta_for_variant, protein_effect_for_indel};
use transcript_model::{
    first_touched_transcript_offset, has_transcript_cds_model, transcript_offset_for_position,
    transcript_sequence_for_gene, variant_overlaps_coding_model,
    variant_touched_transcript_offsets,
};
use transcript_path::{
    apply_frameshift_labeling, frameshift_ptc_protein_pos, merge_transcript_snp_into_groups,
    process_transcript_codon, transcript_oriented_base, TranscriptSnp,
};

fn get_mnv_variants_for_transcript(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
    config: &IndelAnnotationConfig,
    phasing: &FrameshiftPhasing,
) -> Vec<VariantInfo> {
    let Some(ref_cds) = transcript_sequence_for_gene(gene, reference) else {
        return Vec::new();
    };

    let mut variants = Vec::new();
    // See merge_snp_into_groups for the multi-allelic alternative interpretation
    // logic. Each codon_start maps to a list of mutually-exclusive groups.
    let mut codon_to_groups: BTreeMap<usize, Vec<Vec<TranscriptSnp>>> = BTreeMap::new();
    let mut indels: Vec<crate::io::VcfPosition> = Vec::new();

    for variant in snp_list
        .iter()
        .filter(|variant| variant_overlaps_coding_model(gene, variant))
    {
        let substitutions = variant.substitution_components();
        if !substitutions.is_empty() {
            for component in substitutions {
                let Some(offset) = transcript_offset_for_position(gene, component.position) else {
                    continue;
                };
                let Some(transcript_alt_base) =
                    transcript_oriented_base(&component.alt_allele, gene.strand)
                else {
                    continue;
                };
                let codon_start = (offset / 3) * 3;
                let new_snp = TranscriptSnp {
                    snp: Snp {
                        index: component.position,
                        position: component.position,
                        ref_base: component.ref_allele,
                        base: component.alt_allele,
                        original_dp: variant.original_dp,
                        original_freq: variant.original_freq,
                        original_info: variant.original_info.clone(),
                    },
                    transcript_offset: offset,
                    transcript_alt_base,
                };
                let groups = codon_to_groups.entry(codon_start).or_default();
                merge_transcript_snp_into_groups(groups, new_snp);
            }
        } else {
            indels.push(variant.clone());
        }
    }

    let ptc_protein_pos =
        frameshift_ptc_protein_pos(gene, reference, &indels, genetic_code, config);

    let mut codon_starts = codon_to_groups.keys().copied().collect::<Vec<_>>();
    codon_starts.sort_unstable();

    for codon_start in codon_starts {
        let groups = codon_to_groups
            .get(&codon_start)
            .cloned()
            .unwrap_or_default();
        for codon_snps in groups {
            let codon_end = codon_start + 3;
            if codon_end > ref_cds.len() {
                continue;
            }
            let ref_codon = &ref_cds[codon_start..codon_end];
            let mut codon_snps = codon_snps;
            codon_snps.sort_by_key(|s| s.transcript_offset);

            let overlaps_indel = indels.iter().any(|indel| {
                variant_touched_transcript_offsets(gene, indel)
                    .into_iter()
                    .any(|offset| offset >= codon_start && offset < codon_end)
            });

            let codon_positions: Vec<usize> =
                codon_snps.iter().map(|snp| snp.snp.position).collect();
            let mut upstream_shift: isize = 0;
            let mut has_symbolic_sv = false;
            let mut frameshift_linkage = Vec::new();
            for indel in &indels {
                let Some(first_offset) = first_touched_transcript_offset(gene, indel) else {
                    continue;
                };
                if first_offset < codon_start {
                    // Only let sufficiently frequent upstream indels shift the frame
                    // of downstream codons (see IndelAnnotationConfig).
                    if !indel_passes_frameshift_gate(indel, config) {
                        continue;
                    }
                    // Suppress propagation when the BAM shows this indel is in trans
                    // with the codon's SNV: the codon's molecule does not carry it.
                    // Record what the reads said either way, so the label the row
                    // ends up with can be traced back to the evidence for it.
                    if let Some(linkage) = phasing.linkage_for(indel.position, &codon_positions) {
                        frameshift_linkage.push(linkage);
                    }
                    if phasing.indel_in_trans_with(indel.position, &codon_positions) {
                        continue;
                    }
                    if indel.alt_allele.starts_with('<') {
                        has_symbolic_sv = true;
                    } else if let Some(delta) = coding_delta_for_variant(gene, reference, indel) {
                        upstream_shift += delta;
                    } else {
                        has_symbolic_sv = true;
                    }
                }
            }
            let is_frameshifted = has_symbolic_sv || upstream_shift % 3 != 0;

            let mut var_info = process_transcript_codon(
                gene,
                chrom,
                ref_codon,
                codon_start,
                &codon_snps,
                genetic_code,
            );

            if overlaps_indel {
                var_info.change_type = ChangeType::IndelOverlap;
                var_info.aa_changes = vec!["Unknown".to_string()];
                var_info.snp_aa_changes =
                    vec!["Unknown".to_string(); var_info.snp_aa_changes.len()];
                var_info.aa_changes_local = vec!["Unknown".to_string()];
                var_info.snp_aa_changes_local =
                    vec!["Unknown".to_string(); var_info.snp_aa_changes_local.len()];
            } else if is_frameshifted {
                apply_frameshift_labeling(&mut var_info, ptc_protein_pos);
            }

            // A nonsense SNV/MNV (still classified Stop gained after any
            // frameshift/overlap relabelling) begins its premature stop at this
            // codon's CDS offset.
            if var_info.change_type == ChangeType::StopGained {
                var_info.annotations.nmd = nmd::snv_mnv_nmd_prediction(gene, codon_start);
            }
            var_info.annotations.frameshift_linkage = frameshift_linkage;

            variants.push(var_info);
        }
    }

    for indel in indels {
        let (change_type, aa_changes, aa_changes_local, ref_codon, alt_codon) =
            protein_effect_for_indel(gene, reference, &indel, genetic_code);
        let (event_class, event_components) = variant_event_metadata(&indel);
        let nmd = nmd::indel_nmd_prediction(gene, reference, &indel, genetic_code);

        variants.push(VariantInfo {
            chrom: chrom.to_string(),
            gene: gene.name.clone(),
            positions: vec![indel.position],
            ref_bases: vec![indel.ref_allele.clone()],
            base_changes: vec![indel.alt_allele.clone()],
            aa_changes,
            snp_aa_changes: vec!["-".to_string()],
            aa_changes_local,
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
            mnv_phasing_reads: None,
            ref_codon,
            snp_codon: None,
            mnv_codon: alt_codon,
            original_dp: indel.original_dp.map(|value| vec![value]),
            original_freq: indel.original_freq.map(|value| vec![value]),
            original_info: indel.original_info.clone(),
            event_class: Some(event_class),
            event_components,
            annotations: crate::variants::VariantAnnotations {
                nmd,
                ..Default::default()
            },
        });
    }

    variants
}

pub fn get_mnv_variants_for_gene(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
) -> Vec<VariantInfo> {
    get_mnv_variants_for_gene_with_config(
        gene,
        snp_list,
        reference,
        chrom,
        genetic_code,
        &IndelAnnotationConfig::default(),
        &FrameshiftPhasing::default(),
    )
}

/// Like [`get_mnv_variants_for_gene`] but with explicit indel-annotation
/// configuration (e.g. the upstream-indel frequency gate for frameshift
/// propagation) and BAM-derived phasing evidence that suppresses propagation
/// from indels shown to be in trans with a codon. The no-config wrapper above
/// preserves historical behaviour.
pub fn get_mnv_variants_for_gene_with_config(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
    config: &IndelAnnotationConfig,
    phasing: &FrameshiftPhasing,
) -> Vec<VariantInfo> {
    if has_transcript_cds_model(gene) {
        return get_mnv_variants_for_transcript(
            gene,
            snp_list,
            reference,
            chrom,
            genetic_code,
            config,
            phasing,
        );
    }

    let mut variants = Vec::new();

    // Each codon_start maps to a list of mutually-exclusive interpretations.
    // After --split-multiallelic, two ALT alleles at the same position
    // produce alternative interpretations of the same codon: each is
    // processed independently to emit one row per multi-allelic combination.
    let mut codon_to_groups: BTreeMap<usize, Vec<Vec<crate::variants::Snp>>> = BTreeMap::new();
    let mut indels: Vec<crate::io::VcfPosition> = Vec::new();

    for variant in snp_list
        .iter()
        .filter(|variant| variant.overlaps_interval(gene.start, gene.end))
    {
        let substitutions = variant.substitution_components();
        if !substitutions.is_empty() {
            for component in substitutions {
                if let Some((codon_start, _)) = codon_bounds_for_position(gene, component.position)
                {
                    let snp = crate::variants::Snp {
                        index: component.position,
                        position: component.position,
                        ref_base: component.ref_allele,
                        base: component.alt_allele,
                        original_dp: variant.original_dp,
                        original_freq: variant.original_freq,
                        original_info: variant.original_info.clone(),
                    };
                    let groups = codon_to_groups.entry(codon_start).or_default();
                    merge_snp_into_groups(groups, snp);
                }
            }
        } else {
            indels.push(variant.clone());
        }
    }

    if codon_to_groups.is_empty() && indels.is_empty() {
        return variants;
    }

    let ptc_protein_pos =
        frameshift_ptc_protein_pos(gene, reference, &indels, genetic_code, config);

    let mut codon_starts: Vec<usize> = codon_to_groups.keys().copied().collect();
    match gene.strand {
        Strand::Plus => codon_starts.sort_unstable(),
        Strand::Minus => codon_starts.sort_unstable_by(|a, b| b.cmp(a)),
    }

    for codon_start in codon_starts {
        let groups = codon_to_groups
            .get(&codon_start)
            .cloned()
            .unwrap_or_default();
        for codon_snps in groups {
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

            let overlaps_indel = indels
                .iter()
                .any(|indel| indel.overlaps_interval(codon_start, codon_end));

            let codon_positions: Vec<usize> = codon_snps.iter().map(|s| s.position).collect();
            let mut upstream_shift: isize = 0;
            let mut has_symbolic_sv = false;
            let mut frameshift_linkage = Vec::new();

            for indel in &indels {
                let is_upstream = match gene.strand {
                    Strand::Plus => indel.position < codon_start,
                    Strand::Minus => indel.position > codon_end,
                };

                if is_upstream {
                    // Only sufficiently frequent upstream indels shift the frame of
                    // downstream codons (see IndelAnnotationConfig).
                    if !indel_passes_frameshift_gate(indel, config) {
                        continue;
                    }
                    // Suppress propagation when the BAM shows this indel is in trans
                    // with the codon's SNV: the codon's molecule does not carry it.
                    // Record what the reads said either way, so the label the row
                    // ends up with can be traced back to the evidence for it.
                    if let Some(linkage) = phasing.linkage_for(indel.position, &codon_positions) {
                        frameshift_linkage.push(linkage);
                    }
                    if phasing.indel_in_trans_with(indel.position, &codon_positions) {
                        continue;
                    }
                    if indel.alt_allele.starts_with('<') {
                        has_symbolic_sv = true;
                    } else if let Some(delta) = coding_delta_for_variant(gene, reference, indel) {
                        // CDS-clipped coding-length change, consistent with the
                        // transcript model. For indels that span the gene/CDS
                        // boundary this counts only the in-CDS portion, unlike the
                        // raw allele-length difference used previously.
                        upstream_shift += delta;
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
                var_info.snp_aa_changes =
                    vec!["Unknown".to_string(); var_info.snp_aa_changes.len()];
                var_info.aa_changes_local = vec!["Unknown".to_string()];
                var_info.snp_aa_changes_local =
                    vec!["Unknown".to_string(); var_info.snp_aa_changes_local.len()];
            } else if is_frameshifted {
                apply_frameshift_labeling(&mut var_info, ptc_protein_pos);
            }
            var_info.annotations.frameshift_linkage = frameshift_linkage;

            variants.push(var_info);
        }
    }

    for indel in indels {
        let (change_type, aa_changes, aa_changes_local, ref_codon, alt_codon) =
            protein_effect_for_indel(gene, reference, &indel, genetic_code);
        let (event_class, event_components) = variant_event_metadata(&indel);

        variants.push(VariantInfo {
            chrom: chrom.to_string(),
            gene: gene.name.clone(),
            positions: vec![indel.position],
            ref_bases: vec![indel.ref_allele.clone()],
            base_changes: vec![indel.alt_allele.clone()],
            aa_changes,
            snp_aa_changes: vec!["-".to_string()],
            aa_changes_local,
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
            mnv_phasing_reads: None,
            ref_codon,
            snp_codon: None,
            mnv_codon: alt_codon,
            original_dp: indel.original_dp.map(|value| vec![value]),
            original_freq: indel.original_freq.map(|value| vec![value]),
            original_info: indel.original_info.clone(),
            event_class: Some(event_class),
            event_components,
            annotations: crate::variants::VariantAnnotations::default(),
        });
    }

    variants
}

/// Build a `VariantInfo` for a variant in the splice region of a gene's intron.
/// It carries the gene name and splice consequence; everything else mirrors an
/// intergenic entry (no codon annotation), so its read support is counted
/// alongside intergenic SNPs and it renders through the normal genic path.
pub fn build_splice_variant(
    chrom: &str,
    vcf_pos: &crate::io::VcfPosition,
    gene_name: &str,
    splice: crate::variants::SpliceConsequence,
) -> VariantInfo {
    let mut variant = build_intergenic_variant(chrom, vcf_pos);
    variant.gene = gene_name.to_string();
    variant.annotations.splice = Some(splice);
    variant
}

/// Build a `VariantInfo` for a variant inside a real intron of a gene, away from
/// its splice sites. It carries the gene name so the variant is reported against
/// the transcript it actually sits in rather than as intergenic, which is what a
/// position inside a gene is not.
pub fn build_intron_variant(
    chrom: &str,
    vcf_pos: &crate::io::VcfPosition,
    gene_name: &str,
) -> VariantInfo {
    build_non_coding_variant(
        chrom,
        vcf_pos,
        gene_name,
        crate::variants::NonCodingReason::Intron,
    )
}

/// Build a `VariantInfo` for a variant inside a gene that carries no codon: an
/// intron, or a transcript that is never translated. It keeps the gene name and
/// the reason, and has no amino-acid annotation because there is none to make.
pub fn build_non_coding_variant(
    chrom: &str,
    vcf_pos: &crate::io::VcfPosition,
    gene_name: &str,
    reason: crate::variants::NonCodingReason,
) -> VariantInfo {
    let mut variant = build_intergenic_variant(chrom, vcf_pos);
    variant.gene = gene_name.to_string();
    variant.annotations.non_coding = Some(reason);
    variant
}

/// Build a `VariantInfo` for a VCF position that falls outside all annotated
/// genes (intergenic).  The entry preserves the original position, alleles and
/// any original metrics but has no codon / amino-acid annotation.
pub fn build_intergenic_variant(chrom: &str, vcf_pos: &crate::io::VcfPosition) -> VariantInfo {
    let event = vcf_pos.event();
    let variant_type = match event.class {
        crate::variants::AlleleEventClass::Snp => VariantType::Snp,
        crate::variants::AlleleEventClass::Mnv => VariantType::Mnv,
        _ => VariantType::Indel,
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
        mnv_phasing_reads: None,
        original_dp: vcf_pos.original_dp.map(|dp| vec![dp]),
        original_freq: vcf_pos.original_freq.map(|freq| vec![freq]),
        original_info: vcf_pos.original_info.clone(),
        event_class: Some(event.class.as_str().to_string()),
        event_components: event.component_labels(),
        annotations: crate::variants::VariantAnnotations::default(),
    }
}
