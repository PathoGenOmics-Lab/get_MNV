//! Locally-phased indel / complex-haplotype reconstruction.

use crate::io::VcfPosition;
use crate::variants::event::{AlleleComponentKind, AlleleEventClass};
use crate::variants::{Gene, VariantInfo, VariantType};
use std::collections::{BTreeMap, BTreeSet};

use super::gene_path::{codon_bounds_for_position, effective_bounds};
use super::indel_effect::protein_effect_for_indel;
use super::transcript_model::{has_transcript_cds_model, transcript_offset_for_position};

const LOCAL_HAPLOTYPE_JOIN_DISTANCE: usize = 3;

pub(super) fn phased_indel_window(gene: &Gene, variant: &VcfPosition) -> Option<(usize, usize)> {
    let event = variant.event();
    let (eff_start, eff_end) = effective_bounds(gene);
    if eff_start == 0 || eff_end == 0 || eff_start > eff_end {
        return None;
    }

    let mut windows = Vec::new();
    for pos in [
        variant.record_start,
        event.affected_start,
        event.affected_end,
    ] {
        if pos >= eff_start && pos <= eff_end {
            if let Some(bounds) = codon_bounds_for_position(gene, pos) {
                windows.push(bounds);
            }
        }
    }

    windows.push((
        event.affected_start.max(eff_start),
        event.affected_end.min(eff_end),
    ));

    if windows.len() == 1 {
        let start = event.affected_start.saturating_sub(2).max(eff_start);
        let end = event.affected_end.saturating_add(2).min(eff_end);
        windows.push((start, end));
    }

    let start = windows.iter().map(|(start, _)| *start).min()?;
    let end = windows.iter().map(|(_, end)| *end).max()?;
    Some((start, end))
}

pub(super) fn haplotype_windows_linked(a: (usize, usize), b: (usize, usize)) -> bool {
    let a_end = a.1.saturating_add(LOCAL_HAPLOTYPE_JOIN_DISTANCE);
    let b_end = b.1.saturating_add(LOCAL_HAPLOTYPE_JOIN_DISTANCE);
    a.0 <= b_end && b.0 <= a_end
}

pub(super) fn variant_has_indel_component(variant: &VcfPosition) -> bool {
    variant.event().class.has_indel_component()
}

pub(super) fn group_component_flags(group: &[&VcfPosition]) -> (bool, bool, usize) {
    let mut has_indel = false;
    let mut has_substitution = false;
    let mut indel_components = 0usize;

    for variant in group {
        for component in variant.event().components {
            match component.kind {
                AlleleComponentKind::Snp => has_substitution = true,
                AlleleComponentKind::Insertion
                | AlleleComponentKind::Deletion
                | AlleleComponentKind::Delins
                | AlleleComponentKind::Symbolic => {
                    has_indel = true;
                    indel_components += 1;
                }
            }
        }
    }

    (has_indel, has_substitution, indel_components)
}

/// The haplotype a set of variants describes when one molecule carries all of
/// them. `None` unless the set is a genuine local haplotype: at least two
/// variants, at least one of them an indel (pure substitutions in a codon are
/// already reported as MNVs), and an allele the reference can express.
pub fn phased_haplotype_variant(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    chrom: &str,
    group: &[&VcfPosition],
    genetic_code: crate::genetic_code::GeneticCode,
) -> Option<VariantInfo> {
    let (has_indel, _, _) = group_component_flags(group);
    if !has_indel || group.len() < 2 {
        return None;
    }
    phased_variant_from_group(gene, reference, chrom, group, genetic_code)
}

pub(super) fn compound_allele_from_variants(
    reference: &crate::io::Reference<'_>,
    variants: &[&VcfPosition],
) -> Option<VcfPosition> {
    if variants.len() < 2 {
        return None;
    }

    let mut has_indel = false;
    let mut start = usize::MAX;
    let mut end = 0usize;
    let mut snps: BTreeMap<usize, String> = BTreeMap::new();
    let mut insertions: BTreeMap<usize, String> = BTreeMap::new();
    let mut deleted: BTreeSet<usize> = BTreeSet::new();

    for variant in variants {
        let event = variant.event();
        has_indel |= event.class.has_indel_component();
        start = start.min(event.affected_start).min(variant.record_start);
        end = end
            .max(event.affected_end)
            .max(variant.record_start + variant.ref_allele.len().saturating_sub(1));

        for component in event.components {
            match component.kind {
                AlleleComponentKind::Snp => {
                    if component.ref_allele.len() != 1 || component.alt_allele.len() != 1 {
                        return None;
                    }
                    match snps.get(&component.position) {
                        Some(existing) if !existing.eq_ignore_ascii_case(&component.alt_allele) => {
                            return None;
                        }
                        Some(_) => {}
                        None => {
                            snps.insert(component.position, component.alt_allele);
                        }
                    }
                }
                AlleleComponentKind::Insertion => {
                    if let Some(existing) = insertions.get(&component.position) {
                        if !existing.eq_ignore_ascii_case(&component.alt_allele) {
                            return None;
                        }
                    } else {
                        insertions.insert(component.position, component.alt_allele);
                    }
                }
                AlleleComponentKind::Deletion => {
                    let del_end = component.position + component.ref_allele.len().saturating_sub(1);
                    for pos in component.position..=del_end {
                        deleted.insert(pos);
                    }
                }
                AlleleComponentKind::Delins | AlleleComponentKind::Symbolic => return None,
            }
        }
    }

    if !has_indel || start == usize::MAX || end == 0 {
        return None;
    }
    if start == 0 || end > reference.sequence.len() || start > end {
        return None;
    }
    if snps.keys().any(|pos| deleted.contains(pos)) {
        return None;
    }
    if insertions.keys().any(|pos| deleted.contains(pos)) {
        return None;
    }

    let ref_allele = reference.sequence[(start - 1)..end].to_string();
    let mut alt_allele = String::new();
    for pos in start..=end {
        if !deleted.contains(&pos) {
            if let Some(alt_base) = snps.get(&pos) {
                alt_allele.push_str(alt_base);
            } else {
                alt_allele.push(reference.sequence.as_bytes()[pos - 1] as char);
            }
        }
        if let Some(inserted) = insertions.get(&pos) {
            alt_allele.push_str(inserted);
        }
    }

    if ref_allele.eq_ignore_ascii_case(&alt_allele) {
        return None;
    }

    Some(VcfPosition {
        record_start: start,
        ref_allele,
        alt_allele,
        original_dp: None,
        original_freq: None,
        original_info: None,
        // A synthesised compound allele; its phase came from the reads, not
        // from anything the input declared.
        declared_phase: None,
    })
}

pub(super) fn component_kind_sort_key(kind: AlleleComponentKind) -> u8 {
    match kind {
        AlleleComponentKind::Snp => 0,
        AlleleComponentKind::Deletion => 1,
        AlleleComponentKind::Insertion => 2,
        AlleleComponentKind::Delins => 3,
        AlleleComponentKind::Symbolic => 4,
    }
}

pub(super) fn component_labels_from_group(group: &[&VcfPosition]) -> Vec<String> {
    let mut components = Vec::new();
    for variant in group {
        components.extend(variant.event().components);
    }
    components.sort_by(|a, b| {
        a.position
            .cmp(&b.position)
            .then_with(|| component_kind_sort_key(a.kind).cmp(&component_kind_sort_key(b.kind)))
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
    });
    components.dedup();
    components
        .into_iter()
        .map(|component| component.label())
        .collect()
}

pub(super) fn phased_variant_from_group(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    chrom: &str,
    group: &[&VcfPosition],
    genetic_code: crate::genetic_code::GeneticCode,
) -> Option<VariantInfo> {
    let compound = compound_allele_from_variants(reference, group)?;
    let compound_event = compound.event();
    let (has_indel, has_substitution, indel_components) = group_component_flags(group);
    if !has_indel {
        return None;
    }
    let event_class = if group.len() > 1
        && (has_substitution || indel_components > 1 || !compound_event.class.has_indel_component())
    {
        AlleleEventClass::ComplexIndel
    } else {
        compound_event.class
    };

    let event_components = component_labels_from_group(group);
    let (change_type, aa_changes, aa_changes_local, ref_codon, alt_codon) =
        protein_effect_for_indel(gene, reference, &compound, genetic_code);
    let nmd = super::nmd::indel_nmd_prediction(gene, reference, &compound, genetic_code);

    Some(VariantInfo {
        chrom: chrom.to_string(),
        gene: gene.name.clone(),
        positions: vec![compound.record_start],
        ref_bases: vec![compound.ref_allele],
        base_changes: vec![compound.alt_allele],
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
        original_dp: None,
        original_freq: None,
        original_info: None,
        event_class: Some(event_class.as_str().to_string()),
        event_components,
        annotations: crate::variants::VariantAnnotations {
            nmd,
            ..Default::default()
        },
    })
}

/// Variants close enough to have been on the same molecule, grouped into the
/// windows a read must span to tell their combinations apart.
///
/// Proximity only proposes; it never asserts linkage. The caller asks the reads
/// which of these variants actually travel together, so a group here is a
/// question, not an answer. A group is returned only when it could describe a
/// phased haplotype at all: two or more variants, at least one an indel.
pub fn local_haplotype_components<'a>(
    gene: &Gene,
    snp_list: &'a [VcfPosition],
) -> Vec<Vec<&'a VcfPosition>> {
    let local_variants = snp_list
        .iter()
        .filter(|variant| variant.overlaps_interval(gene.start, gene.end))
        .filter(|variant| {
            // In a spliced transcript model the gene span includes introns; only
            // exonic (coding) variants may join a phased coding haplotype, so an
            // intronic indel near an exon boundary is not merged with exonic SNVs.
            // Asked of the bases the record changes, not of the base it starts
            // on: a record padded with reference bases can begin in an intron
            // and change an exonic base, and testing POS dropped it from the
            // coding haplotype it belongs to.
            !has_transcript_cds_model(gene)
                || variant
                    .changed_positions()
                    .iter()
                    .any(|&position| transcript_offset_for_position(gene, position).is_some())
        })
        .filter(|variant| {
            variant_has_indel_component(variant) || !variant.substitution_components().is_empty()
        })
        .filter_map(|variant| phased_indel_window(gene, variant).map(|window| (variant, window)))
        .collect::<Vec<_>>();

    if local_variants.len() < 2
        || !local_variants
            .iter()
            .any(|(variant, _)| variant_has_indel_component(variant))
    {
        return Vec::new();
    }

    let mut components = Vec::new();
    let mut visited = vec![false; local_variants.len()];

    for start_idx in 0..local_variants.len() {
        if visited[start_idx] {
            continue;
        }

        let mut component = Vec::new();
        let mut stack = vec![start_idx];
        visited[start_idx] = true;
        while let Some(idx) = stack.pop() {
            component.push(idx);
            for next_idx in 0..local_variants.len() {
                if visited[next_idx] {
                    continue;
                }
                if component.iter().any(|component_idx| {
                    haplotype_windows_linked(
                        local_variants[*component_idx].1,
                        local_variants[next_idx].1,
                    )
                }) {
                    visited[next_idx] = true;
                    stack.push(next_idx);
                }
            }
        }

        if component.len() < 2
            || !component
                .iter()
                .any(|idx| variant_has_indel_component(local_variants[*idx].0))
        {
            continue;
        }

        component.sort_by_key(|idx| local_variants[*idx].0.record_start);
        components.push(
            component
                .into_iter()
                .map(|idx| local_variants[idx].0)
                .collect::<Vec<_>>(),
        );
    }

    components
}
