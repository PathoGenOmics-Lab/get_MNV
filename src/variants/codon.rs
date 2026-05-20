//! Codon-level variant processing: SNP grouping, MNV detection, and amino
//! acid change calculation.

use super::event::{AlleleComponentKind, AlleleEventClass};
use super::types::*;
use crate::io::VcfPosition;
use crate::utils::{aa_three_letter, determine_change_type, iupac_aa, reverse_complement};
use std::collections::{BTreeMap, BTreeSet};

const LOCAL_HAPLOTYPE_JOIN_DISTANCE: usize = 3;
const MAX_LOCAL_HAPLOTYPE_VARIANTS: usize = 8;

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

fn event_metadata(position: usize, ref_allele: &str, alt_allele: &str) -> (String, Vec<String>) {
    let event = crate::variants::decompose_allele(position, ref_allele, alt_allele);
    (event.class.as_str().to_string(), event.component_labels())
}

fn variant_event_metadata(variant: &crate::io::VcfPosition) -> (String, Vec<String>) {
    event_metadata(variant.position, &variant.ref_allele, &variant.alt_allele)
}

fn translate_cds(seq: &str, genetic_code: crate::genetic_code::GeneticCode) -> String {
    let mut protein = String::with_capacity(seq.len() / 3);
    for codon in seq.as_bytes().chunks_exact(3) {
        let codon = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        protein.push(genetic_code.translate(&codon));
    }
    protein
}

fn aa_segment_three_letter(segment: &[char]) -> String {
    segment
        .iter()
        .map(|aa| aa_three_letter(*aa))
        .collect::<Vec<_>>()
        .join("")
}

fn describe_protein_change(
    ref_protein: &str,
    alt_protein: &str,
    protein_offset: usize,
    frameshift: bool,
) -> (String, String) {
    if ref_protein == alt_protein {
        if frameshift {
            let ref_chars: Vec<char> = ref_protein.chars().collect();
            let idx = ref_chars.len().saturating_sub(1);
            let local_pos = idx + 1;
            let protein_pos = protein_offset + local_pos;
            let ref_aa = ref_chars.get(idx).copied().unwrap_or('X');
            return (
                format!("{}{}fs", aa_three_letter(ref_aa), protein_pos),
                format!("{}{}fs", aa_three_letter(ref_aa), local_pos),
            );
        }
        return ("Synonymous".to_string(), "Synonymous".to_string());
    }

    let ref_chars: Vec<char> = ref_protein.chars().collect();
    let alt_chars: Vec<char> = alt_protein.chars().collect();
    let mut prefix = 0usize;
    while prefix < ref_chars.len()
        && prefix < alt_chars.len()
        && ref_chars[prefix] == alt_chars[prefix]
    {
        prefix += 1;
    }

    let local_pos = prefix + 1;
    let protein_pos = protein_offset + local_pos;
    let ref_aa = ref_chars.get(prefix).copied().unwrap_or('X');
    let alt_aa = alt_chars.get(prefix).copied().unwrap_or('X');

    if frameshift {
        let local = format!(
            "{}{}{}fs",
            aa_three_letter(ref_aa),
            local_pos,
            aa_three_letter(alt_aa)
        );
        let protein = format!(
            "{}{}{}fs",
            aa_three_letter(ref_aa),
            protein_pos,
            aa_three_letter(alt_aa)
        );
        return (protein, local);
    }

    let mut suffix = 0usize;
    while suffix + prefix < ref_chars.len()
        && suffix + prefix < alt_chars.len()
        && ref_chars[ref_chars.len() - 1 - suffix] == alt_chars[alt_chars.len() - 1 - suffix]
    {
        suffix += 1;
    }

    let ref_mid_end = ref_chars.len().saturating_sub(suffix);
    let alt_mid_end = alt_chars.len().saturating_sub(suffix);
    let ref_mid = &ref_chars[prefix..ref_mid_end];
    let alt_mid = &alt_chars[prefix..alt_mid_end];

    let build = |start_pos: usize| -> String {
        if ref_mid.len() == 1 && alt_mid.len() == 1 {
            return format!(
                "{}{}{}",
                aa_three_letter(ref_mid[0]),
                start_pos,
                aa_three_letter(alt_mid[0])
            );
        }

        let end_pos = start_pos + ref_mid.len().saturating_sub(1);
        let ref_start = ref_mid.first().copied().unwrap_or(ref_aa);
        let ref_end = ref_mid.last().copied().unwrap_or(ref_aa);
        let ref_range = if start_pos == end_pos {
            format!("{}{}", aa_three_letter(ref_start), start_pos)
        } else {
            format!(
                "{}{}_{}{}",
                aa_three_letter(ref_start),
                start_pos,
                aa_three_letter(ref_end),
                end_pos
            )
        };

        if ref_mid.is_empty() {
            let anchor = start_pos.saturating_sub(1);
            format!(
                "{anchor}_{start_pos}ins{}",
                aa_segment_three_letter(alt_mid)
            )
        } else if alt_mid.is_empty() {
            format!("{ref_range}del")
        } else {
            format!("{ref_range}delins{}", aa_segment_three_letter(alt_mid))
        }
    };

    (build(protein_pos), build(local_pos))
}

fn complement_base(base: char) -> char {
    match base.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        'N' => 'N',
        other => other,
    }
}

fn has_transcript_cds_model(gene: &Gene) -> bool {
    !gene.cds_segments.is_empty()
}

fn transcript_sequence_for_gene(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
) -> Option<String> {
    if !has_transcript_cds_model(gene) {
        return None;
    }

    let mut cds = String::new();
    for segment in &gene.cds_segments {
        if segment.start == 0
            || segment.end == 0
            || segment.start > segment.end
            || segment.end > reference.sequence.len()
        {
            return None;
        }
        let genomic = &reference.sequence[(segment.start - 1)..segment.end];
        match gene.strand {
            Strand::Plus => cds.push_str(genomic),
            Strand::Minus => cds.push_str(&reverse_complement(genomic)),
        }
    }
    Some(cds)
}

fn transcript_offset_for_position(gene: &Gene, position: usize) -> Option<usize> {
    if !has_transcript_cds_model(gene) {
        return None;
    }

    let mut offset = 0usize;
    for segment in &gene.cds_segments {
        if position >= segment.start && position <= segment.end {
            return Some(match gene.strand {
                Strand::Plus => offset + position.saturating_sub(segment.start),
                Strand::Minus => offset + segment.end.saturating_sub(position),
            });
        }
        offset += segment.end.saturating_sub(segment.start) + 1;
    }
    None
}

fn insertion_anchor_in_segment(segment: &CdsSegment, position: usize) -> bool {
    position >= segment.start && position < segment.end
}

fn variant_overlaps_coding_model(gene: &Gene, variant: &VcfPosition) -> bool {
    if !has_transcript_cds_model(gene) {
        return variant.overlaps_interval(gene.start, gene.end);
    }

    gene.cds_segments
        .iter()
        .any(|segment| variant.overlaps_interval(segment.start, segment.end))
}

fn variant_touched_transcript_offsets(gene: &Gene, variant: &VcfPosition) -> Vec<usize> {
    if !has_transcript_cds_model(gene) {
        return Vec::new();
    }

    let mut offsets = BTreeSet::new();
    for component in variant.event().components {
        match component.kind {
            AlleleComponentKind::Snp => {
                if let Some(offset) = transcript_offset_for_position(gene, component.position) {
                    offsets.insert(offset);
                }
            }
            AlleleComponentKind::Insertion => {
                for segment in &gene.cds_segments {
                    if insertion_anchor_in_segment(segment, component.position) {
                        if let Some(offset) =
                            transcript_offset_for_position(gene, component.position)
                        {
                            offsets.insert(offset);
                        }
                    }
                }
            }
            AlleleComponentKind::Deletion
            | AlleleComponentKind::Delins
            | AlleleComponentKind::Symbolic => {
                let component_end =
                    component.position + component.ref_allele.len().saturating_sub(1);
                for segment in &gene.cds_segments {
                    let overlap_start = component.position.max(segment.start);
                    let overlap_end = component_end.min(segment.end);
                    if overlap_start <= overlap_end {
                        for pos in overlap_start..=overlap_end {
                            if let Some(offset) = transcript_offset_for_position(gene, pos) {
                                offsets.insert(offset);
                            }
                        }
                    }
                }
            }
        }
    }

    offsets.into_iter().collect()
}

fn first_touched_transcript_offset(gene: &Gene, variant: &VcfPosition) -> Option<usize> {
    variant_touched_transcript_offsets(gene, variant)
        .into_iter()
        .min()
}

fn coding_sequence_for_gene(gene: &Gene, reference: &crate::io::Reference<'_>) -> Option<String> {
    if let Some(transcript_cds) = transcript_sequence_for_gene(gene, reference) {
        return Some(transcript_cds);
    }

    let (eff_start, eff_end) = effective_bounds(gene);
    if eff_start == 0 || eff_end == 0 || eff_start > eff_end || eff_end > reference.sequence.len() {
        return None;
    }
    let genomic = &reference.sequence[(eff_start - 1)..eff_end];
    Some(match gene.strand {
        Strand::Plus => genomic.to_string(),
        Strand::Minus => reverse_complement(genomic),
    })
}

fn apply_allele_to_transcript(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &crate::io::VcfPosition,
) -> Option<String> {
    if !has_transcript_cds_model(gene) {
        return None;
    }

    let event = variant.event();
    if matches!(event.class, AlleleEventClass::Symbolic) {
        return None;
    }

    let mut touched_transcript = false;
    let mut cds = String::new();

    for segment in &gene.cds_segments {
        if segment.start == 0
            || segment.end == 0
            || segment.start > segment.end
            || segment.end > reference.sequence.len()
        {
            return None;
        }

        let mut snps: BTreeMap<usize, String> = BTreeMap::new();
        let mut insertions: BTreeMap<usize, String> = BTreeMap::new();
        let mut deleted: BTreeSet<usize> = BTreeSet::new();

        for component in &event.components {
            match component.kind {
                AlleleComponentKind::Snp => {
                    if component.position >= segment.start && component.position <= segment.end {
                        match snps.get(&component.position) {
                            Some(existing)
                                if !existing.eq_ignore_ascii_case(&component.alt_allele) =>
                            {
                                return None;
                            }
                            Some(_) => {}
                            None => {
                                snps.insert(component.position, component.alt_allele.clone());
                            }
                        }
                        touched_transcript = true;
                    }
                }
                AlleleComponentKind::Insertion => {
                    if insertion_anchor_in_segment(segment, component.position) {
                        match insertions.get(&component.position) {
                            Some(existing)
                                if !existing.eq_ignore_ascii_case(&component.alt_allele) =>
                            {
                                return None;
                            }
                            Some(_) => {}
                            None => {
                                insertions.insert(component.position, component.alt_allele.clone());
                            }
                        }
                        touched_transcript = true;
                    }
                }
                AlleleComponentKind::Deletion => {
                    let del_end = component.position + component.ref_allele.len().saturating_sub(1);
                    let overlap_start = component.position.max(segment.start);
                    let overlap_end = del_end.min(segment.end);
                    if overlap_start <= overlap_end {
                        for pos in overlap_start..=overlap_end {
                            deleted.insert(pos);
                        }
                        touched_transcript = true;
                    }
                }
                AlleleComponentKind::Delins | AlleleComponentKind::Symbolic => return None,
            }
        }

        if snps.keys().any(|pos| deleted.contains(pos))
            || insertions.keys().any(|pos| deleted.contains(pos))
        {
            return None;
        }

        let mut genomic = String::with_capacity(segment.end - segment.start + 1);
        for pos in segment.start..=segment.end {
            if !deleted.contains(&pos) {
                if let Some(alt_base) = snps.get(&pos) {
                    genomic.push_str(alt_base);
                } else {
                    genomic.push(reference.sequence.as_bytes()[pos - 1] as char);
                }
            }
            if let Some(inserted) = insertions.get(&pos) {
                genomic.push_str(inserted);
            }
        }

        match gene.strand {
            Strand::Plus => cds.push_str(&genomic),
            Strand::Minus => cds.push_str(&reverse_complement(&genomic)),
        }
    }

    touched_transcript.then_some(cds)
}

fn apply_allele_to_feature(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &crate::io::VcfPosition,
) -> Option<String> {
    if let Some(transcript_alt) = apply_allele_to_transcript(gene, reference, variant) {
        return Some(transcript_alt);
    }

    let (eff_start, eff_end) = effective_bounds(gene);
    if eff_start == 0 || eff_end == 0 || eff_start > eff_end || eff_end > reference.sequence.len() {
        return None;
    }

    let event = variant.event();
    if matches!(event.class, AlleleEventClass::Symbolic) {
        return None;
    }

    let mut snps: BTreeMap<usize, String> = BTreeMap::new();
    let mut insertions: BTreeMap<usize, String> = BTreeMap::new();
    let mut deleted: BTreeSet<usize> = BTreeSet::new();
    let mut touched_feature = false;

    for component in event.components {
        match component.kind {
            AlleleComponentKind::Snp => {
                if component.position >= eff_start && component.position <= eff_end {
                    match snps.get(&component.position) {
                        Some(existing) if !existing.eq_ignore_ascii_case(&component.alt_allele) => {
                            return None;
                        }
                        Some(_) => {}
                        None => {
                            snps.insert(component.position, component.alt_allele);
                        }
                    }
                    touched_feature = true;
                }
            }
            AlleleComponentKind::Insertion => {
                if component.position >= eff_start && component.position < eff_end {
                    match insertions.get(&component.position) {
                        Some(existing) if !existing.eq_ignore_ascii_case(&component.alt_allele) => {
                            return None;
                        }
                        Some(_) => {}
                        None => {
                            insertions.insert(component.position, component.alt_allele);
                        }
                    }
                    touched_feature = true;
                }
            }
            AlleleComponentKind::Deletion => {
                let del_end = component.position + component.ref_allele.len().saturating_sub(1);
                let overlap_start = component.position.max(eff_start);
                let overlap_end = del_end.min(eff_end);
                if overlap_start <= overlap_end {
                    for pos in overlap_start..=overlap_end {
                        deleted.insert(pos);
                    }
                    touched_feature = true;
                }
            }
            AlleleComponentKind::Delins | AlleleComponentKind::Symbolic => return None,
        }
    }

    if !touched_feature {
        return None;
    }

    if snps.keys().any(|pos| deleted.contains(pos))
        || insertions.keys().any(|pos| deleted.contains(pos))
    {
        return None;
    }

    let mut genomic = String::with_capacity(reference.sequence[(eff_start - 1)..eff_end].len());
    for pos in eff_start..=eff_end {
        if !deleted.contains(&pos) {
            if let Some(alt_base) = snps.get(&pos) {
                genomic.push_str(alt_base);
            } else {
                genomic.push(reference.sequence.as_bytes()[pos - 1] as char);
            }
        }
        if let Some(inserted) = insertions.get(&pos) {
            genomic.push_str(inserted);
        }
    }

    Some(match gene.strand {
        Strand::Plus => genomic,
        Strand::Minus => reverse_complement(&genomic),
    })
}

fn protein_effect_for_indel(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
    genetic_code: crate::genetic_code::GeneticCode,
) -> (Vec<String>, Vec<String>, Option<String>, Option<String>) {
    let Some(ref_cds) = coding_sequence_for_gene(gene, reference) else {
        return (
            vec!["Unknown".to_string()],
            vec!["Unknown".to_string()],
            None,
            None,
        );
    };
    let Some(alt_cds) = apply_allele_to_feature(gene, reference, variant) else {
        return (
            vec!["Unknown".to_string()],
            vec!["Unknown".to_string()],
            None,
            None,
        );
    };

    let frameshift = (alt_cds.len() as isize - ref_cds.len() as isize) % 3 != 0;
    let ref_protein = translate_cds(&ref_cds, genetic_code);
    let alt_protein = translate_cds(&alt_cds, genetic_code);
    let (protein_change, local_change) =
        describe_protein_change(&ref_protein, &alt_protein, gene.protein_offset, frameshift);

    let transcript_codon = first_touched_transcript_offset(gene, variant).and_then(|offset| {
        let codon_start = (offset / 3) * 3;
        let codon_end = codon_start + 3;
        let ref_codon = ref_cds.get(codon_start..codon_end).map(str::to_string)?;
        let alt_codon = alt_cds.get(codon_start..codon_end).map(str::to_string);
        Some((ref_codon, alt_codon))
    });

    if let Some((ref_codon, alt_codon)) = transcript_codon {
        return (
            vec![protein_change],
            vec![local_change],
            Some(ref_codon),
            alt_codon,
        );
    }

    let first_codon_start = {
        let (eff_start, eff_end) = effective_bounds(gene);
        let anchor = variant.position.clamp(eff_start, eff_end);
        codon_bounds_for_position(gene, anchor).map(|(start, _)| start)
    };

    let (ref_codon, alt_codon) = if let Some(codon_start) = first_codon_start {
        let (eff_start, eff_end) = effective_bounds(gene);
        let codon_end = codon_start + 2;
        let local = match gene.strand {
            Strand::Plus => codon_start.saturating_sub(eff_start),
            Strand::Minus => eff_end.saturating_sub(codon_end),
        };
        let ref_codon = ref_cds.get(local..local + 3).map(str::to_string);
        let alt_codon = alt_cds.get(local..local + 3).map(str::to_string);
        (ref_codon, alt_codon)
    } else {
        (None, None)
    };

    (
        vec![protein_change],
        vec![local_change],
        ref_codon,
        alt_codon,
    )
}

fn coding_delta_for_variant(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
) -> Option<isize> {
    let ref_cds = coding_sequence_for_gene(gene, reference)?;
    let alt_cds = apply_allele_to_feature(gene, reference, variant)?;
    Some(alt_cds.len() as isize - ref_cds.len() as isize)
}

fn indel_change_type_for_variant(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
) -> ChangeType {
    let frameshift = if variant.alt_allele.starts_with('<') {
        true
    } else if let Some(delta) = coding_delta_for_variant(gene, reference, variant) {
        delta % 3 != 0
    } else {
        (variant.ref_allele.len() as isize - variant.alt_allele.len() as isize) % 3 != 0
    };

    if frameshift {
        ChangeType::FrameshiftIndel
    } else {
        ChangeType::InFrameIndel
    }
}

fn phased_indel_window(gene: &Gene, variant: &VcfPosition) -> Option<(usize, usize)> {
    let event = variant.event();
    let (eff_start, eff_end) = effective_bounds(gene);
    if eff_start == 0 || eff_end == 0 || eff_start > eff_end {
        return None;
    }

    let mut windows = Vec::new();
    for pos in [variant.position, event.affected_start, event.affected_end] {
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

fn haplotype_windows_linked(a: (usize, usize), b: (usize, usize)) -> bool {
    let a_end = a.1.saturating_add(LOCAL_HAPLOTYPE_JOIN_DISTANCE);
    let b_end = b.1.saturating_add(LOCAL_HAPLOTYPE_JOIN_DISTANCE);
    a.0 <= b_end && b.0 <= a_end
}

fn variant_has_indel_component(variant: &VcfPosition) -> bool {
    variant.event().class.has_indel_component()
}

fn group_component_flags(group: &[&VcfPosition]) -> (bool, bool, usize) {
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

fn add_phased_candidate(
    out: &mut Vec<VariantInfo>,
    seen: &mut BTreeSet<(usize, String, String)>,
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
    group: &[&VcfPosition],
) {
    let (has_indel, _, _) = group_component_flags(group);
    if !has_indel || group.len() < 2 {
        return;
    }
    if let Some(candidate) = phased_variant_from_group(gene, reference, chrom, group, genetic_code)
    {
        let key = (
            candidate.positions[0],
            candidate.ref_bases[0].clone(),
            candidate.base_changes[0].clone(),
        );
        if seen.insert(key) {
            out.push(candidate);
        }
    }
}

fn compound_allele_from_variants(
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
        start = start.min(event.affected_start).min(variant.position);
        end = end
            .max(event.affected_end)
            .max(variant.position + variant.ref_allele.len().saturating_sub(1));

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
        position: start,
        ref_allele,
        alt_allele,
        original_dp: None,
        original_freq: None,
        original_info: None,
    })
}

fn component_kind_sort_key(kind: AlleleComponentKind) -> u8 {
    match kind {
        AlleleComponentKind::Snp => 0,
        AlleleComponentKind::Deletion => 1,
        AlleleComponentKind::Insertion => 2,
        AlleleComponentKind::Delins => 3,
        AlleleComponentKind::Symbolic => 4,
    }
}

fn component_labels_from_group(group: &[&VcfPosition]) -> Vec<String> {
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

fn phased_variant_from_group(
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

    let change_type = indel_change_type_for_variant(gene, reference, &compound);
    let event_components = component_labels_from_group(group);
    let (aa_changes, aa_changes_local, ref_codon, alt_codon) =
        protein_effect_for_indel(gene, reference, &compound, genetic_code);

    Some(VariantInfo {
        chrom: chrom.to_string(),
        gene: gene.name.clone(),
        positions: vec![compound.position],
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
        ref_codon,
        snp_codon: None,
        mnv_codon: alt_codon,
        original_dp: None,
        original_freq: None,
        original_info: None,
        event_class: Some(event_class.as_str().to_string()),
        event_components,
    })
}

pub fn build_phased_indel_haplotype_variants(
    gene: &Gene,
    snp_list: &[VcfPosition],
    reference: &crate::io::Reference<'_>,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
) -> Vec<VariantInfo> {
    let local_variants = snp_list
        .iter()
        .filter(|variant| variant.overlaps_interval(gene.start, gene.end))
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

    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
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

        if component.len() <= MAX_LOCAL_HAPLOTYPE_VARIANTS {
            let subset_count = 1usize << component.len();
            for mask in 1usize..subset_count {
                if mask.count_ones() < 2 {
                    continue;
                }
                let group = component
                    .iter()
                    .enumerate()
                    .filter_map(|(bit, idx)| {
                        (mask & (1usize << bit) != 0).then_some(local_variants[*idx].0)
                    })
                    .collect::<Vec<_>>();
                add_phased_candidate(
                    &mut out,
                    &mut seen,
                    gene,
                    reference,
                    chrom,
                    genetic_code,
                    &group,
                );
            }
        } else {
            for i in 0..component.len() {
                for j in (i + 1)..component.len() {
                    let group = [
                        local_variants[component[i]].0,
                        local_variants[component[j]].0,
                    ];
                    add_phased_candidate(
                        &mut out,
                        &mut seen,
                        gene,
                        reference,
                        chrom,
                        genetic_code,
                        &group,
                    );
                }
            }
            let group = component
                .iter()
                .map(|idx| local_variants[*idx].0)
                .collect::<Vec<_>>();
            add_phased_candidate(
                &mut out,
                &mut seen,
                gene,
                reference,
                chrom,
                genetic_code,
                &group,
            );
        }
    }

    out
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
            event_class: None,
            event_components: Vec::new(),
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

#[derive(Debug, Clone)]
struct TranscriptSnp {
    snp: Snp,
    transcript_offset: usize,
    transcript_alt_base: char,
}

fn transcript_oriented_base(base: &str, strand: Strand) -> Option<char> {
    let base = base.chars().next()?;
    Some(match strand {
        Strand::Plus => base.to_ascii_uppercase(),
        Strand::Minus => complement_base(base),
    })
}

fn construct_transcript_codon(ref_codon: &str, target_snps: &[&TranscriptSnp]) -> String {
    let mut codon = ref_codon.chars().collect::<Vec<_>>();
    for snp in target_snps {
        let idx = snp.transcript_offset % 3;
        if let Some(slot) = codon.get_mut(idx) {
            *slot = snp.transcript_alt_base.to_ascii_uppercase();
        }
    }
    codon.into_iter().collect()
}

fn process_transcript_codon(
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
    let combined_aa = iupac_aa(&format!("{orig_aa}{aa_pos}{mut_aa}"));
    let change_type = ChangeType::from_label(&determine_change_type(&combined_aa));

    let snp_changes = codon_snps
        .iter()
        .map(|snp| {
            let single_codon = construct_transcript_codon(ref_codon, &[snp]);
            let single_aa = genetic_code.translate_seq(single_codon.as_bytes());
            iupac_aa(&format!("{orig_aa}{aa_pos}{single_aa}"))
        })
        .collect::<Vec<_>>();

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
    }
}

fn get_mnv_variants_for_transcript(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
    genetic_code: crate::genetic_code::GeneticCode,
) -> Vec<VariantInfo> {
    let Some(ref_cds) = transcript_sequence_for_gene(gene, reference) else {
        return Vec::new();
    };

    let mut variants = Vec::new();
    let mut codon_to_snps: BTreeMap<usize, Vec<TranscriptSnp>> = BTreeMap::new();
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
                let group = codon_to_snps.entry(codon_start).or_default();
                if !group.iter().any(|snp| snp.transcript_offset == offset) {
                    group.push(TranscriptSnp {
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
                    });
                } else {
                    log::debug!(
                        "Skipping duplicate ALT at transcript offset {} in gene '{}' (multi-allelic split)",
                        offset,
                        gene.name
                    );
                }
            }
        } else {
            indels.push(variant.clone());
        }
    }

    let mut codon_starts = codon_to_snps.keys().copied().collect::<Vec<_>>();
    codon_starts.sort_unstable();

    for codon_start in codon_starts {
        let codon_end = codon_start + 3;
        if codon_end > ref_cds.len() {
            continue;
        }
        let ref_codon = &ref_cds[codon_start..codon_end];
        let mut codon_snps = codon_to_snps.get(&codon_start).cloned().unwrap_or_default();
        codon_snps.sort_by_key(|s| s.transcript_offset);

        let overlaps_indel = indels.iter().any(|indel| {
            variant_touched_transcript_offsets(gene, indel)
                .into_iter()
                .any(|offset| offset >= codon_start && offset < codon_end)
        });

        let mut upstream_shift: isize = 0;
        let mut has_symbolic_sv = false;
        for indel in &indels {
            let Some(first_offset) = first_touched_transcript_offset(gene, indel) else {
                continue;
            };
            if first_offset < codon_start {
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
        let change_type = indel_change_type_for_variant(gene, reference, &indel);
        let (aa_changes, aa_changes_local, ref_codon, alt_codon) =
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
            ref_codon,
            snp_codon: None,
            mnv_codon: alt_codon,
            original_dp: indel.original_dp.map(|value| vec![value]),
            original_freq: indel.original_freq.map(|value| vec![value]),
            original_info: indel.original_info.clone(),
            event_class: Some(event_class),
            event_components,
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
    if has_transcript_cds_model(gene) {
        return get_mnv_variants_for_transcript(gene, snp_list, reference, chrom, genetic_code);
    }

    let mut variants = Vec::new();

    let mut codon_to_snps: BTreeMap<usize, Vec<crate::variants::Snp>> = BTreeMap::new();
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
                    let group = codon_to_snps.entry(codon_start).or_default();
                    // After --split-multiallelic, two ALT alleles at the same
                    // position produce separate VCF records.  Only keep the first
                    // ALT per position within each codon group — otherwise
                    // construct_codon would overwrite the same base twice, and the
                    // MNV annotation would be incorrect.
                    if !group.iter().any(|s| s.position == component.position) {
                        group.push(crate::variants::Snp {
                            index: component.position,
                            position: component.position,
                            ref_base: component.ref_allele,
                            base: component.alt_allele,
                            original_dp: variant.original_dp,
                            original_freq: variant.original_freq,
                            original_info: variant.original_info.clone(),
                        });
                    } else {
                        log::debug!(
                            "Skipping duplicate ALT at position {} in gene '{}' (multi-allelic split)",
                            component.position,
                            gene.name
                        );
                    }
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

        let overlaps_indel = indels
            .iter()
            .any(|indel| indel.overlaps_interval(codon_start, codon_end));

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
        let change_type = indel_change_type_for_variant(gene, reference, &indel);
        let (aa_changes, aa_changes_local, ref_codon, alt_codon) =
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
            ref_codon,
            snp_codon: None,
            mnv_codon: alt_codon,
            original_dp: indel.original_dp.map(|value| vec![value]),
            original_freq: indel.original_freq.map(|value| vec![value]),
            original_info: indel.original_info.clone(),
            event_class: Some(event_class),
            event_components,
        });
    }

    variants
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
        original_dp: vcf_pos.original_dp.map(|dp| vec![dp]),
        original_freq: vcf_pos.original_freq.map(|freq| vec![freq]),
        original_info: vcf_pos.original_info.clone(),
        event_class: Some(event.class.as_str().to_string()),
        event_components: event.component_labels(),
    }
}

#[cfg(test)]
mod tests {
    use super::{build_phased_indel_haplotype_variants, codon_bounds_for_position, process_codon};
    use crate::utils::reverse_complement;
    use crate::variants::{CdsSegment, ChangeType, CodonInfo, Gene, Snp, Strand, VariantType};

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
            transcript_id: None,
            cds_segments: Vec::new(),
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
            transcript_id: None,
            cds_segments: Vec::new(),
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
    fn test_vcf_mnv_record_is_decomposed_into_codon_haplotype() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 1,
            end: 9,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAATTT",
        };
        let variants = crate::variants::get_mnv_variants_for_gene(
            &gene,
            &[crate::io::VcfPosition {
                position: 4,
                ref_allele: "AA".to_string(),
                alt_allele: "CC".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            }],
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].positions, vec![4, 5]);
        assert_eq!(variants[0].event_class.as_deref(), Some("mnv"));
        assert_eq!(
            variants[0].event_components,
            vec!["SNV:4:A>C".to_string(), "SNV:5:A>C".to_string()]
        );
    }

    #[test]
    fn test_transcript_model_groups_mnv_across_exon_junction() {
        let gene = Gene {
            name: "tx_cds".to_string(),
            start: 1,
            end: 14,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: Some("tx1".to_string()),
            cds_segments: vec![
                CdsSegment { start: 1, end: 4 },
                CdsSegment { start: 10, end: 14 },
            ],
        };
        let reference = crate::io::Reference {
            sequence: "ATGACCCCCAATTT",
        };
        let variants = crate::variants::get_mnv_variants_for_gene(
            &gene,
            &[
                crate::io::VcfPosition {
                    position: 4,
                    ref_allele: "A".to_string(),
                    alt_allele: "C".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
                crate::io::VcfPosition {
                    position: 10,
                    ref_allele: "A".to_string(),
                    alt_allele: "G".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
            ],
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant_type, VariantType::SnpMnv);
        assert_eq!(variants[0].positions, vec![4, 10]);
        assert_eq!(variants[0].ref_codon.as_deref(), Some("AAA"));
        assert_eq!(variants[0].mnv_codon.as_deref(), Some("CGA"));
        assert_eq!(variants[0].aa_changes, vec!["Lys2Arg".to_string()]);
        assert_eq!(
            variants[0].event_components,
            vec!["SNV:4:A>C".to_string(), "SNV:10:A>G".to_string()]
        );
    }

    #[test]
    fn test_transcript_model_restored_frame_does_not_mark_downstream_snp_frameshift() {
        let gene = Gene {
            name: "tx_cds".to_string(),
            start: 1,
            end: 25,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: Some("tx1".to_string()),
            cds_segments: vec![
                CdsSegment { start: 1, end: 6 },
                CdsSegment { start: 20, end: 25 },
            ],
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAAGGGGGGGGGGGGGTTTCCC",
        };
        let variants = crate::variants::get_mnv_variants_for_gene(
            &gene,
            &[
                crate::io::VcfPosition {
                    position: 3,
                    ref_allele: "G".to_string(),
                    alt_allele: "GA".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
                crate::io::VcfPosition {
                    position: 20,
                    ref_allele: "TT".to_string(),
                    alt_allele: "T".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
                crate::io::VcfPosition {
                    position: 23,
                    ref_allele: "C".to_string(),
                    alt_allele: "A".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
            ],
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        let snp_row = variants
            .iter()
            .find(|variant| variant.positions == vec![23])
            .expect("downstream SNP row");
        assert!(!snp_row.aa_changes[0].contains("(fs)"));
        assert!(!matches!(
            snp_row.change_type,
            ChangeType::FrameshiftSynonymous
                | ChangeType::FrameshiftNonSynonymous
                | ChangeType::FrameshiftStopGained
                | ChangeType::FrameshiftStopLost
                | ChangeType::FrameshiftUnknown
        ));
    }

    #[test]
    fn test_indel_reports_event_components_and_frameshift_protein_effect() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 1,
            end: 9,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAATTT",
        };
        let variants = crate::variants::get_mnv_variants_for_gene(
            &gene,
            &[crate::io::VcfPosition {
                position: 6,
                ref_allele: "A".to_string(),
                alt_allele: "AT".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            }],
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant_type, VariantType::Indel);
        assert_eq!(variants[0].event_class.as_deref(), Some("insertion"));
        assert_eq!(variants[0].event_components, vec!["INS:6:+T".to_string()]);
        assert!(variants[0].aa_changes[0].contains("fs"));
    }

    #[test]
    fn test_deletion_anchored_before_cds_keeps_protein_effect() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 2,
            end: 10,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "CATGAAATTT",
        };
        let variants = crate::variants::get_mnv_variants_for_gene(
            &gene,
            &[crate::io::VcfPosition {
                position: 1,
                ref_allele: "CA".to_string(),
                alt_allele: "C".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            }],
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant_type, VariantType::Indel);
        assert_eq!(variants[0].event_class.as_deref(), Some("deletion"));
        assert_eq!(variants[0].event_components, vec!["DEL:2:A".to_string()]);
        assert_eq!(variants[0].change_type, ChangeType::FrameshiftIndel);
        assert_ne!(variants[0].aa_changes, vec!["Unknown".to_string()]);
        assert!(variants[0].aa_changes[0].contains("fs"));
    }

    #[test]
    fn test_insertion_after_cds_end_is_not_coding() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 1,
            end: 9,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAATTTC",
        };
        let variants = crate::variants::get_mnv_variants_for_gene(
            &gene,
            &[crate::io::VcfPosition {
                position: 9,
                ref_allele: "T".to_string(),
                alt_allele: "TA".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            }],
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert!(variants.is_empty());
    }

    #[test]
    fn test_insertion_after_codon_end_does_not_mask_that_codon_snp() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 1,
            end: 12,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAATTTCCC",
        };
        let variants = crate::variants::get_mnv_variants_for_gene(
            &gene,
            &[
                crate::io::VcfPosition {
                    position: 7,
                    ref_allele: "T".to_string(),
                    alt_allele: "C".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
                crate::io::VcfPosition {
                    position: 9,
                    ref_allele: "T".to_string(),
                    alt_allele: "TA".to_string(),
                    original_dp: None,
                    original_freq: None,
                    original_info: None,
                },
            ],
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );
        let snp_row = variants
            .iter()
            .find(|variant| variant.variant_type == VariantType::Snp)
            .expect("SNP row");

        assert_ne!(snp_row.change_type, ChangeType::IndelOverlap);
        assert_ne!(snp_row.aa_changes, vec!["Unknown".to_string()]);
    }

    #[test]
    fn test_build_phased_indel_haplotype_combines_nearby_snv() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 1,
            end: 9,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAATTT",
        };
        let variants = vec![
            crate::io::VcfPosition {
                position: 4,
                ref_allele: "A".to_string(),
                alt_allele: "C".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            crate::io::VcfPosition {
                position: 6,
                ref_allele: "A".to_string(),
                alt_allele: "AT".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ];

        let phased = build_phased_indel_haplotype_variants(
            &gene,
            &variants,
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert_eq!(phased.len(), 1);
        assert_eq!(phased[0].variant_type, VariantType::Indel);
        assert_eq!(phased[0].event_class.as_deref(), Some("complex_indel"));
        assert_eq!(phased[0].positions, vec![4]);
        assert_eq!(phased[0].ref_bases, vec!["AAA".to_string()]);
        assert_eq!(phased[0].base_changes, vec!["CAAT".to_string()]);
        assert_eq!(
            phased[0].event_components,
            vec!["SNV:4:A>C".to_string(), "INS:6:+T".to_string()]
        );
    }

    #[test]
    fn test_build_phased_indel_haplotype_preserves_deletion_component_coordinate() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 1,
            end: 9,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAATTT",
        };
        let variants = vec![
            crate::io::VcfPosition {
                position: 4,
                ref_allele: "A".to_string(),
                alt_allele: "C".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            crate::io::VcfPosition {
                position: 5,
                ref_allele: "AA".to_string(),
                alt_allele: "A".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ];

        let phased = build_phased_indel_haplotype_variants(
            &gene,
            &variants,
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert_eq!(phased.len(), 1);
        assert_eq!(phased[0].event_class.as_deref(), Some("complex_indel"));
        assert_eq!(phased[0].positions, vec![4]);
        assert_eq!(phased[0].ref_bases, vec!["AAA".to_string()]);
        assert_eq!(phased[0].base_changes, vec!["CA".to_string()]);
        assert_eq!(
            phased[0].event_components,
            vec!["SNV:4:A>C".to_string(), "DEL:6:A".to_string()]
        );
    }

    #[test]
    fn test_build_phased_indel_haplotype_combines_two_indels() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 1,
            end: 9,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAATTT",
        };
        let variants = vec![
            crate::io::VcfPosition {
                position: 4,
                ref_allele: "A".to_string(),
                alt_allele: "AG".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            crate::io::VcfPosition {
                position: 5,
                ref_allele: "AA".to_string(),
                alt_allele: "A".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ];

        let phased = build_phased_indel_haplotype_variants(
            &gene,
            &variants,
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert_eq!(phased.len(), 1);
        assert_eq!(phased[0].event_class.as_deref(), Some("complex_indel"));
        assert_eq!(phased[0].positions, vec![4]);
        assert_eq!(phased[0].ref_bases, vec!["AAA".to_string()]);
        assert_eq!(phased[0].base_changes, vec!["AGA".to_string()]);
        assert_eq!(
            phased[0].event_components,
            vec!["INS:4:+G".to_string(), "DEL:6:A".to_string()]
        );
    }

    #[test]
    fn test_build_phased_indel_haplotype_ignores_distant_snv() {
        let gene = Gene {
            name: "cds".to_string(),
            start: 1,
            end: 12,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
        };
        let reference = crate::io::Reference {
            sequence: "ATGAAATTTCCC",
        };
        let variants = vec![
            crate::io::VcfPosition {
                position: 10,
                ref_allele: "C".to_string(),
                alt_allele: "G".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            crate::io::VcfPosition {
                position: 6,
                ref_allele: "A".to_string(),
                alt_allele: "AT".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ];

        let phased = build_phased_indel_haplotype_variants(
            &gene,
            &variants,
            &reference,
            "chr1",
            crate::genetic_code::GeneticCode::default(),
        );

        assert!(phased.is_empty());
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
            transcript_id: None,
            cds_segments: Vec::new(),
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
            transcript_id: None,
            cds_segments: Vec::new(),
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
            transcript_id: None,
            cds_segments: Vec::new(),
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
            transcript_id: None,
            cds_segments: Vec::new(),
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
            transcript_id: None,
            cds_segments: Vec::new(),
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
            transcript_id: None,
            cds_segments: Vec::new(),
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
            transcript_id: None,
            cds_segments: Vec::new(),
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
