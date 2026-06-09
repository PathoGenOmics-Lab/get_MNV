//! Per-read BAM observation primitives: CIGAR-aware base lookup, observed
//! allele reconstruction over a reference span, component support checks, and
//! region construction.

use crate::error::{AppError, AppResult};
use crate::variants::{AlleleComponent, AlleleComponentKind};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub(super) struct ReadKey {
    qname: Vec<u8>,
    is_first_segment: bool,
    is_last_segment: bool,
    start_pos: i64,
}

pub(super) fn build_read_key(rec: &bam::Record) -> ReadKey {
    let qname = rec
        .name()
        .map(|n| <_ as AsRef<[u8]>>::as_ref(&n).to_vec())
        .unwrap_or_default();
    let flags = rec.flags();
    let start_pos = rec
        .alignment_start()
        .and_then(|p| p.ok())
        .map(|p| {
            let pos: usize = p.into();
            pos as i64 - 1
        })
        .unwrap_or(0);
    ReadKey {
        qname,
        is_first_segment: flags.is_first_segment(),
        is_last_segment: flags.is_last_segment(),
        start_pos,
    }
}

pub(super) fn get_query_pos(rec: &bam::Record, pos: usize) -> Option<usize> {
    // pos is 1-based
    let target_pos = (pos - 1) as i64; // 0-based reference position
    let rec_start: i64 = rec
        .alignment_start()
        .and_then(|p| p.ok())
        .map(|p| {
            let v: usize = p.into();
            v as i64 - 1
        })
        .unwrap_or(0);
    let mut ref_pos = rec_start;
    let mut query_pos = 0usize;

    let cigar = rec.cigar();

    for op_result in cigar.iter() {
        let op: noodles::sam::alignment::record::cigar::Op = match op_result {
            Ok(o) => o,
            Err(_) => return None,
        };
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let len_i64 = len as i64;
                if target_pos >= ref_pos && target_pos < ref_pos + len_i64 {
                    return Some(query_pos + (target_pos - ref_pos) as usize);
                }
                ref_pos += len_i64;
                query_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                query_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                let len_i64 = len as i64;
                if target_pos >= ref_pos && target_pos < ref_pos + len_i64 {
                    return None;
                }
                ref_pos += len_i64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    None
}

/// Phred quality of the read base aligned to 1-based reference `position`, or
/// `None` if the position is not covered by an aligned (match) base in this read
/// (deletion, skip, soft-clip, or outside the alignment).
pub(super) fn anchor_base_quality(rec: &bam::Record, position: usize) -> Option<u8> {
    let qidx = get_query_pos(rec, position)?;
    let qual = rec.quality_scores();
    if qual.iter().next().is_none() {
        // SAM `QUAL=*`: the read carries no per-base qualities (noodles exposes
        // this as an empty slice). Treat the base as passing, mirroring the
        // missing-MAPQ default of 255, instead of substituting quality 0 and
        // silently failing every quality gate.
        return Some(u8::MAX);
    }
    // Defensive: a malformed CIGAR whose query length exceeds the stored
    // sequence could yield an index past the quality scores. Treat that as "no
    // observation" rather than reading a bogus quality.
    if qidx >= rec.sequence().len() {
        return None;
    }
    Some(qual.iter().nth(qidx).unwrap_or(0))
}

#[derive(Debug, Clone)]
pub(super) struct ObservedAllele {
    pub(super) allele: String,
    pub(super) min_quality: u8,
    bases_by_position: HashMap<usize, char>,
    insertions_after: HashMap<usize, String>,
    deleted_positions: HashSet<usize>,
}

pub(super) fn observed_allele_for_ref_span(
    rec: &bam::Record,
    pos: usize,
    ref_len: usize,
) -> Option<ObservedAllele> {
    if pos == 0 || ref_len == 0 {
        return None;
    }
    let target_start = (pos - 1) as i64;
    let target_end = target_start + ref_len as i64;
    let rec_start: i64 = rec
        .alignment_start()
        .and_then(|p| p.ok())
        .map(|p| {
            let v: usize = p.into();
            v as i64 - 1
        })
        .unwrap_or(0);

    let seq = rec.sequence();
    let qual = rec.quality_scores();
    let seq_len = seq.len();
    // SAM `QUAL=*` (no per-base qualities) is exposed as an empty slice; in that
    // case treat every base as max quality so missing quality data does not drop
    // the read from numerator and denominator (mirrors the missing-MAPQ default).
    let qual_present = qual.iter().next().is_some();
    let mut bases: Vec<Option<(char, u8)>> = vec![None; ref_len];
    let mut deleted: Vec<bool> = vec![false; ref_len];
    let mut insertions_after: Vec<Vec<(char, u8)>> = vec![Vec::new(); ref_len];

    let mut ref_pos = rec_start;
    let mut query_pos = 0usize;

    for op_result in rec.cigar().iter() {
        let op: noodles::sam::alignment::record::cigar::Op = op_result.ok()?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for offset in 0..len {
                    let current_ref = ref_pos + offset as i64;
                    if current_ref >= target_start && current_ref < target_end {
                        let idx = (current_ref - target_start) as usize;
                        let qidx = query_pos + offset;
                        if qidx < seq_len {
                            let base = seq.iter().nth(qidx)? as char;
                            let q = if qual_present {
                                qual.iter().nth(qidx).unwrap_or(0)
                            } else {
                                u8::MAX
                            };
                            bases[idx] = Some((base, q));
                        }
                    }
                }
                ref_pos += len as i64;
                query_pos += len;
            }
            Kind::Insertion => {
                let anchor_ref = ref_pos - 1;
                if anchor_ref >= target_start && anchor_ref < target_end {
                    let idx = (anchor_ref - target_start) as usize;
                    for offset in 0..len {
                        let qidx = query_pos + offset;
                        if qidx < seq_len {
                            let base = seq.iter().nth(qidx)? as char;
                            let q = if qual_present {
                                qual.iter().nth(qidx).unwrap_or(0)
                            } else {
                                u8::MAX
                            };
                            insertions_after[idx].push((base, q));
                        }
                    }
                }
                query_pos += len;
            }
            Kind::SoftClip => {
                query_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                for offset in 0..len {
                    let current_ref = ref_pos + offset as i64;
                    if current_ref >= target_start && current_ref < target_end {
                        let idx = (current_ref - target_start) as usize;
                        deleted[idx] = true;
                    }
                }
                ref_pos += len as i64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    let mut allele = String::new();
    let mut min_quality = u8::MAX;
    let mut saw_observed_base = false;
    let mut bases_by_position = HashMap::new();
    let mut observed_insertions_after = HashMap::new();
    let mut deleted_positions = HashSet::new();

    for idx in 0..ref_len {
        let reference_pos = pos + idx;
        if let Some((base, q)) = bases[idx] {
            let base = base.to_ascii_uppercase();
            allele.push(base);
            bases_by_position.insert(reference_pos, base);
            min_quality = min_quality.min(q);
            saw_observed_base = true;
        } else if deleted[idx] {
            // Deleted reference bases contribute no query base to the allele.
            deleted_positions.insert(reference_pos);
        } else {
            return None;
        }

        let mut inserted = String::new();
        for (base, q) in &insertions_after[idx] {
            let base = base.to_ascii_uppercase();
            allele.push(base);
            inserted.push(base);
            min_quality = min_quality.min(*q);
            saw_observed_base = true;
        }
        if !inserted.is_empty() {
            observed_insertions_after.insert(reference_pos, inserted);
        }
    }

    if !saw_observed_base {
        return None;
    }
    Some(ObservedAllele {
        allele,
        min_quality,
        bases_by_position,
        insertions_after: observed_insertions_after,
        deleted_positions,
    })
}

pub(super) fn observed_supports_components(
    observed: &ObservedAllele,
    required_components: &[AlleleComponent],
) -> bool {
    required_components
        .iter()
        .all(|component| match component.kind {
            AlleleComponentKind::Snp => {
                let Some(observed_base) = observed.bases_by_position.get(&component.position)
                else {
                    return false;
                };
                component
                    .alt_allele
                    .chars()
                    .next()
                    .is_some_and(|alt| observed_base.eq_ignore_ascii_case(&alt))
            }
            AlleleComponentKind::Insertion => observed
                .insertions_after
                .get(&component.position)
                .is_some_and(|inserted| inserted.eq_ignore_ascii_case(&component.alt_allele)),
            AlleleComponentKind::Deletion => {
                let del_end = component.position + component.ref_allele.len().saturating_sub(1);
                (component.position..=del_end)
                    .all(|position| observed.deleted_positions.contains(&position))
            }
            AlleleComponentKind::Delins => true,
            AlleleComponentKind::Symbolic => false,
        })
}

pub(super) fn increment_directional_count(
    forward_supported: bool,
    reverse_supported: bool,
    forward_count: &mut usize,
    reverse_count: &mut usize,
) {
    match (forward_supported, reverse_supported) {
        (true, false) => *forward_count += 1,
        (false, true) => *reverse_count += 1,
        (true, true) => *forward_count += 1,
        (false, false) => {}
    }
}

/// Build a coordinate region for a BAM query using the structured API, so a
/// contig name containing ':' (e.g. HLA allele contigs) is not misparsed by the
/// `chrom:start-end` string form.
pub(super) fn build_region(
    chrom: &str,
    start: usize,
    end: usize,
) -> AppResult<noodles::core::Region> {
    use noodles::core::Position;
    let start_pos = Position::try_from(start)
        .map_err(|e| AppError::validation(format!("Invalid region start {chrom}:{start}: {e}")))?;
    let end_pos = Position::try_from(end)
        .map_err(|e| AppError::validation(format!("Invalid region end {chrom}:{end}: {e}")))?;
    Ok(noodles::core::Region::new(chrom, start_pos..=end_pos))
}

#[cfg(test)]
mod tests {
    use super::increment_directional_count;

    #[test]
    fn test_increment_directional_count_single_strand() {
        let mut forward = 0usize;
        let mut reverse = 0usize;

        increment_directional_count(true, false, &mut forward, &mut reverse);
        assert_eq!(forward, 1);
        assert_eq!(reverse, 0);

        increment_directional_count(false, true, &mut forward, &mut reverse);
        assert_eq!(forward, 1);
        assert_eq!(reverse, 1);
    }

    #[test]
    fn test_increment_directional_count_ambiguous_defaults_to_forward() {
        let mut forward = 0usize;
        let mut reverse = 0usize;
        increment_directional_count(true, true, &mut forward, &mut reverse);
        assert_eq!(forward, 1);
        assert_eq!(reverse, 0);
    }
}
