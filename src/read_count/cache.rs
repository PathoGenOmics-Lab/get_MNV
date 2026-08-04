//! Region observation cache and cache-based per-position / MNV haplotype read
//! counting with strand-specific metrics.

use crate::error::{AppError, AppResult};
use noodles::bam;
use noodles::sam::Header;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;

use super::observation::{
    build_molecule_key, build_read_key, build_region, get_query_pos, increment_directional_count,
    MoleculeKey, ReadKey,
};
use super::ReadCountSummary;

#[derive(Debug, Clone)]
pub struct RegionObservationCache {
    index_by_position: HashMap<usize, usize>,
    reads: Vec<CachedReadObservation>,
}

#[derive(Debug, Clone)]
struct CachedReadObservation {
    key: Rc<ReadKey>,
    /// The molecule this record belongs to when mates are treated as one. Kept
    /// alongside the per-record key so the cache serves both counting modes.
    fragment: Rc<MoleculeKey>,
    is_reverse: bool,
    observations: Vec<Option<(char, u8)>>,
}

/// One molecule's view of the requested positions: a single record, or a
/// paired-end fragment's two mates merged.
struct MoleculeObservation {
    identity: Rc<MoleculeKey>,
    has_forward: bool,
    has_reverse: bool,
    observations: Vec<Option<(char, u8)>>,
}

/// Merge what two mates saw at one position.
///
/// Where only one mate reaches the position, it speaks alone. Where both do and
/// agree, the better base quality stands. Where both do and disagree, one of
/// them is wrong and there is no way to tell which, so the molecule is treated
/// as not having observed that position at all.
fn merge_observation(left: Option<(char, u8)>, right: Option<(char, u8)>) -> Option<(char, u8)> {
    match (left, right) {
        (None, other) | (other, None) => other,
        (Some((left_base, left_quality)), Some((right_base, right_quality))) => {
            if left_base.eq_ignore_ascii_case(&right_base) {
                Some((left_base, left_quality.max(right_quality)))
            } else {
                None
            }
        }
    }
}

/// Collapse cached records into one entry per molecule, preserving the order
/// molecules were first seen so counting stays deterministic.
fn molecules(reads: &[CachedReadObservation], pair_aware: bool) -> Vec<MoleculeObservation> {
    let mut order: Vec<Rc<MoleculeKey>> = Vec::new();
    let mut by_identity: HashMap<Rc<MoleculeKey>, usize> = HashMap::new();
    let mut merged: Vec<MoleculeObservation> = Vec::new();

    for read in reads {
        let identity = if pair_aware {
            read.fragment.clone()
        } else {
            Rc::new(MoleculeKey::Record((*read.key).clone()))
        };
        match by_identity.get(&identity) {
            Some(&index) => {
                let molecule: &mut MoleculeObservation = &mut merged[index];
                molecule.has_forward |= !read.is_reverse;
                molecule.has_reverse |= read.is_reverse;
                for (slot, observation) in molecule.observations.iter_mut().zip(&read.observations)
                {
                    *slot = merge_observation(*slot, *observation);
                }
            }
            None => {
                by_identity.insert(identity.clone(), merged.len());
                order.push(identity.clone());
                merged.push(MoleculeObservation {
                    identity,
                    has_forward: !read.is_reverse,
                    has_reverse: read.is_reverse,
                    observations: read.observations.clone(),
                });
            }
        }
    }
    debug_assert_eq!(order.len(), merged.len());
    merged
}

#[derive(Debug, Clone)]
struct MultiReadSupport {
    snp_support: Vec<bool>,
    snp_forward: Vec<bool>,
    snp_reverse: Vec<bool>,
    /// Whether the read observes every requested position at sufficient base
    /// quality. Only such a read can testify about linkage: one that ends
    /// between two positions carries no information about the pair.
    spans_all: bool,
}

impl MultiReadSupport {
    fn new(size: usize) -> Self {
        Self {
            snp_support: vec![false; size],
            snp_forward: vec![false; size],
            snp_reverse: vec![false; size],
            spans_all: false,
        }
    }
}

fn normalize_positions(positions: &[usize]) -> AppResult<Vec<usize>> {
    if positions.is_empty() {
        return Err(AppError::validation(
            "No positions provided for read counting",
        ));
    }
    if positions.contains(&0) {
        return Err(AppError::validation(
            "Invalid position 0 provided for read counting",
        ));
    }
    let mut normalized = positions.to_vec();
    normalized.sort_unstable();
    normalized.dedup();
    Ok(normalized)
}

pub fn build_region_observation_cache(
    bam_reader: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: &Header,
    chrom: &str,
    region_start: usize,
    region_end: usize,
    target_positions: &[usize],
    min_mapq: u8,
) -> AppResult<RegionObservationCache> {
    if region_start == 0 || region_end == 0 || region_start > region_end {
        return Err(AppError::validation(format!(
            "Invalid region bounds for cache: start={region_start}, end={region_end}"
        )));
    }

    let normalized_positions = normalize_positions(target_positions)?;
    if normalized_positions
        .iter()
        .any(|pos| *pos < region_start || *pos > region_end)
    {
        return Err(AppError::validation(format!(
            "Some target positions are outside region {chrom}:{region_start}-{region_end}"
        )));
    }

    let index_by_position = normalized_positions
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, pos)| (pos, idx))
        .collect::<HashMap<_, _>>();

    let mut reads: Vec<CachedReadObservation> = Vec::new();

    // Structured region (1-based, inclusive) — avoids misparsing ':' in contig names.
    let region = build_region(chrom, region_start, region_end)?;
    let mut query = bam_reader.query(header, &region).map_err(|e| {
        AppError::validation(format!(
            "BAM query failed for {chrom}:{region_start}-{region_end}: {e}"
        ))
    })?;

    let mut record = bam::Record::default();
    while query
        .read_record(&mut record)
        .map_err(|e| AppError::validation(format!("BAM read error: {e}")))?
        != 0
    {
        let flags = record.flags();
        if flags.is_duplicate() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let mapq = record
            .mapping_quality()
            .map(|q: noodles::sam::alignment::record::MappingQuality| q.get())
            .unwrap_or(255);
        if mapq < min_mapq {
            continue;
        }

        let seq = record.sequence();
        let qual = record.quality_scores();
        let seq_len = seq.len();
        // SAM `QUAL=*` (no per-base qualities) is exposed as an empty slice;
        // treat such bases as max quality so a quality-less BAM does not lose all
        // read support (mirrors the missing-MAPQ default of 255).
        let qual_present = qual.iter().next().is_some();

        let observations = normalized_positions
            .iter()
            .copied()
            .map(|pos| {
                get_query_pos(&record, pos).and_then(|idx| {
                    if idx < seq_len {
                        let base_byte: u8 = seq.iter().nth(idx)?;
                        let base = base_byte as char;
                        let q: u8 = if qual_present {
                            qual.iter().nth(idx).unwrap_or(0)
                        } else {
                            u8::MAX
                        };
                        Some((base, q))
                    } else {
                        None
                    }
                })
            })
            .collect::<Vec<_>>();

        if observations.iter().all(Option::is_none) {
            continue;
        }

        reads.push(CachedReadObservation {
            key: Rc::new(build_read_key(&record)),
            fragment: Rc::new(build_molecule_key(&record, true)),
            is_reverse: flags.is_reverse_complemented(),
            observations,
        });
    }

    Ok(RegionObservationCache {
        index_by_position,
        reads,
    })
}

fn resolve_requested_indices(
    cache: &RegionObservationCache,
    positions: &[usize],
) -> AppResult<Vec<usize>> {
    positions
        .iter()
        .copied()
        .map(|pos| {
            cache.index_by_position.get(&pos).copied().ok_or_else(|| {
                AppError::validation(format!("Position {pos} missing from observation cache"))
            })
        })
        .collect()
}

pub fn count_reads_from_cache(
    cache: &RegionObservationCache,
    positions: &[usize],
    alt_bases: &[String],
    min_phred_quality: u8,
    pair_aware: bool,
) -> AppResult<ReadCountSummary> {
    if positions.is_empty() {
        return Err(AppError::validation(
            "No positions provided for read counting",
        ));
    }
    if positions.len() != alt_bases.len() {
        return Err(AppError::validation(format!(
            "Positions/ALT mismatch: {} positions vs {} ALT alleles",
            positions.len(),
            alt_bases.len()
        )));
    }

    let requested_indices = resolve_requested_indices(cache, positions)?;
    let alt_chars: Vec<char> = positions
        .iter()
        .copied()
        .zip(alt_bases.iter())
        .map(|(pos, alt)| {
            alt.chars().next().ok_or_else(|| {
                AppError::validation(format!(
                    "Empty ALT allele provided for read counting at position {pos}"
                ))
            })
        })
        .collect::<AppResult<Vec<_>>>()?;

    let mut snp_counts = vec![0; positions.len()];
    let mut snp_only_informative_counts = vec![0; positions.len()];
    let mut total_reads = vec![0; positions.len()];
    let mut total_forward_reads = vec![0; positions.len()];
    let mut total_reverse_reads = vec![0; positions.len()];
    let mut snp_forward_counts = vec![0; positions.len()];
    let mut snp_reverse_counts = vec![0; positions.len()];
    let mut mnv_count = 0usize;
    let mut mnv_forward_count = 0usize;
    let mut mnv_reverse_count = 0usize;
    let mut unique_mnv_total: HashSet<Rc<MoleculeKey>> = HashSet::new();
    let mut unique_mnv_total_forward: HashSet<Rc<MoleculeKey>> = HashSet::new();
    let mut unique_mnv_total_reverse: HashSet<Rc<MoleculeKey>> = HashSet::new();

    let mut unique_snp: Vec<HashSet<Rc<MoleculeKey>>> =
        (0..positions.len()).map(|_| HashSet::new()).collect();
    let mut unique_total: Vec<HashSet<Rc<MoleculeKey>>> =
        (0..positions.len()).map(|_| HashSet::new()).collect();
    let mut unique_total_forward: Vec<HashSet<Rc<MoleculeKey>>> =
        (0..positions.len()).map(|_| HashSet::new()).collect();
    let mut unique_total_reverse: Vec<HashSet<Rc<MoleculeKey>>> =
        (0..positions.len()).map(|_| HashSet::new()).collect();
    let mut per_read_multi_support: HashMap<Rc<MoleculeKey>, MultiReadSupport> = HashMap::new();

    for molecule in &molecules(&cache.reads, pair_aware) {
        let observations = requested_indices
            .iter()
            .map(|idx| molecule.observations[*idx])
            .collect::<Vec<_>>();

        for (idx, observation) in observations.iter().enumerate() {
            if let Some((_, qual)) = observation {
                if *qual < min_phred_quality {
                    continue;
                }
                if unique_total[idx].insert(molecule.identity.clone()) {
                    total_reads[idx] += 1;
                }
                if molecule.has_forward
                    && unique_total_forward[idx].insert(molecule.identity.clone())
                {
                    total_forward_reads[idx] += 1;
                }
                if molecule.has_reverse
                    && unique_total_reverse[idx].insert(molecule.identity.clone())
                {
                    total_reverse_reads[idx] += 1;
                }
            }
        }

        if positions.len() == 1 {
            if observations[0]
                .map(|(base, qual)| {
                    qual >= min_phred_quality && base.eq_ignore_ascii_case(&alt_chars[0])
                })
                .unwrap_or(false)
                && unique_snp[0].insert(molecule.identity.clone())
            {
                snp_counts[0] += 1;
                increment_directional_count(
                    molecule.has_forward,
                    molecule.has_reverse,
                    &mut snp_forward_counts[0],
                    &mut snp_reverse_counts[0],
                );
            }
            continue;
        }

        let covers_all_positions = observations
            .iter()
            .all(|value| matches!(value, Some((_, qual)) if *qual >= min_phred_quality));
        if covers_all_positions {
            unique_mnv_total.insert(molecule.identity.clone());
            if molecule.has_forward {
                unique_mnv_total_forward.insert(molecule.identity.clone());
            }
            if molecule.has_reverse {
                unique_mnv_total_reverse.insert(molecule.identity.clone());
            }
        }

        let support = per_read_multi_support
            .entry(molecule.identity.clone())
            .or_insert_with(|| MultiReadSupport::new(positions.len()));
        support.spans_all |= covers_all_positions;
        for (idx, observation) in observations.iter().enumerate() {
            if observation
                .map(|(base, qual)| {
                    qual >= min_phred_quality && base.eq_ignore_ascii_case(&alt_chars[idx])
                })
                .unwrap_or(false)
            {
                support.snp_support[idx] = true;
                support.snp_forward[idx] |= molecule.has_forward;
                support.snp_reverse[idx] |= molecule.has_reverse;
            }
        }
    }

    if positions.len() > 1 {
        for support in per_read_multi_support.values() {
            if support.snp_support.iter().all(|is_supported| *is_supported) {
                mnv_count += 1;
                let has_forward = support.snp_forward.iter().any(|is_forward| *is_forward);
                let has_reverse = support.snp_reverse.iter().any(|is_reverse| *is_reverse);
                increment_directional_count(
                    has_forward,
                    has_reverse,
                    &mut mnv_forward_count,
                    &mut mnv_reverse_count,
                );
            } else {
                for (idx, is_supported) in support.snp_support.iter().enumerate() {
                    if *is_supported {
                        snp_counts[idx] += 1;
                        if support.spans_all {
                            snp_only_informative_counts[idx] += 1;
                        }
                        increment_directional_count(
                            support.snp_forward[idx],
                            support.snp_reverse[idx],
                            &mut snp_forward_counts[idx],
                            &mut snp_reverse_counts[idx],
                        );
                    }
                }
            }
        }
    }

    let mnv_total_reads = if positions.len() > 1 {
        unique_mnv_total.len()
    } else {
        total_reads.first().copied().unwrap_or(0)
    };
    let mnv_total_forward_reads = if positions.len() > 1 {
        unique_mnv_total_forward.len()
    } else {
        total_forward_reads.first().copied().unwrap_or(0)
    };
    let mnv_total_reverse_reads = if positions.len() > 1 {
        unique_mnv_total_reverse.len()
    } else {
        total_reverse_reads.first().copied().unwrap_or(0)
    };

    Ok(ReadCountSummary {
        snp_counts,
        mnv_count,
        total_reads,
        total_forward_reads,
        total_reverse_reads,
        snp_forward_counts,
        snp_reverse_counts,
        mnv_forward_count,
        mnv_reverse_count,
        mnv_total_reads,
        mnv_total_forward_reads,
        mnv_total_reverse_reads,
        snp_only_informative_counts,
    })
}

#[cfg(test)]
mod tests {
    use super::normalize_positions;

    #[test]
    fn test_normalize_positions_deduplicates_and_sorts() {
        let normalized = normalize_positions(&[5, 2, 5, 1]).expect("should normalize");
        assert_eq!(normalized, vec![1, 2, 5]);
    }
}
