//! BAM read counting: region observation cache, per-position and MNV
//! haplotype read support with strand-specific metrics.

use crate::error::{AppError, AppResult};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::Header;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;

#[derive(Debug, Clone)]
pub struct ReadCountSummary {
    pub snp_counts: Vec<usize>,
    pub mnv_count: usize,
    pub total_reads: Vec<usize>,
    pub total_forward_reads: Vec<usize>,
    pub total_reverse_reads: Vec<usize>,
    pub snp_forward_counts: Vec<usize>,
    pub snp_reverse_counts: Vec<usize>,
    pub mnv_forward_count: usize,
    pub mnv_reverse_count: usize,
    pub mnv_total_reads: usize,
    pub mnv_total_forward_reads: usize,
    pub mnv_total_reverse_reads: usize,
}

#[derive(Debug, Clone)]
pub struct RegionObservationCache {
    index_by_position: HashMap<usize, usize>,
    reads: Vec<CachedReadObservation>,
}

#[derive(Debug, Clone)]
struct CachedReadObservation {
    key: Rc<ReadKey>,
    is_reverse: bool,
    observations: Vec<Option<(char, u8)>>,
}

#[derive(Debug, Clone)]
struct MultiReadSupport {
    snp_support: Vec<bool>,
    snp_forward: Vec<bool>,
    snp_reverse: Vec<bool>,
}

impl MultiReadSupport {
    fn new(size: usize) -> Self {
        Self {
            snp_support: vec![false; size],
            snp_forward: vec![false; size],
            snp_reverse: vec![false; size],
        }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct ReadKey {
    qname: Vec<u8>,
    is_first_segment: bool,
    is_last_segment: bool,
    start_pos: i64,
}

fn build_read_key(rec: &bam::Record) -> ReadKey {
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

fn get_query_pos(rec: &bam::Record, pos: usize) -> Option<usize> {
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

fn increment_directional_count(
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

    // Build region query string: chrom:start-end (1-based, inclusive)
    let region_str = format!("{chrom}:{region_start}-{region_end}");
    let region: noodles::core::Region = region_str
        .parse()
        .map_err(|e| AppError::validation(format!("Invalid region '{region_str}': {e}")))?;

    let mut query = bam_reader
        .query(header, &region)
        .map_err(|e| AppError::validation(format!("BAM query failed for {region_str}: {e}")))?;

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

        let observations = normalized_positions
            .iter()
            .copied()
            .map(|pos| {
                get_query_pos(&record, pos).and_then(|idx| {
                    if idx < seq_len {
                        let base_byte: u8 = seq.iter().nth(idx)?;
                        let base = base_byte as char;
                        let q: u8 = qual.iter().nth(idx).unwrap_or(0);
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
    let mut total_reads = vec![0; positions.len()];
    let mut total_forward_reads = vec![0; positions.len()];
    let mut total_reverse_reads = vec![0; positions.len()];
    let mut snp_forward_counts = vec![0; positions.len()];
    let mut snp_reverse_counts = vec![0; positions.len()];
    let mut mnv_count = 0usize;
    let mut mnv_forward_count = 0usize;
    let mut mnv_reverse_count = 0usize;
    let mut unique_mnv_total: HashSet<Rc<ReadKey>> = HashSet::new();
    let mut unique_mnv_total_forward: HashSet<Rc<ReadKey>> = HashSet::new();
    let mut unique_mnv_total_reverse: HashSet<Rc<ReadKey>> = HashSet::new();

    let mut unique_snp: Vec<HashSet<Rc<ReadKey>>> =
        (0..positions.len()).map(|_| HashSet::new()).collect();
    let mut unique_total: Vec<HashSet<Rc<ReadKey>>> =
        (0..positions.len()).map(|_| HashSet::new()).collect();
    let mut unique_total_forward: Vec<HashSet<Rc<ReadKey>>> =
        (0..positions.len()).map(|_| HashSet::new()).collect();
    let mut unique_total_reverse: Vec<HashSet<Rc<ReadKey>>> =
        (0..positions.len()).map(|_| HashSet::new()).collect();
    let mut per_read_multi_support: HashMap<Rc<ReadKey>, MultiReadSupport> = HashMap::new();

    for cached_read in &cache.reads {
        let observations = requested_indices
            .iter()
            .map(|idx| cached_read.observations[*idx])
            .collect::<Vec<_>>();

        for (idx, observation) in observations.iter().enumerate() {
            if let Some((_, qual)) = observation {
                if *qual >= min_phred_quality && unique_total[idx].insert(cached_read.key.clone()) {
                    total_reads[idx] += 1;
                }
                if *qual >= min_phred_quality {
                    if cached_read.is_reverse {
                        if unique_total_reverse[idx].insert(cached_read.key.clone()) {
                            total_reverse_reads[idx] += 1;
                        }
                    } else if unique_total_forward[idx].insert(cached_read.key.clone()) {
                        total_forward_reads[idx] += 1;
                    }
                }
            }
        }

        if positions.len() == 1 {
            if observations[0]
                .map(|(base, qual)| {
                    qual >= min_phred_quality && base.eq_ignore_ascii_case(&alt_chars[0])
                })
                .unwrap_or(false)
                && unique_snp[0].insert(cached_read.key.clone())
            {
                snp_counts[0] += 1;
                if cached_read.is_reverse {
                    snp_reverse_counts[0] += 1;
                } else {
                    snp_forward_counts[0] += 1;
                }
            }
            continue;
        }

        let covers_all_positions = observations
            .iter()
            .all(|value| matches!(value, Some((_, qual)) if *qual >= min_phred_quality));
        if covers_all_positions {
            unique_mnv_total.insert(cached_read.key.clone());
            if cached_read.is_reverse {
                unique_mnv_total_reverse.insert(cached_read.key.clone());
            } else {
                unique_mnv_total_forward.insert(cached_read.key.clone());
            }
        }

        let support = per_read_multi_support
            .entry(cached_read.key.clone())
            .or_insert_with(|| MultiReadSupport::new(positions.len()));
        for (idx, observation) in observations.iter().enumerate() {
            if observation
                .map(|(base, qual)| {
                    qual >= min_phred_quality && base.eq_ignore_ascii_case(&alt_chars[idx])
                })
                .unwrap_or(false)
            {
                support.snp_support[idx] = true;
                if cached_read.is_reverse {
                    support.snp_reverse[idx] = true;
                } else {
                    support.snp_forward[idx] = true;
                }
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
    })
}

pub fn count_reads_per_position(
    bam_reader: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: &Header,
    chrom: &str,
    positions: &[usize],
    alt_bases: &[String],
    min_phred_quality: u8,
    min_mapq: u8,
) -> AppResult<ReadCountSummary> {
    let min_pos = positions
        .iter()
        .copied()
        .min()
        .ok_or_else(|| AppError::validation("No positions provided for read counting"))?;
    let max_pos = positions
        .iter()
        .copied()
        .max()
        .ok_or_else(|| AppError::validation("No positions provided for read counting"))?;

    let cache = build_region_observation_cache(
        bam_reader, header, chrom, min_pos, max_pos, positions, min_mapq,
    )?;
    count_reads_from_cache(&cache, positions, alt_bases, min_phred_quality)
}

#[cfg(test)]
mod tests {
    use super::{increment_directional_count, normalize_positions};

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

    #[test]
    fn test_normalize_positions_deduplicates_and_sorts() {
        let normalized = normalize_positions(&[5, 2, 5, 1]).expect("should normalize");
        assert_eq!(normalized, vec![1, 2, 5]);
    }
}
