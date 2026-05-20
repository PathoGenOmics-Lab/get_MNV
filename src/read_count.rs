//! BAM read counting: region observation cache, per-position and MNV
//! haplotype read support with strand-specific metrics.

use crate::error::{AppError, AppResult};
use crate::variants::{AlleleComponent, AlleleComponentKind};
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

pub struct IndelReadCountRequest<'a> {
    pub chrom: &'a str,
    pub position: usize,
    pub ref_allele: &'a str,
    pub alt_allele: &'a str,
    pub required_components: &'a [AlleleComponent],
    pub min_phred_quality: u8,
    pub min_mapq: u8,
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

#[derive(Debug, Clone)]
struct ObservedAllele {
    allele: String,
    min_quality: u8,
    bases_by_position: HashMap<usize, char>,
    insertions_after: HashMap<usize, String>,
    deleted_positions: HashSet<usize>,
}

fn observed_allele_for_ref_span(
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
                            let q = qual.iter().nth(qidx).unwrap_or(0);
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
                            let q = qual.iter().nth(qidx).unwrap_or(0);
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

fn observed_supports_components(
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

pub fn count_indel_reads(
    bam_reader: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: &Header,
    request: IndelReadCountRequest<'_>,
) -> AppResult<ReadCountSummary> {
    let IndelReadCountRequest {
        chrom,
        position,
        ref_allele,
        alt_allele,
        required_components,
        min_phred_quality,
        min_mapq,
    } = request;

    if position == 0 || ref_allele.is_empty() {
        return Err(AppError::validation(format!(
            "Invalid indel allele for read counting at {chrom}:{position} REF='{ref_allele}' ALT='{alt_allele}'"
        )));
    }

    let end = position + ref_allele.len().saturating_sub(1);
    let region_str = format!("{chrom}:{position}-{end}");
    let region: noodles::core::Region = region_str
        .parse()
        .map_err(|e| AppError::validation(format!("Invalid region '{region_str}': {e}")))?;
    let mut query = bam_reader
        .query(header, &region)
        .map_err(|e| AppError::validation(format!("BAM query failed for {region_str}: {e}")))?;

    let mut unique_total: HashSet<Rc<ReadKey>> = HashSet::new();
    let mut unique_total_forward: HashSet<Rc<ReadKey>> = HashSet::new();
    let mut unique_total_reverse: HashSet<Rc<ReadKey>> = HashSet::new();
    let mut unique_alt: HashSet<Rc<ReadKey>> = HashSet::new();
    let mut unique_alt_forward: HashSet<Rc<ReadKey>> = HashSet::new();
    let mut unique_alt_reverse: HashSet<Rc<ReadKey>> = HashSet::new();

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

        let Some(observed) = observed_allele_for_ref_span(&record, position, ref_allele.len())
        else {
            continue;
        };
        if observed.min_quality < min_phred_quality {
            continue;
        }

        let key = Rc::new(build_read_key(&record));
        let is_reverse = flags.is_reverse_complemented();
        unique_total.insert(key.clone());
        if is_reverse {
            unique_total_reverse.insert(key.clone());
        } else {
            unique_total_forward.insert(key.clone());
        }

        if observed.allele.eq_ignore_ascii_case(alt_allele)
            && observed_supports_components(&observed, required_components)
        {
            unique_alt.insert(key.clone());
            if is_reverse {
                unique_alt_reverse.insert(key);
            } else {
                unique_alt_forward.insert(key);
            }
        }
    }

    let alt_count = unique_alt.len();
    let alt_forward = unique_alt_forward.len();
    let alt_reverse = unique_alt_reverse.len();
    let total = unique_total.len();
    let total_forward = unique_total_forward.len();
    let total_reverse = unique_total_reverse.len();

    Ok(ReadCountSummary {
        snp_counts: vec![alt_count],
        mnv_count: alt_count,
        total_reads: vec![total],
        total_forward_reads: vec![total_forward],
        total_reverse_reads: vec![total_reverse],
        snp_forward_counts: vec![alt_forward],
        snp_reverse_counts: vec![alt_reverse],
        mnv_forward_count: alt_forward,
        mnv_reverse_count: alt_reverse,
        mnv_total_reads: total,
        mnv_total_forward_reads: total_forward,
        mnv_total_reverse_reads: total_reverse,
    })
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
