//! Region observation cache and cache-based per-position / MNV haplotype read
//! counting with strand-specific metrics.

use crate::error::{AppError, AppResult};
use noodles::bam;
use noodles::sam::Header;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::rc::Rc;

#[cfg(test)]
use super::observation::test_read_key;
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

/// What one molecule saw at one position, and from which strands.
///
/// The strands are per position, not per molecule. A fragment's two mates often
/// land hundreds of bases apart, and the strand that read one variant says
/// nothing about the strand that read another; crediting a position to both
/// arms because the far mate happened to be on the other strand makes every
/// site look strand-balanced and leaves strand-bias detection unable to reject
/// anything.
#[derive(Debug, Clone, Copy)]
struct PositionObservation {
    base: char,
    quality: u8,
    from_forward: bool,
    from_reverse: bool,
}

/// One molecule's view of the requested positions: a single record, or a
/// paired-end fragment's two mates merged.
struct MoleculeObservation {
    identity: Rc<MoleculeKey>,
    observations: Vec<Option<PositionObservation>>,
}

impl MoleculeObservation {
    /// Whether a single strand established every one of these positions. For a
    /// multi-position row that is the strand-bias question: a haplotype pieced
    /// together from a forward mate and a reverse mate was confirmed by neither
    /// strand on its own.
    fn spans_all_on(&self, forward: bool, min_phred_quality: u8, indices: &[usize]) -> bool {
        indices.iter().all(|idx| {
            self.observations[*idx].is_some_and(|seen| {
                seen.quality >= min_phred_quality
                    && if forward {
                        seen.from_forward
                    } else {
                        seen.from_reverse
                    }
            })
        })
    }
}

/// Merge what two mates saw at one position.
///
/// Where only one mate reaches the position, it speaks alone. Where both do and
/// agree, the better base quality stands and the position is credited to both
/// strands. Where both do and disagree, one of them is wrong and there is no
/// way to tell which, so the molecule is treated as not having observed that
/// position at all.
fn merge_observation(
    left: Option<PositionObservation>,
    right: Option<PositionObservation>,
) -> Option<PositionObservation> {
    match (left, right) {
        (None, other) | (other, None) => other,
        (Some(left), Some(right)) => {
            if left.base.eq_ignore_ascii_case(&right.base) {
                Some(PositionObservation {
                    base: left.base,
                    quality: left.quality.max(right.quality),
                    from_forward: left.from_forward || right.from_forward,
                    from_reverse: left.from_reverse || right.from_reverse,
                })
            } else {
                None
            }
        }
    }
}

/// Collapse cached records into one entry per molecule, preserving the order
/// molecules were first seen so counting stays deterministic.
fn molecules(reads: &[CachedReadObservation], pair_aware: bool) -> Vec<MoleculeObservation> {
    let mut by_identity: HashMap<Rc<MoleculeKey>, usize> = HashMap::new();
    let mut merged: Vec<MoleculeObservation> = Vec::new();

    for read in reads {
        let identity = if pair_aware {
            read.fragment.clone()
        } else {
            Rc::new(MoleculeKey::Record((*read.key).clone()))
        };
        // The strand travels with each observation, so it only ever applies to
        // the positions this record actually read.
        let stamped = read
            .observations
            .iter()
            .map(|observation| {
                observation.map(|(base, quality)| PositionObservation {
                    base,
                    quality,
                    from_forward: !read.is_reverse,
                    from_reverse: read.is_reverse,
                })
            })
            .collect::<Vec<_>>();
        match by_identity.get(&identity) {
            Some(&index) => {
                let molecule: &mut MoleculeObservation = &mut merged[index];
                for (slot, observation) in molecule.observations.iter_mut().zip(&stamped) {
                    *slot = merge_observation(*slot, *observation);
                }
            }
            None => {
                by_identity.insert(identity.clone(), merged.len());
                merged.push(MoleculeObservation {
                    identity,
                    observations: stamped,
                });
            }
        }
    }
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

/// One described record for [`cache_from_molecules`]: the fragment it belongs
/// to, whether it is on the reverse strand, and what it read at each requested
/// position.
#[cfg(test)]
pub(super) type DescribedRecord<'a> = (&'a str, bool, Vec<Option<(char, u8)>>);

/// Build a cache directly from described molecules, without a BAM.
///
/// The counting this feeds is where the layer's bugs have lived, and it was
/// reachable only through an indexed BAM, which is why it had almost no unit
/// tests and was excluded from mutation testing. Each entry is
/// `(fragment_name, is_reverse, bases)`, one base per requested position, and
/// two entries sharing a name are the two mates of one molecule.
#[cfg(test)]
pub(super) fn cache_from_molecules(
    positions: &[usize],
    records: &[DescribedRecord<'_>],
) -> RegionObservationCache {
    let index_by_position = positions
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, pos)| (pos, idx))
        .collect::<HashMap<_, _>>();
    let reads = records
        .iter()
        .enumerate()
        .map(
            |(index, (name, is_reverse, observations))| CachedReadObservation {
                key: Rc::new(test_read_key(index)),
                fragment: Rc::new(MoleculeKey::Fragment(name.as_bytes().to_vec())),
                is_reverse: *is_reverse,
                observations: observations.clone(),
            },
        )
        .collect();
    RegionObservationCache {
        index_by_position,
        reads,
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

    // Structured region (1-based, inclusive); avoids misparsing ':' in contig names.
    let region = build_region(chrom, region_start, region_end)?;
    let mut query = bam_reader.query(header, &region).map_err(|e| {
        AppError::validation(format!(
            "BAM query failed for {chrom}:{region_start}-{region_end}: {e}"
        ))
    })?;

    let mut record_index = 0usize;
    let mut record = bam::Record::default();
    while query
        .read_record(&mut record)
        .map_err(|e| AppError::validation(format!("BAM read error: {e}")))?
        != 0
    {
        record_index += 1;
        let flags = record.flags();
        if flags.is_duplicate()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_qc_fail()
        {
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
            fragment: Rc::new(build_molecule_key(&record, true, record_index)),
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
            if let Some(seen) = observation {
                if seen.quality < min_phred_quality {
                    continue;
                }
                if unique_total[idx].insert(molecule.identity.clone()) {
                    total_reads[idx] += 1;
                }
                if seen.from_forward && unique_total_forward[idx].insert(molecule.identity.clone())
                {
                    total_forward_reads[idx] += 1;
                }
                if seen.from_reverse && unique_total_reverse[idx].insert(molecule.identity.clone())
                {
                    total_reverse_reads[idx] += 1;
                }
            }
        }

        if positions.len() == 1 {
            if let Some(seen) = observations[0] {
                if seen.quality >= min_phred_quality
                    && seen.base.eq_ignore_ascii_case(&alt_chars[0])
                    && unique_snp[0].insert(molecule.identity.clone())
                {
                    snp_counts[0] += 1;
                    increment_directional_count(
                        seen.from_forward,
                        seen.from_reverse,
                        &mut snp_forward_counts[0],
                        &mut snp_reverse_counts[0],
                    );
                }
            }
            continue;
        }

        let covers_all_positions = observations
            .iter()
            .all(|value| matches!(value, Some(seen) if seen.quality >= min_phred_quality));
        if covers_all_positions {
            unique_mnv_total.insert(molecule.identity.clone());
            // A haplotype assembled from a forward mate and a reverse mate was
            // established by neither strand alone, so it belongs to neither arm.
            if molecule.spans_all_on(true, min_phred_quality, &requested_indices) {
                unique_mnv_total_forward.insert(molecule.identity.clone());
            }
            if molecule.spans_all_on(false, min_phred_quality, &requested_indices) {
                unique_mnv_total_reverse.insert(molecule.identity.clone());
            }
        }

        let support = per_read_multi_support
            .entry(molecule.identity.clone())
            .or_insert_with(|| MultiReadSupport::new(positions.len()));
        support.spans_all |= covers_all_positions;
        for (idx, observation) in observations.iter().enumerate() {
            if let Some(seen) = observation {
                if seen.quality >= min_phred_quality
                    && seen.base.eq_ignore_ascii_case(&alt_chars[idx])
                {
                    support.snp_support[idx] = true;
                    support.snp_forward[idx] |= seen.from_forward;
                    support.snp_reverse[idx] |= seen.from_reverse;
                }
            }
        }
    }

    // The pattern each spanning molecule showed, which is the joint distribution
    // the linkage statistics are computed from. A molecule that did not observe
    // every position saw no combination and is left out.
    let mut haplotype_patterns: BTreeMap<Vec<bool>, usize> = BTreeMap::new();
    if positions.len() > 1 {
        for support in per_read_multi_support.values() {
            if support.spans_all {
                *haplotype_patterns
                    .entry(support.snp_support.clone())
                    .or_default() += 1;
            }
            if support.snp_support.iter().all(|is_supported| *is_supported) {
                mnv_count += 1;
                // A strand supports the haplotype only if it read every one of
                // its substitutions. `any` would credit both arms for a
                // haplotype pieced together from two mates, which is the same
                // over-claim the depth arms above avoid.
                let has_forward = support.snp_forward.iter().all(|is_forward| *is_forward);
                let has_reverse = support.snp_reverse.iter().all(|is_reverse| *is_reverse);
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
        haplotype_patterns: haplotype_patterns.into_iter().collect(),
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

#[cfg(test)]
mod counting_tests {
    use super::*;

    const POSITIONS: [usize; 2] = [100, 102];
    const GOOD: u8 = 40;
    const POOR: u8 = 5;

    fn alts() -> Vec<String> {
        vec!["T".to_string(), "A".to_string()]
    }

    fn count(records: &[DescribedRecord<'_>], pair_aware: bool) -> ReadCountSummary {
        let cache = cache_from_molecules(&POSITIONS, records);
        count_reads_from_cache(&cache, &POSITIONS, &alts(), 20, pair_aware).expect("counted")
    }

    /// One read carrying both alternates is one haplotype-supporting molecule.
    #[test]
    fn a_read_carrying_every_alternate_supports_the_haplotype() {
        let summary = count(
            &[("r1", false, vec![Some(('T', GOOD)), Some(('A', GOOD))])],
            true,
        );
        assert_eq!(summary.mnv_count, 1);
        assert_eq!(summary.mnv_total_reads, 1);
        assert_eq!(summary.snp_only_informative_counts, vec![0, 0]);
        assert_eq!(summary.haplotype_patterns, vec![(vec![true, true], 1)]);
    }

    /// Two mates of one fragment are one molecule, not two observations.
    #[test]
    fn overlapping_mates_are_one_molecule() {
        let both = vec![Some(('T', GOOD)), Some(('A', GOOD))];
        let summary = count(&[("frag", false, both.clone()), ("frag", true, both)], true);
        assert_eq!(summary.mnv_count, 1);
        assert_eq!(summary.mnv_total_reads, 1);
        assert_eq!(summary.total_reads, vec![1, 1]);
    }

    /// The same records counted per record are two observations.
    #[test]
    fn counting_mates_separately_gives_two_observations() {
        let both = vec![Some(('T', GOOD)), Some(('A', GOOD))];
        let summary = count(
            &[("frag", false, both.clone()), ("frag", true, both)],
            false,
        );
        assert_eq!(summary.mnv_count, 2);
        assert_eq!(summary.total_reads, vec![2, 2]);
    }

    /// A molecule read from both ends is evidence on both strands.
    #[test]
    fn an_overlapping_pair_counts_on_both_strands() {
        let both = vec![Some(('T', GOOD)), Some(('A', GOOD))];
        let summary = count(&[("frag", false, both.clone()), ("frag", true, both)], true);
        assert_eq!(summary.mnv_forward_count, 1);
        assert_eq!(summary.mnv_reverse_count, 1);
        assert_eq!(summary.total_forward_reads, vec![1, 1]);
        assert_eq!(summary.total_reverse_reads, vec![1, 1]);
        // Each strand read the whole codon by itself here, so each arm holds
        // the molecule.
        assert_eq!(summary.mnv_total_forward_reads, 1);
        assert_eq!(summary.mnv_total_reverse_reads, 1);
    }

    /// A base exactly at the floor passes: the gate is "at least", and an
    /// off-by-one there silently drops a whole quality tier.
    #[test]
    fn a_base_exactly_at_the_quality_floor_is_observed() {
        let summary = count(
            &[("r1", false, vec![Some(('T', 20)), Some(('A', 20))])],
            true,
        );
        assert_eq!(summary.total_reads, vec![1, 1]);
        assert_eq!(summary.mnv_count, 1);
        assert_eq!(summary.mnv_total_reads, 1);
        assert_eq!(summary.mnv_total_forward_reads, 1);
    }

    /// The strand belongs to the position that was read, not to the fragment.
    /// A mate covering one position lends nothing to the position the other
    /// mate covered, which is how a variant read only on reverse reads came to
    /// report forward support.
    #[test]
    fn strand_does_not_travel_between_positions_of_a_fragment() {
        let summary = count(
            &[
                ("frag", false, vec![Some(('T', GOOD)), None]),
                ("frag", true, vec![None, Some(('A', GOOD))]),
            ],
            true,
        );
        assert_eq!(summary.total_forward_reads, vec![1, 0]);
        assert_eq!(summary.total_reverse_reads, vec![0, 1]);
        // The molecule does span both positions, so it is a haplotype; but
        // neither strand established it alone.
        assert_eq!(summary.mnv_count, 1);
        assert_eq!(summary.mnv_forward_count, 0);
        assert_eq!(summary.mnv_reverse_count, 0);
        assert_eq!(summary.mnv_total_forward_reads, 0);
        assert_eq!(summary.mnv_total_reverse_reads, 0);
    }

    /// Mates that read different bases leave the molecule with no observation
    /// there: one of them is wrong and there is no telling which.
    #[test]
    fn mates_that_disagree_observe_nothing_at_that_position() {
        let summary = count(
            &[
                ("frag", false, vec![Some(('T', GOOD)), Some(('A', GOOD))]),
                ("frag", true, vec![Some(('G', GOOD)), Some(('A', GOOD))]),
            ],
            true,
        );
        assert_eq!(summary.total_reads, vec![0, 1]);
        assert_eq!(summary.mnv_count, 0);
        assert_eq!(summary.mnv_total_reads, 0);
    }

    /// A read that stops between the positions is not evidence about the pair,
    /// on either side of the ratio.
    #[test]
    fn a_read_covering_one_position_is_not_informative_about_linkage() {
        let summary = count(
            &[
                (
                    "spanning",
                    false,
                    vec![Some(('T', GOOD)), Some(('A', GOOD))],
                ),
                ("partial", true, vec![Some(('T', GOOD)), None]),
            ],
            true,
        );
        // The partial molecule is reverse and gives the first position reverse
        // depth, but it spans nothing, so the haplotype's reverse arm stays
        // empty. The two numbers must not be read from the same place.
        assert_eq!(summary.total_reverse_reads, vec![1, 0]);
        assert_eq!(summary.mnv_total_reverse_reads, 0);
        assert_eq!(summary.mnv_count, 1);
        assert_eq!(summary.mnv_total_reads, 1);
        // The partial read carries the first alternate but never saw the
        // second, so it argues neither for nor against the haplotype.
        assert_eq!(summary.snp_only_informative_counts, vec![0, 0]);
        assert_eq!(summary.snp_counts, vec![1, 0]);
    }

    /// A spanning read carrying one alternate and not the other argues against
    /// linkage, and is the only kind entitled to.
    #[test]
    fn a_spanning_read_with_one_alternate_argues_against_linkage() {
        let summary = count(
            &[
                (
                    "spanning",
                    false,
                    vec![Some(('T', GOOD)), Some(('A', GOOD))],
                ),
                ("solo", false, vec![Some(('T', GOOD)), Some(('C', GOOD))]),
            ],
            true,
        );
        assert_eq!(summary.mnv_count, 1);
        assert_eq!(summary.snp_only_informative_counts, vec![1, 0]);
        assert_eq!(
            summary.haplotype_patterns,
            vec![(vec![true, false], 1), (vec![true, true], 1)]
        );
    }

    /// Bases below the floor are not observed at all, so they leave both the
    /// numerator and the denominator.
    #[test]
    fn a_base_below_the_quality_floor_is_not_observed() {
        let summary = count(
            &[("r1", false, vec![Some(('T', POOR)), Some(('A', GOOD))])],
            true,
        );
        assert_eq!(summary.total_reads, vec![0, 1]);
        // The second position was read cleanly and does carry its alternate,
        // so it counts there; the first was never observed, which is why the
        // molecule spans nothing and joins no combination.
        assert_eq!(summary.snp_counts, vec![0, 1]);
        assert_eq!(summary.snp_only_informative_counts, vec![0, 0]);
        assert_eq!(summary.mnv_total_reads, 0);
        assert!(summary.haplotype_patterns.is_empty());
    }

    /// Where the mates agree, the better base quality carries the position over
    /// the floor.
    #[test]
    fn agreeing_mates_keep_the_better_quality() {
        let summary = count(
            &[
                ("frag", false, vec![Some(('T', POOR)), Some(('A', GOOD))]),
                ("frag", true, vec![Some(('T', GOOD)), Some(('A', GOOD))]),
            ],
            true,
        );
        assert_eq!(summary.total_reads, vec![1, 1]);
        assert_eq!(summary.mnv_count, 1);
    }

    /// Reference-carrying molecules belong in the joint distribution: without
    /// them the linkage statistics cannot tell a haplotype from two common
    /// alleles meeting by chance.
    #[test]
    fn molecules_carrying_neither_alternate_are_still_counted() {
        let summary = count(
            &[
                ("both", false, vec![Some(('T', GOOD)), Some(('A', GOOD))]),
                ("neither", false, vec![Some(('G', GOOD)), Some(('C', GOOD))]),
            ],
            true,
        );
        assert_eq!(summary.mnv_total_reads, 2);
        assert_eq!(
            summary.haplotype_patterns,
            vec![(vec![false, false], 1), (vec![true, true], 1)]
        );
    }

    /// A single-position request counts that position alone and asks no
    /// linkage question.
    #[test]
    fn a_single_position_request_counts_only_that_position() {
        let cache = cache_from_molecules(
            &[100],
            &[
                ("alt", false, vec![Some(('T', GOOD))]),
                ("ref", true, vec![Some(('G', GOOD))]),
            ],
        );
        let summary =
            count_reads_from_cache(&cache, &[100], &["T".to_string()], 20, true).expect("counted");
        assert_eq!(summary.snp_counts, vec![1]);
        assert_eq!(summary.snp_forward_counts, vec![1]);
        assert_eq!(summary.snp_reverse_counts, vec![0]);
        assert_eq!(summary.total_reads, vec![2]);
        assert!(summary.haplotype_patterns.is_empty());
        // With one position the haplotype depths are that position's depths;
        // there is no set of spanning molecules distinct from it.
        assert_eq!(summary.mnv_total_reads, 2);
        assert_eq!(summary.mnv_total_forward_reads, 1);
        assert_eq!(summary.mnv_total_reverse_reads, 1);
    }
}
