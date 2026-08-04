//! BAM read counting: region observation cache, per-position and MNV haplotype
//! read support with strand-specific metrics.
//!
//! The per-read observation primitives, the direct indel counter, and the
//! region cache live in submodules; this file defines the public summary and
//! request types plus the `count_reads_per_position` convenience entry point.

use crate::error::{AppError, AppResult};
use crate::variants::AlleleComponent;
use noodles::bam;
use noodles::sam::Header;

mod cache;
mod indel;
mod observation;
mod phasing;

pub use cache::{build_region_observation_cache, count_reads_from_cache, RegionObservationCache};
pub use indel::count_indel_reads;
pub use phasing::indel_snv_in_trans;

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
    /// Per position: reads that observe *every* requested position and carry
    /// this position's ALT without carrying the full haplotype. These are the
    /// reads that argue against linkage, and they are the only ones entitled
    /// to: a read that stops before a partner position saw no evidence either
    /// way. Always zero for a single-position request.
    pub snp_only_informative_counts: Vec<usize>,
}

pub struct IndelReadCountRequest<'a> {
    pub chrom: &'a str,
    pub position: usize,
    pub ref_allele: &'a str,
    pub alt_allele: &'a str,
    pub required_components: &'a [AlleleComponent],
    pub min_phred_quality: u8,
    pub min_mapq: u8,
    /// When true, locus depth (the EDP/EFREQ denominator) counts every read that
    /// observes the anchor base at sufficient quality, instead of only reads
    /// that fully span the REF allele. Reduces depth under-counting / EFREQ bias
    /// for multi-base deletions. Defaults to false (historical behaviour).
    pub anchor_depth: bool,
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
