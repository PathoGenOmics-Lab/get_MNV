//! Local haplotypes read straight off the molecules.
//!
//! Given a window and the variants called inside it, this reports which
//! combinations the reads actually carry, and how many reads carry each. It
//! replaces asking the BAM "does this combination exist?" once per candidate
//! combination, which needs `2^n` questions and, worse, answers yes to
//! combinations no molecule is: a read carrying A, B and C matches the compound
//! allele of A and B as well, because that allele only spans as far as B.
//!
//! Reading the molecules instead has no combinatorial limit, so a window with
//! more variants than the old cap no longer has to fall back to reporting only
//! pairs and the full set.

use crate::error::{AppError, AppResult};
use crate::variants::AlleleComponent;
use noodles::bam;
use noodles::sam::Header;
use std::collections::{BTreeMap, HashSet};
use std::rc::Rc;

use super::observation::{
    build_read_key, build_region, observed_allele_for_ref_span, observed_supports_components,
    ReadKey,
};

pub struct LocalHaplotypeRequest<'a> {
    pub chrom: &'a str,
    /// First and last reference position a read must cover to be assigned a
    /// haplotype. A read that stops inside the window witnessed only part of
    /// the molecule and is counted separately rather than guessed at.
    pub start: usize,
    pub end: usize,
    /// The components of each candidate variant, in the caller's order. The
    /// `carried` vector of every returned haplotype uses that same order.
    pub variants: &'a [Vec<AlleleComponent>],
    pub min_phred_quality: u8,
    pub min_mapq: u8,
}

#[derive(Debug, Clone)]
pub struct ObservedHaplotype {
    /// `carried[i]` is true when the reads of this haplotype carry variant `i`.
    pub carried: Vec<bool>,
    pub reads: usize,
}

#[derive(Debug, Clone, Default)]
pub struct LocalHaplotypeObservations {
    /// One entry per distinct combination seen, including the all-reference
    /// one, ordered deterministically.
    pub haplotypes: Vec<ObservedHaplotype>,
    /// Reads that covered the whole window and were assigned a combination.
    pub spanning_reads: usize,
    /// Reads that overlapped the window but stopped inside it, or fell below
    /// the base-quality floor. They witnessed no complete molecule here.
    pub partial_reads: usize,
}

pub fn observe_local_haplotypes(
    bam_reader: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: &Header,
    request: LocalHaplotypeRequest<'_>,
) -> AppResult<LocalHaplotypeObservations> {
    let LocalHaplotypeRequest {
        chrom,
        start,
        end,
        variants,
        min_phred_quality,
        min_mapq,
    } = request;

    if start == 0 || end < start {
        return Err(AppError::validation(format!(
            "Invalid local haplotype window {chrom}:{start}-{end}"
        )));
    }
    if variants.is_empty() {
        return Ok(LocalHaplotypeObservations::default());
    }

    let region = build_region(chrom, start, end)?;
    let mut query = bam_reader.query(header, &region).map_err(|e| {
        AppError::validation(format!("BAM query failed for {chrom}:{start}-{end}: {e}"))
    })?;

    let span = end - start + 1;
    let mut by_combination: BTreeMap<Vec<bool>, HashSet<Rc<ReadKey>>> = BTreeMap::new();
    let mut partial: HashSet<Rc<ReadKey>> = HashSet::new();
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

        let key = Rc::new(build_read_key(&record));
        let Some(observed) = observed_allele_for_ref_span(&record, start, span) else {
            partial.insert(key);
            continue;
        };
        if observed.min_quality < min_phred_quality {
            partial.insert(key);
            continue;
        }

        let carried = variants
            .iter()
            .map(|components| observed_supports_components(&observed, components))
            .collect::<Vec<_>>();
        by_combination.entry(carried).or_default().insert(key);
    }

    let spanning_reads = by_combination.values().map(HashSet::len).sum();
    let haplotypes = by_combination
        .into_iter()
        .map(|(carried, reads)| ObservedHaplotype {
            carried,
            reads: reads.len(),
        })
        .collect();

    Ok(LocalHaplotypeObservations {
        haplotypes,
        spanning_reads,
        partial_reads: partial.len(),
    })
}
