//! A codon read straight off the molecules, whatever the variant calls said.
//!
//! Every other read count in get_MNV starts from a VCF record: the reads are
//! asked about variants a caller already reported. A minority population below
//! the caller's threshold is then never asked about, however many reads carry
//! it. This reader starts from a codon instead: each molecule that observes all
//! three bases contributes the triplet it carries, and the triplets are counted.
//! A triplet is one molecule's codon, so two changes seen together here are on
//! one molecule, and two changes on different molecules are two triplets.

use crate::error::{AppError, AppResult};
use noodles::bam;
use noodles::sam::Header;
use std::collections::{BTreeMap, HashMap};

use super::observation::{build_molecule_key, build_region, observed_allele_for_ref_span};

pub struct CodonReadRequest<'a> {
    pub chrom: &'a str,
    /// The codon's three genomic positions, in the order the triplet is to be
    /// written (the caller's order; for a minus-strand gene the caller
    /// complements).
    pub positions: [usize; 3],
    pub min_phred_quality: u8,
    pub min_mapq: u8,
    /// Whether the two mates of a fragment are one molecule.
    pub pair_aware: bool,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CodonTriplets {
    /// Molecules per triplet, in genomic bases in the requested order.
    pub counts: BTreeMap<[u8; 3], usize>,
    /// Molecules that observed all three bases and were assigned a triplet.
    pub spanning: usize,
    /// Molecules that reached the codon but left a base unobserved, below the
    /// quality floor, inside an indel, or with mates that disagree. A partial
    /// codon is not a triplet, and counting it as one would assert bases nobody
    /// read.
    pub partial: usize,
}

/// What one record says about one base: `Some(base)` when it read a single
/// A/C/G/T at or above the quality floor, `None` otherwise.
fn base_at(record: &bam::Record, position: usize, min_phred_quality: u8) -> Option<u8> {
    let observed = observed_allele_for_ref_span(record, position, 1)?;
    if observed.min_quality < min_phred_quality || observed.allele.len() != 1 {
        return None;
    }
    let base = observed.allele.as_bytes()[0].to_ascii_uppercase();
    matches!(base, b'A' | b'C' | b'G' | b'T').then_some(base)
}

/// Whether a record reaches any part of `start..=end`.
fn reaches(record: &bam::Record, start: usize, end: usize) -> bool {
    let Some(Ok(first)) = record.alignment_start() else {
        return false;
    };
    let first: usize = first.into();
    let span: usize = record
        .cigar()
        .iter()
        .filter_map(Result::ok)
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum();
    first <= end && first + span.max(1) > start
}

pub fn read_codon_triplets(
    bam_reader: &mut crate::io::AlignmentReader,
    header: &Header,
    request: CodonReadRequest<'_>,
) -> AppResult<CodonTriplets> {
    let CodonReadRequest {
        chrom,
        positions,
        min_phred_quality,
        min_mapq,
        pair_aware,
    } = request;
    let start = *positions.iter().min().expect("three positions");
    let end = *positions.iter().max().expect("three positions");
    if start == 0 {
        return Err(AppError::validation(format!(
            "Invalid codon position on {chrom}: positions are 1-based"
        )));
    }
    let region = build_region(chrom, start, end)?;
    let mut query = bam_reader.query(header, &region).map_err(|e| {
        AppError::validation(format!(
            "Alignment query failed for {chrom}:{start}-{end}: {e}"
        ))
    })?;

    struct MoleculeView {
        bases: [Option<u8>; 3],
        conflicted: bool,
    }
    let mut views: HashMap<super::observation::MoleculeKey, MoleculeView> = HashMap::new();
    let mut record_index = 0usize;
    let mut record = bam::Record::default();
    while query
        .read_record(&mut record)
        .map_err(|e| AppError::validation(format!("Alignment read error: {e}")))?
        != 0
    {
        record_index += 1;
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_duplicate()
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
        if mapq < min_mapq || !reaches(&record, start, end) {
            continue;
        }
        let read = positions.map(|position| base_at(&record, position, min_phred_quality));
        let view = views
            .entry(build_molecule_key(&record, pair_aware, record_index))
            .or_insert(MoleculeView {
                bases: [None; 3],
                conflicted: false,
            });
        for (slot, base) in view.bases.iter_mut().zip(read) {
            if let Some(base) = base {
                if slot.is_some_and(|seen| seen != base) {
                    view.conflicted = true;
                }
                *slot = Some(base);
            }
        }
    }

    let mut triplets = CodonTriplets::default();
    for view in views.into_values() {
        match view.bases {
            [Some(a), Some(b), Some(c)] if !view.conflicted => {
                *triplets.counts.entry([a, b, c]).or_insert(0) += 1;
                triplets.spanning += 1;
            }
            _ => triplets.partial += 1,
        }
    }
    Ok(triplets)
}
