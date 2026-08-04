//! Read-level cis/trans phasing between an indel and a downstream SNV.
//!
//! Used to suppress frameshift propagation when the BAM shows the two variants
//! sit on different molecules: an upstream frameshift indel only shifts the
//! reading frame of a downstream codon on the molecules that actually carry it.

use crate::error::{AppError, AppResult};
use noodles::bam;
use noodles::sam::Header;
use std::collections::HashMap;

use super::observation::{
    build_molecule_key, build_region, get_query_pos, observed_allele_for_ref_span, MoleculeKey,
};

/// Minimum number of SNV-carrying reads that also span the indel before a pair
/// can be judged; below this there is too little evidence and the pair is left
/// "unknown" (the caller then keeps its frequency-based behaviour).
const MIN_INFORMATIVE_READS: usize = 2;
/// A pair is called trans only when at most this fraction of the informative
/// reads also carry the indel (tolerating a few sequencing errors).
const MAX_CIS_FRACTION: f64 = 0.1;

/// The base observed at genomic `position` on `rec`, if present at sufficient
/// quality. Mirrors the base/quality extraction used by the region cache.
fn observed_base_at(rec: &bam::Record, position: usize, min_phred_quality: u8) -> Option<char> {
    let idx = get_query_pos(rec, position)?;
    let seq = rec.sequence();
    if idx >= seq.len() {
        return None;
    }
    let base = seq.iter().nth(idx)? as char;
    let quals = rec.quality_scores();
    let quality = if quals.iter().next().is_some() {
        quals.iter().nth(idx).unwrap_or(0)
    } else {
        u8::MAX
    };
    (quality >= min_phred_quality).then_some(base)
}

/// What the BAM says about an indel and a downstream SNV sharing molecules.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct IndelSnvLinkage {
    /// Reads carrying the SNV alt that also span the indel locus. Only these
    /// can answer the question.
    pub informative_reads: usize,
    /// Of those, how many also carry the indel.
    pub cis_reads: usize,
}

impl IndelSnvLinkage {
    /// Whether the reads settle the pair as being on different molecules.
    /// Below [`MIN_INFORMATIVE_READS`] there is too little evidence and the
    /// answer is no, which the caller reads as unknown rather than as cis.
    pub fn is_trans(&self) -> bool {
        self.informative_reads >= MIN_INFORMATIVE_READS
            && (self.cis_reads as f64) <= MAX_CIS_FRACTION * (self.informative_reads as f64)
    }

    /// Whether the reads say anything at all.
    pub fn is_informative(&self) -> bool {
        self.informative_reads >= MIN_INFORMATIVE_READS
    }
}

/// How the BAM sees the indel and the downstream SNV: among the reads that
/// carry the SNV alt *and* span the indel locus, how many also carry the indel.
/// Purely an observation; the suppression decision is
/// [`IndelSnvLinkage::is_trans`], which never asserts cis.
#[allow(clippy::too_many_arguments)]
pub fn indel_snv_linkage(
    bam_reader: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: &Header,
    chrom: &str,
    indel_position: usize,
    indel_ref: &str,
    indel_alt: &str,
    snv_position: usize,
    snv_alt: char,
    min_phred_quality: u8,
    min_mapq: u8,
    pair_aware: bool,
) -> AppResult<IndelSnvLinkage> {
    let indel_end = indel_position + indel_ref.len().saturating_sub(1);
    let start = indel_position.min(snv_position);
    let end = indel_end.max(snv_position);
    let region = build_region(chrom, start, end)?;
    let mut query = bam_reader.query(header, &region).map_err(|e| {
        AppError::validation(format!("BAM query failed for {chrom}:{start}-{end}: {e}"))
    })?;

    // Gathered per molecule, not per record: with paired-end reads the mate
    // carrying the SNV is often not the mate that reaches the indel, and those
    // two are the same molecule. Counting records would throw that pairing away
    // and call the question unanswerable.
    #[derive(Default)]
    struct MoleculeEvidence {
        carries_snv: bool,
        spans_indel: bool,
        carries_indel: bool,
    }
    let mut evidence: HashMap<MoleculeKey, MoleculeEvidence> = HashMap::new();
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

        let carries_snv = matches!(
            observed_base_at(&record, snv_position, min_phred_quality),
            Some(base) if base.eq_ignore_ascii_case(&snv_alt)
        );
        let observed_indel_allele =
            observed_allele_for_ref_span(&record, indel_position, indel_ref.len())
                .filter(|observed| observed.min_quality >= min_phred_quality);
        if !carries_snv && observed_indel_allele.is_none() {
            continue;
        }

        let molecule = evidence
            .entry(build_molecule_key(&record, pair_aware))
            .or_default();
        molecule.carries_snv |= carries_snv;
        if let Some(observed) = observed_indel_allele {
            molecule.spans_indel = true;
            molecule.carries_indel |= observed.allele.eq_ignore_ascii_case(indel_alt);
        }
    }

    let informative = evidence
        .values()
        .filter(|molecule| molecule.carries_snv && molecule.spans_indel);
    let mut informative_reads = 0usize;
    let mut cis_reads = 0usize;
    for molecule in informative {
        informative_reads += 1;
        if molecule.carries_indel {
            cis_reads += 1;
        }
    }

    Ok(IndelSnvLinkage {
        informative_reads,
        cis_reads,
    })
}
