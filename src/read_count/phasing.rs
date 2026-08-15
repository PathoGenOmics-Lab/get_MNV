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
    // Same observation span as the exact counting: a read that stops on an
    // insertion's anchor has seen nothing of the junction and must not be
    // recorded as spanning the indel and not carrying it.
    let indel_observation_span = crate::variants::observation_ref_len(
        indel_ref,
        &crate::variants::decompose_allele(indel_position, indel_ref, indel_alt).components,
    );
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
        saw_snv_position: bool,
        spans_indel: bool,
        carries_indel: bool,
        /// The molecule's mates read different things at the SNV or over the
        /// indel. Neither claim survives, so the molecule is not informative.
        conflicted: bool,
    }
    let mut evidence: HashMap<MoleculeKey, MoleculeEvidence> = HashMap::new();
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

        let observed_snv_base = observed_base_at(&record, snv_position, min_phred_quality);
        let carries_snv =
            matches!(observed_snv_base, Some(base) if base.eq_ignore_ascii_case(&snv_alt));
        let saw_the_whole_indel = indel_observation_span == indel_ref.len()
            || observed_allele_for_ref_span(&record, indel_position, indel_observation_span)
                .is_some();
        let observed_indel_allele =
            observed_allele_for_ref_span(&record, indel_position, indel_ref.len()).filter(
                |observed| observed.min_quality >= min_phred_quality && saw_the_whole_indel,
            );
        if observed_snv_base.is_none() && observed_indel_allele.is_none() {
            continue;
        }

        let molecule = evidence
            .entry(build_molecule_key(&record, pair_aware, record_index))
            .or_default();
        if observed_snv_base.is_some() {
            if molecule.saw_snv_position && molecule.carries_snv != carries_snv {
                molecule.conflicted = true;
            }
            molecule.saw_snv_position = true;
            molecule.carries_snv |= carries_snv;
        }
        if let Some(observed) = observed_indel_allele {
            let carries_indel = observed.allele.eq_ignore_ascii_case(indel_alt);
            if molecule.spans_indel && molecule.carries_indel != carries_indel {
                molecule.conflicted = true;
            }
            molecule.spans_indel = true;
            molecule.carries_indel |= carries_indel;
        }
    }

    let informative = evidence
        .values()
        .filter(|molecule| !molecule.conflicted && molecule.carries_snv && molecule.spans_indel);
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

#[cfg(test)]
mod tests {
    use super::*;

    fn linkage(cis_reads: usize, informative_reads: usize) -> IndelSnvLinkage {
        IndelSnvLinkage {
            informative_reads,
            cis_reads,
        }
    }

    /// The rule that decides whether an upstream indel shifts a downstream
    /// codon's reading frame. It had no test of its own, so nothing pinned
    /// either threshold.
    #[test]
    fn no_molecule_carrying_the_indel_is_trans() {
        assert!(linkage(0, 18).is_trans());
        assert!(linkage(0, 18).is_informative());
    }

    #[test]
    fn a_handful_of_cis_molecules_is_tolerated_as_sequencing_error() {
        // Up to a tenth of the informative molecules may carry the indel and
        // the pair is still called trans.
        assert!(linkage(2, 20).is_trans());
        assert!(!linkage(3, 20).is_trans());
    }

    #[test]
    fn molecules_that_mostly_carry_both_are_not_trans() {
        assert!(!linkage(17, 18).is_trans());
        assert!(!linkage(9, 18).is_trans());
    }

    #[test]
    fn one_molecule_settles_nothing() {
        // A single read is not evidence either way, and "not trans" here means
        // unknown rather than cis.
        assert!(!linkage(0, 1).is_trans());
        assert!(!linkage(0, 1).is_informative());
        // Two is the least that can say anything.
        assert!(linkage(0, 2).is_trans());
        assert!(linkage(0, 2).is_informative());
    }

    #[test]
    fn nothing_spanning_both_loci_says_nothing() {
        assert!(!linkage(0, 0).is_trans());
        assert!(!linkage(0, 0).is_informative());
    }
}
