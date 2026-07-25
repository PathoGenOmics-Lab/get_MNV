//! Nonsense-mediated decay (NMD) prediction under the 50-nucleotide rule.
//!
//! A premature termination codon (PTC) triggers NMD when it lies more than 50 nt
//! upstream (5') of the last exon-exon junction of a multi-exon transcript. A PTC
//! in the last exon, or within 50 nt of the last junction, escapes NMD, as does
//! any stop in a single-exon (intronless) transcript, which has no junction.
//!
//! All positions are measured in spliced-CDS coordinates (0 = first base of the
//! start codon), the same coordinate system as [`transcript_offset_for_position`].

use crate::io::VcfPosition;
use crate::variants::{Gene, NmdPrediction};

use super::allele_apply::apply_allele_to_feature;
use super::protein::translate_cds;
use super::transcript_model::{coding_sequence_for_gene, first_touched_transcript_offset};

/// Distance (in coding nucleotides) below which a PTC upstream of the last
/// exon-exon junction still escapes NMD. The canonical rule uses 50-55 nt.
const NMD_LAST_JUNCTION_MARGIN: usize = 50;

/// Total spliced-CDS length and the length of the terminal (3') coding exon,
/// or `None` for a single-exon transcript (fewer than two CDS segments), which
/// has no exon-exon junction and therefore no NMD prediction.
fn cds_total_and_last_exon(gene: &Gene) -> Option<(usize, usize)> {
    if gene.cds_segments.len() < 2 {
        return None;
    }
    let total = gene
        .cds_segments
        .iter()
        .map(|segment| segment.end.saturating_sub(segment.start) + 1)
        .sum();
    // `cds_segments` is in transcript order, so the terminal exon is the last
    // entry on either strand.
    let last = gene.cds_segments.last()?;
    let last_len = last.end.saturating_sub(last.start) + 1;
    Some((total, last_len))
}

/// CDS offset of the last exon-exon junction (the start of the terminal coding
/// exon), or `None` for a single-exon transcript.
fn last_exon_junction_offset(gene: &Gene) -> Option<usize> {
    let (total, last_len) = cds_total_and_last_exon(gene)?;
    Some(total.saturating_sub(last_len))
}

/// Classify a PTC by its distance upstream of the last exon-exon junction.
fn classify(ptc_cds_offset: usize, last_junction_offset: usize) -> NmdPrediction {
    if last_junction_offset > ptc_cds_offset
        && last_junction_offset - ptc_cds_offset > NMD_LAST_JUNCTION_MARGIN
    {
        NmdPrediction::Triggering
    } else {
        NmdPrediction::Escaping
    }
}

/// NMD prediction for a nonsense SNV/MNV whose premature stop codon begins at
/// `ptc_cds_offset` in the reference spliced CDS (no length change, so reference
/// and alternate coordinates coincide). `None` for single-exon transcripts.
pub(super) fn snv_mnv_nmd_prediction(gene: &Gene, ptc_cds_offset: usize) -> Option<NmdPrediction> {
    let last_junction = last_exon_junction_offset(gene)?;
    Some(classify(ptc_cds_offset, last_junction))
}

/// NMD prediction for an indel that introduces a premature stop (a frameshift
/// PTC, or an in-frame indel that gains a stop). `None` when the transcript is
/// single-exon, the CDS cannot be reconstructed, or the alternate protein does
/// not terminate earlier than the reference (no premature stop).
///
/// Coordinates are taken in the *alternate* CDS, whose length changes by the
/// indel's coding delta. The terminal exon keeps its length unless the indel
/// falls inside it, so the last junction sits at `alt_len - terminal_exon_len`
/// when the indel is upstream of it, and at the reference junction otherwise.
pub(super) fn indel_nmd_prediction(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
    genetic_code: crate::genetic_code::GeneticCode,
) -> Option<NmdPrediction> {
    let (_, last_exon_len) = cds_total_and_last_exon(gene)?;
    let ref_last_junction = last_exon_junction_offset(gene)?;

    let ref_cds = coding_sequence_for_gene(gene, reference)?;
    let alt_cds = apply_allele_to_feature(gene, reference, variant)?;
    let ref_protein = translate_cds(&ref_cds, genetic_code);
    let alt_protein = translate_cds(&alt_cds, genetic_code);

    // The premature stop is the alternate protein's first stop, and only when it
    // truncates earlier than the reference's natural stop.
    let ref_stop = ref_protein.find('*').unwrap_or(ref_protein.len());
    let alt_stop = alt_protein.find('*')?;
    if alt_stop >= ref_stop {
        return None;
    }
    let ptc_cds_offset = alt_stop * 3;

    // Place the last junction in alternate-CDS coordinates: the terminal exon's
    // length is unchanged unless the indel lands inside it, so an indel upstream
    // of the junction shifts it by the coding-length delta (captured here as
    // `alt_cds.len() - last_exon_len`); an indel in the terminal exon leaves the
    // junction at its reference offset.
    let indel_upstream_of_last_junction =
        first_touched_transcript_offset(gene, variant).is_some_and(|pos| pos < ref_last_junction);
    let last_junction = if indel_upstream_of_last_junction {
        alt_cds.len().saturating_sub(last_exon_len)
    } else {
        ref_last_junction
    };

    Some(classify(ptc_cds_offset, last_junction))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::{CdsSegment, Strand};

    /// A two-exon plus-strand transcript: exon 1 spans genomic 1..=exon1_len,
    /// exon 2 spans the next `exon2_len` bases. Coding length is exon1+exon2.
    fn two_exon_gene(exon1_len: usize, exon2_len: usize) -> Gene {
        Gene {
            name: "t".to_string(),
            start: 1,
            end: exon1_len + exon2_len,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: Some("tx".to_string()),
            cds_segments: vec![
                CdsSegment {
                    start: 1,
                    end: exon1_len,
                },
                CdsSegment {
                    start: exon1_len + 1,
                    end: exon1_len + exon2_len,
                },
            ],
        }
    }

    #[test]
    fn ptc_far_upstream_of_last_junction_triggers_nmd() {
        // Last junction at CDS offset 120 (exon1 = 120 nt). A stop at offset 30
        // is 90 nt upstream (> 50) -> NMD-triggering.
        let gene = two_exon_gene(120, 60);
        assert_eq!(
            snv_mnv_nmd_prediction(&gene, 30),
            Some(NmdPrediction::Triggering)
        );
    }

    #[test]
    fn ptc_within_50nt_of_last_junction_escapes() {
        // Last junction at offset 120; a stop at offset 90 is exactly 30 nt
        // upstream (<= 50) -> escapes.
        let gene = two_exon_gene(120, 60);
        assert_eq!(
            snv_mnv_nmd_prediction(&gene, 90),
            Some(NmdPrediction::Escaping)
        );
        // Exactly 50 nt upstream still escapes (the rule needs strictly > 50).
        assert_eq!(
            snv_mnv_nmd_prediction(&gene, 70),
            Some(NmdPrediction::Escaping)
        );
        // 51 nt upstream triggers.
        assert_eq!(
            snv_mnv_nmd_prediction(&gene, 69),
            Some(NmdPrediction::Triggering)
        );
    }

    #[test]
    fn ptc_in_last_exon_escapes() {
        // A stop in the terminal exon (offset >= 120) escapes.
        let gene = two_exon_gene(120, 60);
        assert_eq!(
            snv_mnv_nmd_prediction(&gene, 150),
            Some(NmdPrediction::Escaping)
        );
    }

    #[test]
    fn single_exon_transcript_has_no_prediction() {
        let mut gene = two_exon_gene(120, 60);
        gene.cds_segments.truncate(1);
        assert_eq!(snv_mnv_nmd_prediction(&gene, 30), None);
        gene.cds_segments.clear();
        assert_eq!(snv_mnv_nmd_prediction(&gene, 30), None);
    }
}
