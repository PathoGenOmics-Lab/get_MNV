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

/// Total length of the spliced CDS, in transcript coordinates.
fn cds_total(gene: &Gene) -> usize {
    gene.cds_segments
        .iter()
        .map(|segment| segment.end.saturating_sub(segment.start) + 1)
        .sum()
}

/// CDS offset of the last exon-exon junction, i.e. the transcript coordinate at
/// which the segment following the final intron begins. `None` when the
/// transcript is not spliced (a single CDS segment, or segments that abut or
/// overlap in a ribosomal-slippage join): with no junction there is no NMD
/// prediction to make.
///
/// Only segment pairs separated by a real intron count, so a slippage join that
/// follows a genuine intron does not move the junction downstream of it.
fn last_exon_junction_offset(gene: &Gene) -> Option<usize> {
    gene.introns().last().map(|intron| intron.junction_offset)
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
/// indel's coding delta. The terminal region keeps its length unless the indel
/// falls inside it, so the last junction sits at `alt_len - terminal_len` when
/// the indel is upstream of it, and at the reference junction otherwise.
pub(super) fn indel_nmd_prediction(
    gene: &Gene,
    reference: &crate::io::Reference<'_>,
    variant: &VcfPosition,
    genetic_code: crate::genetic_code::GeneticCode,
) -> Option<NmdPrediction> {
    let ref_last_junction = last_exon_junction_offset(gene)?;
    // Everything downstream of the last junction, which may span more than one
    // CDS segment when a slippage join follows the final intron.
    let terminal_len = cds_total(gene).saturating_sub(ref_last_junction);

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

    // Place the last junction in alternate-CDS coordinates: the terminal region's
    // length is unchanged unless the indel lands inside it, so an indel upstream
    // of the junction shifts it by the coding-length delta (captured here as
    // `alt_cds.len() - terminal_len`); an indel downstream of the junction leaves
    // the junction at its reference offset.
    let indel_upstream_of_last_junction =
        first_touched_transcript_offset(gene, variant).is_some_and(|pos| pos < ref_last_junction);
    let last_junction = if indel_upstream_of_last_junction {
        alt_cds.len().saturating_sub(terminal_len)
    } else {
        ref_last_junction
    };

    Some(classify(ptc_cds_offset, last_junction))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{joined_gene, spliced_gene, transcript_gene};
    use crate::variants::Strand;

    /// Length of the intron separating the two exons of the test transcript.
    /// The exons must be genuinely separated, otherwise the transcript is an
    /// unspliced join and has no exon-exon junction at all.
    const INTRON_LEN: usize = 500;

    /// A two-exon plus-strand transcript: exon 1 spans genomic 1..=exon1_len, an
    /// intron follows, then exon 2 spans `exon2_len` bases. Coding length is
    /// exon1+exon2, so CDS offsets are unaffected by the intron.
    ///
    /// Built through `spliced_gene`, which refuses segments that are not
    /// separated by an intron. The hand-written version of this helper produced
    /// abutting exons, so every NMD test below was silently asserting the
    /// behaviour of a transcript with no exon-exon junction at all.
    fn two_exon_gene(exon1_len: usize, exon2_len: usize) -> Gene {
        let exon2_start = exon1_len + INTRON_LEN + 1;
        spliced_gene(
            "t",
            Strand::Plus,
            &[(1, exon1_len), (exon2_start, exon2_start + exon2_len - 1)],
        )
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
        let gene = transcript_gene("t", Strand::Plus, &[(1, 120)]);
        assert_eq!(snv_mnv_nmd_prediction(&gene, 30), None);

        // Written by hand on purpose: no annotation file produces a CDS record
        // with zero segments, so this one cannot come from the parser.
        let empty = Gene {
            cds_segments: Vec::new(),
            ..transcript_gene("t", Strand::Plus, &[(1, 120)])
        };
        assert_eq!(snv_mnv_nmd_prediction(&empty, 30), None);
    }

    #[test]
    fn ribosomal_slippage_join_has_no_junction() {
        // A CDS annotated as a slippage join (SARS-CoV-2 ORF1ab style) has two
        // segments but no intron, so there is no exon-exon junction and no NMD
        // prediction to make.
        // Abutting segments: 1-120 then 121-180.
        let abutting = joined_gene("t", Strand::Plus, &[(1, 120), (121, 180)]);
        assert_eq!(last_exon_junction_offset(&abutting), None);
        assert_eq!(snv_mnv_nmd_prediction(&abutting, 30), None);

        // Overlapping segments, as in `join(266..13468,13468..21555)`.
        let overlapping = joined_gene("t", Strand::Plus, &[(1, 120), (120, 180)]);
        assert_eq!(last_exon_junction_offset(&overlapping), None);
        assert_eq!(snv_mnv_nmd_prediction(&overlapping, 30), None);
    }

    #[test]
    fn junction_is_the_last_real_intron_not_the_last_segment() {
        // exon 1 (1-120), intron, exon 2 (201-260), then a slippage join
        // continuing at 261-320. The last junction is the start of exon 2 at CDS
        // offset 120, not the start of the trailing joined segment at offset 180.
        let gene = transcript_gene("t", Strand::Plus, &[(1, 120), (201, 260), (261, 320)]);
        assert_eq!(last_exon_junction_offset(&gene), Some(120));
        // 90 nt upstream of that junction still triggers NMD.
        assert_eq!(
            snv_mnv_nmd_prediction(&gene, 30),
            Some(NmdPrediction::Triggering)
        );
        // A stop past the junction escapes.
        assert_eq!(
            snv_mnv_nmd_prediction(&gene, 200),
            Some(NmdPrediction::Escaping)
        );
    }
}
