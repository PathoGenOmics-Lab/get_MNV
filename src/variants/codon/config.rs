//! Indel-aware annotation configuration.

use crate::io::VcfPosition;
use crate::variants::types::{FrameshiftLinkage, LinkageVerdict};
use std::collections::HashMap;

/// The reads' answer for one pair, already judged. The cis/trans thresholds
/// live with the read counting that applies them, so nothing here re-derives
/// the rule from the counts.
#[derive(Debug, Clone, Copy)]
pub struct PairLinkage {
    pub verdict: LinkageVerdict,
    pub cis_reads: usize,
    pub informative_reads: usize,
}

/// BAM-derived phasing evidence for frameshift propagation, keyed by
/// `(indel_position, snv_position, snv_alt)`. An upstream indel's frame shift
/// is not propagated to a downstream codon when the reads show the indel is in
/// **trans** with that codon's substitution, since the codon's molecule does
/// not carry it. Empty (the default, and whenever no BAM is available) means
/// nobody was asked, so every pair keeps the frequency-based behaviour. This is
/// a suppression-only signal that never adds propagation the frequency gate
/// would not.
#[derive(Debug, Clone, Default)]
pub struct FrameshiftPhasing {
    pairs: HashMap<(usize, usize, char), PairLinkage>,
}

impl FrameshiftPhasing {
    /// Build from `(indel_position, snv_position, snv_alt)` to the reads'
    /// answer, for every pair they were consulted about. The ALT belongs in the
    /// key: the linkage is queried per allele, and at a multi-allelic site a
    /// position-only key let the last ALT processed decide the other's codon.
    pub fn from_pairs(pairs: HashMap<(usize, usize, char), PairLinkage>) -> Self {
        Self { pairs }
    }

    /// Whether the indel at `indel_position` is confirmed in trans with any of a
    /// codon's SNV positions, so its frame shift must not propagate to that codon.
    pub(super) fn indel_in_trans_with(
        &self,
        indel_position: usize,
        snv_alleles: &[(usize, char)],
    ) -> bool {
        snv_alleles.iter().any(|(position, alt)| {
            self.pairs
                .get(&(indel_position, *position, *alt))
                .is_some_and(|linkage| linkage.verdict == LinkageVerdict::Trans)
        })
    }

    /// What the reads said about this indel and this codon, for the output. The
    /// SNV whose answer drove the decision is the one reported, so the column
    /// explains the label the row carries. `None` when the reads were never
    /// consulted about this pair, which is not the same as their having found
    /// nothing.
    pub(super) fn linkage_for(
        &self,
        indel_position: usize,
        snv_alleles: &[(usize, char)],
    ) -> Option<FrameshiftLinkage> {
        snv_alleles
            .iter()
            .filter_map(|(position, alt)| {
                self.pairs.get(&(indel_position, *position, *alt)).copied()
            })
            .min_by_key(|linkage| (linkage.verdict, linkage.cis_reads))
            .map(|linkage| FrameshiftLinkage {
                indel_position,
                cis_reads: linkage.cis_reads,
                informative_reads: linkage.informative_reads,
                verdict: linkage.verdict,
            })
    }
}

/// Knobs for indel-aware annotation that can change scientific output. The
/// `Default` impl reproduces the historical behaviour exactly, so callers that
/// do not opt in (tests, benchmarks, the public `get_mnv_variants_for_gene`
/// wrapper) see no change.
#[derive(Debug, Clone, Copy)]
pub struct IndelAnnotationConfig {
    /// Minimum allele frequency an *upstream* indel must reach to contribute to
    /// downstream frameshift propagation. `0.0` (default) propagates from every
    /// indel regardless of frequency, matching the original behaviour. Raising
    /// it avoids relabelling high-frequency downstream SNV/MNV codons as
    /// frameshifted because of a low-frequency upstream indel that is almost
    /// certainly on a different molecule (relevant for intra-host data).
    pub frameshift_min_freq: f64,
}

impl Default for IndelAnnotationConfig {
    fn default() -> Self {
        Self {
            frameshift_min_freq: 0.0,
        }
    }
}

/// Whether an upstream indel is allowed to contribute to downstream frameshift
/// propagation under the configured frequency gate. Indels without a known
/// frequency always pass (we cannot filter what we cannot measure).
pub(super) fn indel_passes_frameshift_gate(
    indel: &VcfPosition,
    config: &IndelAnnotationConfig,
) -> bool {
    match indel.original_freq {
        Some(freq) => freq >= config.frameshift_min_freq,
        None => true,
    }
}
