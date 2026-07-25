//! Indel-aware annotation configuration.

use crate::io::VcfPosition;
use std::collections::HashSet;

/// BAM-derived phasing evidence for frameshift propagation: the set of
/// `(indel_position, snv_position)` pairs that reads show are in **trans** (on
/// different molecules). An upstream indel's frame shift is not propagated to a
/// downstream codon when the indel is in trans with that codon's SNV, since the
/// codon's molecule does not carry the indel. Empty (the default, and whenever
/// no BAM is available) means no evidence, so every pair keeps the
/// frequency-based behaviour. This is a suppression-only signal that never adds
/// propagation the frequency gate would not.
#[derive(Debug, Clone, Default)]
pub struct FrameshiftPhasing {
    trans_pairs: HashSet<(usize, usize)>,
}

impl FrameshiftPhasing {
    /// Build from the set of `(indel_position, snv_position)` trans pairs.
    pub fn from_trans_pairs(trans_pairs: HashSet<(usize, usize)>) -> Self {
        Self { trans_pairs }
    }

    /// Whether the indel at `indel_position` is confirmed in trans with any of a
    /// codon's SNV positions, so its frame shift must not propagate to that codon.
    pub(super) fn indel_in_trans_with(
        &self,
        indel_position: usize,
        snv_positions: &[usize],
    ) -> bool {
        !self.trans_pairs.is_empty()
            && snv_positions
                .iter()
                .any(|snv| self.trans_pairs.contains(&(indel_position, *snv)))
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
