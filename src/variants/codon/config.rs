//! Indel-aware annotation configuration.

use crate::io::VcfPosition;

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
