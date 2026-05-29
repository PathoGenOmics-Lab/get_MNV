//! VcfWriter strand-bias, support-filter, and per-allele metric helpers.

use super::*;

impl VcfWriter {
    pub(super) fn should_evaluate_strand_bias(&self) -> bool {
        self.include_strand_bias_info || self.min_strand_bias_p > 0.0
    }

    pub(super) fn snp_strand_bias(&self, variant: &VariantInfo, index: usize) -> Option<f64> {
        if self.should_evaluate_strand_bias() {
            snp_strand_bias_p_value(variant, index)
        } else {
            None
        }
    }

    pub(super) fn mnv_strand_bias(&self, variant: &VariantInfo) -> Option<f64> {
        if self.should_evaluate_strand_bias() {
            mnv_strand_bias_p_value(variant)
        } else {
            None
        }
    }

    pub(super) fn build_support_filters(&self, input: SupportFilterInput) -> Vec<&'static str> {
        let mut filters = Vec::new();
        if input.support_reads < input.min_reads {
            filters.push("LowSupport");
        }
        if input.min_frequency > 0.0
            && allele_frequency(input.support_reads, input.depth) < input.min_frequency
        {
            filters.push("LowFrequency");
        }
        if input.forward_reads < input.min_strand_reads
            || input.reverse_reads < input.min_strand_reads
        {
            filters.push("StrandSupport");
        }
        if let Some(p_value) = input.strand_bias_p {
            if self.min_strand_bias_p > 0.0 && p_value < self.min_strand_bias_p {
                filters.push("StrandBias");
            }
        }
        filters
    }

    pub(super) fn should_emit_record(&self, filters: &[&str]) -> bool {
        filters.is_empty() || self.emit_filtered
    }

    pub(super) fn snp_metrics_at(
        &self,
        variant: &VariantInfo,
        bam_vectors: &SnpBamVectors<'_>,
        index: usize,
    ) -> AppResult<SnpCallMetrics> {
        Ok(SnpCallMetrics {
            support_reads: *get_required(bam_vectors.reads, index, "snp_reads", variant)?,
            forward_reads: *get_required(
                bam_vectors.forward_reads,
                index,
                "snp_forward_reads",
                variant,
            )?,
            reverse_reads: *get_required(
                bam_vectors.reverse_reads,
                index,
                "snp_reverse_reads",
                variant,
            )?,
            depth: *get_required(bam_vectors.total_reads, index, "total_reads", variant)?,
            strand_bias_p: self.snp_strand_bias(variant, index),
        })
    }

    pub(super) fn mnv_metrics(
        &self,
        variant: &VariantInfo,
        total_reads: &[usize],
    ) -> AppResult<MnvCallMetrics> {
        let support_reads = variant
            .mnv_reads
            .ok_or_else(|| format!("Missing MNV read counts for {}", variant_context(variant)))?;
        let forward_reads = variant.mnv_forward_reads.ok_or_else(|| {
            format!(
                "Missing MNV forward read counts for {}",
                variant_context(variant)
            )
        })?;
        let reverse_reads = variant.mnv_reverse_reads.ok_or_else(|| {
            format!(
                "Missing MNV reverse read counts for {}",
                variant_context(variant)
            )
        })?;

        Ok(MnvCallMetrics {
            support_reads,
            forward_reads,
            reverse_reads,
            depth: get_mnv_depth_from_variant(variant, total_reads),
            strand_bias_p: self.mnv_strand_bias(variant),
        })
    }
}
