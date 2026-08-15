//! VcfWriter INFO-entry construction for SNP and MNV records.

use super::*;

impl VcfWriter {
    pub(super) fn build_snp_entry(
        &self,
        variant: &VariantInfo,
        index: usize,
        aa: &str,
        metrics: SnpCallMetrics,
    ) -> AppResult<Option<VcfEntry>> {
        let filters = self.build_support_filters(SupportFilterInput {
            support_reads: metrics.support_reads,
            min_reads: self.min_snp_reads,
            depth: metrics.depth,
            min_frequency: self.min_snp_frequency,
            forward_reads: metrics.forward_reads,
            reverse_reads: metrics.reverse_reads,
            min_strand_reads: self.min_snp_strand_reads,
            strand_bias_p: metrics.strand_bias_p,
        });
        if !self.should_emit_record(&filters) {
            return Ok(None);
        }

        let pos = *get_required(&variant.positions, index, "positions", variant)?;
        let ref_base = get_required(&variant.ref_bases, index, "ref_bases", variant)?;
        let alt_base = get_required(&variant.base_changes, index, "base_changes", variant)?;
        let info = build_info_string(
            variant,
            Some(aa),
            VariantType::Snp.as_str(),
            Some((
                metrics.support_reads,
                metrics.forward_reads,
                metrics.reverse_reads,
            )),
            None,
            Some(index),
            Some(metrics.depth),
            Some(metrics.support_reads),
            if self.include_strand_bias_info {
                metrics.strand_bias_p
            } else {
                None
            },
            None,
            variant.original_info.as_deref(),
        );
        let filter = filter_value(&filters);
        Ok(Some((
            pos,
            vcf_entry_line(&variant.chrom, pos, ref_base, alt_base, &filter, &info),
        )))
    }

    pub(super) fn build_mnv_entry(
        &self,
        variant: &VariantInfo,
        reference_sequence: &str,
        aa: &str,
        metrics: MnvCallMetrics,
    ) -> AppResult<Option<VcfEntry>> {
        let filters = self.build_support_filters(SupportFilterInput {
            support_reads: metrics.support_reads,
            min_reads: self.min_mnv_reads,
            depth: metrics.depth,
            min_frequency: self.min_mnv_frequency,
            forward_reads: metrics.forward_reads,
            reverse_reads: metrics.reverse_reads,
            min_strand_reads: self.min_mnv_strand_reads,
            strand_bias_p: metrics.strand_bias_p,
        });
        if !self.should_emit_record(&filters) {
            return Ok(None);
        }

        let min_pos = *variant
            .positions
            .iter()
            .min()
            .ok_or_else(|| format!("Missing positions for {}", variant_context(variant)))?;
        let max_pos = *variant
            .positions
            .iter()
            .max()
            .ok_or_else(|| format!("Missing positions for {}", variant_context(variant)))?;
        let ref_region = reference_subsequence(reference_sequence, min_pos, max_pos, variant)?;
        let alt_region = build_alt_region(
            reference_sequence,
            &variant.positions,
            &variant.base_changes,
        )?;
        let info = build_info_string(
            variant,
            Some(aa),
            VariantType::Mnv.as_str(),
            None,
            Some((
                metrics.support_reads,
                metrics.forward_reads,
                metrics.reverse_reads,
            )),
            None,
            Some(metrics.depth),
            Some(metrics.support_reads),
            None,
            if self.include_strand_bias_info {
                metrics.strand_bias_p
            } else {
                None
            },
            variant.original_info.as_deref(),
        );
        let filter = filter_value(&filters);
        Ok(Some((
            min_pos,
            vcf_entry_line(
                &variant.chrom,
                min_pos,
                ref_region,
                &alt_region,
                &filter,
                &info,
            ),
        )))
    }
}
