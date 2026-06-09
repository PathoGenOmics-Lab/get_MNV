//! VcfWriter per-record writers (indel, SNP, MNV, SNP+MNV, intergenic).

use super::*;

impl VcfWriter {
    pub(super) fn write_indel(&mut self, variant: &VariantInfo) -> AppResult<()> {
        validate_variant_shape(variant)?;
        let ref_base = get_required(&variant.ref_bases, 0, "ref_bases", variant)?;
        let alt_base = get_required(&variant.base_changes, 0, "base_changes", variant)?;
        let pos = *get_required(&variant.positions, 0, "positions", variant)?;
        let filters = if self.bam_provided {
            self.build_support_filters(SupportFilterInput {
                support_reads: variant.mnv_reads.unwrap_or(0),
                min_reads: self.min_mnv_reads,
                depth: variant.mnv_total_reads.unwrap_or(0),
                min_frequency: self.min_mnv_frequency,
                forward_reads: variant.mnv_forward_reads.unwrap_or(0),
                reverse_reads: variant.mnv_reverse_reads.unwrap_or(0),
                min_strand_reads: self.min_mnv_strand_reads,
                strand_bias_p: None,
            })
        } else {
            Vec::new()
        };
        if !self.should_emit_record(&filters) {
            return Ok(());
        }
        let info = build_info_string(
            variant,
            variant.aa_changes.first().map(String::as_str),
            VariantType::Indel.as_str(),
            None,
            None,
            Some(0),
            None,
            None,
            None,
            None,
            variant.original_info.as_deref(),
        );
        let filter = filter_value(&filters);
        self.write_variant_line(&variant.chrom, pos, ref_base, alt_base, &filter, &info)
    }

    pub(super) fn write_snp(&mut self, variant: &VariantInfo) -> AppResult<()> {
        validate_variant_shape(variant)?;
        if self.bam_provided {
            let bam_vectors = snp_bam_vectors(variant)?;
            for i in 0..variant.positions.len() {
                let metrics = self.snp_metrics_at(variant, &bam_vectors, i)?;
                let aa = variant.aa_changes.join(",");
                if let Some((_, line)) = self.build_snp_entry(variant, i, &aa, metrics)? {
                    writeln!(self.writer, "{line}")?;
                }
            }
            Ok(())
        } else {
            for (i, &pos) in variant.positions.iter().enumerate() {
                let ref_base = get_required(&variant.ref_bases, i, "ref_bases", variant)?;
                let alt_base = get_required(&variant.base_changes, i, "base_changes", variant)?;
                let aa = variant.aa_changes.join(",");
                let info = build_info_string(
                    variant,
                    Some(&aa),
                    VariantType::Snp.as_str(),
                    None,
                    None,
                    Some(i),
                    None,
                    None,
                    None,
                    None,
                    variant.original_info.as_deref(),
                );
                self.write_variant_line(&variant.chrom, pos, ref_base, alt_base, "PASS", &info)?;
            }
            Ok(())
        }
    }

    pub(super) fn write_mnv(
        &mut self,
        variant: &VariantInfo,
        reference_sequence: &str,
    ) -> AppResult<()> {
        validate_variant_shape(variant)?;
        if self.bam_provided {
            let total_reads = variant.total_reads.as_ref().ok_or_else(|| {
                format!("Missing total read depth for {}", variant_context(variant))
            })?;
            let metrics = self.mnv_metrics(variant, total_reads)?;
            let aa = variant.aa_changes.join(",");
            if let Some((_, line)) =
                self.build_mnv_entry(variant, reference_sequence, &aa, metrics)?
            {
                writeln!(self.writer, "{line}")?;
            }
            Ok(())
        } else {
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
            let aa = variant.aa_changes.join(",");
            let info = build_info_string(
                variant,
                Some(&aa),
                VariantType::Mnv.as_str(),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                variant.original_info.as_deref(),
            );
            self.write_variant_line(
                &variant.chrom,
                min_pos,
                ref_region,
                &alt_region,
                "PASS",
                &info,
            )
        }
    }

    pub(super) fn write_snp_mnv(
        &mut self,
        variant: &VariantInfo,
        reference_sequence: &str,
    ) -> AppResult<()> {
        validate_variant_shape(variant)?;
        let mut entries: Vec<VcfEntry> = Vec::new();

        if self.bam_provided {
            let bam_vectors = snp_bam_vectors(variant)?;
            for i in 0..variant.positions.len() {
                let metrics = self.snp_metrics_at(variant, &bam_vectors, i)?;
                if metrics.support_reads == 0 {
                    continue;
                }
                let aa = snp_aa_for_index(variant, i);
                if let Some(entry) = self.build_snp_entry(variant, i, &aa, metrics)? {
                    entries.push(entry);
                }
            }

            let mnv_metrics = self.mnv_metrics(variant, bam_vectors.total_reads)?;
            if mnv_metrics.support_reads > 0 {
                let aa = variant.aa_changes.join(",");
                if let Some(entry) =
                    self.build_mnv_entry(variant, reference_sequence, &aa, mnv_metrics)?
                {
                    entries.push(entry);
                }
            }
        } else {
            for (i, &pos) in variant.positions.iter().enumerate() {
                let ref_base = get_required(&variant.ref_bases, i, "ref_bases", variant)?;
                let alt_base = get_required(&variant.base_changes, i, "base_changes", variant)?;
                let aa = snp_aa_for_index(variant, i);
                let info = build_info_string(
                    variant,
                    Some(&aa),
                    VariantType::Snp.as_str(),
                    None,
                    None,
                    Some(i),
                    None,
                    None,
                    None,
                    None,
                    variant.original_info.as_deref(),
                );
                entries.push((
                    pos,
                    vcf_entry_line(&variant.chrom, pos, ref_base, alt_base, "PASS", &info),
                ));
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

            let aa = variant.aa_changes.join(",");
            let info = build_info_string(
                variant,
                Some(&aa),
                VariantType::Mnv.as_str(),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                variant.original_info.as_deref(),
            );
            entries.push((
                min_pos,
                vcf_entry_line(
                    &variant.chrom,
                    min_pos,
                    ref_region,
                    &alt_region,
                    "PASS",
                    &info,
                ),
            ));
        }

        write_sorted_vcf_entries(&mut self.writer, entries)
    }

    pub(super) fn write_intergenic(&mut self, variant: &VariantInfo) -> AppResult<()> {
        validate_variant_shape(variant)?;
        let ref_base = get_required(&variant.ref_bases, 0, "ref_bases", variant)?;
        let alt_base = get_required(&variant.base_changes, 0, "base_changes", variant)?;
        let pos = *get_required(&variant.positions, 0, "positions", variant)?;

        // Intergenic SNPs are read-counted at their position, so they are
        // filtered by their real support exactly like any SNP. Intergenic
        // indels carry no read counts and are always emitted.
        let snp_support = if self.bam_provided && variant.variant_type == VariantType::Snp {
            let support = variant.snp_reads.as_ref().and_then(|c| c.first()).copied();
            let forward = variant
                .snp_forward_reads
                .as_ref()
                .and_then(|c| c.first())
                .copied();
            let reverse = variant
                .snp_reverse_reads
                .as_ref()
                .and_then(|c| c.first())
                .copied();
            let depth = variant.total_reads.as_ref().and_then(|c| c.first()).copied();
            match (support, forward, reverse, depth) {
                (Some(s), Some(f), Some(r), Some(d)) => Some((s, f, r, d)),
                _ => None,
            }
        } else {
            None
        };

        // Strand-bias p-value, computed like the genic SNP path so the
        // --min-strand-bias-p filter and the SBP INFO field apply to intergenic
        // SNPs too (the strand read counts are populated for them).
        let strand_bias_p = if snp_support.is_some() {
            self.snp_strand_bias(variant, 0)
        } else {
            None
        };
        let filters = match snp_support {
            Some((support, forward, reverse, depth)) => {
                self.build_support_filters(SupportFilterInput {
                    support_reads: support,
                    min_reads: self.min_snp_reads,
                    depth,
                    min_frequency: self.min_snp_frequency,
                    forward_reads: forward,
                    reverse_reads: reverse,
                    min_strand_reads: self.min_snp_strand_reads,
                    strand_bias_p,
                })
            }
            None if self.bam_provided && variant.variant_type != VariantType::Snp => {
                // Intergenic indels carry no recomputed read support, so they
                // cannot satisfy a read-based filter; drop them when the relevant
                // MNV filters are active, matching the TSV output (support = 0).
                self.build_support_filters(SupportFilterInput {
                    support_reads: 0,
                    min_reads: self.min_mnv_reads,
                    depth: 0,
                    min_frequency: self.min_mnv_frequency,
                    forward_reads: 0,
                    reverse_reads: 0,
                    min_strand_reads: self.min_mnv_strand_reads,
                    strand_bias_p: None,
                })
            }
            None => Vec::new(),
        };
        if !self.should_emit_record(&filters) {
            return Ok(());
        }

        let info = build_info_string(
            variant,
            None,
            variant.variant_type.as_str(),
            snp_support.map(|(s, f, r, _)| (s, f, r)),
            None,
            Some(0),
            snp_support.map(|(_, _, _, d)| d),
            snp_support.map(|(s, _, _, _)| s),
            if self.include_strand_bias_info {
                strand_bias_p
            } else {
                None
            },
            None,
            variant.original_info.as_deref(),
        );
        let filter = filter_value(&filters);
        self.write_variant_line(&variant.chrom, pos, ref_base, alt_base, &filter, &info)
    }
}
