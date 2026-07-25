//! Shared output utilities: INFO string building, Fisher exact strand-bias
//! test, VCF entry formatting, variant shape validation, and header generation.

use crate::error::AppResult;
use crate::variants::{VariantInfo, VariantType};
use std::io::Write;

/// INFO keys emitted by get_mnv itself. Original VCF fields with these names
/// are silently skipped when `--keep-original-info` is active to prevent
/// duplicate INFO keys in the output VCF.
fn is_reserved_info_key(key: &str) -> bool {
    matches!(
        key,
        "GENE"
            | "AA"
            | "CT"
            | "TYPE"
            | "SR"
            | "SRF"
            | "SRR"
            | "MR"
            | "MRF"
            | "MRR"
            | "SBP"
            | "MSBP"
            | "ODP"
            | "OFREQ"
            | "DP"
            | "FREQ"
            | "EC"
            | "COMP"
            | "ER"
            | "ERF"
            | "ERR"
            | "EDP"
            | "EFREQ"
            | "SO"
            | "IMPACT"
            | "GD"
            | "MNVSHIFT"
            | "DBS"
            | "MNVPS"
    )
}

/// Map a variant to a Sequence Ontology consequence term and its impact level
/// (`HIGH` / `MODERATE` / `LOW` / `MODIFIER`), following SnpEff/VEP conventions.
pub(crate) fn so_consequence(variant: &VariantInfo) -> (&'static str, &'static str) {
    use crate::variants::ChangeType;
    if variant.gene == "intergenic" {
        return ("intergenic_variant", "MODIFIER");
    }
    match variant.change_type {
        ChangeType::Synonymous => ("synonymous_variant", "LOW"),
        ChangeType::NonSynonymous => ("missense_variant", "MODERATE"),
        ChangeType::StartLost => ("start_lost", "HIGH"),
        ChangeType::StopGained | ChangeType::FrameshiftStopGained => ("stop_gained", "HIGH"),
        ChangeType::StopLost | ChangeType::FrameshiftStopLost => ("stop_lost", "HIGH"),
        ChangeType::FrameshiftIndel
        | ChangeType::FrameshiftSynonymous
        | ChangeType::FrameshiftNonSynonymous
        | ChangeType::FrameshiftUnknown => ("frameshift_variant", "HIGH"),
        // A frameshift after the stop codon does not alter the protein.
        ChangeType::FrameshiftDownstreamOfStop => ("coding_sequence_variant", "MODIFIER"),
        ChangeType::InFrameIndel => {
            let ref_len = variant.ref_bases.first().map_or(0, |b| b.len());
            let alt_len = variant.base_changes.first().map_or(0, |b| b.len());
            if alt_len > ref_len {
                ("inframe_insertion", "MODERATE")
            } else {
                ("inframe_deletion", "MODERATE")
            }
        }
        ChangeType::IndelOverlap | ChangeType::Unknown => ("coding_sequence_variant", "MODIFIER"),
    }
}

#[derive(Default)]
struct InfoBuilder {
    fields: Vec<(String, String)>,
}

impl InfoBuilder {
    fn push(&mut self, key: &str, value: impl ToString) {
        self.fields.push((key.to_string(), value.to_string()));
    }

    /// Push a free-text value (gene name, AA change, change type, event
    /// components), percent-encoding characters that are structurally reserved
    /// in a VCF INFO column so the value cannot corrupt downstream parsing.
    fn push_text(&mut self, key: &str, value: &str) {
        self.fields
            .push((key.to_string(), encode_info_value(value)));
    }

    fn build(self) -> String {
        self.fields
            .into_iter()
            .map(|(key, value)| {
                if value.is_empty() {
                    key
                } else {
                    format!("{key}={value}")
                }
            })
            .collect::<Vec<_>>()
            .join(";")
    }
}

/// Percent-encode the characters that are structurally reserved in a VCF INFO
/// value: the field separator `;`, the key/value separator `=`, the array
/// separator `,`, a literal `%`, and the line-structural tab/newline/CR. Spaces,
/// `:` and `|` are intentionally left intact so controlled labels (e.g.
/// "Stop gained") and `COMP` components (e.g. "SNV:28:G>T") stay readable, while
/// arbitrary GFF-derived gene names can no longer break the INFO column.
fn encode_info_value(value: &str) -> String {
    if !value
        .bytes()
        .any(|b| matches!(b, b'%' | b';' | b'=' | b',' | b'\t' | b'\n' | b'\r'))
    {
        return value.to_string();
    }
    let mut out = String::with_capacity(value.len() + 8);
    for c in value.chars() {
        match c {
            '%' => out.push_str("%25"),
            ';' => out.push_str("%3B"),
            '=' => out.push_str("%3D"),
            ',' => out.push_str("%2C"),
            '\t' => out.push_str("%09"),
            '\n' => out.push_str("%0A"),
            '\r' => out.push_str("%0D"),
            other => out.push(other),
        }
    }
    out
}

pub(crate) fn snp_aa_for_index(variant: &VariantInfo, index: usize) -> String {
    variant
        .snp_aa_changes
        .get(index)
        .cloned()
        .or_else(|| variant.aa_changes.first().cloned())
        .unwrap_or_else(|| "Unknown".to_string())
}

pub(crate) fn format_freq(value: f64) -> String {
    format!("{value:.4}")
}

pub(crate) fn allele_frequency(support_reads: usize, depth: usize) -> f64 {
    if depth > 0 {
        support_reads as f64 / depth as f64
    } else {
        0.0
    }
}

/// MNV phasing (linkage) support from BAM reads: the fraction of the reads
/// carrying the least-supported constituent SNV that also carry the full MNV
/// haplotype. `1.0` means every read with the rarer variant also carries the
/// other (a genuine co-occurring haplotype); low values indicate the SNVs
/// largely fall on different molecules (a same-codon coincidence rather than a
/// real MNV), which matters most for intra-host data. `None` unless this is a
/// multi-SNV MNV / SNP-MNV with per-SNV and MNV read counts available.
pub(crate) fn mnv_phasing_support(variant: &VariantInfo) -> Option<f64> {
    if variant.positions.len() < 2 {
        return None;
    }
    if !matches!(variant.variant_type, VariantType::Mnv | VariantType::SnpMnv) {
        return None;
    }
    let min_snp = variant.snp_reads.as_ref()?.iter().copied().min()?;
    if min_snp == 0 {
        return None;
    }
    let mnv_reads = variant.mnv_reads?;
    // A read in the MNV count carries both variants, so it is also counted in
    // every per-SNV count; the ratio is therefore in [0, 1]. Clamp defensively.
    Some((mnv_reads as f64 / min_snp as f64).min(1.0))
}

fn original_dp_for_index(variant: &VariantInfo, index: usize) -> Option<usize> {
    variant.original_dp.as_ref()?.get(index).copied()
}

fn original_freq_for_index(variant: &VariantInfo, index: usize) -> Option<f64> {
    variant.original_freq.as_ref()?.get(index).copied()
}

fn original_dp_list(variant: &VariantInfo) -> Option<String> {
    let values = variant.original_dp.as_ref()?;
    Some(
        values
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>()
            .join(","),
    )
}

fn original_freq_list(variant: &VariantInfo) -> Option<String> {
    let values = variant.original_freq.as_ref()?;
    Some(
        values
            .iter()
            .map(|v| format_freq(*v))
            .collect::<Vec<_>>()
            .join(","),
    )
}

fn push_original_metrics(builder: &mut InfoBuilder, variant: &VariantInfo, index: Option<usize>) {
    match index {
        Some(i) => {
            if let Some(dp) = original_dp_for_index(variant, i) {
                builder.push("ODP", dp);
            }
            if let Some(freq) = original_freq_for_index(variant, i) {
                builder.push("OFREQ", format_freq(freq));
            }
        }
        None => {
            if let Some(dp_values) = original_dp_list(variant) {
                builder.push("ODP", dp_values);
            }
            if let Some(freq_values) = original_freq_list(variant) {
                builder.push("OFREQ", freq_values);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_info_string(
    variant: &VariantInfo,
    aa: Option<&str>,
    variant_type: &str,
    snp_metrics: Option<(usize, usize, usize)>,
    mnv_metrics: Option<(usize, usize, usize)>,
    original_index: Option<usize>,
    depth: Option<usize>,
    support_reads: Option<usize>,
    snp_strand_bias_p: Option<f64>,
    mnv_strand_bias_p: Option<f64>,
    original_info: Option<&str>,
) -> String {
    let mut builder = InfoBuilder::default();
    builder.push_text("GENE", &variant.gene);
    if let Some(aa_change) = aa {
        builder.push_text("AA", aa_change);
    }
    builder.push_text("CT", &variant.change_type.to_string());
    builder.push_text("TYPE", variant_type);
    let (so_term, impact) = so_consequence(variant);
    builder.push_text("SO", so_term);
    builder.push_text("IMPACT", impact);
    if let Some(gd) = variant.annotations.grantham {
        builder.push("GD", gd);
    }
    if variant.annotations.consequence_shift != crate::variants::ConsequenceShift::NotApplicable {
        builder.push_text("MNVSHIFT", variant.annotations.consequence_shift.as_str());
    }
    if let Some(dbs) = variant.annotations.dbs_class.as_deref() {
        builder.push_text("DBS", dbs);
    }
    if let Some(phasing) = mnv_phasing_support(variant) {
        builder.push("MNVPS", format_freq(phasing));
    }

    if let Some((sr, srf, srr)) = snp_metrics {
        builder.push("SR", sr);
        builder.push("SRF", srf);
        builder.push("SRR", srr);
    }
    if let Some((mr, mrf, mrr)) = mnv_metrics {
        builder.push("MR", mr);
        builder.push("MRF", mrf);
        builder.push("MRR", mrr);
    }
    if let Some(sbp) = snp_strand_bias_p {
        builder.push("SBP", format!("{sbp:.6}"));
    }
    if let Some(msbp) = mnv_strand_bias_p {
        builder.push("MSBP", format!("{msbp:.6}"));
    }

    if let Some(event_class) = variant.event_class.as_deref() {
        builder.push_text("EC", event_class);
    }
    if !variant.event_components.is_empty() {
        builder.push_text("COMP", &variant.event_components.join("|"));
    }

    push_original_metrics(&mut builder, variant, original_index);

    if let (Some(dp), Some(support)) = (depth, support_reads) {
        let freq = allele_frequency(support, dp);
        builder.push("DP", dp);
        builder.push("FREQ", format_freq(freq));
    }

    if variant.variant_type == VariantType::Indel {
        if let (Some(reads), Some(forward), Some(reverse), Some(depth)) = (
            variant.mnv_reads,
            variant.mnv_forward_reads,
            variant.mnv_reverse_reads,
            variant.mnv_total_reads,
        ) {
            builder.push("ER", reads);
            builder.push("ERF", forward);
            builder.push("ERR", reverse);
            builder.push("EDP", depth);
            builder.push("EFREQ", format_freq(allele_frequency(reads, depth)));
        }
    }

    // Append original INFO fields from the input VCF (if --keep-original-info).
    // Skip keys that collide with get_mnv's own INFO fields to avoid duplicate
    // INFO keys (which violate the VCF spec and confuse downstream tools).
    if let Some(orig) = original_info {
        for part in orig.split(';') {
            let key = part.split('=').next().unwrap_or(part);
            if is_reserved_info_key(key) {
                continue;
            }
            if let Some(eq_pos) = part.find('=') {
                builder.push(&part[..eq_pos], &part[eq_pos + 1..]);
            } else {
                // Flag field (no value)
                builder.push(part, "");
            }
        }
    }

    builder.build()
}

pub(crate) fn write_sorted_vcf_entries(
    writer: &mut dyn Write,
    mut entries: Vec<(usize, String)>,
) -> AppResult<()> {
    entries.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    for (_, line) in entries {
        writeln!(writer, "{line}")?;
    }
    Ok(())
}

pub(crate) fn variant_context(variant: &VariantInfo) -> String {
    format!(
        "contig '{}' gene '{}' positions {:?}",
        variant.chrom, variant.gene, variant.positions
    )
}

pub(crate) fn get_mnv_depth_from_variant(variant: &VariantInfo, total_reads: &[usize]) -> usize {
    variant
        .mnv_total_reads
        .unwrap_or_else(|| total_reads.iter().sum::<usize>() / total_reads.len().max(1))
}

fn ensure_non_empty_positions(variant: &VariantInfo) -> AppResult<()> {
    if variant.positions.is_empty() {
        return Err(format!("Variant has no positions ({})", variant_context(variant)).into());
    }
    Ok(())
}

fn ensure_len(field: &str, expected: usize, actual: usize, variant: &VariantInfo) -> AppResult<()> {
    if expected != actual {
        return Err(format!(
            "Length mismatch for {}: expected {}, got {} ({})",
            field,
            expected,
            actual,
            variant_context(variant)
        )
        .into());
    }
    Ok(())
}

pub(crate) fn validate_variant_shape(variant: &VariantInfo) -> AppResult<()> {
    ensure_non_empty_positions(variant)?;
    let positions_len = variant.positions.len();
    ensure_len("ref_bases", positions_len, variant.ref_bases.len(), variant)?;
    ensure_len(
        "base_changes",
        positions_len,
        variant.base_changes.len(),
        variant,
    )?;
    Ok(())
}

pub(crate) fn get_required<'a, T>(
    values: &'a [T],
    index: usize,
    field: &str,
    variant: &VariantInfo,
) -> AppResult<&'a T> {
    values.get(index).ok_or_else(|| {
        format!(
            "Missing '{}' at index {} ({})",
            field,
            index,
            variant_context(variant)
        )
        .into()
    })
}

pub(crate) fn reference_subsequence<'a>(
    reference_sequence: &'a str,
    start_1based: usize,
    end_1based: usize,
    variant: &VariantInfo,
) -> AppResult<&'a str> {
    if start_1based == 0 || end_1based == 0 || start_1based > end_1based {
        return Err(format!(
            "Invalid reference slice bounds {}-{} ({})",
            start_1based,
            end_1based,
            variant_context(variant)
        )
        .into());
    }
    if end_1based > reference_sequence.len() {
        return Err(format!(
            "Reference slice out of bounds {}-{} for sequence length {} ({})",
            start_1based,
            end_1based,
            reference_sequence.len(),
            variant_context(variant)
        )
        .into());
    }
    Ok(&reference_sequence[(start_1based - 1)..end_1based])
}

// Strand-bias statistics moved to output/stats.rs
pub(crate) use super::stats::{mnv_strand_bias_p_value, snp_strand_bias_p_value};

pub(crate) fn filter_value(tags: &[&str]) -> String {
    if tags.is_empty() {
        "PASS".to_string()
    } else {
        tags.join(";")
    }
}

pub(crate) type VcfEntry = (usize, String);

pub(crate) struct SnpBamVectors<'a> {
    pub(crate) reads: &'a [usize],
    pub(crate) forward_reads: &'a [usize],
    pub(crate) reverse_reads: &'a [usize],
    pub(crate) total_reads: &'a [usize],
}

#[derive(Clone, Copy)]
pub(crate) struct SnpCallMetrics {
    pub(crate) support_reads: usize,
    pub(crate) forward_reads: usize,
    pub(crate) reverse_reads: usize,
    pub(crate) depth: usize,
    pub(crate) strand_bias_p: Option<f64>,
}

#[derive(Clone, Copy)]
pub(crate) struct MnvCallMetrics {
    pub(crate) support_reads: usize,
    pub(crate) forward_reads: usize,
    pub(crate) reverse_reads: usize,
    pub(crate) depth: usize,
    pub(crate) strand_bias_p: Option<f64>,
}

pub(crate) fn vcf_entry_line(
    chrom: &str,
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    filter: &str,
    info: &str,
) -> String {
    format!("{chrom}\t{pos}\t.\t{ref_allele}\t{alt_allele}\t.\t{filter}\t{info}")
}

pub(crate) fn snp_bam_vectors(variant: &VariantInfo) -> AppResult<SnpBamVectors<'_>> {
    let reads = variant
        .snp_reads
        .as_ref()
        .ok_or_else(|| format!("Missing SNP read counts for {}", variant_context(variant)))?;
    let total_reads = variant
        .total_reads
        .as_ref()
        .ok_or_else(|| format!("Missing total read depth for {}", variant_context(variant)))?;
    let forward_reads = variant.snp_forward_reads.as_ref().ok_or_else(|| {
        format!(
            "Missing SNP forward read counts for {}",
            variant_context(variant)
        )
    })?;
    let reverse_reads = variant.snp_reverse_reads.as_ref().ok_or_else(|| {
        format!(
            "Missing SNP reverse read counts for {}",
            variant_context(variant)
        )
    })?;

    ensure_len("snp_reads", variant.positions.len(), reads.len(), variant)?;
    ensure_len(
        "snp_forward_reads",
        variant.positions.len(),
        forward_reads.len(),
        variant,
    )?;
    ensure_len(
        "snp_reverse_reads",
        variant.positions.len(),
        reverse_reads.len(),
        variant,
    )?;
    ensure_len(
        "total_reads",
        variant.positions.len(),
        total_reads.len(),
        variant,
    )?;

    Ok(SnpBamVectors {
        reads,
        forward_reads,
        reverse_reads,
        total_reads,
    })
}

pub(crate) fn build_alt_region(
    reference_sequence: &str,
    positions: &[usize],
    base_changes: &[String],
) -> AppResult<String> {
    if positions.is_empty() {
        return Err("positions must not be empty when building ALT region".into());
    }
    if positions.len() != base_changes.len() {
        return Err(format!(
            "positions/base_changes length mismatch when building ALT region: {} vs {}",
            positions.len(),
            base_changes.len()
        )
        .into());
    }
    let min_pos = *positions
        .iter()
        .min()
        .ok_or_else(|| "positions must not be empty when building ALT region".to_string())?;
    let max_pos = *positions
        .iter()
        .max()
        .ok_or_else(|| "positions must not be empty when building ALT region".to_string())?;
    let mut alt_region = String::new();
    let mut current_pos = min_pos;

    let mut pos_changes: Vec<(usize, &String)> =
        positions.iter().copied().zip(base_changes.iter()).collect();
    pos_changes.sort_by_key(|&(p, _)| p);

    for (p, alt) in pos_changes {
        if p > current_pos {
            if p - 1 > reference_sequence.len() || current_pos == 0 {
                return Err("Reference bounds exceeded when building ALT region".into());
            }
            alt_region.push_str(&reference_sequence[(current_pos - 1)..(p - 1)]);
        }
        alt_region.push_str(alt);
        current_pos = p + 1;
    }

    if current_pos <= max_pos {
        if max_pos > reference_sequence.len() || current_pos == 0 {
            return Err("Reference bounds exceeded when finishing ALT region".into());
        }
        alt_region.push_str(&reference_sequence[(current_pos - 1)..max_pos]);
    }

    Ok(alt_region)
}

pub(crate) fn write_info_header(
    writer: &mut dyn Write,
    bam_provided: bool,
    include_strand_bias_info: bool,
    original_info_headers: &[String],
) -> AppResult<()> {
    writeln!(
        writer,
        "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=AA,Number=1,Type=String,Description=\"Amino acid change\">"
    )?;
    writeln!(writer, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type (SNP, MNV, SNP/MNV, INDEL)\">")?;
    writeln!(
        writer,
        "##INFO=<ID=CT,Number=1,Type=String,Description=\"Change type\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=SO,Number=1,Type=String,Description=\"Sequence Ontology consequence term\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Predicted impact (HIGH, MODERATE, LOW, MODIFIER)\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=GD,Number=1,Type=Integer,Description=\"Grantham distance of the (combined) missense change\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=MNVSHIFT,Number=1,Type=String,Description=\"Combined MNV consequence vs. individual SNVs (MNV-gained / MNV-masked / Concordant)\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=DBS,Number=1,Type=String,Description=\"COSMIC-style doublet base substitution class for adjacent 2-SNV MNVs (reverse-complement collapsed)\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=MNVPS,Number=1,Type=Float,Description=\"MNV phasing support: fraction of the least-supported constituent SNV reads that also carry the full MNV haplotype\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=EC,Number=1,Type=String,Description=\"Canonical allele event class derived from REF/ALT\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=COMP,Number=.,Type=String,Description=\"Allele event components derived from REF/ALT\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=ODP,Number=.,Type=Integer,Description=\"Original depth from input VCF\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=OFREQ,Number=.,Type=Float,Description=\"Original allele frequency from input VCF\">"
    )?;
    if bam_provided {
        writeln!(
            writer,
            "##INFO=<ID=SR,Number=.,Type=Integer,Description=\"SNP read counts per position\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=SRF,Number=.,Type=Integer,Description=\"SNP forward read counts per position\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=SRR,Number=.,Type=Integer,Description=\"SNP reverse read counts per position\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=MR,Number=1,Type=Integer,Description=\"MNV read count\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=MRF,Number=1,Type=Integer,Description=\"MNV forward read count\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=MRR,Number=1,Type=Integer,Description=\"MNV reverse read count\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Calculated depth from BAM at the site\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=FREQ,Number=1,Type=Float,Description=\"Calculated allele frequency from BAM at the site\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=ER,Number=1,Type=Integer,Description=\"Exact event read count for indel/complex alleles\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=ERF,Number=1,Type=Integer,Description=\"Exact event forward read count for indel/complex alleles\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=ERR,Number=1,Type=Integer,Description=\"Exact event reverse read count for indel/complex alleles\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=EDP,Number=1,Type=Integer,Description=\"Exact event depth for indel/complex alleles\">"
        )?;
        writeln!(
            writer,
            "##INFO=<ID=EFREQ,Number=1,Type=Float,Description=\"Exact event allele frequency for indel/complex alleles\">"
        )?;
        if include_strand_bias_info {
            writeln!(
                writer,
                "##INFO=<ID=SBP,Number=1,Type=Float,Description=\"SNP strand bias Fisher exact two-tailed p-value\">"
            )?;
            writeln!(
                writer,
                "##INFO=<ID=MSBP,Number=1,Type=Float,Description=\"MNV strand bias Fisher exact two-tailed p-value\">"
            )?;
        }
    }
    for header_line in original_info_headers {
        writeln!(writer, "{header_line}")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::{ChangeType, VariantInfo, VariantType};

    fn make_snp_variant() -> VariantInfo {
        VariantInfo {
            chrom: "chr1".to_string(),
            gene: "geneA".to_string(),
            positions: vec![100],
            ref_bases: vec!["A".to_string()],
            base_changes: vec!["T".to_string()],
            aa_changes: vec!["Ala10Val".to_string()],
            snp_aa_changes: vec!["Ala10Val".to_string()],
            aa_changes_local: vec!["-".to_string()],
            snp_aa_changes_local: vec!["-".to_string()],
            variant_type: VariantType::Snp,
            change_type: ChangeType::NonSynonymous,
            snp_reads: Some(vec![5]),
            snp_forward_reads: Some(vec![3]),
            snp_reverse_reads: Some(vec![2]),
            mnv_reads: None,
            mnv_forward_reads: None,
            mnv_reverse_reads: None,
            mnv_total_reads: None,
            total_reads: Some(vec![20]),
            total_forward_reads: Some(vec![12]),
            total_reverse_reads: Some(vec![8]),
            mnv_total_forward_reads: None,
            mnv_total_reverse_reads: None,
            ref_codon: Some("ACG".to_string()),
            snp_codon: Some("TCG".to_string()),
            mnv_codon: None,
            original_dp: Some(vec![30]),
            original_freq: Some(vec![0.25]),
            original_info: None,
            event_class: Some("snp".to_string()),
            event_components: vec!["SNV:100:A>T".to_string()],
            annotations: crate::variants::VariantAnnotations::default(),
        }
    }

    // ---- format_freq ----

    #[test]
    fn test_format_freq_precision() {
        assert_eq!(format_freq(0.5), "0.5000");
        assert_eq!(format_freq(1.0), "1.0000");
        assert_eq!(format_freq(0.0), "0.0000");
        assert_eq!(format_freq(0.12345), "0.1235");
    }

    // ---- INFO value encoding ----

    #[test]
    fn test_encode_info_value_escapes_structural_chars() {
        assert_eq!(encode_info_value("a;b=c,d%e"), "a%3Bb%3Dc%2Cd%25e");
        assert_eq!(encode_info_value("x\ty\nz\r"), "x%09y%0Az%0D");
        // Spaces, ':' and '|' stay readable.
        assert_eq!(encode_info_value("Stop gained"), "Stop gained");
        assert_eq!(
            encode_info_value("SNV:28:G>T|INS:29:+GCT"),
            "SNV:28:G>T|INS:29:+GCT"
        );
        // No reserved chars → unchanged (fast path).
        assert_eq!(encode_info_value("geneA"), "geneA");
    }

    #[test]
    fn test_build_info_string_special_gene_name_does_not_corrupt_info() {
        let mut v = make_snp_variant();
        v.gene = "weird;gene=x,y".to_string();
        let info = build_info_string(
            &v,
            Some("Ala10Val"),
            "SNP",
            None,
            None,
            Some(0),
            None,
            None,
            None,
            None,
            None,
        );
        // The whole gene value stays inside the GENE field instead of spilling
        // into spurious INFO keys via its ';' and '='.
        let gene_field = info.split(';').next().unwrap();
        assert_eq!(gene_field, "GENE=weird%3Bgene%3Dx%2Cy", "got: {info}");
    }

    // ---- mnv_phasing_support ----

    #[test]
    fn test_mnv_phasing_support_ratio_and_guards() {
        let mut v = make_snp_variant();
        v.variant_type = VariantType::SnpMnv;
        v.positions = vec![100, 101];
        v.ref_bases = vec!["A".to_string(), "C".to_string()];
        v.base_changes = vec!["T".to_string(), "G".to_string()];
        v.snp_reads = Some(vec![8, 10]); // limiting SNV has 8 reads
        v.mnv_reads = Some(6); // 6 of those 8 carry the full haplotype
        assert_eq!(mnv_phasing_support(&v), Some(0.75));

        // Perfect linkage is clamped to 1.0.
        v.mnv_reads = Some(8);
        assert_eq!(mnv_phasing_support(&v), Some(1.0));

        // A single SNP is not an MNV, so phasing is not applicable.
        assert_eq!(mnv_phasing_support(&make_snp_variant()), None);
    }

    // ---- filter_value ----

    #[test]
    fn test_filter_value_pass_when_empty() {
        assert_eq!(filter_value(&[]), "PASS");
    }

    #[test]
    fn test_filter_value_single() {
        assert_eq!(filter_value(&["LowSupport"]), "LowSupport");
    }

    #[test]
    fn test_filter_value_multiple() {
        assert_eq!(
            filter_value(&["LowSupport", "StrandBias"]),
            "LowSupport;StrandBias"
        );
    }

    // ---- build_alt_region ----

    #[test]
    fn test_build_alt_region_single_snp() {
        let ref_seq = "ACGTACGT";
        let result = build_alt_region(ref_seq, &[3], &["T".to_string()]).unwrap();
        assert_eq!(result, "T");
    }

    #[test]
    fn test_build_alt_region_two_snps() {
        let ref_seq = "ACGTACGT";
        // positions 2 and 4 (1-based): C->T, T->A
        let result =
            build_alt_region(ref_seq, &[2, 4], &["T".to_string(), "A".to_string()]).unwrap();
        assert_eq!(result, "TGA");
    }

    #[test]
    fn test_build_alt_region_adjacent_snps() {
        let ref_seq = "ACGTACGT";
        let result =
            build_alt_region(ref_seq, &[1, 2], &["T".to_string(), "G".to_string()]).unwrap();
        assert_eq!(result, "TG");
    }

    #[test]
    fn test_build_alt_region_empty_positions_fails() {
        let result = build_alt_region("ACGT", &[], &[]);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_alt_region_mismatch_lengths_fails() {
        let result = build_alt_region("ACGT", &[1, 2], &["T".to_string()]);
        assert!(result.is_err());
    }

    // ---- build_info_string ----

    #[test]
    fn test_build_info_string_basic() {
        let variant = make_snp_variant();
        let info = build_info_string(
            &variant,
            Some("Ala10Val"),
            "SNP",
            Some((5, 3, 2)),
            None,
            Some(0),
            Some(20),
            Some(5),
            None,
            None,
            None,
        );
        assert!(info.contains("GENE=geneA"));
        assert!(info.contains("AA=Ala10Val"));
        assert!(info.contains("TYPE=SNP"));
        assert!(info.contains("SR=5"));
        assert!(info.contains("SRF=3"));
        assert!(info.contains("SRR=2"));
        assert!(info.contains("DP=20"));
        assert!(info.contains("FREQ=0.2500"));
        assert!(info.contains("ODP=30"));
        assert!(info.contains("OFREQ=0.2500"));
    }

    #[test]
    fn test_build_info_string_with_original_info() {
        let variant = make_snp_variant();
        let info = build_info_string(
            &variant,
            Some("Ala10Val"),
            "SNP",
            None,
            None,
            Some(0),
            None,
            None,
            None,
            None,
            Some("ANN=missense;CLNSIG=pathogenic"),
        );
        assert!(info.contains("ANN=missense"));
        assert!(info.contains("CLNSIG=pathogenic"));
    }

    #[test]
    fn test_build_info_string_skips_reserved_original_keys() {
        let variant = make_snp_variant();
        // DP and FREQ are reserved — should NOT appear from original_info
        let info = build_info_string(
            &variant,
            None,
            "SNP",
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some("DP=99;CUSTOM=yes"),
        );
        assert!(info.contains("CUSTOM=yes"));
        // DP from original should be skipped (it's reserved)
        assert!(!info.contains("DP=99"));
    }

    // ---- validate_variant_shape ----

    #[test]
    fn test_validate_variant_shape_ok() {
        let variant = make_snp_variant();
        assert!(validate_variant_shape(&variant).is_ok());
    }

    #[test]
    fn test_validate_variant_shape_empty_positions() {
        let mut variant = make_snp_variant();
        variant.positions = vec![];
        assert!(validate_variant_shape(&variant).is_err());
    }

    #[test]
    fn test_validate_variant_shape_mismatched_lengths() {
        let mut variant = make_snp_variant();
        variant.ref_bases = vec!["A".to_string(), "C".to_string()];
        assert!(validate_variant_shape(&variant).is_err());
    }

    // ---- snp_aa_for_index ----

    #[test]
    fn test_snp_aa_for_index_valid() {
        let variant = make_snp_variant();
        assert_eq!(snp_aa_for_index(&variant, 0), "Ala10Val");
    }

    #[test]
    fn test_snp_aa_for_index_out_of_range_falls_back() {
        let variant = make_snp_variant();
        // Index 5 doesn't exist, falls back to aa_changes[0]
        assert_eq!(snp_aa_for_index(&variant, 5), "Ala10Val");
    }

    // ---- reference_subsequence ----

    #[test]
    fn test_reference_subsequence_valid() {
        let variant = make_snp_variant();
        let seq = reference_subsequence("ACGTACGT", 2, 5, &variant).unwrap();
        assert_eq!(seq, "CGTA");
    }

    #[test]
    fn test_reference_subsequence_out_of_bounds() {
        let variant = make_snp_variant();
        let result = reference_subsequence("ACGT", 1, 10, &variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_reference_subsequence_zero_start() {
        let variant = make_snp_variant();
        let result = reference_subsequence("ACGT", 0, 2, &variant);
        assert!(result.is_err());
    }
}
