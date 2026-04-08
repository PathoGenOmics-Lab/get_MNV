//! TSV output writer: generates tab-delimited variant tables with optional
//! BAM-derived read support columns.

use crate::error::AppResult;
use crate::variants::{VariantInfo, VariantType};
use csv::WriterBuilder;
use std::fs::File;

use super::common::{
    format_freq, get_mnv_depth_from_variant, snp_bam_vectors, validate_variant_shape,
};

fn is_intergenic(variant: &VariantInfo) -> bool {
    variant.gene == "intergenic"
}

/// Render a "Local …" column. Falls back to the protein-wide column when the
/// local vector is empty (e.g. records produced before this field existed and
/// then deserialized via `#[serde(default)]`).
fn local_aa_or_fallback(local: &[String], protein: &[String]) -> String {
    let chosen = if local.is_empty() { protein } else { local };
    chosen.join(", ")
}

fn build_tsv_row_with_reads(variant: &VariantInfo) -> AppResult<Vec<String>> {
    validate_variant_shape(variant)?;
    if variant.variant_type == VariantType::Indel || is_intergenic(variant) {
        let pos_str = variant
            .positions
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        let ref_base_str = variant.ref_bases.join(", ");
        let base_str = variant.base_changes.join(", ");
        return Ok(vec![
            variant.chrom.clone(),
            variant.gene.clone(),
            pos_str,
            ref_base_str,
            base_str,
            variant.aa_changes.join(", "),
            variant.snp_aa_changes.join(", "),
            local_aa_or_fallback(&variant.aa_changes_local, &variant.aa_changes),
            local_aa_or_fallback(&variant.snp_aa_changes_local, &variant.snp_aa_changes),
            variant.variant_type.to_string(),
            variant.change_type.to_string(),
            variant.ref_codon.clone().unwrap_or_else(|| "-".to_string()),
            variant.snp_codon.clone().unwrap_or_else(|| "-".to_string()),
            variant.mnv_codon.clone().unwrap_or_else(|| "-".to_string()),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
        ]);
    }

    let bam_vectors = snp_bam_vectors(variant)?;
    let snp_counts = bam_vectors.reads;
    let snp_forward_counts = bam_vectors.forward_reads;
    let snp_reverse_counts = bam_vectors.reverse_reads;
    let total_reads = bam_vectors.total_reads;

    let snp_freq: Vec<String> = snp_counts
        .iter()
        .zip(total_reads.iter())
        .map(|(&s, &t)| {
            if t > 0 {
                format_freq(s as f64 / t as f64)
            } else {
                format_freq(0.0)
            }
        })
        .collect();

    let dp_mean: usize = total_reads.iter().sum::<usize>() / total_reads.len().max(1);
    let mnv_depth = get_mnv_depth_from_variant(variant, total_reads);

    let mnv_freq_str = if mnv_depth > 0 {
        format_freq(variant.mnv_reads.unwrap_or(0) as f64 / mnv_depth as f64)
    } else {
        format_freq(0.0)
    };

    let total_str = if variant.variant_type == VariantType::Snp {
        dp_mean.to_string()
    } else {
        mnv_depth.to_string()
    };

    let (ref_cod, snp_cod, mnv_cod) = if variant.variant_type == VariantType::Snp {
        (
            variant.ref_codon.clone().unwrap_or_default(),
            variant.snp_codon.clone().unwrap_or_default(),
            String::new(),
        )
    } else {
        (
            variant.ref_codon.clone().unwrap_or_default(),
            variant.snp_codon.clone().unwrap_or_default(),
            variant.mnv_codon.clone().unwrap_or_default(),
        )
    };

    let pos_str = variant
        .positions
        .iter()
        .map(std::string::ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    let ref_base_str = variant.ref_bases.join(", ");
    let base_str = variant.base_changes.join(", ");
    let aa_str = variant.aa_changes.join(", ");
    let snp_aa_str = variant.snp_aa_changes.join(", ");
    let local_aa_str = local_aa_or_fallback(&variant.aa_changes_local, &variant.aa_changes);
    let local_snp_aa_str =
        local_aa_or_fallback(&variant.snp_aa_changes_local, &variant.snp_aa_changes);
    let snp_reads_str = snp_counts
        .iter()
        .map(std::string::ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    let snp_forward_reads_str = snp_forward_counts
        .iter()
        .map(std::string::ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    let snp_reverse_reads_str = snp_reverse_counts
        .iter()
        .map(std::string::ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    let mnv_reads_str = variant.mnv_reads.unwrap_or(0).to_string();
    let mnv_forward_reads_str = variant.mnv_forward_reads.unwrap_or(0).to_string();
    let mnv_reverse_reads_str = variant.mnv_reverse_reads.unwrap_or(0).to_string();

    Ok(vec![
        variant.chrom.clone(),
        variant.gene.clone(),
        pos_str,
        ref_base_str,
        base_str,
        aa_str,
        snp_aa_str,
        local_aa_str,
        local_snp_aa_str,
        variant.variant_type.to_string(),
        variant.change_type.to_string(),
        ref_cod,
        snp_cod,
        mnv_cod,
        snp_reads_str,
        snp_forward_reads_str,
        snp_reverse_reads_str,
        mnv_reads_str,
        mnv_forward_reads_str,
        mnv_reverse_reads_str,
        total_str,
        snp_freq.join(", "),
        mnv_freq_str,
    ])
}

fn build_tsv_row_without_reads(variant: &VariantInfo) -> AppResult<Vec<String>> {
    validate_variant_shape(variant)?;
    if variant.variant_type == VariantType::Indel || is_intergenic(variant) {
        let pos_str = variant
            .positions
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        let ref_base_str = variant.ref_bases.join(", ");
        let base_str = variant.base_changes.join(", ");
        return Ok(vec![
            variant.chrom.clone(),
            variant.gene.clone(),
            pos_str,
            ref_base_str,
            base_str,
            variant.aa_changes.join(", "),
            variant.snp_aa_changes.join(", "),
            local_aa_or_fallback(&variant.aa_changes_local, &variant.aa_changes),
            local_aa_or_fallback(&variant.snp_aa_changes_local, &variant.snp_aa_changes),
            variant.variant_type.to_string(),
            variant.change_type.to_string(),
            variant.ref_codon.clone().unwrap_or_else(|| "-".to_string()),
            variant.snp_codon.clone().unwrap_or_else(|| "-".to_string()),
            variant.mnv_codon.clone().unwrap_or_else(|| "-".to_string()),
        ]);
    }

    let pos_str = variant
        .positions
        .iter()
        .map(std::string::ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    let ref_base_str = variant.ref_bases.join(", ");
    let base_str = variant.base_changes.join(", ");
    let aa_str = variant.aa_changes.join(", ");
    let snp_aa_str = variant.snp_aa_changes.join(", ");
    let local_aa_str = local_aa_or_fallback(&variant.aa_changes_local, &variant.aa_changes);
    let local_snp_aa_str =
        local_aa_or_fallback(&variant.snp_aa_changes_local, &variant.snp_aa_changes);
    Ok(vec![
        variant.chrom.clone(),
        variant.gene.clone(),
        pos_str,
        ref_base_str,
        base_str,
        aa_str,
        snp_aa_str,
        local_aa_str,
        local_snp_aa_str,
        variant.variant_type.to_string(),
        variant.change_type.to_string(),
        variant.ref_codon.clone().unwrap_or_default(),
        variant.snp_codon.clone().unwrap_or_default(),
        variant.mnv_codon.clone().unwrap_or_default(),
    ])
}

pub struct TsvWriter {
    writer: csv::Writer<File>,
    bam_provided: bool,
}

impl TsvWriter {
    pub fn new(filename: &str, bam_provided: bool) -> AppResult<Self> {
        let out_file = format!("{filename}.MNV.tsv");
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&out_file)?;

        let header = if bam_provided {
            vec![
                "Chromosome",
                "Gene",
                "Positions",
                "Reference Bases",
                "Base Changes",
                "AA Changes",
                "SNP AA Changes",
                "Local AA Changes",
                "Local SNP AA Changes",
                "Variant Type",
                "Change Type",
                "Reference Codon",
                "SNP Codon",
                "MNV Codon",
                "SNP Reads",
                "SNP Forward Reads",
                "SNP Reverse Reads",
                "MNV Reads",
                "MNV Forward Reads",
                "MNV Reverse Reads",
                "Total Reads",
                "SNP Frequencies",
                "MNV Frequencies",
            ]
        } else {
            vec![
                "Chromosome",
                "Gene",
                "Positions",
                "Reference Bases",
                "Base Changes",
                "AA Changes",
                "SNP AA Changes",
                "Local AA Changes",
                "Local SNP AA Changes",
                "Variant Type",
                "Change Type",
                "Reference Codon",
                "SNP Codon",
                "MNV Codon",
            ]
        };
        writer.write_record(&header)?;

        Ok(Self {
            writer,
            bam_provided,
        })
    }

    pub fn write_variants(&mut self, variants: &[VariantInfo]) -> AppResult<()> {
        for variant in variants {
            if self.bam_provided {
                let row = build_tsv_row_with_reads(variant)?;
                self.writer.write_record(&row)?;
            } else {
                let row = build_tsv_row_without_reads(variant)?;
                self.writer.write_record(&row)?;
            }
        }
        self.writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::build_tsv_row_with_reads;
    use crate::variants::{ChangeType, VariantInfo, VariantType};

    fn variant_with_reads() -> VariantInfo {
        VariantInfo {
            chrom: "chrX".to_string(),
            gene: "geneX".to_string(),
            positions: vec![10, 11],
            ref_bases: vec!["A".to_string(), "C".to_string()],
            base_changes: vec!["T".to_string(), "G".to_string()],
            aa_changes: vec!["Ala1Val".to_string()],
            snp_aa_changes: vec!["Ala1Ser".to_string(), "Ala1Thr".to_string()],
            aa_changes_local: vec!["-".to_string()],
            snp_aa_changes_local: vec!["-".to_string()],
            variant_type: VariantType::SnpMnv,
            change_type: ChangeType::NonSynonymous,
            snp_reads: Some(vec![1, 1]),
            snp_forward_reads: Some(vec![1, 0]),
            snp_reverse_reads: Some(vec![0, 1]),
            mnv_reads: Some(1),
            mnv_forward_reads: Some(1),
            mnv_reverse_reads: Some(0),
            mnv_total_reads: Some(4),
            total_reads: Some(vec![10, 2]),
            total_forward_reads: Some(vec![7, 1]),
            total_reverse_reads: Some(vec![3, 1]),
            mnv_total_forward_reads: Some(3),
            mnv_total_reverse_reads: Some(1),
            ref_codon: Some("ACC".to_string()),
            snp_codon: Some("TCC, AGC".to_string()),
            mnv_codon: Some("TGC".to_string()),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }
    }

    #[test]
    fn test_build_tsv_row_uses_mnv_total_depth_for_mnv_metrics() {
        let row = build_tsv_row_with_reads(&variant_with_reads()).expect("row should build");
        // Two new columns (Local AA Changes / Local SNP AA Changes) shifted
        // every later column by +2.
        assert_eq!(row[20], "4");
        assert_eq!(row[22], "0.2500");
    }

    #[test]
    fn test_build_tsv_row_emits_local_aa_columns() {
        let mut variant = variant_with_reads();
        variant.aa_changes = vec!["Phe228Ile".to_string()];
        variant.aa_changes_local = vec!["Val26His".to_string()];
        variant.snp_aa_changes = vec!["Phe228Ile".to_string()];
        variant.snp_aa_changes_local = vec!["Val26His".to_string()];
        let row = build_tsv_row_with_reads(&variant).expect("row should build");
        assert_eq!(row[5], "Phe228Ile", "AA Changes should be protein-wide");
        assert_eq!(row[6], "Phe228Ile", "SNP AA Changes should be protein-wide");
        assert_eq!(row[7], "Val26His", "Local AA Changes should be exon-local");
        assert_eq!(row[8], "Val26His", "Local SNP AA Changes should be exon-local");
    }
}
