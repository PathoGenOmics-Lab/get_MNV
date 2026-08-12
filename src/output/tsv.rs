//! TSV output writer: generates tab-delimited variant tables with optional
//! BAM-derived read support columns.

use crate::error::AppResult;
use crate::variants::{VariantInfo, VariantType};
use csv::WriterBuilder;
use std::fs::File;

use super::common::{
    allele_frequency, format_freq, get_mnv_depth_from_variant, snp_bam_vectors,
    validate_variant_shape, SnpBamVectors,
};

fn is_intergenic(variant: &VariantInfo) -> bool {
    variant.gene == "intergenic"
}

/// Trailing annotation cells appended to every TSV row: SO consequence term,
/// impact, Grantham distance (with category), MNV-vs-SNV consequence shift,
/// COSMIC-style DBS class, BAM-derived MNV phasing support, the 50-nt-rule NMD
/// prediction, and the HGVS genomic / coding descriptors.
///
/// With a BAM the phasing support is followed by the number of reads it was
/// computed from, so a ratio of 1.0 from two reads is not mistaken for the same
/// ratio from two hundred.
fn annotation_cells(variant: &VariantInfo, bam_provided: bool) -> Vec<String> {
    let grantham = match variant.annotations.grantham {
        Some(d) => format!("{d} ({})", crate::utils::grantham_category(d)),
        None => "-".to_string(),
    };
    let (so_term, impact) = super::common::so_consequence(variant);
    let dbs = variant
        .annotations
        .dbs_class
        .clone()
        .unwrap_or_else(|| "-".to_string());
    let phasing = super::common::mnv_phasing_support(variant)
        .map(format_freq)
        .unwrap_or_else(|| "-".to_string());
    let nmd = variant
        .annotations
        .nmd
        .map(|prediction| prediction.to_string())
        .unwrap_or_else(|| "-".to_string());
    let hgvs_g = crate::variants::hgvs::genomic(variant).unwrap_or_else(|| "-".to_string());
    let hgvs_c = variant
        .annotations
        .hgvs_c
        .clone()
        .unwrap_or_else(|| "-".to_string());
    let mut cells = vec![
        so_term,
        impact.to_string(),
        grantham,
        variant.annotations.consequence_shift.to_string(),
        dbs,
        super::common::declared_phase_field(variant),
        phasing,
    ];
    if bam_provided {
        let (d_prime, p_value) = super::common::codon_linkage_fields(variant);
        cells.push(d_prime);
        cells.push(p_value);
        cells.push(
            super::common::mnv_phasing_read_count(variant)
                .map(|reads| reads.to_string())
                .unwrap_or_else(|| "-".to_string()),
        );
        cells.push(super::common::frameshift_linkage_field(variant));
    }
    cells.extend([nmd, hgvs_g, hgvs_c]);
    cells
}

/// Render a "Local …" column. Falls back to the protein-wide column when the
/// local vector is empty (e.g. records produced before this field existed and
/// then deserialized via `#[serde(default)]`).
fn local_aa_or_fallback(local: &[String], protein: &[String]) -> String {
    let chosen = if local.is_empty() { protein } else { local };
    chosen.join(", ")
}

fn mnv_frequency(variant: &VariantInfo, total_reads: &[usize]) -> f64 {
    allele_frequency(
        variant.mnv_reads.unwrap_or(0),
        get_mnv_depth_from_variant(variant, total_reads),
    )
}

fn event_class(variant: &VariantInfo) -> String {
    variant
        .event_class
        .clone()
        .unwrap_or_else(|| variant.variant_type.to_string().to_ascii_lowercase())
}

fn event_components(variant: &VariantInfo) -> String {
    if variant.event_components.is_empty() {
        "-".to_string()
    } else {
        variant.event_components.join(" | ")
    }
}

fn event_read_columns(variant: &VariantInfo) -> Vec<String> {
    if variant.variant_type != VariantType::Indel {
        return vec![
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
        ];
    }

    let reads = variant.mnv_reads.unwrap_or(0);
    let forward = variant.mnv_forward_reads.unwrap_or(0);
    let reverse = variant.mnv_reverse_reads.unwrap_or(0);
    let depth = variant.mnv_total_reads.unwrap_or(0);
    vec![
        reads.to_string(),
        forward.to_string(),
        reverse.to_string(),
        depth.to_string(),
        format_freq(allele_frequency(reads, depth)),
    ]
}

#[derive(Clone, Copy)]
struct TsvFilterConfig {
    min_snp_reads: usize,
    min_snp_frequency: f64,
    min_mnv_reads: usize,
    min_mnv_frequency: f64,
    min_snp_strand_reads: usize,
    min_mnv_strand_reads: usize,
}

impl TsvFilterConfig {
    fn has_active_filters(self) -> bool {
        self.min_snp_reads > 0
            || self.min_snp_frequency > 0.0
            || self.min_mnv_reads > 0
            || self.min_mnv_frequency > 0.0
            || self.min_snp_strand_reads > 0
            || self.min_mnv_strand_reads > 0
    }

    fn has_active_mnv_filters(self) -> bool {
        self.min_mnv_reads > 0 || self.min_mnv_frequency > 0.0 || self.min_mnv_strand_reads > 0
    }
}

pub struct TsvWriterConfig<'a> {
    pub filename: &'a str,
    pub bam_provided: bool,
    pub min_snp_reads: usize,
    pub min_snp_frequency: f64,
    pub min_mnv_reads: usize,
    pub min_mnv_frequency: f64,
    pub min_snp_strand_reads: usize,
    pub min_mnv_strand_reads: usize,
}

fn snp_component_passes_filters(bam_vectors: &SnpBamVectors<'_>, filters: TsvFilterConfig) -> bool {
    bam_vectors
        .reads
        .iter()
        .zip(bam_vectors.total_reads.iter())
        .zip(bam_vectors.forward_reads.iter())
        .zip(bam_vectors.reverse_reads.iter())
        .all(|(((&support, &depth), &forward), &reverse)| {
            support >= filters.min_snp_reads
                && allele_frequency(support, depth) >= filters.min_snp_frequency
                && forward >= filters.min_snp_strand_reads
                && reverse >= filters.min_snp_strand_reads
        })
}

fn mnv_component_passes_filters(
    variant: &VariantInfo,
    total_reads: &[usize],
    filters: TsvFilterConfig,
) -> bool {
    variant.mnv_reads.unwrap_or(0) >= filters.min_mnv_reads
        && mnv_frequency(variant, total_reads) >= filters.min_mnv_frequency
        && variant.mnv_forward_reads.unwrap_or(0) >= filters.min_mnv_strand_reads
        && variant.mnv_reverse_reads.unwrap_or(0) >= filters.min_mnv_strand_reads
}

fn passes_filters(variant: &VariantInfo, filters: TsvFilterConfig) -> AppResult<bool> {
    if !filters.has_active_filters() {
        return Ok(true);
    }
    // Intergenic variants carry no recomputed read support, so they cannot
    // satisfy a read-based filter. Drop them when the filters relevant to their
    // variant type are active (indels/MNV use the MNV filters, SNPs the SNP
    // filters); keep them otherwise. Checked before the indel branch so an
    // intergenic indel is filtered consistently with an intergenic SNP.
    if is_intergenic(variant) {
        // Intergenic SNPs are read-counted at their position, so filter them by
        // their real support exactly like any SNP.
        if variant.variant_type == VariantType::Snp {
            let bam_vectors = snp_bam_vectors(variant)?;
            return Ok(snp_component_passes_filters(&bam_vectors, filters));
        }
        // Intergenic indels carry no recomputed read support, so they cannot
        // satisfy a read-based filter; drop them when the relevant (MNV) filters
        // are active, keep them otherwise.
        return Ok(!filters.has_active_mnv_filters());
    }
    if variant.variant_type == VariantType::Indel {
        if variant.mnv_reads.is_some() || variant.mnv_total_reads.is_some() {
            return Ok(mnv_component_passes_filters(
                variant,
                &[variant.mnv_total_reads.unwrap_or(0)],
                filters,
            ));
        }
        return Ok(true);
    }

    let bam_vectors = snp_bam_vectors(variant)?;
    let snp_passes = || snp_component_passes_filters(&bam_vectors, filters);
    let mnv_passes = || mnv_component_passes_filters(variant, bam_vectors.total_reads, filters);

    Ok(match variant.variant_type {
        VariantType::Snp => snp_passes(),
        VariantType::Mnv => mnv_passes(),
        VariantType::SnpMnv => snp_passes() || mnv_passes(),
        VariantType::Indel => true,
    })
}

fn build_tsv_row_with_reads(variant: &VariantInfo) -> AppResult<Vec<String>> {
    validate_variant_shape(variant)?;
    // Indels (genic or intergenic) have no read counts and use the placeholder
    // layout. Intergenic SNPs are read-counted, so they fall through to the
    // normal SNP path below and report their read support like any SNP.
    if variant.variant_type == VariantType::Indel {
        let pos_str = variant
            .positions
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        let ref_base_str = variant.ref_bases.join(", ");
        let base_str = variant.base_changes.join(", ");
        let mut row = vec![
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
        ];
        row.push(event_class(variant));
        row.push(event_components(variant));
        row.extend(event_read_columns(variant));
        return Ok(row);
    }

    let bam_vectors = snp_bam_vectors(variant)?;
    let snp_counts = bam_vectors.reads;
    let snp_forward_counts = bam_vectors.forward_reads;
    let snp_reverse_counts = bam_vectors.reverse_reads;
    let total_reads = bam_vectors.total_reads;

    let snp_freq: Vec<String> = snp_counts
        .iter()
        .zip(total_reads.iter())
        .map(|(&s, &t)| format_freq(allele_frequency(s, t)))
        .collect();

    let dp_mean: usize = total_reads.iter().sum::<usize>() / total_reads.len().max(1);
    let mnv_depth = get_mnv_depth_from_variant(variant, total_reads);

    let mnv_freq_str = format_freq(allele_frequency(variant.mnv_reads.unwrap_or(0), mnv_depth));

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

    let mut row = vec![
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
    ];
    row.push(event_class(variant));
    row.push(event_components(variant));
    row.extend(event_read_columns(variant));
    Ok(row)
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
        let mut row = vec![
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
        ];
        row.push(event_class(variant));
        row.push(event_components(variant));
        return Ok(row);
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
    let mut row = vec![
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
    ];
    row.push(event_class(variant));
    row.push(event_components(variant));
    Ok(row)
}

pub struct TsvWriter {
    writer: csv::Writer<File>,
    bam_provided: bool,
    filters: TsvFilterConfig,
}

impl TsvWriter {
    pub fn new(cfg: TsvWriterConfig<'_>) -> AppResult<Self> {
        let filename = cfg.filename;
        let bam_provided = cfg.bam_provided;
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
                "Event Class",
                "Event Components",
                "Event Reads",
                "Event Forward Reads",
                "Event Reverse Reads",
                "Event Depth",
                "Event Frequency",
                "SO Term",
                "Impact",
                "Grantham",
                "MNV Consequence Shift",
                "DBS Class",
                "Declared Phase",
                "MNV Phasing Support",
                "Codon LD",
                "Codon LD p",
                "MNV Phasing Reads",
                "Frameshift Phasing",
                "NMD Prediction",
                "HGVS g.",
                "HGVS c.",
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
                "Event Class",
                "Event Components",
                "SO Term",
                "Impact",
                "Grantham",
                "MNV Consequence Shift",
                "DBS Class",
                "Declared Phase",
                "MNV Phasing Support",
                "NMD Prediction",
                "HGVS g.",
                "HGVS c.",
            ]
        };
        writer.write_record(&header)?;

        Ok(Self {
            writer,
            bam_provided,
            filters: TsvFilterConfig {
                min_snp_reads: cfg.min_snp_reads,
                min_snp_frequency: cfg.min_snp_frequency,
                min_mnv_reads: cfg.min_mnv_reads,
                min_mnv_frequency: cfg.min_mnv_frequency,
                min_snp_strand_reads: cfg.min_snp_strand_reads,
                min_mnv_strand_reads: cfg.min_mnv_strand_reads,
            },
        })
    }

    pub fn write_variants(&mut self, variants: &[VariantInfo]) -> AppResult<()> {
        for variant in variants {
            if self.bam_provided {
                if !passes_filters(variant, self.filters)? {
                    continue;
                }
                let mut row = build_tsv_row_with_reads(variant)?;
                row.extend(annotation_cells(variant, true));
                self.writer.write_record(&row)?;
            } else {
                let mut row = build_tsv_row_without_reads(variant)?;
                row.extend(annotation_cells(variant, false));
                self.writer.write_record(&row)?;
            }
        }
        self.writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{build_tsv_row_with_reads, passes_filters, TsvFilterConfig};
    use crate::variants::{ChangeType, VariantInfo, VariantType};

    fn filters(
        min_snp_reads: usize,
        min_snp_frequency: f64,
        min_mnv_reads: usize,
        min_mnv_frequency: f64,
        min_snp_strand_reads: usize,
        min_mnv_strand_reads: usize,
    ) -> TsvFilterConfig {
        TsvFilterConfig {
            min_snp_reads,
            min_snp_frequency,
            min_mnv_reads,
            min_mnv_frequency,
            min_snp_strand_reads,
            min_mnv_strand_reads,
        }
    }

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
            mnv_phasing_reads: None,
            ref_codon: Some("ACC".to_string()),
            snp_codon: Some("TCC, AGC".to_string()),
            mnv_codon: Some("TGC".to_string()),
            original_dp: None,
            original_freq: None,
            original_info: None,
            event_class: Some("mnv".to_string()),
            event_components: vec!["SNV:10:A>T".to_string(), "SNV:11:C>G".to_string()],
            annotations: crate::variants::VariantAnnotations::default(),
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
        assert_eq!(
            row[8], "Val26His",
            "Local SNP AA Changes should be exon-local"
        );
    }

    fn intergenic_indel() -> VariantInfo {
        // As produced by build_intergenic_variant: no recomputed read support.
        let mut v = variant_with_reads();
        v.gene = "intergenic".to_string();
        v.variant_type = VariantType::Indel;
        v.positions = vec![10];
        v.ref_bases = vec!["AC".to_string()];
        v.base_changes = vec!["A".to_string()];
        v.snp_reads = None;
        v.snp_forward_reads = None;
        v.snp_reverse_reads = None;
        v.mnv_reads = None;
        v.mnv_forward_reads = None;
        v.mnv_reverse_reads = None;
        v.mnv_total_reads = None;
        v.total_reads = None;
        v.total_forward_reads = None;
        v.total_reverse_reads = None;
        v.mnv_total_forward_reads = None;
        v.mnv_total_reverse_reads = None;
        v.event_class = Some("deletion".to_string());
        v
    }

    #[test]
    fn test_intergenic_indel_respects_mnv_filters() {
        // Regression: an intergenic indel carries no recomputed read support, so
        // it must be dropped when the relevant (MNV) filters are active, the same
        // way an intergenic SNP is dropped under active SNP filters - instead of
        // being kept unconditionally by the indel branch.
        let variant = intergenic_indel();
        // No filters -> kept.
        assert!(passes_filters(&variant, filters(0, 0.0, 0, 0.0, 0, 0)).unwrap());
        // Active MNV frequency / read / strand filters -> dropped.
        assert!(!passes_filters(&variant, filters(0, 0.0, 0, 0.20, 0, 0)).unwrap());
        assert!(!passes_filters(&variant, filters(0, 0.0, 2, 0.0, 0, 0)).unwrap());
        assert!(!passes_filters(&variant, filters(0, 0.0, 0, 0.0, 0, 1)).unwrap());
        // A SNP-only filter does not affect an intergenic indel.
        assert!(passes_filters(&variant, filters(5, 0.0, 0, 0.0, 0, 0)).unwrap());
    }

    #[test]
    fn test_tsv_frequency_filters_use_bam_derived_values() {
        let variant = variant_with_reads();

        assert!(passes_filters(&variant, filters(0, 0.09, 0, 0.20, 0, 0)).unwrap());
        assert!(passes_filters(&variant, filters(0, 0.20, 0, 0.20, 0, 0)).unwrap());
        assert!(!passes_filters(&variant, filters(0, 0.20, 0, 0.30, 0, 0)).unwrap());
    }

    #[test]
    fn test_tsv_frequency_filters_keep_snp_mnv_when_mnv_component_passes() {
        let variant = variant_with_reads();

        assert!(
            passes_filters(&variant, filters(0, 0.20, 0, 0.20, 0, 0)).unwrap(),
            "a strong MNV haplotype should keep a mixed SNP/MNV row even when SNP frequencies fail"
        );
        assert!(
            passes_filters(&variant, filters(0, 0.20, 0, 0.0, 0, 0)).unwrap(),
            "an SNP-only frequency filter should not remove the MNV component of a mixed row"
        );
        assert!(
            passes_filters(&variant, filters(0, 0.09, 0, 0.30, 0, 0)).unwrap(),
            "a passing SNP component should keep a mixed row even when the MNV frequency fails"
        );
        assert!(
            !passes_filters(&variant, filters(0, 0.20, 0, 0.30, 0, 0)).unwrap(),
            "a mixed row should be removed only when every relevant component fails its own active filter"
        );
    }

    #[test]
    fn test_tsv_strand_filters_keep_snp_mnv_when_mnv_component_passes() {
        let mut variant = variant_with_reads();
        variant.mnv_reads = Some(4);
        variant.mnv_forward_reads = Some(2);
        variant.mnv_reverse_reads = Some(2);

        assert!(
            passes_filters(&variant, filters(0, 0.0, 0, 0.0, 1, 2)).unwrap(),
            "a mixed SNP/MNV row should be kept when the MNV strand support passes"
        );
        assert!(
            !passes_filters(&variant, filters(0, 0.0, 0, 0.0, 1, 3)).unwrap(),
            "a mixed SNP/MNV row should be removed when neither SNP nor MNV strand support passes"
        );
    }

    fn intergenic_snp(support: usize, forward: usize, reverse: usize, depth: usize) -> VariantInfo {
        VariantInfo {
            chrom: "chr1".to_string(),
            gene: "intergenic".to_string(),
            positions: vec![35],
            ref_bases: vec!["G".to_string()],
            base_changes: vec!["A".to_string()],
            aa_changes: vec!["-".to_string()],
            snp_aa_changes: vec!["-".to_string()],
            aa_changes_local: vec!["-".to_string()],
            snp_aa_changes_local: vec!["-".to_string()],
            variant_type: VariantType::Snp,
            change_type: ChangeType::Unknown,
            snp_reads: Some(vec![support]),
            snp_forward_reads: Some(vec![forward]),
            snp_reverse_reads: Some(vec![reverse]),
            mnv_reads: Some(support),
            mnv_forward_reads: Some(forward),
            mnv_reverse_reads: Some(reverse),
            mnv_total_reads: Some(depth),
            total_reads: Some(vec![depth]),
            total_forward_reads: Some(vec![forward]),
            total_reverse_reads: Some(vec![reverse]),
            mnv_total_forward_reads: Some(forward),
            mnv_total_reverse_reads: Some(reverse),
            mnv_phasing_reads: None,
            ref_codon: None,
            snp_codon: None,
            mnv_codon: None,
            original_dp: None,
            original_freq: None,
            original_info: None,
            event_class: None,
            event_components: Vec::new(),
            annotations: crate::variants::VariantAnnotations::default(),
        }
    }

    #[test]
    fn test_intergenic_snp_filtered_by_real_read_support() {
        // Intergenic SNPs are read-counted, so thresholds apply to them like any
        // SNP: kept when their real support passes, dropped when it does not.
        let variant = intergenic_snp(10, 5, 5, 10);
        assert!(passes_filters(&variant, filters(1, 0.0, 0, 0.0, 0, 0)).unwrap());
        assert!(passes_filters(&variant, filters(10, 0.5, 0, 0.0, 5, 0)).unwrap());
        assert!(!passes_filters(&variant, filters(11, 0.0, 0, 0.0, 0, 0)).unwrap());
        assert!(!passes_filters(&variant, filters(0, 0.0, 0, 0.0, 6, 0)).unwrap());

        let weak = intergenic_snp(1, 1, 0, 1);
        assert!(!passes_filters(&weak, filters(5, 0.0, 0, 0.0, 0, 0)).unwrap());
    }
}
