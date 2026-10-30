//! Gene annotation parsing: GFF/GFF3 and TSV formats, format detection, and
//! gene-to-SNP interval filtering.
//!
//! Parsing primitives, the TSV loader, the GFF CDS-model builder, and the GFF
//! loaders live in submodules; this file keeps the gene/variant overlap
//! predicates and the top-level `load_genes` / `preload_gff_genes` entry points.

use super::variant_index::VariantIndex;
use super::vcf::VcfPosition;
use crate::error::AppResult;
use crate::variants::{AlleleComponentKind, Gene};
use std::collections::HashMap;

mod cds_model;
mod format;
mod gff_load;
mod gff_parse;
mod tsv;

#[cfg(test)]
mod tests;

pub use format::{detect_annotation_format, AnnotationFormat};

#[cfg(test)]
pub(crate) use cds_model::parse_gff_gene_records_from_reader;
pub(crate) use cds_model::{
    assign_cds_protein_offsets, build_transcript_cds_records, parse_gff_gene_records,
};
pub(crate) use gff_load::{gff_has_cds_features, load_genes_from_gff};

// Re-exports that exist only so the unit tests keep a stable `annotation::*`
// path to internal helpers; not referenced by non-test code.
#[cfg(test)]
pub(crate) use cds_model::{count_multi_exon_cds_transcripts, GffGeneRecord};
#[cfg(test)]
pub(crate) use gff_parse::{
    decode_percent_encoded, gene_name_from_gff, parse_gff_attributes, split_attribute_fields,
};

use cds_model::warn_if_multi_exon_cds_detected;
use gff_load::gff_has_non_zero_phase_cds;
use tsv::load_genes_from_tsv;
pub(crate) use tsv::{load_tsv_genes_with, read_tsv_gene_rows, TsvGeneRow};
pub use tsv::{tsv_named_contigs, tsv_rows_without_contig};

/// Return true when any variant event overlaps the interval. Insertions only
/// overlap when their inserted bases fall between bases inside the interval;
/// deletions use the deleted reference span.
pub(crate) fn has_snp_in_interval(variants: &VariantIndex, start: usize, end: usize) -> bool {
    variants
        .near(start, end)
        .iter()
        .any(|variant| variant.overlaps_interval(start, end))
}

pub fn gene_overlaps_variant(gene: &Gene, variant: &VcfPosition) -> bool {
    if gene.cds_segments.is_empty() {
        return variant.overlaps_interval(gene.start, gene.end);
    }
    if gene
        .cds_segments
        .iter()
        .any(|segment| variant.overlaps_interval(segment.start, segment.end))
    {
        return true;
    }
    // An insertion anchored at an internal exon-exon junction is coding even
    // though overlaps_interval (half-open) excludes the exon's last base; keep it
    // out of the intergenic set so it is annotated on the spliced CDS instead.
    variant.event().components.iter().any(|component| {
        matches!(component.kind, AlleleComponentKind::Insertion)
            && gene.insertion_at_internal_junction(component.position)
    })
}

fn gene_has_variant(gene: &Gene, variants: &VariantIndex) -> bool {
    variants.near_gene(gene).iter().any(|variant| {
        gene_overlaps_variant(gene, variant)
            // Keep genes whose only nearby variant is intronic, whether at a
            // splice site or deeper in the intron, so it can be annotated
            // against its gene instead of being dropped to intergenic.
            // `gene_overlaps_variant` (used for intergenic coverage) is
            // deliberately left CDS-only.
            // Asked of the bases the record changes, not of the base it starts
            // on: a padded record such as `40 GCTGCTG>GCTGCTA` begins on an
            // exonic base it leaves alone and destroys a splice donor at 46, and
            // testing POS dropped the gene here, so the fallback downstream had
            // no gene left to name and called the variant intergenic.
            || variant.changed_positions().iter().any(|&position| {
                crate::variants::splice::splice_consequence_for_position(gene, position).is_some()
                    || crate::variants::splice::is_intronic_position(gene, position)
            })
    })
}

pub fn load_genes(
    genes_file: &str,
    snp_list: &[VcfPosition],
    expected_chrom: Option<&str>,
    feature_types: &[String],
) -> AppResult<Vec<Gene>> {
    match detect_annotation_format(genes_file)? {
        AnnotationFormat::Tsv => load_genes_from_tsv(genes_file, Some(snp_list), expected_chrom),
        AnnotationFormat::Gff => {
            load_genes_from_gff(genes_file, snp_list, expected_chrom, feature_types)
        }
    }
}

/// Every gene of a TSV annotation on `chrom`, whether or not a variant falls in
/// it: a catalogue entry upstream of a gene is placed from the gene's start,
/// and the gene itself holds no variant. A GFF is pre-loaded whole already.
pub fn load_all_tsv_genes(genes_file: &str, chrom: &str) -> AppResult<Vec<Gene>> {
    load_genes_from_tsv(genes_file, None, Some(chrom))
}

/// Pre-load all gene entries from a GFF/GFF3 file, grouped by contig.
/// No SNP filtering is applied; use `filter_genes_with_snps` afterwards.
pub fn preload_gff_genes(
    genes_file: &str,
    feature_types: &[String],
) -> AppResult<HashMap<String, Vec<Gene>>> {
    log::info!("Pre-loading GFF/GFF3 annotation file: {genes_file}");
    warn_if_cds_phase_ignored(genes_file, feature_types)?;
    let mut records = parse_gff_gene_records(genes_file, feature_types)?;
    assign_cds_protein_offsets(&mut records);
    warn_if_multi_exon_cds_detected(&records);
    let records = build_transcript_cds_records(records);
    let total_entries = records.len();

    let mut genes_by_contig: HashMap<String, Vec<Gene>> = HashMap::new();
    for rec in records {
        genes_by_contig
            .entry(rec.contig)
            .or_default()
            .push(rec.gene);
    }

    let contig_count = genes_by_contig.len();
    log::info!("GFF/GFF3 pre-loaded: {total_entries} gene entries across {contig_count} contigs");
    Ok(genes_by_contig)
}

/// Emit a single warning line if the GFF contains CDS rows with non-zero
/// phase but the user did not select `CDS` in `--gff-features`. In that
/// configuration the phase-aware codon grouping used by current releases is
/// effectively disabled and codon boundaries may be wrong for eukaryotic
/// annotations.
fn warn_if_cds_phase_ignored(genes_file: &str, feature_types: &[String]) -> AppResult<()> {
    let selects_cds = feature_types.iter().any(|f| f == "CDS");
    if selects_cds {
        return Ok(());
    }
    if gff_has_non_zero_phase_cds(genes_file)? {
        log::warn!(
            "GFF/GFF3 file '{genes_file}' contains CDS rows with non-zero phase, \
             but --gff-features does not include 'CDS' (currently: {}). \
             Codon grouping will ignore the phase information and may be incorrect \
             for eukaryotic annotations. Re-run with '--gff-features CDS' to use \
             phase-aware codon boundaries.",
            feature_types.join(",")
        );
    }
    Ok(())
}

/// Filter genes that overlap with at least one SNP position. The list need
/// not be sorted, but the search is quickest when it is.
pub fn filter_genes_with_snps(genes: &[Gene], snp_list: &[VcfPosition]) -> Vec<Gene> {
    filter_genes_near(genes, &VariantIndex::new(snp_list))
}

/// [`filter_genes_with_snps`] for a contig whose variants are indexed already.
pub(crate) fn filter_genes_near(genes: &[Gene], variants: &VariantIndex) -> Vec<Gene> {
    genes
        .iter()
        .filter(|gene| gene_has_variant(gene, variants))
        .cloned()
        .collect()
}
