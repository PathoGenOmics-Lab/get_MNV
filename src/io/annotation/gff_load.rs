//! GFF/GFF3 CDS probes and the GFF gene loader.

use crate::error::AppResult;
use crate::io::VcfPosition;
use crate::variants::Gene;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

use super::cds_model::{
    assign_cds_protein_offsets, build_transcript_cds_records, parse_gff_gene_records,
    warn_if_multi_exon_cds_detected,
};
use super::{gene_has_variant, warn_if_cds_phase_ignored};

/// Scan a GFF file once and report whether it contains any CDS row with a
/// non-zero phase. Used to warn users that they should pass
/// `--gff-features CDS` when working with eukaryotic annotations.
pub(crate) fn gff_has_non_zero_phase_cds(genes_file: &str) -> AppResult<bool> {
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let entry = line?;
        let trimmed = entry.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() != 9 {
            continue;
        }
        if fields[2] != "CDS" {
            continue;
        }
        match fields[7].trim() {
            "1" | "2" => return Ok(true),
            _ => {}
        }
    }
    Ok(false)
}

/// Whether the GFF/GFF3 file contains any `CDS` feature row. Used to auto-select
/// the phase/splice-aware CDS model when the user did not pass `--gff-features`,
/// so a eukaryotic annotation is not silently analysed over whole-gene spans.
pub(crate) fn gff_has_cds_features(genes_file: &str) -> AppResult<bool> {
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let entry = line?;
        let trimmed = entry.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() == 9 && fields[2] == "CDS" {
            return Ok(true);
        }
    }
    Ok(false)
}

pub(crate) fn load_genes_from_gff(
    genes_file: &str,
    snp_list: &[VcfPosition],
    expected_chrom: Option<&str>,
    feature_types: &[String],
) -> AppResult<Vec<Gene>> {
    log::info!("Loading GFF/GFF3 annotation file: {genes_file}");
    warn_if_cds_phase_ignored(genes_file, feature_types)?;
    let mut records = parse_gff_gene_records(genes_file, feature_types)?;
    assign_cds_protein_offsets(&mut records);
    warn_if_multi_exon_cds_detected(&records);
    let records = build_transcript_cds_records(records);

    let mut genes: Vec<Gene> = Vec::new();
    let mut parsed_entries = 0usize;
    let mut genes_with_snps = 0usize;
    let mut genes_without_snps = 0usize;
    let mut genes_other_contigs = 0usize;
    let mut other_contigs: HashSet<String> = HashSet::new();

    for rec in records {
        if let Some(chrom) = expected_chrom {
            if rec.contig != chrom {
                genes_other_contigs += 1;
                other_contigs.insert(rec.contig);
                continue;
            }
        }

        parsed_entries += 1;

        if gene_has_variant(&rec.gene, snp_list) {
            genes_with_snps += 1;
            genes.push(rec.gene);
        } else {
            genes_without_snps += 1;
        }
    }

    log::info!(
        "GFF/GFF3 gene entries parsed: {} | mapped to SNPs: {} | without SNPs: {} | skipped other contigs: {}{}",
        parsed_entries,
        genes_with_snps,
        genes_without_snps,
        genes_other_contigs,
        if other_contigs.is_empty() {
            String::new()
        } else {
            let mut values = other_contigs.into_iter().collect::<Vec<_>>();
            values.sort();
            format!(" ({})", values.join(", "))
        }
    );

    Ok(genes)
}
