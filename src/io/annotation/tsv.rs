//! Loading genes from the simple 4 to 7-column TSV annotation format.
//!
//! The seventh column, when present, names the contig a gene lives on. A row
//! without it applies to whichever contig is being annotated, which is the only
//! reading a single-contig file ever needed; a multi-contig run refuses such
//! rows rather than lay one contig's genes over another's variants.

use crate::error::AppResult;
use crate::io::VcfPosition;
use crate::variants::Gene;
use std::fs::File;
use std::io::{BufRead, BufReader};

use super::gff_parse::{parse_gff_phase, parse_interval, parse_strand};
use super::has_snp_in_interval;

/// The contig a TSV row names in its seventh column, if it names one. `.` and
/// an empty cell name none, as they do in the phase column.
fn row_contig<'a>(fields: &[&'a str]) -> Option<&'a str> {
    fields
        .get(6)
        .map(|value| value.trim())
        .filter(|value| !value.is_empty() && *value != ".")
}

/// How many annotation rows name no contig. A multi-contig run needs this to
/// be zero, because such a row cannot be placed on one contig rather than
/// another.
pub fn tsv_rows_without_contig(genes_file: &str) -> AppResult<usize> {
    let reader = BufReader::new(File::open(genes_file)?);
    let mut without = 0usize;
    for line in reader.lines() {
        let entry = line?;
        let trimmed = entry.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if row_contig(&fields).is_none() {
            without += 1;
        }
    }
    Ok(without)
}

/// Genes of the TSV annotation on `expected_chrom`, keeping only those with a
/// variant from `snp_list` inside, or every one when `snp_list` is `None`.
pub(super) fn load_genes_from_tsv(
    genes_file: &str,
    snp_list: Option<&[VcfPosition]>,
    expected_chrom: Option<&str>,
) -> AppResult<Vec<Gene>> {
    log::info!("Loading TSV gene file: {genes_file}");
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    let mut genes: Vec<Gene> = Vec::new();
    let mut parsed_entries = 0usize;
    let mut genes_with_snps = 0usize;
    let mut genes_without_snps = 0usize;
    let mut on_other_contigs = 0usize;

    for (line_idx, line) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let entry =
            line.map_err(|e| format!("Failed to read line {line_number} in genes file: {e}"))?;
        let trimmed = entry.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        parsed_entries += 1;

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 4 {
            return Err(format!(
                "Invalid genes entry at line {}: expected 4 tab-separated fields, got {}. Content: '{}'",
                line_number,
                fields.len(),
                entry
            ).into());
        }

        let (start, end) = parse_interval(fields[1], fields[2], line_number)?;
        let strand = parse_strand(fields[3], line_number)?;
        // Optional 5th column: phase (0|1|2|.). Defaults to 0 (prokaryote-style)
        // when omitted, preserving the historical 4-column TSV format.
        let phase = if fields.len() >= 5 {
            parse_gff_phase(fields[4], line_number)?
        } else {
            0
        };
        // Optional 6th column: biotype. Without it every feature is assumed to be
        // translated, which turns a non-coding RNA gene into an invented protein.
        let biotype = if fields.len() >= 6 {
            crate::variants::Biotype::parse(fields[5]).map_err(|err| {
                format!("Invalid genes entry at line {line_number}: {err}. Content: '{entry}'")
            })?
        } else {
            crate::variants::Biotype::default()
        };

        if let (Some(contig), Some(expected)) = (row_contig(&fields), expected_chrom) {
            if contig != expected {
                on_other_contigs += 1;
                continue;
            }
        }

        if snp_list.is_none_or(|snps| has_snp_in_interval(snps, start, end)) {
            genes_with_snps += 1;
            genes.push(crate::variants::Gene {
                name: fields[0].to_string(),
                start,
                end,
                strand,
                phase,
                protein_offset: 0,
                transcript_id: None,
                cds_segments: Vec::new(),
                biotype,
            });
        } else {
            genes_without_snps += 1;
        }
    }

    log::info!(
        "TSV gene entries parsed: {parsed_entries} | mapped to SNPs: {genes_with_snps} | without SNPs: {genes_without_snps} | on other contigs: {on_other_contigs}"
    );

    Ok(genes)
}
