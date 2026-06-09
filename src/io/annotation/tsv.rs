//! Loading genes from the simple 4/5-column TSV annotation format.

use crate::error::AppResult;
use crate::io::VcfPosition;
use crate::variants::Gene;
use std::fs::File;
use std::io::{BufRead, BufReader};

use super::gff_parse::{parse_gff_phase, parse_interval, parse_strand};
use super::has_snp_in_interval;

pub(super) fn load_genes_from_tsv(
    genes_file: &str,
    snp_list: &[VcfPosition],
) -> AppResult<Vec<Gene>> {
    log::info!("Loading TSV gene file: {genes_file}");
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    let mut genes: Vec<Gene> = Vec::new();
    let mut parsed_entries = 0usize;
    let mut genes_with_snps = 0usize;
    let mut genes_without_snps = 0usize;

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

        if has_snp_in_interval(snp_list, start, end) {
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
            });
        } else {
            genes_without_snps += 1;
        }
    }

    log::info!(
        "TSV gene entries parsed: {parsed_entries} | mapped to SNPs: {genes_with_snps} | without SNPs: {genes_without_snps}"
    );

    Ok(genes)
}
