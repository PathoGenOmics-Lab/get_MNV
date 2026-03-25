//! Gene annotation parsing: GFF/GFF3 and TSV formats, format detection, and
//! gene-to-SNP interval filtering.

use super::vcf::VcfPosition;
use crate::error::AppResult;
use crate::variants::Gene;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnnotationFormat {
    Tsv,
    Gff,
}

pub fn detect_annotation_format(genes_file: &str) -> AppResult<AnnotationFormat> {
    if let Some(ext) = Path::new(genes_file).extension().and_then(|e| e.to_str()) {
        if ext.eq_ignore_ascii_case("gff") || ext.eq_ignore_ascii_case("gff3") {
            return Ok(AnnotationFormat::Gff);
        }
        if ext.eq_ignore_ascii_case("txt") || ext.eq_ignore_ascii_case("tsv") {
            return Ok(AnnotationFormat::Tsv);
        }
    }

    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let entry = line?;
        let trimmed = entry.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() >= 9 {
            return Ok(AnnotationFormat::Gff);
        }
        return Ok(AnnotationFormat::Tsv);
    }

    Ok(AnnotationFormat::Tsv)
}

pub(crate) fn decode_percent_encoded(value: &str) -> String {
    let mut buf: Vec<u8> = Vec::with_capacity(value.len());
    let mut bytes = value.as_bytes().iter().copied();
    while let Some(current) = bytes.next() {
        if current == b'%' {
            let hi = bytes.next();
            let lo = bytes.next();
            if let (Some(hi), Some(lo)) = (hi, lo) {
                let hex_pair = [hi, lo];
                if let Ok(hex) = std::str::from_utf8(&hex_pair) {
                    if let Ok(decoded) = u8::from_str_radix(hex, 16) {
                        buf.push(decoded);
                        continue;
                    }
                }
                buf.extend_from_slice(&[b'%', hi, lo]);
                continue;
            }
            buf.push(b'%');
            if let Some(hi) = hi {
                buf.push(hi);
            }
            if let Some(lo) = lo {
                buf.push(lo);
            }
            continue;
        }
        buf.push(current);
    }
    String::from_utf8(buf).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned())
}

pub(crate) fn split_attribute_fields(attributes: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut start = 0usize;
    let mut in_quotes = false;
    let mut escaped = false;
    for (idx, ch) in attributes.char_indices() {
        match ch {
            '"' if !escaped => in_quotes = !in_quotes,
            ';' if !in_quotes => {
                let field = attributes[start..idx].trim();
                if !field.is_empty() {
                    fields.push(field);
                }
                start = idx + 1;
            }
            _ => {}
        }
        escaped = ch == '\\' && !escaped;
    }
    let trailing = attributes[start..].trim();
    if !trailing.is_empty() {
        fields.push(trailing);
    }
    fields
}

pub(crate) fn parse_gff_attributes(attributes: &str) -> HashMap<String, String> {
    let mut parsed = HashMap::new();
    for field in split_attribute_fields(attributes) {
        let (key_raw, value_raw) = if let Some((key, value)) = field.split_once('=') {
            (key.trim(), value.trim())
        } else if let Some((key, value)) = field.split_once(' ') {
            (key.trim(), value.trim().trim_matches('"'))
        } else {
            continue;
        };
        if key_raw.is_empty() || value_raw.is_empty() {
            continue;
        }
        parsed.insert(
            key_raw.to_string(),
            decode_percent_encoded(value_raw.trim_matches('"').trim()),
        );
    }
    parsed
}

/// Requires `snp_list` to be sorted by position (guaranteed by `load_vcf_positions_by_contig`).
pub(crate) fn has_snp_in_interval(snp_list: &[VcfPosition], start: usize, end: usize) -> bool {
    let idx = snp_list.partition_point(|snp| snp.position < start);
    idx < snp_list.len() && snp_list[idx].position <= end
}

fn parse_strand(raw: &str, line_number: usize) -> AppResult<crate::variants::Strand> {
    raw.parse::<crate::variants::Strand>().map_err(|_| {
        format!("Invalid strand at line {line_number} ('{raw}'). Expected '+' or '-'").into()
    })
}

fn parse_interval(start_raw: &str, end_raw: &str, line_number: usize) -> AppResult<(usize, usize)> {
    let start = start_raw.parse::<usize>().map_err(|e| {
        format!("Invalid start coordinate at line {line_number} ('{start_raw}'): {e}")
    })?;
    let end = end_raw
        .parse::<usize>()
        .map_err(|e| format!("Invalid end coordinate at line {line_number} ('{end_raw}'): {e}"))?;

    if start > end {
        return Err(format!(
            "Invalid gene interval at line {line_number}: start ({start}) is greater than end ({end})"
        )
        .into());
    }
    if start == 0 || end == 0 {
        return Err(format!(
            "Invalid gene interval at line {line_number}: coordinates must be 1-based positive integers (start={start}, end={end})"
        )
        .into());
    }

    Ok((start, end))
}

pub(crate) fn gene_name_from_gff(attrs: &HashMap<String, String>) -> String {
    let primary = attrs
        .get("gene")
        .or_else(|| attrs.get("Name"))
        .or_else(|| attrs.get("locus_tag"))
        .or_else(|| attrs.get("ID"))
        .map(|value| value.trim_start_matches("gene-").to_string())
        .unwrap_or_else(|| "unknown_gene".to_string());

    let suffix = attrs
        .get("locus_tag")
        .map(std::string::ToString::to_string)
        .unwrap_or_else(|| primary.clone());

    format!("{primary}_{suffix}")
}

fn load_genes_from_tsv(genes_file: &str, snp_list: &[VcfPosition]) -> AppResult<Vec<Gene>> {
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

        if has_snp_in_interval(snp_list, start, end) {
            genes_with_snps += 1;
            genes.push(crate::variants::Gene {
                name: fields[0].to_string(),
                start,
                end,
                strand,
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

/// A raw parsed gene record from a GFF line, before any filtering.
pub(crate) struct GffGeneRecord {
    pub(crate) contig: String,
    pub(crate) gene: Gene,
}

/// Parse a GFF/GFF3 file, yielding one GffGeneRecord per matching feature type.
pub(crate) fn parse_gff_gene_records(
    genes_file: &str,
    feature_types: &[String],
) -> AppResult<Vec<GffGeneRecord>> {
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for (line_idx, line) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let entry = line.map_err(|e| {
            format!("Failed to read line {line_number} in GFF/GFF3 annotation file: {e}")
        })?;
        let trimmed = entry.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() != 9 {
            return Err(format!(
                "Invalid GFF/GFF3 entry at line {}: expected 9 tab-separated fields, got {}. Content: '{}'",
                line_number,
                fields.len(),
                entry
            )
            .into());
        }

        if !feature_types.iter().any(|ft| ft == fields[2]) {
            continue;
        }

        let (start, end) = parse_interval(fields[3], fields[4], line_number)?;
        let strand = parse_strand(fields[6], line_number)?;
        let attrs = parse_gff_attributes(fields[8]);
        let gene_name = gene_name_from_gff(&attrs);

        records.push(GffGeneRecord {
            contig: fields[0].to_string(),
            gene: Gene {
                name: gene_name,
                start,
                end,
                strand,
            },
        });
    }

    Ok(records)
}

pub(crate) fn load_genes_from_gff(
    genes_file: &str,
    snp_list: &[VcfPosition],
    expected_chrom: Option<&str>,
    feature_types: &[String],
) -> AppResult<Vec<Gene>> {
    log::info!("Loading GFF/GFF3 annotation file: {genes_file}");
    let records = parse_gff_gene_records(genes_file, feature_types)?;

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

        if has_snp_in_interval(snp_list, rec.gene.start, rec.gene.end) {
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

pub fn load_genes(
    genes_file: &str,
    snp_list: &[VcfPosition],
    expected_chrom: Option<&str>,
    feature_types: &[String],
) -> AppResult<Vec<Gene>> {
    match detect_annotation_format(genes_file)? {
        AnnotationFormat::Tsv => load_genes_from_tsv(genes_file, snp_list),
        AnnotationFormat::Gff => {
            load_genes_from_gff(genes_file, snp_list, expected_chrom, feature_types)
        }
    }
}

/// Pre-load all gene entries from a GFF/GFF3 file, grouped by contig.
/// No SNP filtering is applied; use `filter_genes_with_snps` afterwards.
pub fn preload_gff_genes(
    genes_file: &str,
    feature_types: &[String],
) -> AppResult<HashMap<String, Vec<Gene>>> {
    log::info!("Pre-loading GFF/GFF3 annotation file: {genes_file}");
    let records = parse_gff_gene_records(genes_file, feature_types)?;
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

/// Filter genes that overlap with at least one SNP position.
/// Requires `snp_list` to be sorted by position.
pub fn filter_genes_with_snps(genes: &[Gene], snp_list: &[VcfPosition]) -> Vec<Gene> {
    genes
        .iter()
        .filter(|gene| has_snp_in_interval(snp_list, gene.start, gene.end))
        .cloned()
        .collect()
}
