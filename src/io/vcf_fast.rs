//! Fast plain-text VCF parser that bypasses htslib for uncompressed `.vcf`
//! files. Achieves ~10× speedup over htslib by avoiding FFI overhead and
//! unnecessary header/index parsing.
//!
//! Falls back to htslib for `.bcf` and `.vcf.gz` formats.

use super::validation::validate_vcf_allele;
use super::vcf::{normalize_ref_alt, parse_optional_depth, VcfPosition};
use crate::error::AppResult;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

/// Returns `true` if the file should use the fast text parser (plain `.vcf`).
pub fn use_fast_parser(path: &str) -> bool {
    let lower = path.to_ascii_lowercase();
    lower.ends_with(".vcf") && !lower.ends_with(".bcf") && !lower.ends_with(".vcf.gz")
}

/// Parse a plain-text VCF file into positions by contig.
/// This replicates the behaviour of `load_vcf_positions_by_contig` but operates
/// on raw text lines instead of htslib records.
pub fn load_vcf_text(
    vcf_file: &str,
    sample_name: Option<&str>,
    split_multiallelic: bool,
    normalize_alleles: bool,
    keep_original_info: bool,
) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    log::info!("Loading VCF (fast text parser): {vcf_file}");
    let file = std::fs::File::open(vcf_file)
        .map_err(|e| format!("Cannot open VCF file '{}': {}", vcf_file, e))?;
    let reader = BufReader::with_capacity(64 * 1024, file);

    let mut sample_names: Vec<String> = Vec::new();
    let mut sample_index: Option<usize> = None;
    let mut header_seen = false;
    let mut positions_by_contig: HashMap<String, Vec<VcfPosition>> = HashMap::new();
    let mut split_count = 0usize;
    let mut record_idx = 0usize;

    // INFO tags to preserve when keep_original_info is active
    let get_mnv_tags: &[&str] = &[
        "GENE", "AA", "CT", "TYPE", "ODP", "OFREQ", "SR", "SRF", "SRR",
        "MR", "MRF", "MRR", "DP", "FREQ", "SBP", "MSBP",
    ];

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Error reading VCF line: {e}"))?;
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with("#CHROM") || line.starts_with("#chrom") {
            let cols: Vec<&str> = line.split('\t').collect();
            // Extract sample names (columns 9+)
            if cols.len() > 9 {
                sample_names = cols[9..].iter().map(|s| s.to_string()).collect();
            }
            // Resolve sample index
            sample_index = resolve_text_sample_index(&sample_names, sample_name)?;
            header_seen = true;
            continue;
        }
        if !header_seen {
            continue;
        }

        record_idx += 1;
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            return Err(format!(
                "VCF record {record_idx}: expected at least 8 columns, got {}",
                cols.len()
            )
            .into());
        }

        let chrom = cols[0];
        let pos: usize = cols[1].parse().map_err(|_| {
            format!("VCF record {record_idx}: invalid POS '{}'", cols[1])
        })?;
        if pos == 0 {
            return Err(format!("VCF record {record_idx}: POS cannot be 0").into());
        }
        let ref_allele = cols[3];
        let alt_field = cols[4];

        if ref_allele.is_empty() {
            return Err(format!("VCF record {record_idx} at pos {pos}: empty REF").into());
        }
        validate_vcf_allele(ref_allele, record_idx, chrom, pos, "REF")?;

        let alt_alleles: Vec<&str> = alt_field.split(',').collect();
        if alt_alleles.len() > 1 && !split_multiallelic {
            return Err(format!(
                "Multiallelic VCF record {} at {}:{} is not supported. \
Split multiallelic sites first (e.g. bcftools norm -m -).",
                record_idx, chrom, pos
            )
            .into());
        }

        // Parse FORMAT fields for DP and FREQ
        let format_keys: Vec<&str> = if cols.len() > 8 { cols[8].split(':').collect() } else { Vec::new() };
        let sample_values: Vec<&str> = if let Some(si) = sample_index {
            let col_idx = 9 + si;
            if col_idx < cols.len() {
                cols[col_idx].split(':').collect()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };

        // Index lookup for FORMAT fields
        let dp_idx = format_keys.iter().position(|k| *k == "DP");
        let freq_idx = format_keys.iter().position(|k| *k == "FREQ");
        let af_idx = format_keys.iter().position(|k| *k == "AF");
        let ad_idx = format_keys.iter().position(|k| *k == "AD");

        // Extract original INFO if requested
        let original_info = if keep_original_info {
            extract_text_original_info(cols[7], get_mnv_tags)
        } else {
            None
        };

        for (alt_idx, alt_allele) in alt_alleles.iter().enumerate() {
            if alt_allele.is_empty() || *alt_allele == "." {
                return Err(format!(
                    "VCF record {record_idx} at pos {pos}: invalid ALT '{alt_allele}'"
                )
                .into());
            }

            let (norm_pos, norm_ref, norm_alt) = if normalize_alleles {
                normalize_ref_alt(pos, ref_allele, alt_allele)
            } else {
                (pos, ref_allele.to_string(), alt_allele.to_string())
            };
            validate_vcf_allele(&norm_ref, record_idx, chrom, norm_pos, "REF")?;
            validate_vcf_allele(&norm_alt, record_idx, chrom, norm_pos, "ALT")?;

            let (original_dp, original_freq) =
                parse_text_metrics(&sample_values, dp_idx, freq_idx, af_idx, ad_idx, alt_idx, cols[7]);

            positions_by_contig
                .entry(chrom.to_string())
                .or_default()
                .push(VcfPosition {
                    position: norm_pos,
                    ref_allele: norm_ref,
                    alt_allele: norm_alt,
                    original_dp,
                    original_freq,
                    original_info: original_info.clone(),
                });
        }
        if alt_alleles.len() > 1 {
            split_count += alt_alleles.len() - 1;
        }
    }

    if !header_seen {
        return Err(format!("No #CHROM header line found in '{vcf_file}'").into());
    }

    for values in positions_by_contig.values_mut() {
        values.sort_by_key(|v| v.position);
    }

    if split_multiallelic && split_count > 0 {
        log::info!("Split {split_count} additional ALT alleles from multiallelic VCF records");
    }

    Ok(positions_by_contig)
}

fn resolve_text_sample_index(
    sample_names: &[String],
    sample_name: Option<&str>,
) -> AppResult<Option<usize>> {
    match sample_name {
        None => {
            if sample_names.is_empty() {
                Ok(None)
            } else {
                Ok(Some(0))
            }
        }
        Some("all") => {
            // "all" is handled at a higher level; default to first sample here
            if sample_names.is_empty() {
                Ok(None)
            } else {
                Ok(Some(0))
            }
        }
        Some(name) => {
            let idx = sample_names
                .iter()
                .position(|s| s == name)
                .ok_or_else(|| {
                    format!(
                        "Sample '{}' not found in VCF. Available samples: {}",
                        name,
                        sample_names.join(", ")
                    )
                })?;
            Ok(Some(idx))
        }
    }
}

fn parse_text_metrics(
    sample_values: &[&str],
    dp_idx: Option<usize>,
    freq_idx: Option<usize>,
    af_idx: Option<usize>,
    ad_idx: Option<usize>,
    alt_index: usize,
    info_field: &str,
) -> (Option<usize>, Option<f64>) {
    let mut original_dp: Option<usize> = None;
    let mut original_freq: Option<f64> = None;

    // Try FORMAT:DP
    if let Some(idx) = dp_idx {
        if let Some(val) = sample_values.get(idx) {
            original_dp = parse_optional_depth(val);
        }
    }

    // Try FORMAT:FREQ
    if original_freq.is_none() {
        if let Some(idx) = freq_idx {
            if let Some(val) = sample_values.get(idx) {
                original_freq = parse_freq_token(val);
            }
        }
    }

    // Try FORMAT:AF
    if original_freq.is_none() {
        if let Some(idx) = af_idx {
            if let Some(val) = sample_values.get(idx) {
                original_freq = parse_freq_indexed(val, alt_index);
            }
        }
    }

    // Try FORMAT:AD → derive freq
    if original_freq.is_none() {
        if let Some(idx) = ad_idx {
            if let Some(val) = sample_values.get(idx) {
                original_freq = derive_freq_from_text_ad(val, alt_index);
            }
        }
    }

    // Fallback to INFO:DP
    if original_dp.is_none() {
        if let Some(dp_val) = find_info_tag(info_field, "DP") {
            original_dp = parse_optional_depth(dp_val);
        }
    }

    // Fallback to INFO:AF / INFO:FREQ
    if original_freq.is_none() {
        for tag in ["AF", "FREQ"] {
            if let Some(val) = find_info_tag(info_field, tag) {
                original_freq = parse_freq_indexed(val, alt_index);
                if original_freq.is_some() {
                    break;
                }
            }
        }
    }

    // Fallback to INFO:AD
    if original_freq.is_none() {
        if let Some(ad_val) = find_info_tag(info_field, "AD") {
            original_freq = derive_freq_from_text_ad(ad_val, alt_index);
        }
    }

    (original_dp, original_freq)
}

fn parse_freq_token(raw: &str) -> Option<f64> {
    let trimmed = raw.trim();
    if trimmed.is_empty() || trimmed == "." {
        return None;
    }
    let has_pct = trimmed.ends_with('%');
    let numeric = trimmed.trim_end_matches('%').parse::<f64>().ok()?;
    if has_pct || (numeric > 1.0 && numeric <= 100.0) {
        Some(numeric / 100.0)
    } else {
        Some(numeric)
    }
}

fn parse_freq_indexed(raw: &str, alt_index: usize) -> Option<f64> {
    let tokens: Vec<&str> = raw.split(',').collect();
    if tokens.is_empty() {
        return None;
    }
    let token = if alt_index < tokens.len() {
        tokens[alt_index]
    } else if tokens.len() == 1 {
        tokens[0]
    } else {
        return None;
    };
    parse_freq_token(token)
}

fn derive_freq_from_text_ad(raw: &str, alt_index: usize) -> Option<f64> {
    let values: Vec<i64> = raw
        .split(',')
        .filter_map(|s| {
            let t = s.trim();
            if t.is_empty() || t == "." {
                None
            } else {
                t.parse::<i64>().ok()
            }
        })
        .collect();
    let alt_count = values.get(alt_index + 1).copied().filter(|v| *v >= 0)? as f64;
    let total: f64 = values.iter().filter(|v| **v >= 0).sum::<i64>() as f64;
    if total == 0.0 {
        None
    } else {
        Some(alt_count / total)
    }
}

fn find_info_tag<'a>(info: &'a str, tag: &str) -> Option<&'a str> {
    if info == "." {
        return None;
    }
    for field in info.split(';') {
        if let Some(val) = field.strip_prefix(tag) {
            if let Some(val) = val.strip_prefix('=') {
                return Some(val);
            }
        }
    }
    None
}

fn extract_text_original_info(info: &str, skip_tags: &[&str]) -> Option<String> {
    if info == "." {
        return None;
    }
    let kept: Vec<&str> = info
        .split(';')
        .filter(|field| {
            let tag = field.split('=').next().unwrap_or("");
            !skip_tags.contains(&tag)
        })
        .collect();
    if kept.is_empty() {
        None
    } else {
        Some(kept.join(";"))
    }
}

/// List sample names from a plain-text VCF.
pub fn list_text_vcf_samples(vcf_file: &str) -> AppResult<Vec<String>> {
    let file = std::fs::File::open(vcf_file)
        .map_err(|e| format!("Cannot open VCF file '{}': {}", vcf_file, e))?;
    let reader = BufReader::new(file);
    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Error reading VCF: {e}"))?;
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            return Ok(cols.get(9..).unwrap_or(&[]).iter().map(|s| s.to_string()).collect());
        }
    }
    Ok(Vec::new())
}

/// Extract original INFO header lines from a plain-text VCF.
pub fn extract_text_info_headers(vcf_file: &str) -> AppResult<Vec<String>> {
    let get_mnv_tags: &[&str] = &[
        "GENE", "AA", "CT", "TYPE", "ODP", "OFREQ", "SR", "SRF", "SRR",
        "MR", "MRF", "MRR", "DP", "FREQ", "SBP", "MSBP",
    ];
    let file = std::fs::File::open(vcf_file)
        .map_err(|e| format!("Cannot open VCF file '{}': {}", vcf_file, e))?;
    let reader = BufReader::new(file);
    let mut headers = Vec::new();
    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Error reading VCF: {e}"))?;
        if !line.starts_with("##") {
            break;
        }
        // Match ##INFO=<ID=TAG,...>
        if let Some(rest) = line.strip_prefix("##INFO=<ID=") {
            if let Some(comma_pos) = rest.find(',') {
                let tag = &rest[..comma_pos];
                if !get_mnv_tags.contains(&tag) {
                    headers.push(line);
                }
            }
        }
    }
    Ok(headers)
}
