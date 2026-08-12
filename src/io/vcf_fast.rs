//! Fast plain-text VCF parser for uncompressed `.vcf` files. This avoids
//! unnecessary BGZF/header work and keeps common VCF parsing on the fastest
//! path.
//!
//! `.vcf.gz` files are handled by the BGZF-aware parser. BCF input is not
//! supported; convert BCF to VCF before running get_MNV.

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
/// on raw text lines instead of the BGZF-aware parser.
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
        "GENE", "AA", "CT", "TYPE", "ODP", "OFREQ", "SR", "SRF", "SRR", "MR", "MRF", "MRR", "DP",
        "FREQ", "SBP", "MSBP", "EC", "COMP", "ER", "ERF", "ERR", "EDP", "EFREQ", "SO", "IMPACT",
        "GD", "MNVSHIFT", "DBS", "MNVPS", "MNVPR", "FSPH", "DPHASE", "LD", "LDP", "NMD", "HGVSG",
        "HGVSC",
    ];

    // Per-allele (Number=A/R/G) INFO fields, so they can be subset to a single
    // ALT when a multiallelic record is split (avoids invalid copied arrays).
    let per_allele_info = if keep_original_info {
        crate::io::vcf::per_allele_info_numbers(&extract_text_info_headers(vcf_file)?)
    } else {
        HashMap::new()
    };

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
        let pos: usize = cols[1]
            .parse()
            .map_err(|_| format!("VCF record {record_idx}: invalid POS '{}'", cols[1]))?;
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
        let format_keys: Vec<&str> = if cols.len() > 8 {
            cols[8].split(':').collect()
        } else {
            Vec::new()
        };
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
        let ao_idx = format_keys.iter().position(|k| *k == "AO");
        let ro_idx = format_keys.iter().position(|k| *k == "RO");
        let gt_idx = format_keys.iter().position(|k| *k == "GT");
        let ps_idx = format_keys.iter().position(|k| *k == "PS");
        let genotype = gt_idx.and_then(|idx| sample_values.get(idx).copied());
        let phase_set = ps_idx.and_then(|idx| sample_values.get(idx).copied());

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

            let (original_dp, original_freq) = parse_text_metrics(
                &sample_values,
                dp_idx,
                freq_idx,
                af_idx,
                ad_idx,
                ao_idx,
                ro_idx,
                alt_idx,
                cols[7],
            );
            let original_info = if keep_original_info {
                extract_text_original_info(cols[7], get_mnv_tags, &per_allele_info, alt_idx)
            } else {
                None
            };

            positions_by_contig
                .entry(chrom.to_string())
                .or_default()
                .push(VcfPosition {
                    position: norm_pos,
                    ref_allele: norm_ref,
                    alt_allele: norm_alt,
                    original_dp,
                    original_freq,
                    original_info,
                    declared_phase: crate::io::vcf::parse_declared_phase(
                        genotype, phase_set, alt_idx,
                    ),
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
            let idx = sample_names.iter().position(|s| s == name).ok_or_else(|| {
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

#[allow(clippy::too_many_arguments)]
fn parse_text_metrics(
    sample_values: &[&str],
    dp_idx: Option<usize>,
    freq_idx: Option<usize>,
    af_idx: Option<usize>,
    ad_idx: Option<usize>,
    ao_idx: Option<usize>,
    ro_idx: Option<usize>,
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

    // Try FORMAT:AO/RO → derive freq (FreeBayes-style)
    if original_freq.is_none() {
        if let Some(idx) = ao_idx {
            if let Some(ao_val) = sample_values.get(idx) {
                let ro_val = ro_idx.and_then(|i| sample_values.get(i)).copied();
                original_freq = derive_freq_from_ao_ro(ao_val, ro_val, original_dp, alt_index);
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

    // Fallback to INFO:AO/RO (FreeBayes-style)
    if original_freq.is_none() {
        if let Some(ao_val) = find_info_tag(info_field, "AO") {
            let ro_val = find_info_tag(info_field, "RO");
            original_freq = derive_freq_from_ao_ro(ao_val, ro_val, original_dp, alt_index);
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

/// Derive an alternate-allele frequency from FreeBayes-style AO (alt observation
/// counts, one per ALT) and RO (reference observation count). The denominator is
/// the original depth when known, otherwise `sum(AO) + RO`. Mirrors the BGZF
/// parser so a plain `.vcf` and its `.vcf.gz` equivalent derive the same OFREQ.
fn derive_freq_from_ao_ro(
    ao_raw: &str,
    ro_raw: Option<&str>,
    original_dp: Option<usize>,
    alt_index: usize,
) -> Option<f64> {
    let ao_values: Vec<i64> = ao_raw
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    let alt_count = *ao_values.get(alt_index)?;
    let ro = ro_raw
        .and_then(|s| s.trim().parse::<i64>().ok())
        .unwrap_or(0);
    let total = if let Some(dp) = original_dp {
        dp as i64
    } else {
        let ao_sum: i64 = ao_values.iter().filter(|v| **v >= 0).sum();
        ao_sum + ro
    };
    if total > 0 && alt_count >= 0 {
        Some(alt_count as f64 / total as f64)
    } else {
        None
    }
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

fn extract_text_original_info(
    info: &str,
    skip_tags: &[&str],
    per_allele: &HashMap<String, char>,
    alt_idx: usize,
) -> Option<String> {
    super::vcf::filter_and_subset_original_info(
        info,
        |k| skip_tags.contains(&k),
        per_allele,
        alt_idx,
    )
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
            return Ok(cols
                .get(9..)
                .unwrap_or(&[])
                .iter()
                .map(|s| s.to_string())
                .collect());
        }
    }
    Ok(Vec::new())
}

/// Extract original INFO header lines from a plain-text VCF.
pub fn extract_text_info_headers(vcf_file: &str) -> AppResult<Vec<String>> {
    let get_mnv_tags: &[&str] = &[
        "GENE", "AA", "CT", "TYPE", "ODP", "OFREQ", "SR", "SRF", "SRR", "MR", "MRF", "MRR", "DP",
        "FREQ", "SBP", "MSBP", "EC", "COMP", "ER", "ERF", "ERR", "EDP", "EFREQ", "SO", "IMPACT",
        "GD", "MNVSHIFT", "DBS", "MNVPS", "MNVPR", "FSPH", "DPHASE", "LD", "LDP", "NMD", "HGVSG",
        "HGVSC",
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_use_fast_parser_vcf() {
        assert!(use_fast_parser("sample.vcf"));
        assert!(use_fast_parser("/path/to/DATA.VCF"));
    }

    #[test]
    fn test_use_fast_parser_not_vcf_gz() {
        assert!(!use_fast_parser("sample.vcf.gz"));
        assert!(!use_fast_parser("sample.bcf"));
    }

    #[test]
    fn test_use_fast_parser_edge_cases() {
        assert!(!use_fast_parser("sample.vcf.gz"));
        assert!(use_fast_parser("sample.VCF"));
        assert!(!use_fast_parser("sample.BCF"));
    }

    #[test]
    fn test_load_vcf_text_example() {
        // Use the bundled example VCF
        let vcf = concat!(env!("CARGO_MANIFEST_DIR"), "/example/G35894.var.snp.vcf");
        let result = load_vcf_text(vcf, None, false, false, false).unwrap();
        // Should have at least one contig
        assert!(!result.is_empty());
        // Total positions should be substantial
        let total: usize = result.values().map(|v| v.len()).sum();
        assert!(total > 700, "Expected >700 positions, got {total}");
    }

    #[test]
    fn test_list_samples_example() {
        let vcf = concat!(env!("CARGO_MANIFEST_DIR"), "/example/G35894.var.snp.vcf");
        let samples = list_text_vcf_samples(vcf).unwrap();
        // Example VCF has at least one sample
        assert!(!samples.is_empty());
    }

    #[test]
    fn test_extract_info_headers_example() {
        let vcf = concat!(env!("CARGO_MANIFEST_DIR"), "/example/G35894.var.snp.vcf");
        let headers = extract_text_info_headers(vcf).unwrap();
        // VCF should have some INFO headers
        for h in &headers {
            assert!(h.starts_with("##INFO="), "Expected INFO header, got: {h}");
        }
    }

    #[test]
    fn test_parse_text_metrics_derives_freq_from_info_ao_ro() {
        // FreeBayes INFO carrying only AO/RO must derive the same OFREQ as the
        // BGZF parser (regression: the fast parser ignored AO/RO entirely, so a
        // plain .vcf and its .vcf.gz produced different OFREQ).
        let (dp, freq) =
            parse_text_metrics(&[], None, None, None, None, None, None, 0, "AO=30;RO=70");
        assert_eq!(dp, None);
        assert_eq!(freq, Some(0.3));
    }

    #[test]
    fn test_parse_text_metrics_derives_freq_from_format_ao_ro() {
        // FORMAT AO at column 0, RO at column 1.
        let sample = ["30", "70"];
        let (_dp, freq) =
            parse_text_metrics(&sample, None, None, None, None, Some(0), Some(1), 0, ".");
        assert_eq!(freq, Some(0.3));
    }

    #[test]
    fn test_parse_freq_token_percent() {
        assert!((parse_freq_token("50%").unwrap() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_parse_freq_token_decimal() {
        assert!((parse_freq_token("0.25").unwrap() - 0.25).abs() < 1e-6);
    }

    #[test]
    fn test_parse_freq_token_missing() {
        assert!(parse_freq_token(".").is_none());
        assert!(parse_freq_token("").is_none());
    }

    #[test]
    fn test_parse_freq_indexed_multiallelic() {
        assert!((parse_freq_indexed("0.1,0.2,0.3", 1).unwrap() - 0.2).abs() < 1e-6);
    }

    #[test]
    fn test_parse_freq_indexed_out_of_bounds_single() {
        // If only one value and index is beyond, use that single value
        assert!(parse_freq_indexed("0.5", 3).is_some());
    }

    #[test]
    fn test_derive_freq_from_ad() {
        // REF=10, ALT1=5, ALT2=3 → alt_index=0 → 5/18
        let freq = derive_freq_from_text_ad("10,5,3", 0).unwrap();
        assert!((freq - 5.0 / 18.0).abs() < 1e-6);
    }

    #[test]
    fn test_derive_freq_from_ad_zero_total() {
        assert!(derive_freq_from_text_ad("0,0", 0).is_none());
    }

    #[test]
    fn test_derive_freq_from_ad_missing_values() {
        assert!(derive_freq_from_text_ad(".,.", 0).is_none());
    }

    #[test]
    fn test_find_info_tag() {
        assert_eq!(find_info_tag("DP=100;AF=0.5;MQ=60", "AF"), Some("0.5"));
        assert_eq!(find_info_tag("DP=100;AF=0.5", "MQ"), None);
        assert_eq!(find_info_tag(".", "DP"), None);
    }

    #[test]
    fn test_extract_original_info_filters_tags() {
        let get_mnv_tags = &["GENE", "AA", "DP"];
        let result = extract_text_original_info(
            "GENE=rpoB;CUSTOM=yes;DP=100",
            get_mnv_tags,
            &std::collections::HashMap::new(),
            0,
        );
        assert_eq!(result.as_deref(), Some("CUSTOM=yes"));
    }

    #[test]
    fn test_extract_original_info_all_filtered() {
        let get_mnv_tags = &["GENE"];
        assert!(extract_text_original_info(
            "GENE=x",
            get_mnv_tags,
            &std::collections::HashMap::new(),
            0
        )
        .is_none());
    }

    #[test]
    fn test_extract_original_info_dot() {
        assert!(
            extract_text_original_info(".", &["GENE"], &std::collections::HashMap::new(), 0)
                .is_none()
        );
    }
}
