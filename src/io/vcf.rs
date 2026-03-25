//! VCF loading, metrics extraction, allele normalisation, and original INFO
//! field preservation.

use super::validation::validate_vcf_allele;
use crate::error::AppResult;
use rust_htslib::bcf;
use rust_htslib::bcf::header::{HeaderRecord, TagType};
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::Read;
use std::collections::{HashMap, HashSet};

const GET_MNV_INFO_TAGS: &[&[u8]] = &[
    b"GENE", b"AA", b"CT", b"TYPE", b"ODP", b"OFREQ", b"SR", b"SRF", b"SRR", b"MR", b"MRF", b"MRR",
    b"DP", b"FREQ", b"SBP", b"MSBP",
];

#[derive(Debug, Clone)]
pub struct VcfPosition {
    pub position: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub original_dp: Option<usize>,
    pub original_freq: Option<f64>,
    /// Pre-serialised original INFO fields (semicolon-separated) to carry over
    /// when `--keep-original-info` is active.
    pub original_info: Option<String>,
}

fn parse_optional_freq_token(raw_token: &str) -> Option<f64> {
    let token = raw_token.trim();
    if token.is_empty() || token == "." {
        return None;
    }
    let has_percent = token.ends_with('%');
    let numeric = token.trim_end_matches('%').parse::<f64>().ok()?;
    if has_percent || (numeric > 1.0 && numeric <= 100.0) {
        Some(numeric / 100.0)
    } else {
        Some(numeric)
    }
}

fn parse_optional_freq_index(raw: &str, alt_index: usize) -> Option<f64> {
    let values = raw.split(',').map(str::trim).collect::<Vec<_>>();
    if values.is_empty() {
        return None;
    }
    let token = if alt_index < values.len() {
        values[alt_index]
    } else if values.len() == 1 {
        values[0]
    } else {
        return None;
    };
    parse_optional_freq_token(token)
}

pub(crate) fn parse_optional_depth(raw: &str) -> Option<usize> {
    let raw_value = raw.trim();
    let first = raw_value.split(',').next()?;
    if first.is_empty() || first == "." {
        return None;
    }
    if let Ok(value) = first.parse::<usize>() {
        return Some(value);
    }
    first
        .parse::<f64>()
        .ok()
        .filter(|value| *value >= 0.0)
        .map(|value| value.round() as usize)
}

fn first_non_missing_i32(values: &[i32]) -> Option<usize> {
    values
        .iter()
        .copied()
        .find(|value| !value.is_missing() && *value >= 0)
        .map(|value| value as usize)
}

fn first_non_missing_f32(values: &[f32]) -> Option<f64> {
    values
        .iter()
        .copied()
        .find(|value| !value.is_missing())
        .map(f64::from)
}

fn non_missing_i32_at(values: &[i32], index: usize) -> Option<usize> {
    if let Some(value) = values.get(index).copied() {
        if !value.is_missing() && value >= 0 {
            return Some(value as usize);
        }
    }
    if values.len() == 1 {
        return first_non_missing_i32(values);
    }
    None
}

fn non_missing_f32_at(values: &[f32], index: usize) -> Option<f64> {
    if let Some(value) = values.get(index).copied() {
        if !value.is_missing() {
            return Some(f64::from(value));
        }
    }
    if values.len() == 1 {
        return first_non_missing_f32(values);
    }
    None
}

fn derive_freq_from_ad(values: &[i32], alt_index: usize) -> Option<f64> {
    let alt_count = non_missing_i32_at(values, alt_index + 1)?;
    let total = values
        .iter()
        .copied()
        .filter(|value| !value.is_missing() && *value >= 0)
        .map(|value| value as usize)
        .sum::<usize>();
    if total == 0 {
        None
    } else {
        Some(alt_count as f64 / total as f64)
    }
}

fn derive_freq_from_ao_ro(
    ao_values: &[i32],
    ro_values: Option<&[i32]>,
    alt_index: usize,
    dp: Option<usize>,
) -> Option<f64> {
    let alt_count = non_missing_i32_at(ao_values, alt_index)?;
    let total = if let Some(depth) = dp {
        depth
    } else {
        let ao_sum = ao_values
            .iter()
            .copied()
            .filter(|value| !value.is_missing() && *value >= 0)
            .map(|value| value as usize)
            .sum::<usize>();
        let ro = ro_values.and_then(first_non_missing_i32).unwrap_or(0);
        ao_sum + ro
    };
    if total == 0 {
        None
    } else {
        Some(alt_count as f64 / total as f64)
    }
}

fn parse_sample_names(header: &bcf::header::HeaderView) -> Vec<String> {
    if header.sample_count() == 0 {
        return Vec::new();
    }
    header
        .samples()
        .iter()
        .map(|raw| String::from_utf8_lossy(raw).to_string())
        .collect()
}

pub fn list_vcf_samples(vcf_file: &str) -> AppResult<Vec<String>> {
    let vcf = bcf::Reader::from_path(vcf_file)?;
    Ok(parse_sample_names(vcf.header()))
}

fn resolve_sample_index(
    header: &bcf::header::HeaderView,
    sample_name: Option<&str>,
) -> AppResult<Option<usize>> {
    let samples = parse_sample_names(header);

    if let Some(name) = sample_name {
        if samples.is_empty() {
            return Err(format!("Requested sample '{name}' but VCF has no sample columns").into());
        }
        if let Some(index) = samples.iter().position(|value| value == name) {
            return Ok(Some(index));
        }
        return Err(format!(
            "Sample '{}' not found in VCF header. Available samples: {}",
            name,
            samples.join(", ")
        )
        .into());
    }

    if samples.is_empty() {
        return Ok(None);
    }
    if samples.len() > 1 {
        log::warn!(
            "Multi-sample VCF detected ({} samples). Using first sample '{}' for original FORMAT metrics. Use --sample to select another sample.",
            samples.len(),
            samples[0]
        );
    }
    Ok(Some(0))
}

fn parse_original_metrics(
    rec: &bcf::Record,
    sample_index: Option<usize>,
    alt_index: usize,
) -> (Option<usize>, Option<f64>) {
    let mut original_dp: Option<usize> = None;
    let mut original_freq: Option<f64> = None;

    if let Some(sample_idx) = sample_index {
        if let Ok(format_dp) = rec.format(b"DP").integer() {
            if let Some(sample_values) = format_dp.get(sample_idx) {
                original_dp = first_non_missing_i32(sample_values);
            }
        }
        if original_dp.is_none() {
            if let Ok(format_dp_str) = rec.format(b"DP").string() {
                if let Some(sample_values) = format_dp_str.get(sample_idx) {
                    if let Ok(raw) = std::str::from_utf8(sample_values) {
                        original_dp = parse_optional_depth(raw);
                    }
                }
            }
        }

        for tag in [b"FREQ".as_slice(), b"AF".as_slice()] {
            if original_freq.is_none() {
                if let Ok(format_freq) = rec.format(tag).float() {
                    if let Some(sample_values) = format_freq.get(sample_idx) {
                        original_freq = non_missing_f32_at(sample_values, alt_index);
                    }
                }
            }
            if original_freq.is_none() {
                if let Ok(format_freq_str) = rec.format(tag).string() {
                    if let Some(sample_values) = format_freq_str.get(sample_idx) {
                        if let Ok(raw) = std::str::from_utf8(sample_values) {
                            original_freq = parse_optional_freq_index(raw, alt_index);
                        }
                    }
                }
            }
        }

        if original_freq.is_none() {
            if let Ok(format_ad) = rec.format(b"AD").integer() {
                if let Some(sample_values) = format_ad.get(sample_idx) {
                    original_freq = derive_freq_from_ad(sample_values, alt_index);
                }
            }
        }
        if original_freq.is_none() {
            if let Ok(format_ao) = rec.format(b"AO").integer() {
                if let Some(ao_values) = format_ao.get(sample_idx) {
                    let ro_values =
                        rec.format(b"RO").integer().ok().and_then(|matrix| {
                            matrix.get(sample_idx).map(|values| values.to_vec())
                        });
                    original_freq = derive_freq_from_ao_ro(
                        ao_values,
                        ro_values.as_deref(),
                        alt_index,
                        original_dp,
                    );
                }
            }
        }
    }

    if original_dp.is_none() {
        if let Ok(Some(info_dp)) = rec.info(b"DP").integer() {
            original_dp = first_non_missing_i32(&info_dp);
        }
    }
    if original_dp.is_none() {
        if let Ok(Some(info_dp_str)) = rec.info(b"DP").string() {
            if let Some(raw) = info_dp_str.first() {
                if let Ok(value) = std::str::from_utf8(raw) {
                    original_dp = parse_optional_depth(value);
                }
            }
        }
    }

    for tag in [b"AF".as_slice(), b"FREQ".as_slice()] {
        if original_freq.is_none() {
            if let Ok(Some(info_freq)) = rec.info(tag).float() {
                original_freq = non_missing_f32_at(&info_freq, alt_index);
            }
        }
        if original_freq.is_none() {
            if let Ok(Some(info_freq_str)) = rec.info(tag).string() {
                if let Some(raw) = info_freq_str.first() {
                    if let Ok(value) = std::str::from_utf8(raw) {
                        original_freq = parse_optional_freq_index(value, alt_index);
                    }
                }
            }
        }
    }

    if original_freq.is_none() {
        if let Ok(Some(info_ad)) = rec.info(b"AD").integer() {
            original_freq = derive_freq_from_ad(&info_ad, alt_index);
        }
    }
    if original_freq.is_none() {
        if let Ok(Some(info_ao)) = rec.info(b"AO").integer() {
            let info_ro = rec.info(b"RO").integer().ok().flatten();
            original_freq = derive_freq_from_ao_ro(
                &info_ao,
                info_ro.map(|values| &**values),
                alt_index,
                original_dp,
            );
        }
    }

    (original_dp, original_freq)
}

fn normalize_ref_alt(pos: usize, ref_allele: &str, alt_allele: &str) -> (usize, String, String) {
    let is_symbolic = alt_allele.starts_with('<') && alt_allele.ends_with('>');
    if is_symbolic {
        return (pos, ref_allele.to_string(), alt_allele.to_string());
    }

    let ref_chars: Vec<char> = ref_allele.chars().collect();
    let alt_chars: Vec<char> = alt_allele.chars().collect();
    let mut start = 0usize;
    let mut ref_end = ref_chars.len();
    let mut alt_end = alt_chars.len();

    while ref_end - start > 1
        && alt_end - start > 1
        && ref_chars[ref_end - 1] == alt_chars[alt_end - 1]
    {
        ref_end -= 1;
        alt_end -= 1;
    }
    while ref_end - start > 1 && alt_end - start > 1 && ref_chars[start] == alt_chars[start] {
        start += 1;
    }

    (
        pos + start,
        ref_chars[start..ref_end].iter().collect(),
        alt_chars[start..alt_end].iter().collect(),
    )
}

/// Collect the IDs and types of original INFO fields (those not defined by
/// get_mnv itself) from a VCF header.
fn original_info_tags(header: &bcf::header::HeaderView) -> Vec<(String, TagType)> {
    let own: HashSet<&[u8]> = GET_MNV_INFO_TAGS.iter().copied().collect();
    let mut tags = Vec::new();
    for rec in header.header_records() {
        if let HeaderRecord::Info { values, .. } = rec {
            if let Some(id) = values.get("ID") {
                if !own.contains(id.as_bytes()) {
                    let tag_type = header
                        .info_type(id.as_bytes())
                        .map(|(t, _)| t)
                        .unwrap_or(TagType::String);
                    tags.push((id.clone(), tag_type));
                }
            }
        }
    }
    tags
}

/// Reconstruct the `##INFO=<...>` header line for a given tag from the
/// structured `HeaderRecord`.
fn info_header_line(values: &linear_map::LinearMap<String, String>) -> String {
    let mut parts = Vec::new();
    // Emit in canonical order: ID, Number, Type, Description, then rest
    for key in &["ID", "Number", "Type", "Description"] {
        if let Some(val) = values.get(*key) {
            if *key == "Description" {
                parts.push(format!("{key}=\"{val}\""));
            } else {
                parts.push(format!("{key}={val}"));
            }
        }
    }
    for (key, val) in values.iter() {
        if !["ID", "Number", "Type", "Description"].contains(&key.as_str()) {
            parts.push(format!("{key}={val}"));
        }
    }
    format!("##INFO=<{}>", parts.join(","))
}

/// Collect the original `##INFO=` header lines from a VCF whose ID is **not**
/// one of the tags get_mnv writes itself. Returns the full `##INFO=<...>`
/// lines ready to be re-emitted in the output VCF header.
pub fn extract_original_info_headers(vcf_file: &str) -> AppResult<Vec<String>> {
    let vcf = bcf::Reader::from_path(vcf_file)?;
    let own: HashSet<&[u8]> = GET_MNV_INFO_TAGS.iter().copied().collect();
    let mut lines = Vec::new();
    for rec in vcf.header().header_records() {
        if let HeaderRecord::Info { values, .. } = rec {
            if let Some(id) = values.get("ID") {
                if !own.contains(id.as_bytes()) {
                    lines.push(info_header_line(&values));
                }
            }
        }
    }
    Ok(lines)
}

/// Build a pre-serialised string of original INFO fields from a single VCF
/// record, excluding get_mnv's own tags. `tags` is the pre-computed list of
/// (ID, TagType) pairs obtained from `original_info_tags()`.
fn extract_original_info(rec: &bcf::Record, tags: &[(String, TagType)]) -> Option<String> {
    let mut parts: Vec<String> = Vec::new();

    for (id, tag_type) in tags {
        let tag = id.as_bytes();
        match tag_type {
            TagType::Integer => {
                if let Ok(Some(values)) = rec.info(tag).integer() {
                    let formatted: Vec<String> = values
                        .iter()
                        .filter(|&&v| v != i32::missing())
                        .map(std::string::ToString::to_string)
                        .collect();
                    if !formatted.is_empty() {
                        parts.push(format!("{}={}", id, formatted.join(",")));
                    }
                }
            }
            TagType::Float => {
                if let Ok(Some(values)) = rec.info(tag).float() {
                    let formatted: Vec<String> = values
                        .iter()
                        .filter(|&&v| !v.is_nan() && v != f32::missing())
                        .map(|v| format!("{v}"))
                        .collect();
                    if !formatted.is_empty() {
                        parts.push(format!("{}={}", id, formatted.join(",")));
                    }
                }
            }
            TagType::Flag => {
                if let Ok(true) = rec.info(tag).flag() {
                    parts.push(id.clone());
                }
            }
            // TagType::String and anything else
            _ => {
                if let Ok(Some(values)) = rec.info(tag).string() {
                    let formatted: Vec<String> = values
                        .iter()
                        .filter_map(|v| {
                            let s = std::str::from_utf8(v).ok()?;
                            if s == "." {
                                None
                            } else {
                                Some(s.to_string())
                            }
                        })
                        .collect();
                    if !formatted.is_empty() {
                        parts.push(format!("{}={}", id, formatted.join(",")));
                    }
                }
            }
        }
    }

    if parts.is_empty() {
        None
    } else {
        Some(parts.join(";"))
    }
}

pub fn load_vcf_positions_by_contig(
    vcf_file: &str,
    sample_name: Option<&str>,
    split_multiallelic: bool,
    normalize_alleles: bool,
    keep_original_info: bool,
) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    log::info!("Loading VCF: {vcf_file}");
    let mut vcf = bcf::Reader::from_path(vcf_file)?;
    let sample_index = resolve_sample_index(vcf.header(), sample_name)?;
    // Pre-compute the list of original INFO tags before the records loop
    // to avoid a simultaneous mutable+immutable borrow on `vcf`.
    let orig_info_tags = if keep_original_info {
        original_info_tags(vcf.header())
    } else {
        Vec::new()
    };
    let mut positions_by_contig: HashMap<String, Vec<VcfPosition>> = HashMap::new();
    let mut split_count = 0usize;

    for (record_idx, rec_result) in vcf.records().enumerate() {
        let rec = rec_result
            .map_err(|e| format!("Failed to parse VCF record {}: {}", record_idx + 1, e))?;
        let raw_pos = rec.pos();
        if raw_pos < 0 {
            return Err(format!(
                "Invalid VCF position at record {}: {}",
                record_idx + 1,
                raw_pos
            )
            .into());
        }
        let pos = raw_pos as usize + 1;
        let alleles = rec.alleles();
        let rid = rec.rid().ok_or_else(|| {
            format!(
                "VCF record {} at position {} has undefined contig in header. \
Add matching ##contig lines to the VCF header.",
                record_idx + 1,
                pos
            )
        })?;
        let chrom_bytes = rec.header().rid2name(rid).map_err(|e| {
            format!(
                "Failed to resolve contig name for VCF record {} (RID {}): {}",
                record_idx + 1,
                rid,
                e
            )
        })?;
        let chrom = std::str::from_utf8(chrom_bytes).map_err(|e| {
            format!(
                "Invalid UTF-8 in contig name for VCF record {}: {}",
                record_idx + 1,
                e
            )
        })?;

        if alleles.len() < 2 {
            return Err(format!(
                "Invalid VCF record {} at position {}: expected at least REF and one ALT allele",
                record_idx + 1,
                pos
            )
            .into());
        }
        if alleles.len() > 2 && !split_multiallelic {
            return Err(format!(
                "Multiallelic VCF record {} at {}:{} is not supported. Split multiallelic sites first (e.g. bcftools norm -m -).",
                record_idx + 1,
                chrom,
                pos
            )
            .into());
        }

        let ref_allele = std::str::from_utf8(alleles[0]).map_err(|e| {
            format!(
                "Invalid UTF-8 in REF allele at VCF record {} (pos {}): {}",
                record_idx + 1,
                pos,
                e
            )
        })?;
        if ref_allele.is_empty() {
            return Err(format!(
                "Invalid VCF allele at record {} (pos {}): empty REF",
                record_idx + 1,
                pos
            )
            .into());
        }
        validate_vcf_allele(ref_allele, record_idx + 1, chrom, pos, "REF")?;

        for (alt_idx, alt_raw) in alleles.iter().enumerate().skip(1) {
            let alt_allele = std::str::from_utf8(alt_raw).map_err(|e| {
                format!(
                    "Invalid UTF-8 in ALT allele {} at VCF record {} (pos {}): {}",
                    alt_idx,
                    record_idx + 1,
                    pos,
                    e
                )
            })?;
            if alt_allele.is_empty() || alt_allele == "." {
                return Err(format!(
                    "Invalid VCF allele at record {} (pos {}): REF='{}' ALT='{}'",
                    record_idx + 1,
                    pos,
                    ref_allele,
                    alt_allele
                )
                .into());
            }
            let (normalized_pos, normalized_ref, normalized_alt) = if normalize_alleles {
                normalize_ref_alt(pos, ref_allele, alt_allele)
            } else {
                (pos, ref_allele.to_string(), alt_allele.to_string())
            };
            validate_vcf_allele(
                &normalized_ref,
                record_idx + 1,
                chrom,
                normalized_pos,
                "REF",
            )?;
            validate_vcf_allele(
                &normalized_alt,
                record_idx + 1,
                chrom,
                normalized_pos,
                "ALT",
            )?;
            let (original_dp, original_freq) =
                parse_original_metrics(&rec, sample_index, alt_idx - 1);
            let original_info = if keep_original_info {
                extract_original_info(&rec, &orig_info_tags)
            } else {
                None
            };

            positions_by_contig
                .entry(chrom.to_string())
                .or_default()
                .push(VcfPosition {
                    position: normalized_pos,
                    ref_allele: normalized_ref,
                    alt_allele: normalized_alt,
                    original_dp,
                    original_freq,
                    original_info,
                });
        }
        if alleles.len() > 2 {
            split_count += alleles.len() - 2;
        }
    }

    for values in positions_by_contig.values_mut() {
        values.sort_by_key(|v| v.position);
    }

    if split_multiallelic && split_count > 0 {
        log::info!("Split {split_count} additional ALT alleles from multiallelic VCF records");
    }

    Ok(positions_by_contig)
}
