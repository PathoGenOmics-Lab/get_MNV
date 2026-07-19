//! VCF loading for plain and BGZF-compressed VCF files, metrics extraction,
//! allele normalisation, and original INFO field preservation.
//!
//! Plain `.vcf` files use the fast text parser in `vcf_fast.rs`.

use super::validation::validate_vcf_allele;
use crate::error::AppResult;
use crate::variants::{
    decompose_allele, substitution_components, AlleleComponent, AlleleComponentKind, AlleleEvent,
};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};

const GET_MNV_INFO_TAGS: &[&str] = &[
    "GENE", "AA", "CT", "TYPE", "ODP", "OFREQ", "SR", "SRF", "SRR", "MR", "MRF", "MRR", "DP",
    "FREQ", "SBP", "MSBP", "EC", "COMP", "ER", "ERF", "ERR", "EDP", "EFREQ", "SO", "IMPACT", "GD",
    "MNVSHIFT",
];

#[derive(Debug, Clone)]
pub struct VcfPosition {
    pub position: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub original_dp: Option<usize>,
    pub original_freq: Option<f64>,
    pub original_info: Option<String>,
}

impl VcfPosition {
    pub fn event(&self) -> AlleleEvent {
        decompose_allele(self.position, &self.ref_allele, &self.alt_allele)
    }

    pub fn substitution_components(&self) -> Vec<AlleleComponent> {
        substitution_components(self.position, &self.ref_allele, &self.alt_allele)
    }

    pub fn is_single_nucleotide_substitution(&self) -> bool {
        let components = self.substitution_components();
        components.len() == 1
            && self.ref_allele.chars().count() == 1
            && self.alt_allele.chars().count() == 1
    }

    pub fn overlaps_interval(&self, start: usize, end: usize) -> bool {
        if start == 0 || end == 0 || start > end {
            return false;
        }
        let event = self.event();
        event
            .components
            .iter()
            .any(|component| match component.kind {
                AlleleComponentKind::Snp => {
                    component.position >= start && component.position <= end
                }
                AlleleComponentKind::Insertion => {
                    component.position >= start && component.position < end
                }
                AlleleComponentKind::Deletion
                | AlleleComponentKind::Delins
                | AlleleComponentKind::Symbolic => {
                    let component_end = component
                        .position
                        .saturating_add(component.ref_allele.len().saturating_sub(1));
                    component.position <= end && component_end >= start
                }
            })
    }
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

pub(crate) fn normalize_ref_alt(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> (usize, String, String) {
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

    let norm_ref: String = ref_chars[start..ref_end].iter().collect();
    let norm_alt: String = alt_chars[start..alt_end].iter().collect();

    if norm_ref.is_empty() || norm_alt.is_empty() {
        return (pos, ref_allele.to_string(), alt_allele.to_string());
    }

    (pos + start, norm_ref, norm_alt)
}

/// Parse a VCF line into fields. Returns None if the line is a header or empty.
fn parse_vcf_line(line: &str) -> Option<Vec<&str>> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    Some(line.split('\t').collect())
}

/// Extract INFO field value by key from a raw INFO string.
fn get_info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    for field in info.split(';') {
        if let Some((k, v)) = field.split_once('=') {
            if k == key {
                return Some(v);
            }
        } else if field == key {
            return Some("");
        }
    }
    None
}

/// Extract FORMAT field value for a given sample index.
fn get_format_value<'a>(format_keys: &[&str], sample_field: &'a str, key: &str) -> Option<&'a str> {
    let idx = format_keys.iter().position(|k| *k == key)?;
    let values: Vec<&str> = sample_field.split(':').collect();
    values
        .get(idx)
        .copied()
        .filter(|v| !v.is_empty() && *v != ".")
}

fn parse_original_metrics_from_fields(
    info: &str,
    format_keys: &[&str],
    sample_field: Option<&str>,
    alt_index: usize,
) -> (Option<usize>, Option<f64>) {
    let mut original_dp: Option<usize> = None;
    let mut original_freq: Option<f64> = None;

    // Try FORMAT fields first (sample-specific)
    if let Some(sample) = sample_field {
        // DP
        if let Some(dp_str) = get_format_value(format_keys, sample, "DP") {
            original_dp = parse_optional_depth(dp_str);
        }
        // FREQ / AF
        for tag in &["FREQ", "AF"] {
            if original_freq.is_none() {
                if let Some(val) = get_format_value(format_keys, sample, tag) {
                    original_freq = parse_optional_freq_index(val, alt_index);
                }
            }
        }
        // AD → derive freq
        if original_freq.is_none() {
            if let Some(ad_str) = get_format_value(format_keys, sample, "AD") {
                let values: Vec<i64> = ad_str
                    .split(',')
                    .filter_map(|s| s.trim().parse().ok())
                    .collect();
                if let Some(&alt_count) = values.get(alt_index + 1) {
                    let total: i64 = values.iter().filter(|v| **v >= 0).sum();
                    if total > 0 && alt_count >= 0 {
                        original_freq = Some(alt_count as f64 / total as f64);
                    }
                }
            }
        }
        // AO/RO → derive freq
        if original_freq.is_none() {
            if let Some(ao_str) = get_format_value(format_keys, sample, "AO") {
                let ao_values: Vec<i64> = ao_str
                    .split(',')
                    .filter_map(|s| s.trim().parse().ok())
                    .collect();
                if let Some(&alt_count) = ao_values.get(alt_index) {
                    let ro = get_format_value(format_keys, sample, "RO")
                        .and_then(|s| s.trim().parse::<i64>().ok())
                        .unwrap_or(0);
                    let total = if let Some(dp) = original_dp {
                        dp as i64
                    } else {
                        let ao_sum: i64 = ao_values.iter().filter(|v| **v >= 0).sum();
                        ao_sum + ro
                    };
                    if total > 0 && alt_count >= 0 {
                        original_freq = Some(alt_count as f64 / total as f64);
                    }
                }
            }
        }
    }

    // Fall back to INFO fields
    if original_dp.is_none() {
        if let Some(dp_str) = get_info_value(info, "DP") {
            original_dp = parse_optional_depth(dp_str);
        }
    }
    for tag in &["AF", "FREQ"] {
        if original_freq.is_none() {
            if let Some(val) = get_info_value(info, tag) {
                original_freq = parse_optional_freq_index(val, alt_index);
            }
        }
    }
    if original_freq.is_none() {
        if let Some(ad_str) = get_info_value(info, "AD") {
            let values: Vec<i64> = ad_str
                .split(',')
                .filter_map(|s| s.trim().parse().ok())
                .collect();
            if let Some(&alt_count) = values.get(alt_index + 1) {
                let total: i64 = values.iter().filter(|v| **v >= 0).sum();
                if total > 0 && alt_count >= 0 {
                    original_freq = Some(alt_count as f64 / total as f64);
                }
            }
        }
    }
    if original_freq.is_none() {
        if let Some(ao_str) = get_info_value(info, "AO") {
            let ao_values: Vec<i64> = ao_str
                .split(',')
                .filter_map(|s| s.trim().parse().ok())
                .collect();
            if let Some(&alt_count) = ao_values.get(alt_index) {
                let ro = get_info_value(info, "RO")
                    .and_then(|s| s.trim().parse::<i64>().ok())
                    .unwrap_or(0);
                let total = if let Some(dp) = original_dp {
                    dp as i64
                } else {
                    let ao_sum: i64 = ao_values.iter().filter(|v| **v >= 0).sum();
                    ao_sum + ro
                };
                if total > 0 && alt_count >= 0 {
                    original_freq = Some(alt_count as f64 / total as f64);
                }
            }
        }
    }

    (original_dp, original_freq)
}

/// Read a `key=...` attribute value out of a `##INFO` header line, up to the
/// next `,` or `>`.
fn header_attr(line: &str, key: &str) -> Option<String> {
    let start = line.find(key)? + key.len();
    let rest = &line[start..];
    let end = rest.find([',', '>']).unwrap_or(rest.len());
    Some(rest[..end].to_string())
}

/// Parse `Number=` from `##INFO` header lines into a map of field ID → its
/// per-allele cardinality code (`A`, `R` or `G`). Fixed-Number fields (1, 0,
/// ., a constant) are omitted because they need no per-allele subsetting.
pub(crate) fn per_allele_info_numbers(header_lines: &[String]) -> HashMap<String, char> {
    let mut map = HashMap::new();
    for line in header_lines {
        if let (Some(id), Some(number)) = (header_attr(line, "ID="), header_attr(line, "Number=")) {
            match number.as_str() {
                "A" => drop(map.insert(id, 'A')),
                "R" => drop(map.insert(id, 'R')),
                "G" => drop(map.insert(id, 'G')),
                _ => {}
            }
        }
    }
    map
}

/// Subset a per-allele INFO field value to a single ALT so it is valid on a
/// single-ALT output record. `Number=A` keeps the value at `alt_idx`;
/// `Number=R` keeps `[REF, ALT@alt_idx]`. `Number=G` and arrays too short for
/// `alt_idx` are dropped (cannot be safely attributed to one allele).
fn subset_per_allele_field(key: &str, value: &str, kind: char, alt_idx: usize) -> Option<String> {
    let elems: Vec<&str> = value.split(',').collect();
    match kind {
        'A' => elems.get(alt_idx).map(|v| format!("{key}={v}")),
        'R' => match (elems.first(), elems.get(alt_idx + 1)) {
            (Some(r), Some(a)) => Some(format!("{key}={r},{a}")),
            _ => None,
        },
        _ => None,
    }
}

/// Filter get_mnv's own INFO tags out of an original record's INFO string and
/// subset any per-allele field to `alt_idx`, so the result is valid when copied
/// onto a single-ALT output record (`--keep-original-info` + multiallelic
/// input previously copied whole `Number=A/R` arrays, producing invalid VCF).
pub(crate) fn filter_and_subset_original_info<F: Fn(&str) -> bool>(
    info: &str,
    is_own_tag: F,
    per_allele: &HashMap<String, char>,
    alt_idx: usize,
) -> Option<String> {
    if info.is_empty() || info == "." {
        return None;
    }
    let parts: Vec<String> = info
        .split(';')
        .filter_map(|field| {
            let (key, value) = match field.split_once('=') {
                Some((k, v)) => (k, Some(v)),
                None => (field, None),
            };
            if is_own_tag(key) {
                return None;
            }
            match (per_allele.get(key), value) {
                (Some(&kind), Some(v)) => subset_per_allele_field(key, v, kind, alt_idx),
                _ => Some(field.to_string()),
            }
        })
        .collect();
    if parts.is_empty() {
        None
    } else {
        Some(parts.join(";"))
    }
}

pub fn list_vcf_samples(vcf_file: &str) -> AppResult<Vec<String>> {
    if vcf_file.ends_with(".bcf") {
        return Err(
            "BCF input is not supported. Convert to VCF first: bcftools view input.bcf > input.vcf"
                .into(),
        );
    }
    let reader: Box<dyn BufRead> = if vcf_file.ends_with(".gz") {
        let file = std::fs::File::open(vcf_file)?;
        let bgzf = noodles::bgzf::io::Reader::new(file);
        Box::new(BufReader::new(bgzf))
    } else {
        let file = std::fs::File::open(vcf_file)?;
        Box::new(BufReader::new(file))
    };

    // Parse header to find #CHROM line
    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("#CHROM") || line.starts_with("#chrom") {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() > 9 {
                return Ok(fields[9..].iter().map(|s| s.to_string()).collect());
            }
            return Ok(Vec::new());
        }
        if !line.starts_with('#') {
            break;
        }
    }
    Ok(Vec::new())
}

fn resolve_sample_index(samples: &[String], sample_name: Option<&str>) -> AppResult<Option<usize>> {
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

pub fn extract_original_info_headers(vcf_file: &str) -> AppResult<Vec<String>> {
    let own: HashSet<&str> = GET_MNV_INFO_TAGS.iter().copied().collect();
    let mut lines = Vec::new();

    if vcf_file.ends_with(".bcf") {
        return Err(
            "BCF input is not supported. Convert to VCF first: bcftools view input.bcf > input.vcf"
                .into(),
        );
    }
    let reader: Box<dyn BufRead> = if vcf_file.ends_with(".gz") {
        let file = std::fs::File::open(vcf_file)?;
        let bgzf = noodles::bgzf::io::Reader::new(file);
        Box::new(BufReader::new(bgzf))
    } else {
        let file = std::fs::File::open(vcf_file)?;
        Box::new(BufReader::new(file))
    };

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("##INFO=") {
            // Extract ID from the header line
            if let Some(id_start) = line.find("ID=") {
                let rest = &line[id_start + 3..];
                if let Some(id_end) = rest.find([',', '>']) {
                    let id = &rest[..id_end];
                    if !own.contains(id) {
                        lines.push(line.clone());
                    }
                }
            }
        }
        if !line.starts_with('#') {
            break;
        }
    }
    Ok(lines)
}

pub fn load_vcf_positions_by_contig(
    vcf_file: &str,
    sample_name: Option<&str>,
    split_multiallelic: bool,
    normalize_alleles: bool,
    keep_original_info: bool,
) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    log::info!("Loading VCF: {vcf_file}");

    let samples = list_vcf_samples(vcf_file)?;
    let sample_index = resolve_sample_index(&samples, sample_name)?;
    let own_tags: HashSet<&str> = GET_MNV_INFO_TAGS.iter().copied().collect();
    let per_allele_info = if keep_original_info {
        per_allele_info_numbers(&extract_original_info_headers(vcf_file)?)
    } else {
        HashMap::new()
    };

    if vcf_file.ends_with(".bcf") {
        return Err(
            "BCF input is not supported. Convert to VCF first: bcftools view input.bcf > input.vcf"
                .into(),
        );
    }
    let reader: Box<dyn BufRead> = if vcf_file.ends_with(".gz") {
        let file = std::fs::File::open(vcf_file)?;
        let bgzf = noodles::bgzf::io::Reader::new(file);
        Box::new(BufReader::new(bgzf))
    } else {
        let file = std::fs::File::open(vcf_file)?;
        Box::new(BufReader::new(file))
    };

    let mut positions_by_contig: HashMap<String, Vec<VcfPosition>> = HashMap::new();
    let mut split_count = 0usize;
    let mut record_idx = 0usize;
    let mut header_seen = false;

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("#CHROM") || line.starts_with("#chrom") {
            header_seen = true;
        }
        let fields = match parse_vcf_line(&line) {
            Some(f) => f,
            None => continue,
        };
        // Skip any data line appearing before the #CHROM header, matching the
        // plain-text fast parser. A well-formed VCF has no records before it.
        if !header_seen {
            continue;
        }
        record_idx += 1;

        if fields.len() < 8 {
            return Err(format!("VCF record {} has fewer than 8 fields", record_idx).into());
        }

        let chrom = fields[0];
        let pos: usize = fields[1].parse().map_err(|_| {
            format!(
                "Invalid VCF position at record {}: {}",
                record_idx, fields[1]
            )
        })?;
        if pos == 0 {
            return Err(format!(
                "Invalid VCF position 0 at record {record_idx} (VCF coordinates are 1-based)"
            )
            .into());
        }
        let ref_allele = fields[3];
        let alt_field = fields[4];
        let info = fields[7];

        if ref_allele.is_empty() {
            return Err(format!(
                "Invalid VCF allele at record {} (pos {}): empty REF",
                record_idx, pos
            )
            .into());
        }
        validate_vcf_allele(ref_allele, record_idx, chrom, pos, "REF")?;

        let alts: Vec<&str> = alt_field.split(',').collect();
        if alts.len() > 1 && !split_multiallelic {
            return Err(format!(
                "Multiallelic VCF record {} at {}:{} is not supported. Split multiallelic sites first (e.g. bcftools norm -m -).",
                record_idx, chrom, pos
            ).into());
        }

        // Get FORMAT and sample fields
        let format_keys: Vec<&str> = if fields.len() > 8 {
            fields[8].split(':').collect()
        } else {
            Vec::new()
        };
        let sample_field = sample_index.and_then(|idx| fields.get(9 + idx).copied());

        for (alt_idx, alt_allele) in alts.iter().enumerate() {
            if alt_allele.is_empty() || *alt_allele == "." {
                return Err(format!(
                    "Invalid VCF allele at record {} (pos {}): REF='{}' ALT='{}'",
                    record_idx, pos, ref_allele, alt_allele
                )
                .into());
            }
            let (normalized_pos, normalized_ref, normalized_alt) = if normalize_alleles {
                normalize_ref_alt(pos, ref_allele, alt_allele)
            } else {
                (pos, ref_allele.to_string(), alt_allele.to_string())
            };
            validate_vcf_allele(&normalized_ref, record_idx, chrom, normalized_pos, "REF")?;
            validate_vcf_allele(&normalized_alt, record_idx, chrom, normalized_pos, "ALT")?;

            let (original_dp, original_freq) =
                parse_original_metrics_from_fields(info, &format_keys, sample_field, alt_idx);
            let original_info = if keep_original_info {
                filter_and_subset_original_info(
                    info,
                    |k| own_tags.contains(k),
                    &per_allele_info,
                    alt_idx,
                )
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
        if alts.len() > 1 {
            split_count += alts.len() - 1;
        }
    }

    if !header_seen {
        return Err("No #CHROM header line found in VCF".into());
    }

    for values in positions_by_contig.values_mut() {
        values.sort_by_key(|v| v.position);
    }

    if split_multiallelic && split_count > 0 {
        log::info!("Split {split_count} additional ALT alleles from multiallelic VCF records");
    }

    Ok(positions_by_contig)
}

// BCF input is not supported in the pure-Rust build.
// Use `bcftools view input.bcf > input.vcf` to convert.

#[cfg(test)]
mod tests {
    use super::VcfPosition;

    fn variant(position: usize, ref_allele: &str, alt_allele: &str) -> VcfPosition {
        VcfPosition {
            position,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }
    }

    #[test]
    fn test_list_vcf_samples_accepts_lowercase_chrom_header() {
        use std::io::Write;
        let path =
            std::env::temp_dir().join(format!("get_mnv_lc_samples_{}.vcf", std::process::id()));
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "#chrom\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A"
        )
        .unwrap();
        drop(f);
        let samples = super::list_vcf_samples(path.to_str().unwrap()).unwrap();
        assert_eq!(samples, vec!["SAMPLE_A".to_string()]);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_load_vcf_accepts_lowercase_chrom_and_skips_preheader_data() {
        use std::io::Write;
        let path = std::env::temp_dir().join(format!("get_mnv_lc_load_{}.vcf", std::process::id()));
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        // A stray data-shaped line before the header must be ignored, and a
        // lowercase #chrom header must be accepted (parity with the fast parser).
        writeln!(f, "chr1\t1\t.\tA\tG\t.\tPASS\t.").unwrap();
        writeln!(f, "#chrom\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(f, "chr1\t5\t.\tA\tT\t.\tPASS\t.").unwrap();
        drop(f);
        let by_contig =
            super::load_vcf_positions_by_contig(path.to_str().unwrap(), None, false, false, false)
                .unwrap();
        assert_eq!(by_contig.get("chr1").map(std::vec::Vec::len), Some(1));
        assert_eq!(by_contig["chr1"][0].position, 5);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_per_allele_info_numbers_parses_number() {
        let headers = vec![
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"alt count\">".to_string(),
            "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"allele depth\">".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">".to_string(),
        ];
        let map = super::per_allele_info_numbers(&headers);
        assert_eq!(map.get("AC"), Some(&'A'));
        assert_eq!(map.get("AD"), Some(&'R'));
        assert_eq!(map.get("DP"), None); // fixed-Number field omitted
    }

    #[test]
    fn test_filter_and_subset_original_info_multiallelic() {
        let headers = vec![
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">".to_string(),
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">".to_string(),
            "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"x\">".to_string(),
        ];
        let map = super::per_allele_info_numbers(&headers);
        let info = "AC=5,7;AF=0.5,0.7;AD=10,5,7;DP=22;FOO=bar";
        // Split allele index 1: per-allele fields subset to that ALT.
        let out = super::filter_and_subset_original_info(info, |k| k == "GENE", &map, 1).unwrap();
        assert_eq!(out, "AC=7;AF=0.7;AD=10,7;DP=22;FOO=bar");
    }

    #[test]
    fn test_filter_and_subset_original_info_biallelic_unchanged() {
        let headers = vec![
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">".to_string(),
            "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"x\">".to_string(),
        ];
        let map = super::per_allele_info_numbers(&headers);
        let out =
            super::filter_and_subset_original_info("AC=5;AD=10,5", |_| false, &map, 0).unwrap();
        assert_eq!(out, "AC=5;AD=10,5");
    }

    #[test]
    fn insertion_overlaps_only_between_interval_bases() {
        let insertion = variant(10, "A", "AT");
        assert!(insertion.overlaps_interval(9, 11));
        assert!(insertion.overlaps_interval(10, 11));
        assert!(!insertion.overlaps_interval(11, 12));
        assert!(!insertion.overlaps_interval(10, 10));
    }

    #[test]
    fn anchored_deletion_overlaps_deleted_reference_span() {
        let deletion = variant(10, "AT", "A");
        assert!(deletion.overlaps_interval(11, 11));
        assert!(deletion.overlaps_interval(10, 11));
        assert!(!deletion.overlaps_interval(10, 10));
    }
}
