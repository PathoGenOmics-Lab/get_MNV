//! Fast plain-text VCF parser for uncompressed `.vcf` files. This avoids
//! unnecessary BGZF/header work and keeps common VCF parsing on the fastest
//! path.
//!
//! `.vcf.gz` and `.bcf` files are handled by the parser in `vcf.rs`, which
//! decompresses or decodes them to the same VCF text first.

use super::validation::validate_vcf_allele;
use super::vcf::{
    freq_from_ad, freq_from_ao_ro, normalize_ref_alt, parse_optional_depth,
    parse_optional_freq_index, VcfPosition,
};
use super::vcf_text::{ReadLine, Splitting, VcfLine, VcfText};
use crate::error::AppResult;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::sync::Arc;

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
    references: Option<&crate::io::ReferenceMap>,
) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    log::info!("Loading VCF (fast text parser): {vcf_file}");
    let file = std::fs::File::open(vcf_file)
        .map_err(|e| format!("Cannot open VCF file '{}': {}", vcf_file, e))?;
    let reader = BufReader::with_capacity(64 * 1024, file);

    // Per-allele (Number=A/R/G) INFO fields, so they can be subset to a single
    // ALT when a multiallelic record is split (avoids invalid copied arrays).
    let per_allele_info = if keep_original_info {
        crate::io::vcf::per_allele_info_numbers(&extract_text_info_headers(vcf_file)?)
    } else {
        HashMap::new()
    };
    let mut records = TextRecords::new(
        sample_name,
        split_multiallelic,
        normalize_alleles,
        keep_original_info,
        references,
        per_allele_info,
    );
    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Error reading VCF line: {e}"))?;
        records.read(&ReadLine {
            text: &line,
            splitting: Splitting::AsWritten,
        })?;
    }
    records.finish(vcf_file)
}

/// [`load_vcf_text`] for a file a cohort holds in memory: the same lines, the
/// same checks and the same log lines, without reading the file again.
pub(crate) fn load_held_vcf_text(
    text: &VcfText,
    sample_name: Option<&str>,
    split_multiallelic: bool,
    normalize_alleles: bool,
    keep_original_info: bool,
    references: Option<&crate::io::ReferenceMap>,
) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    log::info!("Loading VCF (fast text parser): {}", text.path());
    let per_allele_info = if keep_original_info {
        crate::io::vcf::per_allele_info_numbers(&text_info_headers_from(
            text.replay().map(|line| line.map(str::to_string)),
        )?)
    } else {
        HashMap::new()
    };
    let mut records = TextRecords::new(
        sample_name,
        split_multiallelic,
        normalize_alleles,
        keep_original_info,
        references,
        per_allele_info,
    );
    for line in text.lines() {
        records.read(line)?;
    }
    if let Some(e) = text.unread() {
        return Err(format!("Error reading VCF line: {e}").into());
    }
    records.finish(text.path())
}

/// The plain-text reader's state as it goes through the file's lines.
struct TextRecords<'a> {
    sample_name: Option<&'a str>,
    split_multiallelic: bool,
    normalize_alleles: bool,
    keep_original_info: bool,
    references: Option<&'a crate::io::ReferenceMap>,
    per_allele_info: HashMap<String, char>,
    sample_names: Arc<Vec<String>>,
    sample_index: Option<usize>,
    header_seen: bool,
    positions_by_contig: HashMap<String, Vec<VcfPosition>>,
    split_count: usize,
    not_carried: usize,
    overlapped: usize,
    spelled_out: usize,
    record_idx: usize,
}

impl<'a> TextRecords<'a> {
    fn new(
        sample_name: Option<&'a str>,
        split_multiallelic: bool,
        normalize_alleles: bool,
        keep_original_info: bool,
        references: Option<&'a crate::io::ReferenceMap>,
        per_allele_info: HashMap<String, char>,
    ) -> Self {
        TextRecords {
            sample_name,
            split_multiallelic,
            normalize_alleles,
            keep_original_info,
            references,
            per_allele_info,
            sample_names: Arc::new(Vec::new()),
            sample_index: None,
            header_seen: false,
            positions_by_contig: HashMap::new(),
            split_count: 0,
            not_carried: 0,
            overlapped: 0,
            spelled_out: 0,
            record_idx: 0,
        }
    }

    fn read(&mut self, line: &impl VcfLine) -> AppResult<()> {
        let text = line.text();
        if text.starts_with("##") {
            return Ok(());
        }
        if super::vcf::is_column_header(text) {
            // Extract sample names (columns 9+)
            if let Some(names) = line.header_samples() {
                self.sample_names = names;
            }
            // Resolve sample index
            self.sample_index = resolve_text_sample_index(&self.sample_names, self.sample_name)?;
            self.header_seen = true;
            return Ok(());
        }
        // A blank line, often the last one, is not a record: the reader used
        // for `.vcf.gz` skips it, and this one stopped the run on it.
        if !self.header_seen || text.trim().is_empty() {
            return Ok(());
        }

        self.record_idx += 1;
        let record_idx = self.record_idx;
        let record = line.record(self.sample_index.map(|si| 9 + si));
        let cols = record.fixed;
        if record.count < 8 {
            return Err(format!(
                "VCF record {record_idx}: expected at least 8 columns, got {}",
                record.count
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
        if alt_alleles.len() > 1 && !self.split_multiallelic {
            return Err(format!(
                "Multiallelic VCF record {} at {}:{} is not supported. \
Split multiallelic sites first (e.g. bcftools norm -m -).",
                record_idx, chrom, pos
            )
            .into());
        }

        // Parse FORMAT fields for DP and FREQ
        let format_keys: Vec<&str> = if record.count > 8 {
            cols[8].split(':').collect()
        } else {
            Vec::new()
        };
        let sample_values: Vec<&str> = record
            .sample
            .map(|column| column.split(':').collect())
            .unwrap_or_default();

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

        // Register the contig as soon as a record for it is read, before any
        // allele is kept or skipped. An input whose every allele is skipped, a
        // cohort sample that carries none of them, then still reports the contig
        // it covered, which is what separates "this run has nothing to annotate"
        // from "this file holds no records at all". Only the second is an error.
        self.positions_by_contig
            .entry(chrom.to_string())
            .or_default();

        // The INFO keys get_MNV writes itself, dropped from what keep_original_info
        // keeps. One list for both readers: a copy here missed the keys added since.
        let get_mnv_tags = super::vcf::GET_MNV_INFO_TAGS;

        let mut pushed_here = 0usize;
        for (alt_idx, alt_allele) in alt_alleles.iter().enumerate() {
            if alt_allele.is_empty() || *alt_allele == "." {
                return Err(format!(
                    "VCF record {record_idx} at pos {pos}: invalid ALT '{alt_allele}'"
                )
                .into());
            }

            if !crate::io::vcf::sample_carries_alt(genotype, alt_idx) {
                self.not_carried += 1;
                continue;
            }
            if *alt_allele == "*" {
                self.overlapped += 1;
                continue;
            }
            let sequence = self
                .references
                .and_then(|refs| refs.get(chrom))
                .map(String::as_str);
            let (start, explicit_ref, explicit_alt) = match crate::io::vcf::spell_out_symbolic(
                pos, ref_allele, alt_allele, cols[7], alt_idx, sequence,
            ) {
                Some(spelled) => {
                    self.spelled_out += 1;
                    spelled
                }
                None => (pos, ref_allele.to_string(), alt_allele.to_string()),
            };
            let (norm_pos, norm_ref, norm_alt) = if self.normalize_alleles {
                normalize_ref_alt(start, &explicit_ref, &explicit_alt)
            } else {
                (start, explicit_ref, explicit_alt)
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
            let original_info = if self.keep_original_info {
                crate::io::vcf::retarget_input_ann(
                    extract_text_original_info(
                        cols[7],
                        get_mnv_tags,
                        &self.per_allele_info,
                        alt_idx,
                    ),
                    alt_allele,
                    norm_pos,
                    &norm_ref,
                    &norm_alt,
                )
            } else {
                None
            };

            self.positions_by_contig
                .entry(chrom.to_string())
                .or_default()
                .push(VcfPosition {
                    record_start: norm_pos,
                    ref_allele: norm_ref,
                    alt_allele: norm_alt,
                    original_dp,
                    original_freq,
                    original_info,
                    declared_phase: crate::io::vcf::parse_declared_phase(
                        genotype, phase_set, alt_idx,
                    ),
                });
            pushed_here += 1;
        }
        self.split_count += pushed_here.saturating_sub(1);
        Ok(())
    }

    fn finish(self, vcf_file: &str) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
        let TextRecords {
            split_multiallelic,
            header_seen,
            mut positions_by_contig,
            split_count,
            not_carried,
            overlapped,
            spelled_out,
            ..
        } = self;
        if !header_seen {
            return Err(format!("No #CHROM header line found in '{vcf_file}'").into());
        }

        for values in positions_by_contig.values_mut() {
            values.sort_by_key(|v| v.record_start);
        }

        if split_multiallelic && split_count > 0 {
            log::info!("Split {split_count} additional ALT alleles from multiallelic VCF records");
        }
        if not_carried > 0 {
            log::info!(
                "Skipped {not_carried} ALT alleles the selected sample's genotype does not carry"
            );
        }
        crate::io::vcf::log_overlapped(overlapped);
        crate::io::vcf::log_spelled_out(spelled_out);

        Ok(positions_by_contig)
    }
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

    // Try FORMAT:FREQ, one value per ALT like every other per-allele tag. This
    // parsed the whole field as a single number, so a multiallelic record's
    // `30%,10%` parsed as nothing at all and the frequency was silently
    // dropped: the reader used for `.vcf.gz` picks the element for this ALT and
    // the one used for plain `.vcf` did not, so the same bytes in the two
    // containers disagreed. An indel with no known frequency propagates its
    // frame shift by design, so a synonymous substitution downstream came out
    // frameshift_variant at HIGH from the plain file and synonymous_variant at
    // LOW from the gzipped one. The AF branch below already did this.
    if original_freq.is_none() {
        if let Some(idx) = freq_idx {
            if let Some(val) = sample_values.get(idx) {
                original_freq = parse_optional_freq_index(val, alt_index);
            }
        }
    }

    // Try FORMAT:AF
    if original_freq.is_none() {
        if let Some(idx) = af_idx {
            if let Some(val) = sample_values.get(idx) {
                original_freq = parse_optional_freq_index(val, alt_index);
            }
        }
    }

    // Try FORMAT:AD → derive freq
    if original_freq.is_none() {
        if let Some(idx) = ad_idx {
            if let Some(val) = sample_values.get(idx) {
                original_freq = freq_from_ad(val, alt_index);
            }
        }
    }

    // Try FORMAT:AO/RO → derive freq (FreeBayes-style)
    if original_freq.is_none() {
        if let Some(idx) = ao_idx {
            if let Some(ao_val) = sample_values.get(idx) {
                let ro_val = ro_idx.and_then(|i| sample_values.get(i)).copied();
                original_freq = freq_from_ao_ro(ao_val, ro_val, original_dp, alt_index);
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
                original_freq = parse_optional_freq_index(val, alt_index);
                if original_freq.is_some() {
                    break;
                }
            }
        }
    }

    // Fallback to INFO:AD
    if original_freq.is_none() {
        if let Some(ad_val) = find_info_tag(info_field, "AD") {
            original_freq = freq_from_ad(ad_val, alt_index);
        }
    }

    // Fallback to INFO:AO/RO (FreeBayes-style)
    if original_freq.is_none() {
        if let Some(ao_val) = find_info_tag(info_field, "AO") {
            let ro_val = find_info_tag(info_field, "RO");
            original_freq = freq_from_ao_ro(ao_val, ro_val, original_dp, alt_index);
        }
    }

    (original_dp, original_freq)
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
        if super::vcf::is_column_header(&line) {
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
    let file = std::fs::File::open(vcf_file)
        .map_err(|e| format!("Cannot open VCF file '{}': {}", vcf_file, e))?;
    text_info_headers_from(BufReader::new(file).lines())
}

/// The original INFO header lines of [`extract_text_info_headers`] from the
/// file's lines, wherever they come from.
pub(crate) fn text_info_headers_from(
    lines: impl Iterator<Item = std::io::Result<String>>,
) -> AppResult<Vec<String>> {
    let mut headers = Vec::new();
    for line_result in lines {
        let line = line_result.map_err(|e| format!("Error reading VCF: {e}"))?;
        if !line.starts_with("##") {
            break;
        }
        // Match ##INFO=<ID=TAG,...>
        if let Some(rest) = line.strip_prefix("##INFO=<ID=") {
            if let Some(comma_pos) = rest.find(',') {
                let tag = &rest[..comma_pos];
                if !super::vcf::header_written_by_get_mnv(tag) {
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
        let result = load_vcf_text(vcf, None, false, false, false, None).unwrap();
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
        assert!((crate::io::vcf::parse_optional_freq_token("50%").unwrap() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_parse_freq_token_decimal() {
        assert!((crate::io::vcf::parse_optional_freq_token("0.25").unwrap() - 0.25).abs() < 1e-6);
    }

    #[test]
    fn test_parse_freq_token_missing() {
        assert!(crate::io::vcf::parse_optional_freq_token(".").is_none());
        assert!(crate::io::vcf::parse_optional_freq_token("").is_none());
    }

    #[test]
    fn test_parse_freq_indexed_multiallelic() {
        assert!((parse_optional_freq_index("0.1,0.2,0.3", 1).unwrap() - 0.2).abs() < 1e-6);
    }

    #[test]
    fn test_parse_freq_indexed_out_of_bounds_single() {
        // If only one value and index is beyond, use that single value
        assert!(parse_optional_freq_index("0.5", 3).is_some());
    }

    #[test]
    fn test_derive_freq_from_ad() {
        // REF=10, ALT1=5, ALT2=3 → alt_index=0 → 5/18
        let freq = freq_from_ad("10,5,3", 0).unwrap();
        assert!((freq - 5.0 / 18.0).abs() < 1e-6);
    }

    #[test]
    fn test_derive_freq_from_ad_zero_total() {
        assert!(freq_from_ad("0,0", 0).is_none());
    }

    #[test]
    fn test_derive_freq_from_ad_missing_values() {
        assert!(freq_from_ad(".,.", 0).is_none());
    }

    /// A `.` keeps its place: dropped before indexing, the first ALT of
    /// `.,5,3` took the second ALT's 3 reads.
    #[test]
    fn a_missing_count_does_not_move_the_others_up() {
        assert!(
            freq_from_ad(".,5,3", 0).is_none(),
            "REF unknown, total unknown"
        );
        assert!(
            freq_from_ad("10,.,3", 0).is_none(),
            "this ALT's own count unknown"
        );
        let merged = freq_from_ad("10,5,.", 0).unwrap();
        assert!((merged - 5.0 / 15.0).abs() < 1e-9, "{merged}");
        let second = freq_from_ad("10,.,3", 1).unwrap();
        assert!((second - 3.0 / 13.0).abs() < 1e-9, "{second}");
        assert!(freq_from_ao_ro(".,4", Some("6"), None, 0).is_none());
        let ao = freq_from_ao_ro(".,4", Some("6"), None, 1).unwrap();
        assert!((ao - 0.4).abs() < 1e-9, "{ao}");
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
