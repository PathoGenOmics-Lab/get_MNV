//! iVar `variants` TSV parser.
//!
//! iVar reports variant calls as a tab-separated table with columns such as
//! REGION, POS, REF, ALT, TOTAL_DP, ALT_FREQ and PASS.  get_MNV internally
//! works with VCF-like positions, so this parser maps passing single-nucleotide
//! iVar rows onto `VcfPosition`.

use super::validation::validate_vcf_allele;
use super::vcf::{parse_optional_depth, VcfPosition};
use crate::error::AppResult;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
struct IvarColumns {
    region: usize,
    pos: usize,
    ref_allele: usize,
    alt_allele: usize,
    total_dp: Option<usize>,
    ref_dp: Option<usize>,
    alt_dp: Option<usize>,
    alt_freq: Option<usize>,
    pass: Option<usize>,
}

fn header_key(value: &str) -> String {
    value.trim().trim_start_matches('#').to_ascii_uppercase()
}

fn column_index(headers: &[String], name: &str) -> Option<usize> {
    headers.iter().position(|header| header == name)
}

fn parse_header(line: &str) -> Option<IvarColumns> {
    let headers = line.split('\t').map(header_key).collect::<Vec<_>>();

    Some(IvarColumns {
        region: column_index(&headers, "REGION")?,
        pos: column_index(&headers, "POS")?,
        ref_allele: column_index(&headers, "REF")?,
        alt_allele: column_index(&headers, "ALT")?,
        total_dp: column_index(&headers, "TOTAL_DP"),
        ref_dp: column_index(&headers, "REF_DP"),
        alt_dp: column_index(&headers, "ALT_DP"),
        alt_freq: column_index(&headers, "ALT_FREQ"),
        pass: column_index(&headers, "PASS"),
    })
}

fn get_field<'a>(
    fields: &'a [&str],
    idx: usize,
    name: &str,
    record_idx: usize,
) -> AppResult<&'a str> {
    fields
        .get(idx)
        .copied()
        .ok_or_else(|| format!("iVar record {record_idx}: missing {name} column").into())
}

fn optional_field<'a>(fields: &'a [&str], idx: Option<usize>) -> Option<&'a str> {
    fields.get(idx?).copied()
}

fn parse_freq(raw: &str) -> Option<f64> {
    let token = raw.trim();
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

fn parse_count(raw: &str) -> Option<usize> {
    parse_optional_depth(raw)
}

fn parse_original_metrics(cols: &IvarColumns, fields: &[&str]) -> (Option<usize>, Option<f64>) {
    let total_dp = optional_field(fields, cols.total_dp)
        .and_then(parse_count)
        .or_else(|| {
            let ref_dp = optional_field(fields, cols.ref_dp).and_then(parse_count)?;
            let alt_dp = optional_field(fields, cols.alt_dp).and_then(parse_count)?;
            Some(ref_dp + alt_dp)
        });

    let alt_freq = optional_field(fields, cols.alt_freq)
        .and_then(parse_freq)
        .or_else(|| {
            let alt_dp = optional_field(fields, cols.alt_dp).and_then(parse_count)?;
            let total = total_dp?;
            if total > 0 {
                Some(alt_dp as f64 / total as f64)
            } else {
                None
            }
        });

    (total_dp, alt_freq)
}

fn row_passes(raw: Option<&str>) -> bool {
    match raw.map(str::trim) {
        None | Some("") | Some(".") => true,
        Some(value) => matches!(
            value.to_ascii_uppercase().as_str(),
            "TRUE" | "PASS" | "1" | "YES"
        ),
    }
}

fn is_single_nucleotide(ref_allele: &str, alt_allele: &str) -> bool {
    ref_allele.chars().count() == 1
        && alt_allele.chars().count() == 1
        && !alt_allele.starts_with('+')
        && !alt_allele.starts_with('-')
}

fn open_reader(path: &str) -> AppResult<BufReader<std::fs::File>> {
    let file = std::fs::File::open(path)
        .map_err(|e| format!("Cannot open iVar TSV file '{}': {}", path, e))?;
    Ok(BufReader::with_capacity(64 * 1024, file))
}

/// Return true when the first data-like line has the iVar TSV header shape.
pub fn looks_like_ivar_tsv(path: &str) -> AppResult<bool> {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") || lower.ends_with(".bcf") {
        return Ok(false);
    }

    let reader = open_reader(path)?;
    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Error reading iVar TSV '{}': {}", path, e))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        return Ok(parse_header(trimmed).is_some());
    }
    Ok(false)
}

/// Parse an iVar variants TSV file into VCF-like positions grouped by contig.
pub fn load_ivar_tsv(path: &str) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    log::info!("Loading iVar TSV: {path}");
    let reader = open_reader(path)?;

    let mut columns: Option<IvarColumns> = None;
    let mut positions_by_contig: HashMap<String, Vec<VcfPosition>> = HashMap::new();
    let mut record_idx = 0usize;
    let mut kept = 0usize;
    let mut skipped_failed = 0usize;
    let mut skipped_non_snv = 0usize;

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Error reading iVar TSV line: {e}"))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if columns.is_none() {
            columns = parse_header(trimmed);
            if columns.is_some() {
                continue;
            }
            return Err(format!(
                "No iVar TSV header found in '{path}'. Expected columns REGION, POS, REF and ALT."
            )
            .into());
        }

        record_idx += 1;
        let cols = columns.as_ref().expect("columns are set");
        let fields = line.split('\t').collect::<Vec<_>>();

        if !row_passes(optional_field(&fields, cols.pass)) {
            skipped_failed += 1;
            continue;
        }

        let chrom = get_field(&fields, cols.region, "REGION", record_idx)?.trim();
        let pos_raw = get_field(&fields, cols.pos, "POS", record_idx)?.trim();
        let pos: usize = pos_raw
            .parse()
            .map_err(|_| format!("iVar record {record_idx}: invalid POS '{pos_raw}'"))?;
        if pos == 0 {
            return Err(format!("iVar record {record_idx}: POS cannot be 0").into());
        }

        let ref_allele = get_field(&fields, cols.ref_allele, "REF", record_idx)?
            .trim()
            .to_ascii_uppercase();
        let alt_allele = get_field(&fields, cols.alt_allele, "ALT", record_idx)?
            .trim()
            .to_ascii_uppercase();

        if !is_single_nucleotide(&ref_allele, &alt_allele) || ref_allele == alt_allele {
            skipped_non_snv += 1;
            continue;
        }

        validate_vcf_allele(&ref_allele, record_idx, chrom, pos, "REF")?;
        validate_vcf_allele(&alt_allele, record_idx, chrom, pos, "ALT")?;

        let (original_dp, original_freq) = parse_original_metrics(cols, &fields);
        positions_by_contig
            .entry(chrom.to_string())
            .or_default()
            .push(VcfPosition {
                position: pos,
                ref_allele,
                alt_allele,
                original_dp,
                original_freq,
                original_info: None,
            });
        kept += 1;
    }

    if columns.is_none() {
        return Err(format!(
            "No iVar TSV header found in '{path}'. Expected columns REGION, POS, REF and ALT."
        )
        .into());
    }

    for values in positions_by_contig.values_mut() {
        values.sort_by_key(|v| v.position);
    }

    log::info!(
        "Loaded iVar TSV records: kept={}, skipped_failed={}, skipped_non_snv={}",
        kept,
        skipped_failed,
        skipped_non_snv
    );

    Ok(positions_by_contig)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp_ivar(content: &str) -> String {
        let path = std::env::temp_dir().join(format!(
            "get_mnv_ivar_test_{}.tsv",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut file = std::fs::File::create(&path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        path.to_string_lossy().into_owned()
    }

    #[test]
    fn test_looks_like_ivar_tsv() {
        let path = write_temp_ivar("REGION\tPOS\tREF\tALT\tTOTAL_DP\tALT_FREQ\tPASS\n");
        assert!(looks_like_ivar_tsv(&path).unwrap());
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_load_ivar_tsv_keeps_passing_snvs_only() {
        let path = write_temp_ivar(
            "REGION\tPOS\tREF\tALT\tREF_DP\tALT_DP\tALT_FREQ\tTOTAL_DP\tPASS\n\
chr1\t2\tC\tT\t5\t5\t0.5\t10\tTRUE\n\
chr1\t3\tG\t+A\t5\t5\t0.5\t10\tTRUE\n\
chr1\t4\tT\tA\t5\t5\t0.5\t10\tFALSE\n\
chr1\t5\tA\tA\t5\t5\t0.5\t10\tTRUE\n",
        );
        let parsed = load_ivar_tsv(&path).unwrap();
        let chr1 = parsed.get("chr1").unwrap();
        assert_eq!(chr1.len(), 1);
        assert_eq!(chr1[0].position, 2);
        assert_eq!(chr1[0].ref_allele, "C");
        assert_eq!(chr1[0].alt_allele, "T");
        assert_eq!(chr1[0].original_dp, Some(10));
        assert_eq!(chr1[0].original_freq, Some(0.5));
        let _ = std::fs::remove_file(path);
    }
}
