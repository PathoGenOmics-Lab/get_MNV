use crate::error::AppResult;
use crate::variants::Gene;
use bio::io::fasta;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::Read;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Error as IoError, ErrorKind};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct VcfPosition {
    pub position: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub original_dp: Option<usize>,
    pub original_freq: Option<f64>,
}

pub struct Reference {
    pub sequence: String,
}

pub type ReferenceMap = HashMap<String, String>;

fn is_iupac_dna_base(base: char) -> bool {
    matches!(
        base,
        'A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N'
    )
}

fn validate_iupac_sequence(value: &str, context: &str) -> AppResult<()> {
    if value.is_empty() {
        return Err(format!("Empty sequence found for {}", context).into());
    }
    if let Some(invalid_base) = value
        .chars()
        .map(|base| base.to_ascii_uppercase())
        .find(|base| !is_iupac_dna_base(*base))
    {
        return Err(format!(
            "Invalid base '{}' in {}. Allowed IUPAC DNA bases: A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N",
            invalid_base, context
        )
        .into());
    }
    Ok(())
}

fn validate_vcf_allele(
    allele: &str,
    record_idx: usize,
    chrom: &str,
    pos: usize,
    allele_kind: &str,
) -> AppResult<()> {
    if allele == "*" {
        return Err(format!(
            "Unsupported '{}' allele '*' at VCF record {} ({}:{}). Use normalized VCF without spanning-deletion alleles.",
            allele_kind, record_idx, chrom, pos
        )
        .into());
    }

    let is_symbolic = allele.starts_with('<') && allele.ends_with('>');
    if is_symbolic {
        if allele_kind == "REF" {
            return Err(format!(
                "Unsupported symbolic REF allele '{}' at VCF record {} ({}:{})",
                allele, record_idx, chrom, pos
            )
            .into());
        }
        return Ok(());
    }

    validate_iupac_sequence(
        allele,
        &format!(
            "{} allele at VCF record {} ({}:{})",
            allele_kind, record_idx, chrom, pos
        ),
    )
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

fn parse_optional_depth(raw: &str) -> Option<usize> {
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
            return Err(
                format!("Requested sample '{}' but VCF has no sample columns", name).into(),
            );
        }
        if let Some(index) = samples.iter().position(|value| value == name) {
            return Ok(Some(index));
        }
        return Err(format!(
            "Sample '{}' not found in VCF header. Available samples: {}",
            name,
            samples.join(",")
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

pub fn load_references(fasta_file: &str) -> AppResult<ReferenceMap> {
    log::info!("Loading reference FASTA: {}", fasta_file);
    let reader = fasta::Reader::new(File::open(fasta_file)?);
    let mut references: ReferenceMap = HashMap::new();
    for (record_idx, record_result) in reader.records().enumerate() {
        let record = record_result.map_err(|e| {
            format!(
                "Failed reading FASTA record {} from '{}': {}",
                record_idx + 1,
                fasta_file,
                e
            )
        })?;
        if record.id().is_empty() {
            return Err(format!(
                "Invalid FASTA record {} in '{}': empty contig name",
                record_idx + 1,
                fasta_file
            )
            .into());
        }
        let sequence = String::from_utf8(record.seq().to_vec()).map_err(|e| {
            format!(
                "Invalid UTF-8 sequence for FASTA contig '{}' in '{}': {}",
                record.id(),
                fasta_file,
                e
            )
        })?;
        let sequence_upper = sequence.to_uppercase();
        validate_iupac_sequence(
            &sequence_upper,
            &format!("FASTA contig '{}' in '{}'", record.id(), fasta_file),
        )?;
        if references
            .insert(record.id().to_string(), sequence_upper)
            .is_some()
        {
            return Err(format!(
                "Duplicate FASTA contig '{}' found in '{}'",
                record.id(),
                fasta_file
            )
            .into());
        }
    }
    if references.is_empty() {
        return Err(format!("No FASTA records found in '{}'", fasta_file).into());
    }
    Ok(references)
}

pub fn reference_for_chrom(references: &ReferenceMap, chrom: &str) -> AppResult<Reference> {
    let sequence = references.get(chrom).ok_or_else(|| {
        format!(
            "Reference contig '{}' is missing in FASTA. Available contigs: {}",
            chrom,
            {
                let mut names = references.keys().cloned().collect::<Vec<_>>();
                names.sort();
                names.join(",")
            }
        )
    })?;
    Ok(Reference {
        sequence: sequence.clone(),
    })
}

fn normalize_ref_alt(pos: usize, ref_allele: &str, alt_allele: &str) -> (usize, String, String) {
    let is_symbolic = alt_allele.starts_with('<') && alt_allele.ends_with('>');
    if is_symbolic {
        return (pos, ref_allele.to_string(), alt_allele.to_string());
    }

    let mut normalized_pos = pos;
    let mut ref_chars = ref_allele.chars().collect::<Vec<_>>();
    let mut alt_chars = alt_allele.chars().collect::<Vec<_>>();

    while ref_chars.len() > 1 && alt_chars.len() > 1 && ref_chars.last() == alt_chars.last() {
        ref_chars.pop();
        alt_chars.pop();
    }
    while ref_chars.len() > 1 && alt_chars.len() > 1 && ref_chars.first() == alt_chars.first() {
        ref_chars.remove(0);
        alt_chars.remove(0);
        normalized_pos += 1;
    }

    (
        normalized_pos,
        ref_chars.iter().collect::<String>(),
        alt_chars.iter().collect::<String>(),
    )
}

pub fn load_vcf_positions_by_contig(
    vcf_file: &str,
    sample_name: Option<&str>,
    split_multiallelic: bool,
    normalize_alleles: bool,
) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
    log::info!("Loading VCF: {}", vcf_file);
    let mut vcf = bcf::Reader::from_path(vcf_file)?;
    let sample_index = resolve_sample_index(vcf.header(), sample_name)?;
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

            positions_by_contig
                .entry(chrom.to_string())
                .or_default()
                .push(VcfPosition {
                    position: normalized_pos,
                    ref_allele: normalized_ref,
                    alt_allele: normalized_alt,
                    original_dp,
                    original_freq,
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
        log::info!(
            "Split {} additional ALT alleles from multiallelic VCF records",
            split_count
        );
    }

    Ok(positions_by_contig)
}

pub fn validate_vcf_reference_alleles(
    chrom: &str,
    snp_list: &[VcfPosition],
    reference: &Reference,
) -> AppResult<()> {
    for site in snp_list {
        if site.position == 0 {
            return Err(format!(
                "Invalid VCF position 0 found at contig '{}' (positions must be 1-based)",
                chrom
            )
            .into());
        }
        let ref_len = site.ref_allele.len();
        let start = site.position - 1;
        let end = start + ref_len;
        if end > reference.sequence.len() {
            return Err(format!(
                "VCF REF out of FASTA bounds at {}:{} (REF='{}', FASTA length={})",
                chrom,
                site.position,
                site.ref_allele,
                reference.sequence.len()
            )
            .into());
        }
        let fasta_ref = &reference.sequence[start..end];
        if !fasta_ref.eq_ignore_ascii_case(&site.ref_allele) {
            return Err(format!(
                "VCF REF/FASTA mismatch at {}:{}: VCF REF='{}' FASTA='{}'",
                chrom, site.position, site.ref_allele, fasta_ref
            )
            .into());
        }
    }
    Ok(())
}

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

fn decode_percent_encoded(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    let mut bytes = value.as_bytes().iter().copied();
    while let Some(current) = bytes.next() {
        if current == b'%' {
            let hi = bytes.next();
            let lo = bytes.next();
            if let (Some(hi), Some(lo)) = (hi, lo) {
                let maybe = [hi, lo];
                if let Ok(hex) = std::str::from_utf8(&maybe) {
                    if let Ok(decoded) = u8::from_str_radix(hex, 16) {
                        out.push(decoded as char);
                        continue;
                    }
                }
                out.push('%');
                out.push(hi as char);
                out.push(lo as char);
                continue;
            }
            out.push('%');
            if let Some(hi) = hi {
                out.push(hi as char);
            }
            if let Some(lo) = lo {
                out.push(lo as char);
            }
            continue;
        }
        out.push(current as char);
    }
    out
}

fn split_attribute_fields(attributes: &str) -> Vec<&str> {
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

fn parse_gff_attributes(attributes: &str) -> HashMap<String, String> {
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

fn has_snp_in_interval(snp_list: &[VcfPosition], start: usize, end: usize) -> bool {
    snp_list
        .iter()
        .any(|snp| snp.position >= start && snp.position <= end)
}

fn parse_strand(raw: &str, line_number: usize) -> AppResult<crate::variants::Strand> {
    raw.parse::<crate::variants::Strand>().map_err(|_| {
        format!(
            "Invalid strand at line {} ('{}'). Expected '+' or '-'",
            line_number, raw
        )
        .into()
    })
}

fn parse_interval(start_raw: &str, end_raw: &str, line_number: usize) -> AppResult<(usize, usize)> {
    let start = start_raw.parse::<usize>().map_err(|e| {
        format!(
            "Invalid start coordinate at line {} ('{}'): {}",
            line_number, start_raw, e
        )
    })?;
    let end = end_raw.parse::<usize>().map_err(|e| {
        format!(
            "Invalid end coordinate at line {} ('{}'): {}",
            line_number, end_raw, e
        )
    })?;

    if start > end {
        return Err(format!(
            "Invalid gene interval at line {}: start ({}) is greater than end ({})",
            line_number, start, end
        )
        .into());
    }
    if start == 0 || end == 0 {
        return Err(format!(
            "Invalid gene interval at line {}: coordinates must be 1-based positive integers (start={}, end={})",
            line_number, start, end
        )
        .into());
    }

    Ok((start, end))
}

fn gene_name_from_gff(attrs: &HashMap<String, String>) -> String {
    let primary = attrs
        .get("gene")
        .or_else(|| attrs.get("Name"))
        .or_else(|| attrs.get("locus_tag"))
        .or_else(|| attrs.get("ID"))
        .map(|value| value.trim_start_matches("gene-").to_string())
        .unwrap_or_else(|| "unknown_gene".to_string());

    let suffix = attrs
        .get("locus_tag")
        .map(|value| value.to_string())
        .unwrap_or_else(|| primary.clone());

    format!("{}_{}", primary, suffix)
}

fn load_genes_from_tsv(genes_file: &str, snp_list: &[VcfPosition]) -> AppResult<Vec<Gene>> {
    log::info!("Loading TSV gene file: {}", genes_file);
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    let mut genes: Vec<Gene> = Vec::new();
    let mut parsed_entries = 0usize;
    let mut genes_with_snps = 0usize;
    let mut genes_without_snps = 0usize;

    for (line_idx, line) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let entry =
            line.map_err(|e| format!("Failed to read line {} in genes file: {}", line_number, e))?;
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
        "TSV gene entries parsed: {} | mapped to SNPs: {} | without SNPs: {}",
        parsed_entries,
        genes_with_snps,
        genes_without_snps
    );

    Ok(genes)
}

fn load_genes_from_gff(
    genes_file: &str,
    snp_list: &[VcfPosition],
    expected_chrom: Option<&str>,
) -> AppResult<Vec<Gene>> {
    log::info!("Loading GFF/GFF3 annotation file: {}", genes_file);
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    let mut genes: Vec<Gene> = Vec::new();
    let mut parsed_entries = 0usize;
    let mut genes_with_snps = 0usize;
    let mut genes_without_snps = 0usize;
    let mut genes_other_contigs = 0usize;
    let mut other_contigs: HashSet<String> = HashSet::new();

    for (line_idx, line) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let entry = line.map_err(|e| {
            format!(
                "Failed to read line {} in GFF/GFF3 annotation file: {}",
                line_number, e
            )
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

        if fields[2] != "gene" {
            continue;
        }

        if let Some(chrom) = expected_chrom {
            if fields[0] != chrom {
                genes_other_contigs += 1;
                other_contigs.insert(fields[0].to_string());
                continue;
            }
        }

        parsed_entries += 1;

        let (start, end) = parse_interval(fields[3], fields[4], line_number)?;
        let strand = parse_strand(fields[6], line_number)?;
        let attrs = parse_gff_attributes(fields[8]);
        let gene_name = gene_name_from_gff(&attrs);

        if has_snp_in_interval(snp_list, start, end) {
            genes_with_snps += 1;
            genes.push(Gene {
                name: gene_name,
                start,
                end,
                strand,
            });
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
            format!(" ({})", values.join(","))
        }
    );

    Ok(genes)
}

pub fn load_genes(
    genes_file: &str,
    snp_list: &[VcfPosition],
    expected_chrom: Option<&str>,
) -> AppResult<Vec<Gene>> {
    match detect_annotation_format(genes_file)? {
        AnnotationFormat::Tsv => load_genes_from_tsv(genes_file, snp_list),
        AnnotationFormat::Gff => load_genes_from_gff(genes_file, snp_list, expected_chrom),
    }
}

pub fn get_base_name(file_path: &str) -> AppResult<String> {
    let path = Path::new(file_path);
    let stem = path.file_stem().ok_or_else(|| {
        IoError::new(
            ErrorKind::InvalidInput,
            format!("Invalid input VCF path '{}': missing file stem", file_path),
        )
    })?;
    let stem_utf8 = stem.to_str().ok_or_else(|| {
        IoError::new(
            ErrorKind::InvalidInput,
            format!(
                "Invalid input VCF path '{}': file name is not valid UTF-8",
                file_path
            ),
        )
    })?;
    Ok(stem_utf8.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_temp_path(prefix: &str, suffix: &str) -> std::path::PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time before UNIX_EPOCH")
            .as_nanos();
        std::env::temp_dir().join(format!("{}_{}_{}", prefix, nanos, suffix))
    }

    #[test]
    fn test_parse_annotation_format_by_extension() {
        let gff_format = detect_annotation_format("/tmp/example_annotations.gff3")
            .expect("should detect GFF by extension");
        let tsv_format =
            detect_annotation_format("/tmp/example_annotations.tsv").expect("should detect TSV");
        assert_eq!(gff_format, AnnotationFormat::Gff);
        assert_eq!(tsv_format, AnnotationFormat::Tsv);
    }

    #[test]
    fn test_parse_annotation_format_by_content() {
        let path = unique_temp_path("get_mnv_annotations", "noext");
        fs::write(&path, "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=gene-abc\n")
            .expect("failed to write temp file");
        let format = detect_annotation_format(path.to_string_lossy().as_ref())
            .expect("should detect GFF by content");
        assert_eq!(format, AnnotationFormat::Gff);
        let _ = fs::remove_file(path);
    }

    #[test]
    fn test_gene_name_from_gff_attributes() {
        let mut attrs = HashMap::new();
        attrs.insert("gene".to_string(), "dnaN".to_string());
        attrs.insert("locus_tag".to_string(), "Rv0002".to_string());
        assert_eq!(gene_name_from_gff(&attrs), "dnaN_Rv0002");
    }

    #[test]
    fn test_load_genes_from_gff_filters_contigs() {
        let gff = "\
##gff-version 3
chrA\tsrc\tgene\t10\t30\t.\t+\t.\tID=gene-one;gene=one;locus_tag=L1
chrB\tsrc\tgene\t10\t30\t.\t+\t.\tID=gene-two;gene=two;locus_tag=L2
";
        let path = unique_temp_path("get_mnv_gff", "gff3");
        fs::write(&path, gff).expect("failed to write temp gff");

        let snp_list = vec![
            VcfPosition {
                position: 15,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                original_dp: None,
                original_freq: None,
            },
            VcfPosition {
                position: 20,
                ref_allele: "C".to_string(),
                alt_allele: "G".to_string(),
                original_dp: None,
                original_freq: None,
            },
        ];

        let genes = load_genes_from_gff(path.to_string_lossy().as_ref(), &snp_list, Some("chrA"))
            .expect("should load genes");
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "one_L1");
        assert_eq!(genes[0].start, 10);
        assert_eq!(genes[0].end, 30);
        let _ = fs::remove_file(path);
    }

    #[test]
    fn test_parse_gff_attributes_handles_encoded_and_gtf_styles() {
        let attrs = parse_gff_attributes(
            "ID=gene-abc;Name=dnaN%2Fbeta;locus_tag=Rv0002;gene_name \"dnaN\";note=\"ATPase;essential\"",
        );
        assert_eq!(attrs.get("ID"), Some(&"gene-abc".to_string()));
        assert_eq!(attrs.get("Name"), Some(&"dnaN/beta".to_string()));
        assert_eq!(attrs.get("locus_tag"), Some(&"Rv0002".to_string()));
        assert_eq!(attrs.get("gene_name"), Some(&"dnaN".to_string()));
        assert_eq!(attrs.get("note"), Some(&"ATPase;essential".to_string()));
    }

    #[test]
    fn test_load_references_multiple_contigs() {
        let fasta_content = ">chr1\nACTG\n>chr2\nTTAA\n";
        let path = unique_temp_path("get_mnv_fasta_multi", "fas");
        fs::write(&path, fasta_content).expect("failed to write temp fasta");
        let refs =
            load_references(path.to_string_lossy().as_ref()).expect("should load references");
        assert_eq!(refs.get("chr1"), Some(&"ACTG".to_string()));
        assert_eq!(refs.get("chr2"), Some(&"TTAA".to_string()));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn test_validate_vcf_reference_alleles_detects_mismatch() {
        let reference = Reference {
            sequence: "ACTG".to_string(),
        };
        let snp_list = vec![VcfPosition {
            position: 2,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
        }];
        let error = validate_vcf_reference_alleles("chr1", &snp_list, &reference)
            .expect_err("expected mismatch");
        assert!(error.to_string().contains("VCF REF/FASTA mismatch"));
    }

    #[test]
    fn test_validate_iupac_sequence_rejects_invalid_base() {
        let error = validate_iupac_sequence("ACXT", "test context")
            .expect_err("should reject invalid IUPAC base");
        assert!(error.to_string().contains("Invalid base"));
    }

    #[test]
    fn test_parse_optional_depth_handles_integer_and_float_strings() {
        assert_eq!(parse_optional_depth("22"), Some(22));
        assert_eq!(parse_optional_depth("22.0"), Some(22));
        assert_eq!(parse_optional_depth("."), None);
        assert_eq!(parse_optional_depth(""), None);
    }

    fn next_u64(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        *seed
    }

    fn random_char(seed: &mut u64) -> char {
        const ALPHABET: &[u8] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_=;:%\" ";
        let idx = (next_u64(seed) as usize) % ALPHABET.len();
        ALPHABET[idx] as char
    }

    #[test]
    fn test_parse_gff_attributes_property_randomized_input() {
        let mut seed = 0xC0FFEE_u64;
        for _ in 0..500 {
            let len = (next_u64(&mut seed) as usize % 120) + 1;
            let mut raw = String::with_capacity(len);
            for _ in 0..len {
                raw.push(random_char(&mut seed));
            }
            let parsed = parse_gff_attributes(&raw);
            for (key, value) in parsed {
                assert!(!key.trim().is_empty());
                assert!(!value.trim().is_empty());
            }
        }
    }

    #[test]
    fn test_load_vcf_positions_multiallelic_split_mode() {
        let path = unique_temp_path("get_mnv_vcf_multi_split", "vcf");
        fs::write(
            &path,
            "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t2\t.\tT\tC,G\t.\tPASS\tDP=10;AF=0.5\n",
        )
        .expect("failed to write temp vcf");

        let err = load_vcf_positions_by_contig(path.to_string_lossy().as_ref(), None, false, false)
            .expect_err("multiallelic should fail without split mode");
        assert!(err.to_string().contains("Multiallelic"));

        let parsed =
            load_vcf_positions_by_contig(path.to_string_lossy().as_ref(), None, true, false)
                .expect("split mode should parse");
        let positions = parsed.get("chr1").expect("missing chr1");
        assert_eq!(positions.len(), 2);
        assert_eq!(positions[0].position, 2);
        assert_eq!(positions[0].alt_allele, "C");
        assert_eq!(positions[1].alt_allele, "G");

        let _ = fs::remove_file(path);
    }
}
