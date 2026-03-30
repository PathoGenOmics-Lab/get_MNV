//! Input validation: IUPAC sequence checks, VCF allele validation, and file
//! path utilities.

use crate::error::AppResult;
use std::io::{Error as IoError, ErrorKind};
use std::path::Path;

fn is_iupac_dna_base(base: char) -> bool {
    matches!(
        base,
        'A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N'
    )
}

pub(crate) fn validate_iupac_sequence(value: &str, context: &str) -> AppResult<()> {
    if value.is_empty() {
        return Err(format!("Empty sequence found for {context}").into());
    }
    if let Some(invalid_base) = value
        .chars()
        .map(|base| base.to_ascii_uppercase())
        .find(|base| !is_iupac_dna_base(*base))
    {
        return Err(format!(
            "Invalid base '{invalid_base}' in {context}. Allowed IUPAC DNA bases: A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N"
        )
        .into());
    }
    Ok(())
}

pub(crate) fn validate_vcf_allele(
    allele: &str,
    record_idx: usize,
    chrom: &str,
    pos: usize,
    allele_kind: &str,
) -> AppResult<()> {
    if allele == "*" {
        return Err(format!(
            "Unsupported '{allele_kind}' allele '*' at VCF record {record_idx} ({chrom}:{pos}). Use normalized VCF without spanning-deletion alleles."
        )
        .into());
    }

    let is_symbolic = allele.starts_with('<') && allele.ends_with('>');
    if is_symbolic {
        if allele_kind == "REF" {
            return Err(format!(
                "Unsupported symbolic REF allele '{allele}' at VCF record {record_idx} ({chrom}:{pos})"
            )
            .into());
        }
        return Ok(());
    }

    validate_iupac_sequence(
        allele,
        &format!("{allele_kind} allele at VCF record {record_idx} ({chrom}:{pos})"),
    )
}

/// Extract the file stem from a path, stripping `.vcf`, `.vcf.gz`, and `.gz`
/// suffixes so that output files get clean names (e.g. `sample.MNV.tsv`
/// instead of `sample.vcf.MNV.tsv`).
///
/// # Examples
/// ```
/// use get_mnv::io::get_base_name;
/// assert_eq!(get_base_name("sample.vcf").unwrap(), "sample");
/// assert_eq!(get_base_name("/data/run.vcf.gz").unwrap(), "run");
/// assert_eq!(get_base_name("plain.gz").unwrap(), "plain");
/// assert_eq!(get_base_name("noext").unwrap(), "noext");
/// ```
pub fn get_base_name(file_path: &str) -> AppResult<String> {
    let path = Path::new(file_path);
    let file_name = path.file_name().ok_or_else(|| {
        IoError::new(
            ErrorKind::InvalidInput,
            format!("Invalid input VCF path '{file_path}': missing file name"),
        )
    })?;
    let name = file_name.to_str().ok_or_else(|| {
        IoError::new(
            ErrorKind::InvalidInput,
            format!("Invalid input VCF path '{file_path}': file name is not valid UTF-8"),
        )
    })?;

    // Strip compound extensions: .vcf.gz, .vcf, .gz
    let stem = if let Some(base) = name.strip_suffix(".vcf.gz") {
        base
    } else if let Some(base) = name.strip_suffix(".vcf") {
        base
    } else if let Some(base) = name.strip_suffix(".gz") {
        base
    } else {
        // Fallback: use Path::file_stem for any other extension
        path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(name)
    };

    if stem.is_empty() {
        return Err(IoError::new(
            ErrorKind::InvalidInput,
            format!("Invalid input VCF path '{file_path}': empty file stem after stripping extensions"),
        ).into());
    }

    Ok(stem.to_string())
}
