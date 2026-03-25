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

pub(crate) fn validate_vcf_allele(
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

/// Extract the file stem from a path, stripping `.vcf` / `.vcf.gz` suffixes.
///
/// # Examples
/// ```
/// use get_mnv::io::get_base_name;
/// assert_eq!(get_base_name("sample.vcf").unwrap(), "sample");
/// assert_eq!(get_base_name("/data/run.vcf.gz").unwrap(), "run.vcf");
/// ```
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
