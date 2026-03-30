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

#[cfg(test)]
mod tests {
    use super::*;

    // ---- validate_iupac_sequence ----

    #[test]
    fn test_iupac_valid() {
        assert!(validate_iupac_sequence("ACGTRYWSKMBDHVN", "test").is_ok());
    }

    #[test]
    fn test_iupac_lowercase_valid() {
        assert!(validate_iupac_sequence("acgt", "test").is_ok());
    }

    #[test]
    fn test_iupac_invalid_base() {
        let err = validate_iupac_sequence("ACGX", "test");
        assert!(err.is_err());
        assert!(err.unwrap_err().to_string().contains("Invalid base 'X'"));
    }

    #[test]
    fn test_iupac_empty_rejected() {
        assert!(validate_iupac_sequence("", "test").is_err());
    }

    // ---- validate_vcf_allele ----

    #[test]
    fn test_vcf_allele_valid_snp() {
        assert!(validate_vcf_allele("A", 1, "chr1", 100, "ALT").is_ok());
    }

    #[test]
    fn test_vcf_allele_symbolic_alt_ok() {
        assert!(validate_vcf_allele("<DEL>", 1, "chr1", 100, "ALT").is_ok());
    }

    #[test]
    fn test_vcf_allele_symbolic_ref_rejected() {
        assert!(validate_vcf_allele("<DEL>", 1, "chr1", 100, "REF").is_err());
    }

    #[test]
    fn test_vcf_allele_star_rejected() {
        assert!(validate_vcf_allele("*", 1, "chr1", 100, "ALT").is_err());
    }

    // ---- get_base_name ----

    #[test]
    fn test_base_name_vcf() {
        assert_eq!(get_base_name("sample.vcf").unwrap(), "sample");
    }

    #[test]
    fn test_base_name_vcf_gz() {
        assert_eq!(get_base_name("/data/run.vcf.gz").unwrap(), "run");
    }

    #[test]
    fn test_base_name_plain_gz() {
        assert_eq!(get_base_name("data.gz").unwrap(), "data");
    }

    #[test]
    fn test_base_name_no_extension() {
        assert_eq!(get_base_name("noext").unwrap(), "noext");
    }

    #[test]
    fn test_base_name_other_extension() {
        assert_eq!(get_base_name("file.bam").unwrap(), "file");
    }
}
