//! FASTA reference loading and validation.
//!
//! Uses a hand-rolled parser instead of bio-rs to eliminate intermediate
//! allocations: sequences are built directly into a `Vec<u8>` that is
//! uppercased and IUPAC-validated in-place, then converted to `String`
//! without a second copy.

use super::vcf::VcfPosition;
use crate::error::AppResult;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

/// Lookup table: IUPAC DNA bytes → true. 256-entry table for O(1) check.
const fn build_iupac_table() -> [bool; 256] {
    let mut table = [false; 256];
    let chars = b"ACGTRYSSWKMBDHVN";
    let mut i = 0;
    while i < chars.len() {
        table[chars[i] as usize] = true;
        i += 1;
    }
    table
}
static IUPAC_TABLE: [bool; 256] = build_iupac_table();

#[inline]
fn is_iupac_byte(b: u8) -> bool {
    IUPAC_TABLE[b as usize]
}

pub struct Reference<'a> {
    pub sequence: &'a str,
}

pub type ReferenceMap = HashMap<String, String>;

pub fn load_references(fasta_file: &str) -> AppResult<ReferenceMap> {
    log::info!("Loading reference FASTA: {fasta_file}");
    let file = std::fs::File::open(fasta_file)
        .map_err(|e| format!("Cannot open FASTA file '{}': {}", fasta_file, e))?;
    let reader = BufReader::with_capacity(128 * 1024, file);
    let mut references: ReferenceMap = HashMap::new();
    let mut current_id: Option<String> = None;
    let mut seq_buf: Vec<u8> = Vec::with_capacity(4_500_000); // pre-size for typical genome
    let mut record_idx = 0usize;

    for line_result in reader.lines() {
        let line =
            line_result.map_err(|e| format!("Failed reading FASTA '{}': {}", fasta_file, e))?;
        if line.starts_with('>') {
            // Flush previous record
            if let Some(id) = current_id.take() {
                let seq_str = finalize_sequence(&mut seq_buf, &id, fasta_file)?;
                if references.insert(id.clone(), seq_str).is_some() {
                    return Err(format!(
                        "Duplicate FASTA contig '{}' found in '{}'",
                        id, fasta_file
                    )
                    .into());
                }
            }
            record_idx += 1;
            // Parse header: >id description...
            let header = line.strip_prefix('>').unwrap_or(&line).trim();
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            if id.is_empty() {
                return Err(format!(
                    "Invalid FASTA record {} in '{}': empty contig name",
                    record_idx, fasta_file
                )
                .into());
            }
            current_id = Some(id);
            seq_buf.clear();
        } else if current_id.is_some() {
            // Append sequence bytes, skipping whitespace
            for &b in line.as_bytes() {
                if !b.is_ascii_whitespace() {
                    seq_buf.push(b);
                }
            }
        }
    }

    // Flush last record
    if let Some(id) = current_id.take() {
        let seq_str = finalize_sequence(&mut seq_buf, &id, fasta_file)?;
        if references.insert(id.clone(), seq_str).is_some() {
            return Err(
                format!("Duplicate FASTA contig '{}' found in '{}'", id, fasta_file).into(),
            );
        }
    }

    if references.is_empty() {
        return Err(format!("No FASTA records found in '{fasta_file}'").into());
    }
    Ok(references)
}

/// Uppercase, validate, and convert sequence buffer to String in-place.
fn finalize_sequence(buf: &mut Vec<u8>, contig_id: &str, fasta_file: &str) -> AppResult<String> {
    buf.make_ascii_uppercase();
    if let Some(pos) = buf.iter().position(|&b| !is_iupac_byte(b)) {
        return Err(format!(
            "Invalid base '{}' at position {} in FASTA contig '{}' in '{}'. \
Allowed IUPAC DNA bases: A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N",
            buf[pos] as char,
            pos + 1,
            contig_id,
            fasta_file
        )
        .into());
    }
    // SAFETY: buf contains only ASCII uppercase IUPAC chars after validation
    Ok(unsafe { String::from_utf8_unchecked(std::mem::take(buf)) })
}

pub fn reference_for_chrom<'a>(
    references: &'a ReferenceMap,
    chrom: &str,
) -> AppResult<Reference<'a>> {
    let sequence = references.get(chrom).ok_or_else(|| {
        format!(
            "Reference contig '{}' is missing in FASTA. Available contigs: {}",
            chrom,
            {
                let mut names = references.keys().cloned().collect::<Vec<_>>();
                names.sort();
                names.join(", ")
            }
        )
    })?;
    Ok(Reference { sequence })
}

pub fn validate_vcf_reference_alleles(
    chrom: &str,
    snp_list: &[VcfPosition],
    reference: &Reference<'_>,
) -> AppResult<()> {
    for site in snp_list {
        if site.position == 0 {
            return Err(format!(
                "Invalid VCF position 0 found at contig '{chrom}' (positions must be 1-based)"
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);

    fn write_temp_fasta(content: &str) -> String {
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!("test_fasta_{}_{n}.fa", std::process::id()));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path.to_string_lossy().to_string()
    }

    #[test]
    fn test_load_single_contig() {
        let path = write_temp_fasta(">chr1\nACGT\nTTGG\n");
        let refs = load_references(&path).unwrap();
        assert_eq!(refs.len(), 1);
        assert_eq!(refs["chr1"], "ACGTTTGG");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_load_multi_contig() {
        let path = write_temp_fasta(">chr1\nACGT\n>chr2\nNNNN\n");
        let refs = load_references(&path).unwrap();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs["chr1"], "ACGT");
        assert_eq!(refs["chr2"], "NNNN");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_load_lowercase_uppercased() {
        let path = write_temp_fasta(">c\nacgt\n");
        let refs = load_references(&path).unwrap();
        assert_eq!(refs["c"], "ACGT");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_load_iupac_ambiguity() {
        let path = write_temp_fasta(">c\nACGTRYSWKMBDHVN\n");
        let refs = load_references(&path).unwrap();
        assert_eq!(refs["c"], "ACGTRYSWKMBDHVN");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_load_invalid_base_rejected() {
        let path = write_temp_fasta(">c\nACGTX\n");
        assert!(load_references(&path).is_err());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_load_empty_file_rejected() {
        let path = write_temp_fasta("");
        assert!(load_references(&path).is_err());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_load_duplicate_contig_rejected() {
        let path = write_temp_fasta(">c\nACGT\n>c\nTTTT\n");
        assert!(load_references(&path).is_err());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_load_header_description_stripped() {
        let path = write_temp_fasta(">chr1 some description here\nACGT\n");
        let refs = load_references(&path).unwrap();
        assert!(refs.contains_key("chr1"));
        assert!(!refs.contains_key("chr1 some description here"));
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_reference_for_chrom_found() {
        let mut refs = ReferenceMap::new();
        refs.insert("chr1".to_string(), "ACGT".to_string());
        let r = reference_for_chrom(&refs, "chr1").unwrap();
        assert_eq!(r.sequence, "ACGT");
    }

    #[test]
    fn test_reference_for_chrom_missing() {
        let refs = ReferenceMap::new();
        assert!(reference_for_chrom(&refs, "chr1").is_err());
    }

    #[test]
    fn test_validate_vcf_ref_alleles_ok() {
        let refs = ReferenceMap::from([("c".to_string(), "ACGT".to_string())]);
        let r = reference_for_chrom(&refs, "c").unwrap();
        let snps = vec![VcfPosition {
            position: 1,
            ref_allele: "A".into(),
            alt_allele: "T".into(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }];
        assert!(validate_vcf_reference_alleles("c", &snps, &r).is_ok());
    }

    #[test]
    fn test_validate_vcf_ref_mismatch() {
        let refs = ReferenceMap::from([("c".to_string(), "ACGT".to_string())]);
        let r = reference_for_chrom(&refs, "c").unwrap();
        let snps = vec![VcfPosition {
            position: 1,
            ref_allele: "T".into(),
            alt_allele: "A".into(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }];
        assert!(validate_vcf_reference_alleles("c", &snps, &r).is_err());
    }

    #[test]
    fn test_validate_vcf_ref_out_of_bounds() {
        let refs = ReferenceMap::from([("c".to_string(), "AC".to_string())]);
        let r = reference_for_chrom(&refs, "c").unwrap();
        let snps = vec![VcfPosition {
            position: 2,
            ref_allele: "CG".into(),
            alt_allele: "TT".into(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }];
        assert!(validate_vcf_reference_alleles("c", &snps, &r).is_err());
    }
}
