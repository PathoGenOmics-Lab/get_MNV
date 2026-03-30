//! FASTA reference loading and validation.

use super::vcf::VcfPosition;
use crate::error::AppResult;
use bio::io::fasta;
use std::collections::HashMap;
use std::fs::File;

#[inline]
fn is_iupac_byte(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N')
}

pub struct Reference<'a> {
    pub sequence: &'a str,
}

pub type ReferenceMap = HashMap<String, String>;

pub fn load_references(fasta_file: &str) -> AppResult<ReferenceMap> {
    log::info!("Loading reference FASTA: {fasta_file}");
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
        // In-place uppercase: single allocation, no intermediate String
        let mut seq_bytes = record.seq().to_vec();
        seq_bytes.make_ascii_uppercase();
        // Validate IUPAC in-place on the byte slice (avoids second UTF-8 check)
        if let Some(pos) = seq_bytes.iter().position(|&b| !is_iupac_byte(b)) {
            return Err(format!(
                "Invalid base '{}' at position {} in FASTA contig '{}' in '{}'. Allowed IUPAC DNA bases: A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N",
                seq_bytes[pos] as char,
                pos + 1,
                record.id(),
                fasta_file
            ).into());
        }
        // SAFETY: seq_bytes contains only ASCII uppercase IUPAC chars after validation
        let sequence_upper = unsafe { String::from_utf8_unchecked(seq_bytes) };
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
        return Err(format!("No FASTA records found in '{fasta_file}'").into());
    }
    Ok(references)
}

pub fn reference_for_chrom<'a>(references: &'a ReferenceMap, chrom: &str) -> AppResult<Reference<'a>> {
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
    Ok(Reference {
        sequence,
    })
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
