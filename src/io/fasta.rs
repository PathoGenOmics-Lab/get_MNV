//! FASTA reference loading and validation.

use super::validation::validate_iupac_sequence;
use super::vcf::VcfPosition;
use crate::error::AppResult;
use bio::io::fasta;
use std::collections::HashMap;
use std::fs::File;

pub struct Reference {
    pub sequence: String,
}

pub type ReferenceMap = HashMap<String, String>;

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
                names.join(", ")
            }
        )
    })?;
    Ok(Reference {
        sequence: sequence.clone(),
    })
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
