//! Reproducibility: SHA-256 checksums, input metadata, and run manifest.

use super::summary::{InputChecksums, RunInputs, RunSummary};
use crate::cli::Args;
use crate::error::{AppResult, IoResultExt};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{BufReader, Read as IoRead, Write};
use std::time::{SystemTime, UNIX_EPOCH};

pub(crate) fn compute_sha256(path: &str) -> AppResult<String> {
    let file = File::open(path).with_path(path)?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 8192];
    loop {
        let read_bytes = reader.read(&mut buffer).with_path(path)?;
        if read_bytes == 0 {
            break;
        }
        hasher.update(&buffer[..read_bytes]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

pub(crate) fn build_input_metadata(args: &Args) -> AppResult<RunInputs> {
    Ok(RunInputs {
        vcf: args.vcf_file.clone(),
        fasta: args.fasta_file.clone(),
        annotation: args.genes_file.clone(),
        bam: args.bam_file.clone(),
        checksums: InputChecksums {
            vcf_sha256: compute_sha256(&args.vcf_file)?,
            fasta_sha256: compute_sha256(&args.fasta_file)?,
            annotation_sha256: compute_sha256(&args.genes_file)?,
            bam_sha256: match args.bam_file.as_ref() {
                Some(path) => Some(compute_sha256(path)?),
                None => None,
            },
        },
    })
}

pub(crate) fn build_run_manifest_value(
    summary: &RunSummary,
    command_line: &str,
) -> AppResult<Value> {
    let timestamp_unix = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs();

    let output_checksums = json!({
        "output_tsv_sha256": summary.output_tsv.as_deref().map(compute_sha256).transpose()?,
        "output_vcf_sha256": summary.output_vcf.as_deref().map(compute_sha256).transpose()?,
        "output_bcf_sha256": summary.output_bcf.as_deref().map(compute_sha256).transpose()?,
    });

    Ok(json!({
        "schema_version": "1.0.0",
        "timestamp_unix": timestamp_unix,
        "tool_version": env!("CARGO_PKG_VERSION"),
        "command_line": command_line,
        "summary": summary,
        "output_checksums": output_checksums
    }))
}

pub(crate) fn write_json_value(path: &str, value: &Value) -> AppResult<()> {
    let mut file = File::create(path).with_path(path)?;
    serde_json::to_writer_pretty(&mut file, value)?;
    file.write_all(b"\n").with_path(path)?;
    Ok(())
}
