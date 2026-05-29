//! Annotation format detection (GFF/GFF3 vs TSV).

use crate::error::AppResult;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

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
