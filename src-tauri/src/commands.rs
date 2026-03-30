//! Tauri commands for the get_MNV desktop application.

use std::collections::{BTreeSet, HashMap};
use std::io::{BufRead, BufReader};

use rust_htslib::bam::{self, Read as BamRead};
use serde::{Deserialize, Serialize};
use tauri::Emitter;

// ---------------------------------------------------------------------------
// get_core_version
// ---------------------------------------------------------------------------

#[tauri::command]
pub fn get_core_version() -> String {
    get_mnv::VERSION.to_string()
}

// ---------------------------------------------------------------------------
// run_analysis
// ---------------------------------------------------------------------------

/// Analysis configuration sent from the frontend (camelCase).
#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AnalysisConfig {
    pub vcf_file: String,
    pub bam_file: Option<String>,
    pub fasta_file: String,
    pub genes_file: String,
    pub sample: Option<String>,
    pub chrom: Option<String>,
    pub normalize_alleles: Option<bool>,
    pub min_quality: Option<u8>,
    pub min_mapq: Option<u8>,
    pub threads: Option<usize>,
    pub min_snp_reads: Option<usize>,
    pub min_mnv_reads: Option<usize>,
    pub min_snp_strand_reads: Option<usize>,
    pub min_mnv_strand_reads: Option<usize>,
    pub min_strand_bias_p: Option<f64>,
    pub dry_run: Option<bool>,
    pub strict: Option<bool>,
    pub split_multiallelic: Option<bool>,
    pub emit_filtered: Option<bool>,
    pub vcf_gz: Option<bool>,
    pub index_vcf_gz: Option<bool>,
    pub strand_bias_info: Option<bool>,
    pub keep_original_info: Option<bool>,
    pub bcf: Option<bool>,
    pub summary_json: Option<String>,
    pub error_json: Option<String>,
    pub run_manifest: Option<String>,
    pub gff_features: Option<Vec<String>>,
    pub exclude_intergenic: Option<bool>,
    pub translation_table: Option<u8>,
    pub output_tsv: Option<bool>,
    pub output_vcf: Option<bool>,
    pub convert: Option<bool>,
    pub both: Option<bool>,
    pub output_dir: Option<String>,
    pub output_prefix: Option<String>,
}

impl AnalysisConfig {
    /// Convert the GUI config into the core library's `Args` struct.
    fn into_args(self) -> get_mnv::cli::Args {
        // Route genes_file to the correct field based on file extension.
        let is_gff = self.genes_file.ends_with(".gff")
            || self.genes_file.ends_with(".gff3")
            || self.genes_file.ends_with(".GFF")
            || self.genes_file.ends_with(".GFF3");
        let (gff_file, genes_file_tsv) = if is_gff {
            (Some(self.genes_file), None)
        } else {
            (None, Some(self.genes_file))
        };

        get_mnv::cli::Args {
            vcf_file: self.vcf_file,
            bam_file: self.bam_file,
            fasta_file: self.fasta_file,
            genes_file_tsv,
            gff_file,
            gff_features_raw: self
                .gff_features
                .map(|v| v.join(",")),
            sample: self.sample,
            chrom: self.chrom,
            normalize_alleles: self.normalize_alleles.unwrap_or(false),
            min_quality: self.min_quality.unwrap_or(20),
            min_mapq: self.min_mapq.unwrap_or(0),
            threads: self.threads,
            min_snp_reads: self.min_snp_reads.unwrap_or(0),
            min_mnv_reads: self.min_mnv_reads.unwrap_or(0),
            min_snp_strand_reads: self.min_snp_strand_reads.unwrap_or(0),
            min_mnv_strand_reads: self.min_mnv_strand_reads.unwrap_or(0),
            min_strand_bias_p: self.min_strand_bias_p.unwrap_or(0.0),
            dry_run: self.dry_run.unwrap_or(false),
            strict: self.strict.unwrap_or(false),
            split_multiallelic: self.split_multiallelic.unwrap_or(false),
            emit_filtered: self.emit_filtered.unwrap_or(false),
            vcf_gz: self.vcf_gz.unwrap_or(false),
            index_vcf_gz: self.index_vcf_gz.unwrap_or(self.vcf_gz.unwrap_or(false)),
            strand_bias_info: self.strand_bias_info.unwrap_or(false),
            keep_original_info: self.keep_original_info.unwrap_or(false),
            bcf: self.bcf.unwrap_or(false),
            summary_json: self.summary_json,
            error_json: self.error_json,
            run_manifest: self.run_manifest,
            exclude_intergenic: self.exclude_intergenic.unwrap_or(false),
            translation_table: self.translation_table.unwrap_or(11),
            // Map GUI's outputTsv/outputVcf booleans to CLI's convert/both flags.
            // Default: TSV only (convert=false, both=false).
            // outputTsv + outputVcf → both=true
            // outputVcf only → convert=true
            convert: {
                let tsv = self.output_tsv.unwrap_or(true);
                let vcf = self.output_vcf.unwrap_or(false);
                !tsv && vcf
            },
            both: {
                let tsv = self.output_tsv.unwrap_or(true);
                let vcf = self.output_vcf.unwrap_or(false);
                tsv && vcf
            },
            output_dir: self.output_dir,
            output_prefix: self.output_prefix,
        }
    }
}

#[tauri::command]
pub fn run_analysis(
    app_handle: tauri::AppHandle,
    config: AnalysisConfig,
) -> Result<serde_json::Value, String> {
    let args = config.into_args();

    let progress_callback = move |evt: get_mnv::pipeline::ProgressEvent| {
        let _ = app_handle.emit("analysis-progress", &evt);
    };

    let summary = get_mnv::pipeline::run_with_progress(&args, &progress_callback)
        .map_err(|e| format!("{}", e))?;

    serde_json::to_value(&summary).map_err(|e| format!("Failed to serialize summary: {}", e))
}

// ---------------------------------------------------------------------------
// get_gff_features
// ---------------------------------------------------------------------------

#[tauri::command]
pub fn get_gff_features(path: String) -> Result<Vec<String>, String> {
    let file = std::fs::File::open(&path).map_err(|e| format!("Cannot open GFF file: {}", e))?;
    let reader = BufReader::new(file);
    let mut features = BTreeSet::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() >= 3 {
            features.insert(cols[2].to_string());
        }
    }

    Ok(features.into_iter().collect())
}

// ---------------------------------------------------------------------------
// read_tsv_file
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct TsvData {
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

#[tauri::command]
pub fn read_tsv_file(path: String) -> Result<TsvData, String> {
    let file = std::fs::File::open(&path).map_err(|e| format!("Cannot open TSV file: {}", e))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header_line = lines
        .next()
        .ok_or_else(|| "TSV file is empty".to_string())?
        .map_err(|e| format!("Read error: {}", e))?;

    let headers: Vec<String> = header_line.split('\t').map(|s| s.to_string()).collect();
    let mut rows = Vec::new();

    for line in lines {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        if line.trim().is_empty() {
            continue;
        }
        rows.push(line.split('\t').map(|s| s.to_string()).collect());
    }

    Ok(TsvData { headers, rows })
}

// ---------------------------------------------------------------------------
// get_bam_view
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BamViewRequest {
    pub bam_path: String,
    pub fasta_path: String,
    pub chrom: String,
    pub positions: Vec<u64>,
    pub ref_bases: Vec<String>,
    pub alt_bases: Vec<String>,
    pub min_mapq: Option<u8>,
    pub min_base_quality: Option<u8>,
    pub max_reads: Option<usize>,
    pub window_start: u64,
    pub window_end: u64,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BamViewResponse {
    pub chrom: String,
    pub display_start: u64,
    pub display_end: u64,
    pub reference: String,
    pub sites: Vec<BamVariantSite>,
    pub reads: Vec<BamReadView>,
    pub counts: BamSupportCounts,
    pub total_reads: usize,
    pub truncated: bool,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BamVariantSite {
    pub position: u64,
    pub reference_base: String,
    pub alt_base: String,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BamReadView {
    pub name: String,
    pub strand: String,
    pub support: String,
    pub start: u64,
    pub end: u64,
    pub mapq: u8,
    pub bases: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
#[serde(rename_all = "camelCase")]
pub struct BamSupportCounts {
    pub total: usize,
    pub mnv: usize,
    pub partial: usize,
    pub reference: usize,
    pub other: usize,
}

/// Extract the base a read contributes at a given reference position, using
/// the CIGAR alignment. Returns `None` if the read does not cover that
/// position or the base quality is below threshold.
fn base_at_ref_pos(
    record: &bam::Record,
    ref_pos: u64,
    min_bq: u8,
) -> Option<String> {
    let read_start = record.pos() as u64;
    let seq = record.seq();
    let quals = record.qual();
    let cigar = record.cigar();

    let mut r_pos = read_start; // current reference position
    let mut q_pos: usize = 0; // current query (read) position

    for op in cigar.iter() {
        let len = op.len() as u64;
        match op.char() as u8 {
            b'M' | b'=' | b'X' => {
                // Alignment match / sequence match / mismatch
                if ref_pos >= r_pos && ref_pos < r_pos + len {
                    let offset = (ref_pos - r_pos) as usize;
                    let qi = q_pos + offset;
                    if qi < seq.len() && quals[qi] >= min_bq {
                        let base = seq[qi] as char;
                        return Some(base.to_string());
                    } else {
                        return None;
                    }
                }
                r_pos += len;
                q_pos += len as usize;
            }
            b'I' | b'S' => {
                // Insertion / soft clip: consumes query only
                q_pos += len as usize;
            }
            b'D' | b'N' => {
                // Deletion / reference skip: consumes reference only
                if ref_pos >= r_pos && ref_pos < r_pos + len {
                    return Some("-".to_string()); // deletion at this position
                }
                r_pos += len;
            }
            b'H' | b'P' => {
                // Hard clip / padding: consumes nothing
            }
            _ => {
                // Unknown CIGAR op — skip
            }
        }
    }
    None
}

#[tauri::command]
pub fn get_bam_view(request: BamViewRequest) -> Result<BamViewResponse, String> {
    let min_mapq = request.min_mapq.unwrap_or(0);
    let min_bq = request.min_base_quality.unwrap_or(0);
    let max_reads = request.max_reads.unwrap_or(500);

    // Build variant site lookup: position → (ref_base, alt_base)
    let mut site_map: HashMap<u64, (String, String)> = HashMap::new();
    let mut sites: Vec<BamVariantSite> = Vec::new();
    for i in 0..request.positions.len() {
        let pos = request.positions[i];
        let rb = request.ref_bases[i].clone();
        let ab = request.alt_bases[i].clone();
        site_map.insert(pos, (rb.clone(), ab.clone()));
        sites.push(BamVariantSite {
            position: pos,
            reference_base: rb,
            alt_base: ab,
        });
    }

    // Read reference sequence for the window
    let fasta_reader = bio::io::fasta::IndexedReader::from_file(&request.fasta_path)
        .map_err(|e| format!("Cannot open FASTA: {}", e))?;
    let mut fasta = fasta_reader;
    let window_len = (request.window_end - request.window_start) as usize;
    let mut ref_seq = vec![0u8; window_len];
    fasta
        .fetch(&request.chrom, request.window_start, request.window_end)
        .map_err(|e| format!("FASTA fetch error: {}", e))?;
    fasta
        .read(&mut ref_seq)
        .map_err(|e| format!("FASTA read error: {}", e))?;
    let reference: String = ref_seq.iter().map(|&b| b as char).collect();

    // Open BAM
    let mut bam_reader = bam::IndexedReader::from_path(&request.bam_path)
        .map_err(|e| format!("Cannot open BAM: {}", e))?;

    // Get target id for the chromosome
    let tid = {
        let header = bam_reader.header().clone();
        let tid = header
            .tid(request.chrom.as_bytes())
            .ok_or_else(|| format!("Chromosome '{}' not found in BAM header", request.chrom))?;
        tid
    };

    bam_reader
        .fetch(bam::FetchDefinition::Region(
            tid as i32,
            request.window_start as i64,
            request.window_end as i64,
        ))
        .map_err(|e| format!("BAM fetch error: {}", e))?;

    let mut all_reads: Vec<BamReadView> = Vec::new();
    let mut counts = BamSupportCounts::default();
    let mut total_reads: usize = 0;

    for result in bam_reader.records() {
        let record = result.map_err(|e| format!("BAM record error: {}", e))?;

        // Skip unmapped, secondary, supplementary, QC-fail, duplicate
        if record.is_unmapped()
            || record.is_secondary()
            || record.is_supplementary()
            || record.is_quality_check_failed()
            || record.is_duplicate()
        {
            continue;
        }

        if record.mapq() < min_mapq {
            continue;
        }

        total_reads += 1;

        // Extract bases at each window position
        let mut bases: Vec<String> = Vec::with_capacity(window_len);
        for offset in 0..window_len {
            let rp = request.window_start + offset as u64;
            match base_at_ref_pos(&record, rp, min_bq) {
                Some(b) => bases.push(b),
                None => bases.push(" ".to_string()), // no coverage
            }
        }

        // Classify read support at variant sites
        let mut matches_alt = 0usize;
        let mut matches_ref = 0usize;
        let mut covered_sites = 0usize;

        for &pos in &request.positions {
            if pos < request.window_start || pos >= request.window_end {
                continue;
            }
            let offset = (pos - request.window_start) as usize;
            if offset < bases.len() {
                let base = &bases[offset];
                if base == " " || base == "-" {
                    continue;
                }
                covered_sites += 1;
                if let Some((rb, ab)) = site_map.get(&pos) {
                    if base.eq_ignore_ascii_case(ab) {
                        matches_alt += 1;
                    } else if base.eq_ignore_ascii_case(rb) {
                        matches_ref += 1;
                    }
                }
            }
        }

        let n_sites = request.positions.len();
        let support = if covered_sites == 0 {
            "other".to_string()
        } else if matches_alt == n_sites {
            "mnv".to_string()
        } else if matches_alt > 0 && matches_alt < n_sites {
            "partial".to_string()
        } else if matches_ref == n_sites {
            "reference".to_string()
        } else {
            "other".to_string()
        };

        // Update counts
        match support.as_str() {
            "mnv" => counts.mnv += 1,
            "partial" => counts.partial += 1,
            "reference" => counts.reference += 1,
            _ => counts.other += 1,
        }

        let strand = if record.is_reverse() {
            "-".to_string()
        } else {
            "+".to_string()
        };

        let read_name = std::str::from_utf8(record.qname())
            .unwrap_or("?")
            .to_string();

        all_reads.push(BamReadView {
            name: read_name,
            strand,
            support,
            start: record.pos() as u64,
            end: record.cigar().end_pos() as u64,
            mapq: record.mapq(),
            bases,
        });
    }

    counts.total = total_reads;

    // Sort: MNV first, then partial, reference, other
    all_reads.sort_by(|a, b| {
        let order = |s: &str| -> u8 {
            match s {
                "mnv" => 0,
                "partial" => 1,
                "reference" => 2,
                _ => 3,
            }
        };
        order(&a.support).cmp(&order(&b.support))
    });

    let truncated = all_reads.len() > max_reads;
    all_reads.truncate(max_reads);

    Ok(BamViewResponse {
        chrom: request.chrom,
        display_start: request.window_start,
        display_end: request.window_end,
        reference,
        sites,
        reads: all_reads,
        counts,
        total_reads,
        truncated,
    })
}

// ---------------------------------------------------------------------------
// write_text_file
// ---------------------------------------------------------------------------

#[tauri::command]
pub fn write_text_file(path: String, content: String) -> Result<(), String> {
    std::fs::write(&path, &content).map_err(|e| format!("Failed to write file '{}': {}", path, e))
}
