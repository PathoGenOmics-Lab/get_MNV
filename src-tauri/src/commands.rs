//! Tauri commands for the get_MNV desktop application.

use std::collections::{BTreeSet, HashMap};
use std::io::{BufRead, BufReader};

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
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
#[allow(dead_code)] // `convert` and `both` are computed in into_args(), not read directly
pub struct AnalysisConfig {
    pub vcf_file: String,
    pub input_format: Option<get_mnv::cli::VariantInputFormat>,
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
    pub min_snp_frequency: Option<f64>,
    pub min_mnv_reads: Option<usize>,
    pub min_mnv_frequency: Option<f64>,
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
        let input_format = self
            .input_format
            .unwrap_or(get_mnv::cli::VariantInputFormat::Auto);
        let variant_file = self.vcf_file;
        let output_dir = self
            .output_dir
            .filter(|dir| !dir.trim().is_empty())
            .or_else(|| default_output_dir_for_variant_file(&variant_file));

        // Route genes_file to the correct field based on file extension.
        let genes_lower = self.genes_file.to_lowercase();
        let is_gff = genes_lower.ends_with(".gff")
            || genes_lower.ends_with(".gff3")
            || genes_lower.ends_with(".gtf")
            || genes_lower.ends_with(".gff.gz")
            || genes_lower.ends_with(".gff3.gz")
            || genes_lower.ends_with(".gtf.gz");
        let (gff_file, genes_file_tsv) = if is_gff {
            (Some(self.genes_file), None)
        } else {
            (None, Some(self.genes_file))
        };
        let (vcf_file, tsv_file) = if input_format == get_mnv::cli::VariantInputFormat::Tsv {
            (None, Some(variant_file))
        } else {
            (Some(variant_file), None)
        };

        get_mnv::cli::Args {
            vcf_file,
            tsv_file,
            input_format,
            bam_file: self.bam_file,
            fasta_file: self.fasta_file,
            genes_file_tsv,
            gff_file,
            gff_features_raw: self.gff_features.map(|v| v.join(",")),
            sample: self.sample,
            chrom: self.chrom,
            normalize_alleles: self.normalize_alleles.unwrap_or(false),
            min_quality: self.min_quality.unwrap_or(20),
            min_mapq: self.min_mapq.unwrap_or(0),
            threads: self.threads,
            min_snp_reads: self.min_snp_reads.unwrap_or(0),
            min_snp_frequency: self.min_snp_frequency.unwrap_or(0.0),
            min_mnv_reads: self.min_mnv_reads.unwrap_or(0),
            min_mnv_frequency: self.min_mnv_frequency.unwrap_or(0.0),
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
            output_dir,
            output_prefix: self.output_prefix,
        }
    }
}

fn default_output_dir_for_variant_file(variant_file: &str) -> Option<String> {
    let parent = std::path::Path::new(variant_file).parent()?;
    if parent.as_os_str().is_empty() {
        None
    } else {
        Some(parent.to_string_lossy().into_owned())
    }
}

fn expected_output_paths(config: &AnalysisConfig) -> Result<Vec<String>, String> {
    let args = config.clone().into_args();
    if args.dry_run {
        return Ok(Vec::new());
    }

    let base_name = get_mnv::io::get_base_name(args.variant_file()).map_err(|e| e.to_string())?;
    let stem_name = args.output_prefix.clone().unwrap_or(base_name);
    let output_stem = match &args.output_dir {
        Some(dir) => std::path::Path::new(dir)
            .join(&stem_name)
            .to_string_lossy()
            .into_owned(),
        None => stem_name,
    };

    let mut paths = Vec::new();
    if args.both || !args.convert {
        paths.push(format!("{output_stem}.MNV.tsv"));
    }
    if args.both || args.convert {
        paths.push(if args.vcf_gz {
            format!("{output_stem}.MNV.vcf.gz")
        } else {
            format!("{output_stem}.MNV.vcf")
        });
    }
    if args.bcf {
        paths.push(format!("{output_stem}.MNV.bcf"));
    }

    Ok(paths)
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
// detect_variant_input_format
// ---------------------------------------------------------------------------

#[tauri::command]
pub fn detect_variant_input_format(path: String) -> Result<String, String> {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") || lower.ends_with(".bcf") {
        return Ok("vcf".to_string());
    }

    match get_mnv::io::ivar::looks_like_ivar_tsv(&path) {
        Ok(true) => Ok("tsv".to_string()),
        Ok(false) => Ok("unknown".to_string()),
        Err(e) => Err(e.to_string()),
    }
}

// ---------------------------------------------------------------------------
// check_output_conflicts
// ---------------------------------------------------------------------------

#[tauri::command]
pub fn resolve_output_paths(config: AnalysisConfig) -> Result<Vec<String>, String> {
    expected_output_paths(&config)
}

/// Check if any output files already exist. Returns paths that would be overwritten.
#[tauri::command]
pub fn check_output_conflicts(config: AnalysisConfig) -> Result<Vec<String>, String> {
    Ok(expected_output_paths(&config)?
        .into_iter()
        .filter(|p| std::path::Path::new(p).exists())
        .collect())
}

// ---------------------------------------------------------------------------
// ensure_fasta_index
// ---------------------------------------------------------------------------

/// Create a .fai index for a FASTA file if it doesn't already exist.
/// Parses the FASTA to compute name, length, byte offset, line bases, line width.
/// Returns the path to the .fai file.
#[tauri::command]
pub fn ensure_fasta_index(fasta_path: String) -> Result<String, String> {
    let fai_path = format!("{}.fai", fasta_path);
    if std::path::Path::new(&fai_path).exists() {
        return Ok(fai_path);
    }

    // Cannot index gzipped FASTA — byte offsets would be meaningless
    if fasta_path.ends_with(".gz") || fasta_path.ends_with(".bgz") {
        return Err(
            "Cannot create .fai index for gzipped FASTA. Please decompress first.".to_string(),
        );
    }

    let bytes = std::fs::read(&fasta_path).map_err(|e| format!("Cannot read FASTA: {}", e))?;

    let mut entries: Vec<String> = Vec::new();
    let mut name = String::new();
    let mut seq_len: usize = 0;
    let mut seq_offset: usize = 0;
    let mut line_bases: usize = 0;
    let mut line_width: usize = 0;
    let mut first_seq_line = true;

    let mut pos: usize = 0;
    while pos < bytes.len() {
        // Find end of line
        let line_start = pos;
        while pos < bytes.len() && bytes[pos] != b'\n' {
            pos += 1;
        }
        let line_end = pos; // exclusive, before \n
        let has_newline = pos < bytes.len();
        if has_newline {
            pos += 1; // skip \n
        }

        if line_start < bytes.len() && bytes[line_start] == b'>' {
            // Flush previous sequence
            if !name.is_empty() {
                entries.push(format!(
                    "{}\t{}\t{}\t{}\t{}",
                    name, seq_len, seq_offset, line_bases, line_width
                ));
            }
            // Parse header
            let header = &bytes[line_start + 1..line_end];
            name = std::str::from_utf8(header)
                .unwrap_or("")
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            seq_len = 0;
            seq_offset = pos; // first base starts after header \n
            first_seq_line = true;
        } else if !name.is_empty() && line_end > line_start {
            // Sequence line
            let raw_len = line_end - line_start;
            // Trim trailing \r if present
            let base_len = if raw_len > 0 && bytes[line_end - 1] == b'\r' {
                raw_len - 1
            } else {
                raw_len
            };
            seq_len += base_len;
            if first_seq_line && has_newline {
                line_bases = base_len;
                line_width = pos - line_start; // includes \n
                first_seq_line = false;
            } else if first_seq_line {
                // Last line without trailing newline
                line_bases = base_len;
                line_width = base_len + 1; // assume standard line width
                first_seq_line = false;
            }
        }
    }
    // Flush last sequence
    if !name.is_empty() {
        entries.push(format!(
            "{}\t{}\t{}\t{}\t{}",
            name, seq_len, seq_offset, line_bases, line_width
        ));
    }

    let fai_content = entries.join("\n") + "\n";
    std::fs::write(&fai_path, &fai_content).map_err(|e| format!("Cannot write .fai: {}", e))?;

    Ok(fai_path)
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

    // Strip trailing \r (Windows CRLF)
    let header_clean = header_line.trim_end_matches('\r');
    let headers: Vec<String> = header_clean.split('\t').map(|s| s.to_string()).collect();
    let mut rows = Vec::new();

    for line in lines {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        let clean = line.trim_end_matches('\r');
        if clean.trim().is_empty() {
            continue;
        }
        rows.push(clean.split('\t').map(|s| s.to_string()).collect());
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
    /// Per-position depth from ALL reads (not just the displayed subset).
    pub coverage: Vec<u32>,
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

/// Extract the base a read contributes at a given reference position (0-based).
fn base_at_ref_pos(record: &bam::Record, ref_pos: u64, min_bq: u8) -> Option<String> {
    let read_start = record
        .alignment_start()
        .and_then(|p| p.ok())
        .map(|p| {
            let v: usize = p.into();
            v as u64 - 1
        })
        .unwrap_or(0);

    let seq = record.sequence();
    let quals = record.quality_scores();
    let cigar = record.cigar();

    let mut r_pos = read_start;
    let mut q_pos: usize = 0;

    for op_result in cigar.iter() {
        let op = match op_result {
            Ok(o) => o,
            Err(_) => return None,
        };
        let len = op.len() as u64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if ref_pos >= r_pos && ref_pos < r_pos + len {
                    let offset = (ref_pos - r_pos) as usize;
                    let qi = q_pos + offset;
                    if qi < seq.len() {
                        let bq: u8 = quals.iter().nth(qi).unwrap_or(0);
                        if bq >= min_bq {
                            let base: u8 = seq.iter().nth(qi)?;
                            return Some((base as char).to_string());
                        }
                    }
                    return None;
                }
                r_pos += len;
                q_pos += len as usize;
            }
            Kind::Insertion | Kind::SoftClip => {
                q_pos += len as usize;
            }
            Kind::Deletion | Kind::Skip => {
                if ref_pos >= r_pos && ref_pos < r_pos + len {
                    return Some("-".to_string());
                }
                r_pos += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    None
}

/// Read reference sequence from an indexed FASTA file.
fn read_fasta_region(
    fasta_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
) -> Result<String, String> {
    let fai_path = format!("{}.fai", fasta_path);
    let fai_file = std::fs::File::open(&fai_path)
        .map_err(|e| format!("Cannot open FASTA index '{}': {}", fai_path, e))?;
    let reader = BufReader::new(fai_file);

    let mut offset: u64 = 0;
    let mut line_bases: u64 = 0;
    let mut line_width: u64 = 0;
    let mut seq_len: u64 = 0;
    let mut found = false;

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Error reading .fai: {}", e))?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 5 && fields[0] == chrom {
            seq_len = fields[1]
                .parse()
                .map_err(|_| "Bad .fai length".to_string())?;
            offset = fields[2]
                .parse()
                .map_err(|_| "Bad .fai offset".to_string())?;
            line_bases = fields[3]
                .parse()
                .map_err(|_| "Bad .fai line_bases".to_string())?;
            line_width = fields[4]
                .parse()
                .map_err(|_| "Bad .fai line_width".to_string())?;
            found = true;
            break;
        }
    }
    if !found {
        return Err(format!("Chromosome '{}' not found in FASTA index", chrom));
    }

    if start >= seq_len {
        return Err(format!(
            "Requested FASTA region {}:{}-{} starts beyond contig length {}",
            chrom,
            start + 1,
            end,
            seq_len
        ));
    }

    let actual_end = end.min(seq_len);
    let window_len = (actual_end - start) as usize;

    if line_bases == 0 || line_width == 0 {
        return Err(format!(
            "Invalid FASTA index for '{}': line_bases={}, line_width={}",
            chrom, line_bases, line_width
        ));
    }

    let start_line = start / line_bases;
    let start_byte = offset + start_line * line_width + (start % line_bases);
    let end_line = actual_end / line_bases;
    let end_byte = offset + end_line * line_width + (actual_end % line_bases);
    let bytes_to_read = (end_byte - start_byte) as usize;

    use std::io::{Read, Seek, SeekFrom};
    let mut fasta_file =
        std::fs::File::open(fasta_path).map_err(|e| format!("Cannot open FASTA: {}", e))?;
    fasta_file
        .seek(SeekFrom::Start(start_byte))
        .map_err(|e| format!("FASTA seek error: {}", e))?;
    let mut buf = vec![0u8; bytes_to_read + 100];
    let n = fasta_file
        .read(&mut buf)
        .map_err(|e| format!("FASTA read error: {}", e))?;

    let mut result: String = buf[..n]
        .iter()
        .filter(|&&b| b != b'\n' && b != b'\r')
        .map(|&b| b.to_ascii_uppercase() as char)
        .take(window_len)
        .collect();
    while result.len() < window_len {
        result.push('N');
    }
    Ok(result)
}

#[tauri::command]
pub fn get_bam_view(request: BamViewRequest) -> Result<BamViewResponse, String> {
    let min_mapq = request.min_mapq.unwrap_or(0);
    let min_bq = request.min_base_quality.unwrap_or(0);
    let max_reads = request.max_reads.unwrap_or(500);

    if request.positions.is_empty() {
        return Err("BAM viewer requires at least one variant position".to_string());
    }

    if request.positions.len() != request.ref_bases.len()
        || request.positions.len() != request.alt_bases.len()
    {
        return Err(format!(
            "BAM viewer request has inconsistent site arrays: positions={}, refBases={}, altBases={}",
            request.positions.len(),
            request.ref_bases.len(),
            request.alt_bases.len()
        ));
    }

    if request.positions.contains(&0) {
        return Err("BAM viewer positions must be 1-based positive coordinates".to_string());
    }

    if request.window_start == 0 || request.window_end < request.window_start {
        return Err(format!(
            "Invalid BAM viewer window: {}-{}",
            request.window_start, request.window_end
        ));
    }

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

    // Frontend and TSV positions are 1-based inclusive. Internally, BAM/FASTA
    // helpers use 0-based reference offsets, so convert only at IO boundaries.
    let window_len = (request.window_end - request.window_start + 1) as usize;
    let zero_based_start = request.window_start - 1;
    let zero_based_end_exclusive = request.window_end;
    let reference = read_fasta_region(
        &request.fasta_path,
        &request.chrom,
        zero_based_start,
        zero_based_end_exclusive,
    )?;

    let mut bam_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&request.bam_path)
        .map_err(|e| format!("Cannot open BAM: {}", e))?;
    let header = bam_reader
        .read_header()
        .map_err(|e| format!("Cannot read BAM header: {}", e))?;

    let region_str = format!(
        "{}:{}-{}",
        request.chrom, request.window_start, request.window_end
    );
    let region: noodles::core::Region = region_str
        .parse()
        .map_err(|e| format!("Invalid region '{}': {}", region_str, e))?;
    let mut query = bam_reader
        .query(&header, &region)
        .map_err(|e| format!("BAM query error: {}", e))?;

    let mut all_reads: Vec<BamReadView> = Vec::new();
    let mut counts = BamSupportCounts::default();
    let mut total_reads: usize = 0;

    let mut record = bam::Record::default();
    while query
        .read_record(&mut record)
        .map_err(|e| format!("BAM read error: {}", e))?
        != 0
    {
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_qc_fail()
            || flags.is_duplicate()
        {
            continue;
        }

        let mapq = record
            .mapping_quality()
            .map(|q: noodles::sam::alignment::record::MappingQuality| q.get())
            .unwrap_or(255);
        if mapq < min_mapq {
            continue;
        }

        total_reads += 1;

        let mut bases: Vec<String> = Vec::with_capacity(window_len);
        for off in 0..window_len {
            let one_based_pos = request.window_start + off as u64;
            let zero_based_pos = one_based_pos - 1;
            match base_at_ref_pos(&record, zero_based_pos, min_bq) {
                Some(b) => bases.push(b),
                None => bases.push(" ".to_string()),
            }
        }

        let mut matches_alt = 0usize;
        let mut matches_ref = 0usize;
        let mut covered_sites = 0usize;
        for &pos in &request.positions {
            if pos < request.window_start || pos > request.window_end {
                continue;
            }
            let off = (pos - request.window_start) as usize;
            if off < bases.len() {
                let base = &bases[off];
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
            "other"
        } else if matches_alt == n_sites {
            "mnv"
        } else if matches_alt > 0 && matches_alt < n_sites {
            "partial"
        } else if matches_ref == n_sites {
            "reference"
        } else {
            "other"
        }
        .to_string();

        match support.as_str() {
            "mnv" => counts.mnv += 1,
            "partial" => counts.partial += 1,
            "reference" => counts.reference += 1,
            _ => counts.other += 1,
        }

        let strand = if flags.is_reverse_complemented() {
            "-"
        } else {
            "+"
        }
        .to_string();
        let read_name = record
            .name()
            .map(|n| String::from_utf8_lossy(<_ as AsRef<[u8]>>::as_ref(&n)).to_string())
            .unwrap_or_else(|| "?".to_string());

        let read_start = record
            .alignment_start()
            .and_then(|p| p.ok())
            .map(|p| {
                let v: usize = p.into();
                v as u64 - 1
            })
            .unwrap_or(0);
        let aligned_len: u64 = record
            .cigar()
            .iter()
            .filter_map(|r| r.ok())
            .map(|op| match op.kind() {
                Kind::Match
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Skip => op.len() as u64,
                _ => 0,
            })
            .sum();

        all_reads.push(BamReadView {
            name: read_name,
            strand,
            support,
            start: read_start + 1,
            end: read_start + aligned_len,
            mapq,
            bases,
        });
    }

    counts.total = total_reads;

    let coverage: Vec<u32> = {
        let mut depths = vec![0u32; window_len];
        for read in &all_reads {
            for (i, base) in read.bases.iter().enumerate() {
                if base != " " {
                    depths[i] += 1;
                }
            }
        }
        depths
    };

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
        coverage,
    })
}

// ---------------------------------------------------------------------------
// write_text_file
// ---------------------------------------------------------------------------

#[tauri::command]
pub fn write_text_file(path: String, content: String) -> Result<(), String> {
    std::fs::write(&path, &content).map_err(|e| format!("Failed to write file '{}': {}", path, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unique_temp_path(extension: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "get_mnv_gui_test_{}.{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            extension
        ))
    }

    fn write_indexed_fasta(content: &str) -> String {
        use std::io::Write;

        let path = unique_temp_path("fasta");
        let mut file = std::fs::File::create(&path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        let path_str = path.to_string_lossy().into_owned();
        ensure_fasta_index(path_str.clone()).unwrap();
        path_str
    }

    fn minimal_config(variant_file: &str) -> AnalysisConfig {
        AnalysisConfig {
            vcf_file: variant_file.to_string(),
            input_format: None,
            bam_file: None,
            fasta_file: "/tmp/ref.fasta".to_string(),
            genes_file: "/tmp/ref.gff".to_string(),
            sample: None,
            chrom: None,
            normalize_alleles: None,
            min_quality: None,
            min_mapq: None,
            threads: None,
            min_snp_reads: None,
            min_snp_frequency: None,
            min_mnv_reads: None,
            min_mnv_frequency: None,
            min_snp_strand_reads: None,
            min_mnv_strand_reads: None,
            min_strand_bias_p: None,
            dry_run: None,
            strict: None,
            split_multiallelic: None,
            emit_filtered: None,
            vcf_gz: None,
            index_vcf_gz: None,
            strand_bias_info: None,
            keep_original_info: None,
            bcf: None,
            summary_json: None,
            error_json: None,
            run_manifest: None,
            gff_features: None,
            exclude_intergenic: None,
            translation_table: None,
            output_tsv: None,
            output_vcf: None,
            convert: None,
            both: None,
            output_dir: None,
            output_prefix: None,
        }
    }

    #[test]
    fn test_gui_defaults_output_dir_to_variant_parent() {
        let args = minimal_config("/tmp/sample/sample_variants.tsv").into_args();
        assert_eq!(args.output_dir.as_deref(), Some("/tmp/sample"));
    }

    #[test]
    fn test_gui_preserves_explicit_output_dir() {
        let mut config = minimal_config("/tmp/sample/sample_variants.tsv");
        config.output_dir = Some("/tmp/output".to_string());
        let args = config.into_args();
        assert_eq!(args.output_dir.as_deref(), Some("/tmp/output"));
    }

    #[test]
    fn test_resolve_output_paths_matches_core_suffixes() {
        let mut config = minimal_config("/tmp/sample/sample_variants.tsv");
        config.output_tsv = Some(true);
        config.output_vcf = Some(true);
        config.vcf_gz = Some(true);

        let paths = resolve_output_paths(config).unwrap();

        assert_eq!(
            paths,
            vec![
                "/tmp/sample/sample_variants.MNV.tsv".to_string(),
                "/tmp/sample/sample_variants.MNV.vcf.gz".to_string()
            ]
        );
    }

    #[test]
    fn test_resolve_output_paths_honors_output_prefix() {
        let mut config = minimal_config("/tmp/sample.vcf.gz");
        config.output_dir = Some("/tmp/out".to_string());
        config.output_prefix = Some("custom".to_string());
        config.output_tsv = Some(false);
        config.output_vcf = Some(true);

        let paths = resolve_output_paths(config).unwrap();

        assert_eq!(paths, vec!["/tmp/out/custom.MNV.vcf".to_string()]);
    }

    #[test]
    fn test_gui_config_deserializes_lowercase_input_format() {
        let value = serde_json::json!({
            "vcfFile": "/tmp/sample.vcf",
            "inputFormat": "auto",
            "fastaFile": "/tmp/ref.fasta",
            "genesFile": "/tmp/ref.gff"
        });

        let config: AnalysisConfig = serde_json::from_value(value).unwrap();
        assert_eq!(
            config.input_format,
            Some(get_mnv::cli::VariantInputFormat::Auto)
        );
    }

    #[test]
    fn test_gui_config_deserializes_tsv_input_format_and_ivar_alias() {
        let tsv_value = serde_json::json!({
            "vcfFile": "/tmp/sample_variants.tsv",
            "inputFormat": "tsv",
            "fastaFile": "/tmp/ref.fasta",
            "genesFile": "/tmp/ref.gff"
        });
        let alias_value = serde_json::json!({
            "vcfFile": "/tmp/sample_variants.tsv",
            "inputFormat": "ivar",
            "fastaFile": "/tmp/ref.fasta",
            "genesFile": "/tmp/ref.gff"
        });

        let tsv_config: AnalysisConfig = serde_json::from_value(tsv_value).unwrap();
        let alias_config: AnalysisConfig = serde_json::from_value(alias_value).unwrap();

        assert_eq!(
            tsv_config.input_format,
            Some(get_mnv::cli::VariantInputFormat::Tsv)
        );
        assert_eq!(
            alias_config.input_format,
            Some(get_mnv::cli::VariantInputFormat::Tsv)
        );
    }

    #[test]
    fn test_gui_config_forwards_variant_filters() {
        let mut config = minimal_config("/tmp/sample_variants.tsv");
        config.input_format = Some(get_mnv::cli::VariantInputFormat::Tsv);
        config.min_mapq = Some(20);
        config.min_snp_reads = Some(2);
        config.min_mnv_reads = Some(3);
        config.min_snp_frequency = Some(0.05);
        config.min_mnv_frequency = Some(0.10);
        config.min_snp_strand_reads = Some(2);
        config.min_mnv_strand_reads = Some(2);
        config.strand_bias_info = Some(true);
        config.output_tsv = Some(true);
        config.output_vcf = Some(true);

        let args = config.into_args();

        assert_eq!(args.input_format, get_mnv::cli::VariantInputFormat::Tsv);
        assert_eq!(args.tsv_file.as_deref(), Some("/tmp/sample_variants.tsv"));
        assert!(args.vcf_file.is_none());
        assert_eq!(args.min_mapq, 20);
        assert_eq!(args.min_snp_reads, 2);
        assert_eq!(args.min_mnv_reads, 3);
        assert_eq!(args.min_snp_frequency, 0.05);
        assert_eq!(args.min_mnv_frequency, 0.10);
        assert_eq!(args.min_snp_strand_reads, 2);
        assert_eq!(args.min_mnv_strand_reads, 2);
        assert!(args.strand_bias_info);
        assert!(args.both);
        assert!(!args.convert);
    }

    #[test]
    fn test_detect_variant_input_format_recognizes_vcf_extension() {
        assert_eq!(
            detect_variant_input_format("/tmp/sample.vcf.gz".to_string()).unwrap(),
            "vcf"
        );
    }

    #[test]
    fn test_detect_variant_input_format_recognizes_ivar_tsv() {
        use std::io::Write;

        let path = std::env::temp_dir().join(format!(
            "get_mnv_gui_ivar_detect_{}.tsv",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut file = std::fs::File::create(&path).unwrap();
        writeln!(file, "REGION\tPOS\tREF\tALT\tTOTAL_DP\tALT_FREQ\tPASS").unwrap();

        assert_eq!(
            detect_variant_input_format(path.to_string_lossy().into_owned()).unwrap(),
            "tsv"
        );
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_detect_variant_input_format_leaves_annotation_tsv_unknown() {
        use std::io::Write;

        let path = std::env::temp_dir().join(format!(
            "get_mnv_gui_gene_detect_{}.tsv",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut file = std::fs::File::create(&path).unwrap();
        writeln!(file, "gene\tstart\tend\tstrand").unwrap();

        assert_eq!(
            detect_variant_input_format(path.to_string_lossy().into_owned()).unwrap(),
            "unknown"
        );
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_read_fasta_region_uses_zero_based_offsets_for_tsv_coordinates() {
        let fasta = write_indexed_fasta(">ref_amplicon\nATGGACTTAA\n");

        assert_eq!(
            read_fasta_region(&fasta, "ref_amplicon", 0, 3).unwrap(),
            "ATG"
        );

        // TSV positions 4-5 are 1-based and should resolve to GA.
        assert_eq!(
            read_fasta_region(&fasta, "ref_amplicon", 3, 5).unwrap(),
            "GA"
        );

        let _ = std::fs::remove_file(&fasta);
        let _ = std::fs::remove_file(format!("{fasta}.fai"));
    }

    #[test]
    fn test_bam_view_rejects_mismatched_site_arrays() {
        let err = get_bam_view(BamViewRequest {
            bam_path: "/tmp/missing.bam".to_string(),
            fasta_path: "/tmp/missing.fasta".to_string(),
            chrom: "ref_amplicon".to_string(),
            positions: vec![2704, 2705],
            ref_bases: vec!["G".to_string()],
            alt_bases: vec!["A".to_string(), "G".to_string()],
            min_mapq: Some(0),
            min_base_quality: Some(20),
            max_reads: Some(10),
            window_start: 2624,
            window_end: 2785,
        })
        .unwrap_err();

        assert!(err.contains("inconsistent site arrays"));
    }

    #[test]
    fn test_read_fasta_region_rejects_start_past_contig() {
        let fasta = write_indexed_fasta(">ref_amplicon\nACGT\n");
        let err = read_fasta_region(&fasta, "ref_amplicon", 10_000, 10_010).unwrap_err();

        assert!(err.contains("starts beyond contig length"));

        let _ = std::fs::remove_file(&fasta);
        let _ = std::fs::remove_file(format!("{fasta}.fai"));
    }
}
