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
    pub frameshift_min_freq: Option<f64>,
    pub indel_anchor_depth: Option<bool>,
    pub phased_indel_min_reads: Option<usize>,
    pub count_mates_separately: Option<bool>,
    pub phased_indel_min_freq: Option<f64>,
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
            // Indel-annotation knobs (exposed in the desktop UI); fall back to the
            // CLI defaults when the frontend does not send a value. These must
            // track `src/cli.rs`: the app runs the same engine, so a default
            // that drifts here makes the desktop build answer differently from
            // the command line on the same inputs.
            frameshift_min_freq: self.frameshift_min_freq.unwrap_or(0.5),
            indel_anchor_depth: self.indel_anchor_depth.unwrap_or(true),
            phased_indel_min_reads: self.phased_indel_min_reads.unwrap_or(2),
            phased_indel_min_freq: self.phased_indel_min_freq.unwrap_or(0.0),
            count_mates_separately: self.count_mates_separately.unwrap_or(false),
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
            // The desktop app has its own results viewer, so it does not ask the
            // engine for the standalone HTML report.
            report: None,
            report_from: Vec::new(),
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

/// Where the run's outputs go, directory and all, without the suffixes.
fn output_stem(args: &get_mnv::cli::Args) -> Result<String, String> {
    let base_name = get_mnv::io::get_base_name(args.variant_file()).map_err(|e| e.to_string())?;
    let stem_name = args.output_prefix.clone().unwrap_or(base_name);
    Ok(match &args.output_dir {
        Some(dir) => std::path::Path::new(dir)
            .join(&stem_name)
            .to_string_lossy()
            .into_owned(),
        None => stem_name,
    })
}

/// The suffixes this configuration writes, which is what decides both the
/// preview and what an existing file has to look like to be in the way.
fn output_suffixes(args: &get_mnv::cli::Args) -> Vec<&'static str> {
    let mut suffixes = Vec::new();
    if args.both || !args.convert {
        suffixes.push(".MNV.tsv");
    }
    if args.both || args.convert {
        suffixes.push(if args.vcf_gz {
            ".MNV.vcf.gz"
        } else {
            ".MNV.vcf"
        });
    }
    if args.bcf {
        suffixes.push(".MNV.bcf");
    }
    suffixes
}

fn expected_output_paths(config: &AnalysisConfig) -> Result<Vec<String>, String> {
    let args = config.clone().into_args();
    if args.dry_run {
        return Ok(Vec::new());
    }
    let stem = output_stem(&args)?;
    // With --sample all the run writes one file set per sample, each carrying the
    // sample in its name, so a single stem would name a file that is never
    // written. Listing per-sample names needs the sample list, which this preview
    // does not read, so say the outputs are per sample instead of naming one that
    // will not exist. The chosen directory belongs in that pattern like anywhere
    // else: leaving it out told the user their cohort files land somewhere else.
    if args.sample.as_deref() == Some("all") {
        return Ok(vec![format!("{stem}.sample_<SAMPLE>.MNV.*")]);
    }
    Ok(output_suffixes(&args)
        .into_iter()
        .map(|suffix| format!("{stem}{suffix}"))
        .collect())
}

/// The per-sample files a `--sample all` run would overwrite. The preview names
/// a pattern, because the sample list is not read here, and a pattern is not a
/// path: asking whether it exists always answered no, so a cohort run reported
/// nothing in the way while its own outputs from the last run sat right there.
/// The files can be looked for without knowing which samples the VCF holds.
fn existing_per_sample_outputs(args: &get_mnv::cli::Args) -> Result<Vec<String>, String> {
    let stem = output_stem(args)?;
    let stem_path = std::path::Path::new(&stem);
    let dir = match stem_path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent.to_path_buf(),
        _ => std::path::PathBuf::from("."),
    };
    let prefix = format!(
        "{}.sample_",
        stem_path.file_name().unwrap_or_default().to_string_lossy()
    );
    let suffixes = output_suffixes(args);
    let entries = match std::fs::read_dir(&dir) {
        Ok(entries) => entries,
        Err(_) => return Ok(Vec::new()),
    };
    let mut found: Vec<String> = entries
        .filter_map(|entry| entry.ok())
        .filter(|entry| {
            let name = entry.file_name().to_string_lossy().into_owned();
            name.starts_with(&prefix) && suffixes.iter().any(|suffix| name.ends_with(suffix))
        })
        .map(|entry| entry.path().to_string_lossy().into_owned())
        .collect();
    found.sort();
    Ok(found)
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
    let args = config.clone().into_args();
    if args.dry_run {
        return Ok(Vec::new());
    }
    if args.sample.as_deref() == Some("all") {
        return existing_per_sample_outputs(&args);
    }
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
    pub columns: Vec<BamViewColumn>,
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
pub struct BamViewColumn {
    pub key: String,
    pub position: u64,
    pub kind: String,
    pub insertion_index: Option<u64>,
    pub label: String,
    pub reference_base: String,
    pub is_variant: bool,
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

#[derive(Debug, Clone)]
struct RawBamReadView {
    name: String,
    strand: String,
    support: String,
    start: u64,
    end: u64,
    mapq: u8,
    reference_bases: HashMap<u64, String>,
    insertions_after: HashMap<u64, Vec<String>>,
}

fn bam_support_rank(support: &str) -> u8 {
    match support {
        "mnv" => 0,
        "partial" => 1,
        "reference" => 2,
        _ => 3,
    }
}

fn sort_bam_reads_for_display(reads: &mut [BamReadView]) {
    reads.sort_by(|a, b| {
        bam_support_rank(&a.support)
            .cmp(&bam_support_rank(&b.support))
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
            .then_with(|| a.name.cmp(&b.name))
    });
}

fn select_bam_reads_for_display(
    mut reads: Vec<BamReadView>,
    max_reads: usize,
) -> (Vec<BamReadView>, bool) {
    let truncated = reads.len() > max_reads;
    sort_bam_reads_for_display(&mut reads);

    if !truncated {
        return (reads, false);
    }

    if max_reads == 0 {
        return (Vec::new(), true);
    }

    let mut buckets: [Vec<BamReadView>; 4] = std::array::from_fn(|_| Vec::new());
    for read in reads {
        let rank = bam_support_rank(&read.support) as usize;
        buckets[rank].push(read);
    }

    let mut selected = Vec::with_capacity(max_reads);
    let mut offsets = [0usize; 4];
    while selected.len() < max_reads {
        let mut added = false;
        for bucket_idx in 0..buckets.len() {
            if selected.len() >= max_reads {
                break;
            }
            let offset = offsets[bucket_idx];
            if let Some(read) = buckets[bucket_idx].get(offset).cloned() {
                selected.push(read);
                offsets[bucket_idx] += 1;
                added = true;
            }
        }
        if !added {
            break;
        }
    }

    sort_bam_reads_for_display(&mut selected);
    (selected, true)
}

fn expected_insertion_after(
    position: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<(u64, usize)> {
    let ref_chars: Vec<char> = ref_allele.chars().collect();
    let alt_chars: Vec<char> = alt_allele.chars().collect();
    if alt_chars.len() <= ref_chars.len() {
        return None;
    }

    let min_len = ref_chars.len().min(alt_chars.len());
    let mut prefix = 0usize;
    while prefix < min_len && ref_chars[prefix].eq_ignore_ascii_case(&alt_chars[prefix]) {
        prefix += 1;
    }
    let inserted_len = alt_chars.len().saturating_sub(ref_chars.len());
    let anchor = if prefix == 0 {
        position
    } else {
        position + prefix as u64 - 1
    };
    Some((anchor, inserted_len))
}

/// An inserted base the quality floor rejected. A read whose insertion cannot be
/// read still carries one, which is not the same as carrying none, so the slot
/// is kept and the allele it belongs to is withdrawn rather than shortened.
const UNREADABLE_BASE: &str = "?";

fn observed_allele_from_layout(
    reference_bases: &HashMap<u64, String>,
    insertions_after: &HashMap<u64, Vec<String>>,
    position: u64,
    ref_allele: &str,
) -> Option<String> {
    let ref_len = ref_allele.chars().count() as u64;
    if ref_len == 0 {
        return None;
    }
    let mut allele = String::new();
    let mut observed_any = false;
    for offset in 0..ref_len {
        let pos = position + offset;
        let base = reference_bases.get(&pos)?;
        if base != "-" {
            allele.push_str(base);
            observed_any = true;
        }
        if let Some(inserted) = insertions_after.get(&pos) {
            for base in inserted {
                // An allele half of which could not be read is not an
                // observation of a shorter allele: the read has nothing to say.
                if base == UNREADABLE_BASE {
                    return None;
                }
                allele.push_str(base);
                observed_any = true;
            }
        }
    }
    observed_any.then_some(allele)
}

fn read_alignment_layout(
    record: &bam::Record,
    window_start: u64,
    window_end: u64,
    min_bq: u8,
) -> (HashMap<u64, String>, HashMap<u64, Vec<String>>) {
    let read_start = record
        .alignment_start()
        .and_then(|p| p.ok())
        .map(|p| {
            let v: usize = p.into();
            v as u64
        })
        .unwrap_or(1);

    let seq = record.sequence();
    let quals = record.quality_scores();
    let cigar = record.cigar();

    let mut r_pos = read_start;
    let mut q_pos: usize = 0;
    let mut reference_bases: HashMap<u64, String> = HashMap::new();
    let mut insertions_after: HashMap<u64, Vec<String>> = HashMap::new();

    for op_result in cigar.iter() {
        let op = match op_result {
            Ok(o) => o,
            Err(_) => return (reference_bases, insertions_after),
        };
        let len = op.len() as u64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for offset in 0..len {
                    let pos = r_pos + offset;
                    if pos < window_start || pos > window_end {
                        continue;
                    }
                    let qi = q_pos + offset as usize;
                    if qi < seq.len() {
                        let bq: u8 = quals.iter().nth(qi).unwrap_or(0);
                        if bq >= min_bq {
                            if let Some(base) = seq.iter().nth(qi) {
                                reference_bases
                                    .insert(pos, (base as char).to_ascii_uppercase().to_string());
                            }
                        }
                    }
                }
                r_pos += len;
                q_pos += len as usize;
            }
            Kind::Insertion => {
                let anchor = r_pos.saturating_sub(1);
                if anchor >= window_start && anchor <= window_end {
                    let inserted = insertions_after.entry(anchor).or_default();
                    for offset in 0..len {
                        let qi = q_pos + offset as usize;
                        if qi >= seq.len() {
                            continue;
                        }
                        let bq: u8 = quals.iter().nth(qi).unwrap_or(0);
                        // Dropping a base the floor rejected shortened the
                        // insertion this read carries, and a read whose whole
                        // insertion was rejected came out holding none, which
                        // counts as evidence for the reference allele. The slot
                        // is kept and marked unreadable instead, so the read
                        // withdraws from the judgement the way a rejected base
                        // in a match already does.
                        let base = match seq.iter().nth(qi) {
                            Some(base) if bq >= min_bq => {
                                (base as char).to_ascii_uppercase().to_string()
                            }
                            _ => UNREADABLE_BASE.to_string(),
                        };
                        inserted.push(base);
                    }
                }
                q_pos += len as usize;
            }
            Kind::SoftClip => {
                q_pos += len as usize;
            }
            // A deletion is the read saying those reference bases are not in it,
            // which is evidence about a deletion allele.
            Kind::Deletion => {
                for offset in 0..len {
                    let pos = r_pos + offset;
                    if pos >= window_start && pos <= window_end {
                        reference_bases.insert(pos, "-".to_string());
                    }
                }
                r_pos += len;
            }
            // A skip is the read not being there at all: an aligner writes it for
            // an intron, and get_MNV builds spliced CDS models, so these reach the
            // viewer. Recording it as a gap made a read that observed nothing
            // count as carrying a deletion. Leaving the positions unset is what
            // the viewer already draws for a base a read does not reach.
            Kind::Skip => {
                r_pos += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    (reference_bases, insertions_after)
}

fn classify_layout_support(
    reference_bases: &HashMap<u64, String>,
    insertions_after: &HashMap<u64, Vec<String>>,
    sites: &[BamVariantSite],
) -> String {
    let mut matches_alt = 0usize;
    let mut matches_ref = 0usize;
    let mut covered_sites = 0usize;

    for site in sites {
        let Some(observed) = observed_allele_from_layout(
            reference_bases,
            insertions_after,
            site.position,
            &site.reference_base,
        ) else {
            continue;
        };
        covered_sites += 1;
        if observed.eq_ignore_ascii_case(&site.alt_base) {
            matches_alt += 1;
        } else if observed.eq_ignore_ascii_case(&site.reference_base) {
            matches_ref += 1;
        }
    }

    let n_sites = sites.len();
    if covered_sites == 0 {
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
    .to_string()
}

fn build_bam_view_columns(
    display_start: u64,
    display_end: u64,
    reference: &str,
    sites: &[BamVariantSite],
    insertion_widths: &HashMap<u64, usize>,
) -> Vec<BamViewColumn> {
    debug_assert_eq!(
        display_end,
        display_start + reference.chars().count().saturating_sub(1) as u64
    );
    let mut variant_positions = BTreeSet::new();
    let mut variant_insertion_anchors = BTreeSet::new();
    for site in sites {
        let ref_len = site.reference_base.chars().count().max(1) as u64;
        for pos in site.position..site.position + ref_len {
            variant_positions.insert(pos);
        }
        if let Some((anchor, _)) =
            expected_insertion_after(site.position, &site.reference_base, &site.alt_base)
        {
            variant_insertion_anchors.insert(anchor);
        }
    }
    let mut columns = Vec::new();
    for (idx, base) in reference.chars().enumerate() {
        let pos = display_start + idx as u64;
        let is_variant = variant_positions.contains(&pos);
        columns.push(BamViewColumn {
            key: format!("ref:{pos}"),
            position: pos,
            kind: "ref".to_string(),
            insertion_index: None,
            label: pos.to_string(),
            reference_base: base.to_string(),
            is_variant,
        });

        let width = insertion_widths.get(&pos).copied().unwrap_or(0);
        for insertion_idx in 0..width {
            columns.push(BamViewColumn {
                key: format!("ins:{pos}:{}", insertion_idx + 1),
                position: pos,
                kind: "ins".to_string(),
                insertion_index: Some((insertion_idx + 1) as u64),
                label: format!("+{}", insertion_idx + 1),
                reference_base: "+".to_string(),
                is_variant: variant_insertion_anchors.contains(&pos),
            });
        }
    }
    columns
}

fn expand_read_bases(read: &RawBamReadView, columns: &[BamViewColumn]) -> Vec<String> {
    columns
        .iter()
        .map(|column| {
            if column.kind == "ins" {
                let idx = column.insertion_index.unwrap_or(1).saturating_sub(1) as usize;
                read.insertions_after
                    .get(&column.position)
                    .and_then(|bases| bases.get(idx))
                    .cloned()
                    .unwrap_or_else(|| " ".to_string())
            } else {
                read.reference_bases
                    .get(&column.position)
                    .cloned()
                    .unwrap_or_else(|| " ".to_string())
            }
        })
        .collect()
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
    // Read to the end of the span instead of keeping whatever one call returns:
    // a short read is not a short file, and padding made the two look alike.
    let mut buf = Vec::with_capacity(bytes_to_read + 100);
    fasta_file
        .by_ref()
        .take((bytes_to_read + 100) as u64)
        .read_to_end(&mut buf)
        .map_err(|e| format!("FASTA read error: {}", e))?;

    let result: String = buf
        .iter()
        // The span is over-read by a hundred bytes to cover line breaks. A
        // record header can never be part of a sequence, so stopping at one
        // keeps the next contig out of the answer.
        .take_while(|&&b| b != b'>')
        .filter(|&&b| b != b'\n' && b != b'\r')
        .map(|&b| b.to_ascii_uppercase() as char)
        .take(window_len)
        .collect();
    // Padding answered with bases the file does not hold. The index says the
    // contig is seq_len long; when the file cannot supply that, the two
    // disagree, and saying so is the only honest answer.
    if result.len() < window_len {
        return Err(format!(
            "FASTA '{}' is shorter than its index claims: {}:{}-{} needs {} base(s) \
             and the file supplies {}. Rebuild the index, for example with \
             `samtools faidx '{}'`",
            fasta_path,
            chrom,
            start + 1,
            actual_end,
            window_len,
            result.len(),
            fasta_path
        ));
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

    let mut sites: Vec<BamVariantSite> = Vec::new();
    let mut insertion_widths: HashMap<u64, usize> = HashMap::new();
    for i in 0..request.positions.len() {
        let pos = request.positions[i];
        let rb = request.ref_bases[i].clone();
        let ab = request.alt_bases[i].clone();
        if let Some((anchor, inserted_len)) = expected_insertion_after(pos, &rb, &ab) {
            let width = insertion_widths.entry(anchor).or_default();
            *width = (*width).max(inserted_len);
        }
        sites.push(BamVariantSite {
            position: pos,
            reference_base: rb,
            alt_base: ab,
        });
    }

    // Frontend and TSV positions are 1-based inclusive. Internally, BAM/FASTA
    // helpers use 0-based reference offsets, so convert only at IO boundaries.
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

    let mut raw_reads: Vec<RawBamReadView> = Vec::new();
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
        let (reference_bases, insertions_after) =
            read_alignment_layout(&record, request.window_start, request.window_end, min_bq);
        for (&anchor, inserted) in &insertions_after {
            let width = insertion_widths.entry(anchor).or_default();
            *width = (*width).max(inserted.len());
        }

        let support = classify_layout_support(&reference_bases, &insertions_after, &sites);

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

        raw_reads.push(RawBamReadView {
            name: read_name,
            strand,
            support,
            start: read_start + 1,
            end: read_start + aligned_len,
            mapq,
            reference_bases,
            insertions_after,
        });
    }

    counts.total = total_reads;

    let columns = build_bam_view_columns(
        request.window_start,
        request.window_end,
        &reference,
        &sites,
        &insertion_widths,
    );

    let mut all_reads: Vec<BamReadView> = raw_reads
        .iter()
        .map(|read| BamReadView {
            name: read.name.clone(),
            strand: read.strand.clone(),
            support: read.support.clone(),
            start: read.start,
            end: read.end,
            mapq: read.mapq,
            bases: expand_read_bases(read, &columns),
        })
        .collect();

    let coverage: Vec<u32> = {
        let mut depths = vec![0u32; columns.len()];
        for read in &all_reads {
            for (i, base) in read.bases.iter().enumerate() {
                if base != " " {
                    depths[i] += 1;
                }
            }
        }
        depths
    };

    let (display_reads, truncated) =
        select_bam_reads_for_display(std::mem::take(&mut all_reads), max_reads);

    Ok(BamViewResponse {
        chrom: request.chrom,
        display_start: request.window_start,
        display_end: request.window_end,
        reference,
        columns,
        sites,
        reads: display_reads,
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

    /// A clock is not an identity. This asked for the time in nanoseconds and
    /// used it as a name, but the clock does not tick that finely: over sixteen
    /// thousand calls from eight threads, eleven thousand of them named a file
    /// another thread had just named. Tests run in parallel, so two of them
    /// shared a FASTA, and one truncating it failed the other.
    fn unique_temp_path(extension: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        std::env::temp_dir().join(format!(
            "get_mnv_gui_test_{}_{}.{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed),
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

    /// A FASTA that no longer matches its index used to answer with `N` for
    /// every base it could not supply, so a truncated file read as a contig
    /// that genuinely holds unknown bases there.
    #[test]
    fn test_read_fasta_region_reports_a_file_shorter_than_its_index() {
        use std::io::Write;
        let path = write_indexed_fasta(">chr1\nACGTACGTACGTACGTACGT\n");
        let mut file = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        file.set_len(">chr1\n".len() as u64 + 10).unwrap();
        file.flush().unwrap();

        // The bases the file still holds come back unchanged.
        assert_eq!(
            read_fasta_region(&path, "chr1", 0, 10).unwrap(),
            "ACGTACGTAC"
        );

        let err = read_fasta_region(&path, "chr1", 0, 20)
            .expect_err("a window past the end of a truncated file must not be padded");
        assert!(err.contains("shorter than its index"), "{err}");
        assert!(err.contains("needs 20 base(s)"), "{err}");
        assert!(err.contains("supplies 10"), "{err}");
    }

    /// The read window is over-read to cover line breaks, and the slack must
    /// never let the following record's sequence stand in for this one's.
    #[test]
    fn test_read_fasta_region_stops_at_the_next_record() {
        let path = write_indexed_fasta(">chr1\nACGTACGTAC\n>chr2\nTTTTTTTTTT\n");
        assert_eq!(
            read_fasta_region(&path, "chr1", 0, 10).unwrap(),
            "ACGTACGTAC"
        );
        assert_eq!(read_fasta_region(&path, "chr1", 5, 10).unwrap(), "CGTAC");
        assert_eq!(
            read_fasta_region(&path, "chr2", 0, 10).unwrap(),
            "TTTTTTTTTT"
        );
    }

    /// A cohort run, its output directory, and the files a previous run left.
    fn cohort_config_in_a_fresh_dir(existing: &[&str]) -> (AnalysisConfig, std::path::PathBuf) {
        let dir = unique_temp_path("outdir");
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("cohort.vcf");
        std::fs::write(&vcf, "x").unwrap();
        for name in existing {
            std::fs::write(dir.join(name), "x").unwrap();
        }
        let mut config = minimal_config(vcf.to_str().unwrap());
        config.output_dir = Some(dir.to_string_lossy().into_owned());
        config.sample = Some("all".to_string());
        (config, dir)
    }

    /// The preview names a pattern for a cohort, and a pattern never exists as a
    /// file, so the conflict check answered that nothing was in the way while
    /// the last run's per-sample outputs sat in the very directory chosen.
    #[test]
    fn test_cohort_conflicts_name_the_per_sample_files_already_there() {
        let (config, dir) =
            cohort_config_in_a_fresh_dir(&["cohort.sample_S1.MNV.tsv", "cohort.sample_S2.MNV.tsv"]);
        let conflicts = check_output_conflicts(config).unwrap();
        assert_eq!(
            conflicts,
            vec![
                dir.join("cohort.sample_S1.MNV.tsv")
                    .to_string_lossy()
                    .into_owned(),
                dir.join("cohort.sample_S2.MNV.tsv")
                    .to_string_lossy()
                    .into_owned(),
            ]
        );
    }

    /// Only the formats this run writes are in its way.
    #[test]
    fn test_cohort_conflicts_ignore_files_the_run_will_not_write() {
        let (config, _dir) = cohort_config_in_a_fresh_dir(&[
            "cohort.sample_S1.MNV.vcf",
            "cohort.sample_S1.MNV.bcf",
            "other.sample_S1.MNV.tsv",
            "cohort.MNV.tsv",
        ]);
        assert!(
            check_output_conflicts(config).unwrap().is_empty(),
            "a default cohort run writes TSVs under its own stem and nothing else"
        );
    }

    /// The pattern the preview shows has to point where the files land.
    #[test]
    fn test_cohort_preview_carries_the_output_directory() {
        let (config, dir) = cohort_config_in_a_fresh_dir(&[]);
        assert_eq!(
            resolve_output_paths(config).unwrap(),
            vec![dir
                .join("cohort.sample_<SAMPLE>.MNV.*")
                .to_string_lossy()
                .into_owned()]
        );
    }

    /// A one-read BAM, written and read back, so the layout builder can be
    /// exercised on a real record rather than on maps written by hand.
    fn layout_of(
        cigar_ops: &str,
        start: u64,
        seq: &str,
        window: (u64, u64),
    ) -> (HashMap<u64, String>, HashMap<u64, Vec<String>>) {
        layout_of_q(cigar_ops, start, seq, window, &vec![40u8; seq.len()], 20)
    }

    fn layout_of_q(
        cigar_ops: &str,
        start: u64,
        seq: &str,
        window: (u64, u64),
        quals: &[u8],
        min_bq: u8,
    ) -> (HashMap<u64, String>, HashMap<u64, Vec<String>>) {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use noodles::sam::header::record::value::Map;
        use std::num::NonZeroUsize;

        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(10000).unwrap()),
            )
            .build();

        let ops: Vec<noodles::sam::alignment::record::cigar::Op> = cigar_ops
            .split(',')
            .map(|token| {
                let (len, kind) = token.split_at(token.len() - 1);
                let kind = match kind {
                    "M" => Kind::Match,
                    "N" => Kind::Skip,
                    "D" => Kind::Deletion,
                    "I" => Kind::Insertion,
                    other => panic!("unsupported cigar op {other}"),
                };
                noodles::sam::alignment::record::cigar::Op::new(kind, len.parse().unwrap())
            })
            .collect();

        let record = RecordBuf::builder()
            .set_name("r1")
            .set_flags(noodles::sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles::core::Position::try_from(start as usize).unwrap())
            .set_mapping_quality(MappingQuality::try_from(60).unwrap())
            .set_cigar(Cigar::from(ops))
            .set_sequence(Sequence::from(seq.as_bytes().to_vec()))
            .set_quality_scores(QualityScores::from(quals.to_vec()))
            .build();

        let path = unique_temp_path("bam");
        {
            let mut writer = noodles::bam::io::Writer::new(std::fs::File::create(&path).unwrap());
            writer.write_header(&header).unwrap();
            writer.write_alignment_record(&header, &record).unwrap();
        }
        let mut reader = noodles::bam::io::Reader::new(std::fs::File::open(&path).unwrap());
        reader.read_header().unwrap();
        let mut got = noodles::bam::Record::default();
        reader.read_record(&mut got).unwrap();
        read_alignment_layout(&got, window.0, window.1, min_bq)
    }

    /// A deletion is the read saying those reference bases are not in it. A skip
    /// is the read not being there: an aligner writes N for an intron, and
    /// get_MNV builds spliced CDS models, so both reach this viewer. Marking
    /// them the same made a read that never observed a deletion count as
    /// carrying it.
    #[test]
    fn test_a_skipped_region_is_not_read_as_a_deletion() {
        let deleted = layout_of("1M,3D,2M", 11, "ACG", (1, 100));
        let skipped = layout_of("1M,3N,2M", 11, "ACG", (1, 100));

        // The read that deletes says so at each base it deletes.
        assert_eq!(deleted.0.get(&12), Some(&"-".to_string()));
        assert_eq!(deleted.0.get(&14), Some(&"-".to_string()));
        // The read that skips says nothing there, which is what an empty cell in
        // the viewer means, and the bases it does align are untouched.
        assert_eq!(skipped.0.get(&12), None);
        assert_eq!(skipped.0.get(&14), None);
        assert_eq!(skipped.0.get(&11), Some(&"A".to_string()));
        assert_eq!(skipped.0.get(&15), Some(&"C".to_string()));

        let sites = vec![BamVariantSite {
            position: 11,
            reference_base: "AGGG".to_string(),
            alt_base: "A".to_string(),
        }];
        assert_eq!(
            classify_layout_support(&deleted.0, &deleted.1, &sites),
            "mnv",
            "the read that deletes those bases carries the deletion"
        );
        assert_eq!(
            classify_layout_support(&skipped.0, &skipped.1, &sites),
            "other",
            "the read that skips them observed nothing about the deletion"
        );
    }

    /// A read whose insertion the quality floor rejects still carries an
    /// insertion. Dropping those bases left it holding none, which reads as the
    /// reference allele, so a read that could say nothing was counted as
    /// evidence against the variant.
    #[test]
    fn test_an_unreadable_insertion_is_not_evidence_for_the_reference() {
        let sites = vec![BamVariantSite {
            position: 11,
            reference_base: "A".to_string(),
            alt_base: "ATT".to_string(),
        }];
        let readable = layout_of_q("1M,2I,2M", 11, "ATTCG", (1, 100), &[40, 40, 40, 40, 40], 20);
        let unreadable = layout_of_q("1M,2I,2M", 11, "ATTCG", (1, 100), &[40, 5, 5, 40, 40], 20);
        let half = layout_of_q("1M,2I,2M", 11, "ATTCG", (1, 100), &[40, 40, 5, 40, 40], 20);
        let without = layout_of_q("3M", 11, "ACG", (1, 100), &[40, 40, 40], 20);

        // The insertion is kept whole either way, so the column keeps its width.
        assert_eq!(
            unreadable.1.get(&11),
            Some(&vec![
                UNREADABLE_BASE.to_string(),
                UNREADABLE_BASE.to_string()
            ])
        );
        assert_eq!(
            half.1.get(&11),
            Some(&vec!["T".to_string(), UNREADABLE_BASE.to_string()])
        );

        // Three different things, three different answers.
        assert_eq!(
            classify_layout_support(&readable.0, &readable.1, &sites),
            "mnv"
        );
        assert_eq!(
            classify_layout_support(&without.0, &without.1, &sites),
            "reference"
        );
        for (name, layout) in [("all bases", &unreadable), ("one base", &half)] {
            assert_eq!(
                classify_layout_support(&layout.0, &layout.1, &sites),
                "other",
                "a read with {name} of its insertion rejected says nothing about it"
            );
        }
    }

    fn bam_read_for_test(name: &str, support: &str, start: u64) -> BamReadView {
        BamReadView {
            name: name.to_string(),
            strand: "+".to_string(),
            support: support.to_string(),
            start,
            end: start + 10,
            mapq: 60,
            bases: vec!["A".to_string()],
        }
    }

    pub(super) fn minimal_config(variant_file: &str) -> AnalysisConfig {
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
            frameshift_min_freq: None,
            indel_anchor_depth: None,
            phased_indel_min_reads: None,
            count_mates_separately: None,
            phased_indel_min_freq: None,
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
    fn test_gui_config_forwards_indel_knobs() {
        let mut config = minimal_config("/tmp/sample.vcf");
        config.frameshift_min_freq = Some(0.5);
        config.indel_anchor_depth = Some(true);
        config.phased_indel_min_reads = Some(4);
        config.phased_indel_min_freq = Some(0.25);

        let args = config.into_args();

        assert_eq!(args.frameshift_min_freq, 0.5);
        assert!(args.indel_anchor_depth);
        assert_eq!(args.phased_indel_min_reads, 4);
        assert_eq!(args.phased_indel_min_freq, 0.25);
    }

    #[test]
    fn test_gui_config_indel_knobs_default_to_cli_defaults() {
        // Superseded by `default_parity_tests`, which asks the CLI's own parser
        // what the defaults are instead of writing them down here. This test
        // carried its own copy of them, so when the CLI moved it went on
        // asserting the old values and the drift passed unnoticed.
        let args = minimal_config("/tmp/sample.vcf").into_args();
        assert_eq!(args.phased_indel_min_freq, 0.0);
    }

    #[test]
    fn test_bam_view_balanced_truncation_keeps_partial_reads() {
        let mut reads = Vec::new();
        for i in 0..100 {
            reads.push(bam_read_for_test(&format!("mnv_{i}"), "mnv", i));
        }
        for i in 0..3 {
            reads.push(bam_read_for_test(
                &format!("partial_{i}"),
                "partial",
                1_000 + i,
            ));
        }

        let (display_reads, truncated) = select_bam_reads_for_display(reads, 10);

        assert!(truncated);
        assert_eq!(display_reads.len(), 10);
        assert!(
            display_reads.iter().any(|read| read.support == "partial"),
            "truncated BAM viewer reads should still include partial-support examples"
        );
    }

    #[test]
    fn test_bam_view_untruncated_reads_keep_support_order() {
        let reads = vec![
            bam_read_for_test("ref", "reference", 30),
            bam_read_for_test("partial", "partial", 20),
            bam_read_for_test("mnv", "mnv", 10),
            bam_read_for_test("other", "other", 40),
        ];

        let (display_reads, truncated) = select_bam_reads_for_display(reads, 10);
        let supports: Vec<&str> = display_reads
            .iter()
            .map(|read| read.support.as_str())
            .collect();

        assert!(!truncated);
        assert_eq!(supports, vec!["mnv", "partial", "reference", "other"]);
    }

    #[test]
    fn test_bam_view_columns_include_insertion_slots() {
        let sites = vec![BamVariantSite {
            position: 11,
            reference_base: "C".to_string(),
            alt_base: "CTT".to_string(),
        }];
        let mut insertion_widths = HashMap::new();
        let (anchor, width) = expected_insertion_after(11, "C", "CTT").expect("expected insertion");
        insertion_widths.insert(anchor, width);

        let columns = build_bam_view_columns(10, 12, "ACG", &sites, &insertion_widths);
        let keys: Vec<&str> = columns.iter().map(|column| column.key.as_str()).collect();

        assert_eq!(
            keys,
            vec!["ref:10", "ref:11", "ins:11:1", "ins:11:2", "ref:12"]
        );
        assert_eq!(columns[2].label, "+1");
        assert!(columns[2].is_variant);
    }

    #[test]
    fn test_bam_view_expands_reads_into_insertion_slots() {
        let columns =
            build_bam_view_columns(10, 12, "ACG", &[], &HashMap::from([(11_u64, 2_usize)]));
        let read = RawBamReadView {
            name: "read_1".to_string(),
            strand: "+".to_string(),
            support: "mnv".to_string(),
            start: 10,
            end: 12,
            mapq: 60,
            reference_bases: HashMap::from([
                (10, "A".to_string()),
                (11, "C".to_string()),
                (12, "G".to_string()),
            ]),
            insertions_after: HashMap::from([(11, vec!["T".to_string(), "T".to_string()])]),
        };

        assert_eq!(
            expand_read_bases(&read, &columns),
            vec![
                "A".to_string(),
                "C".to_string(),
                "T".to_string(),
                "T".to_string(),
                "G".to_string()
            ]
        );
    }

    #[test]
    fn test_bam_view_classifies_insertion_support_from_layout() {
        let sites = vec![BamVariantSite {
            position: 11,
            reference_base: "C".to_string(),
            alt_base: "CTT".to_string(),
        }];
        let reference_bases = HashMap::from([(11, "C".to_string())]);
        let alt_insertions = HashMap::from([(11, vec!["T".to_string(), "T".to_string()])]);
        let reference_insertions: HashMap<u64, Vec<String>> = HashMap::new();

        assert_eq!(
            classify_layout_support(&reference_bases, &alt_insertions, &sites),
            "mnv"
        );
        assert_eq!(
            classify_layout_support(&reference_bases, &reference_insertions, &sites),
            "reference"
        );
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

    /// Regenerates the browser demo's read-pileup fixture from the example data
    /// shipped in `example/`, so the demo shows a real `get_bam_view` answer
    /// rather than a hand-written imitation of one.
    ///
    /// Ignored by default because it writes into the working tree. Run it after
    /// changing `BamViewResponse`:
    ///
    /// ```text
    /// cargo test -p get-mnv-gui --bins -- --ignored regenerate_demo_bam_view
    /// ```
    ///
    /// The locus is the `Rv2036` codon that `docs/getting-started.md` walks
    /// through: two substitutions in one codon, carried by all 24 reads. The
    /// request mirrors what `BamViewer.tsx` sends, including its 80 bp padding
    /// and the form's MAPQ and base-quality floors.
    #[test]
    #[ignore = "writes frontend/demo/sample_bam_view.json"]
    fn regenerate_demo_bam_view_fixture() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/..");
        let response = get_bam_view(BamViewRequest {
            bam_path: format!("{root}/example/G35894.demo.bam"),
            fasta_path: format!("{root}/example/MTB_ancestor.fas"),
            chrom: "MTB_anc".to_string(),
            positions: vec![2_282_376, 2_282_377],
            ref_bases: vec!["T".to_string(), "T".to_string()],
            alt_bases: vec!["C".to_string(), "C".to_string()],
            min_mapq: Some(20),
            min_base_quality: Some(20),
            max_reads: Some(80),
            window_start: 2_282_296,
            window_end: 2_282_457,
        })
        .expect("get_bam_view failed on the example data");

        assert!(
            response.total_reads > 0,
            "the example BAM produced no reads at the Rv2036 locus"
        );
        // Pretty-printed: it is checked in, so a regeneration should produce a
        // diff a reviewer can read.
        let path = format!("{root}/frontend/demo/sample_bam_view.json");
        std::fs::write(&path, serde_json::to_string_pretty(&response).unwrap())
            .unwrap_or_else(|e| panic!("cannot write {path}: {e}"));
        println!("wrote {path} ({} reads)", response.total_reads);
    }
}

#[cfg(test)]
mod default_parity_tests {
    use super::tests::minimal_config;
    use clap::Parser;

    /// Every knob the front end leaves unset must fall back to the CLI's own
    /// default. These drifted apart once before, and the app quietly answered
    /// differently from the CLI on the same inputs: a frameshift gate of 0.0
    /// against 0.5, the legacy indel depth denominator, and a phased-haplotype
    /// floor of one read against two.
    ///
    /// This covers the fallback layer only. The React front end sends a value
    /// for nearly every field, so the defaults a user actually gets are the ones
    /// in `frontend/src/types.ts`; those are checked by
    /// `shipped_frontend_defaults_match_the_cli_except_where_intended`.
    #[test]
    fn unset_gui_fields_fall_back_to_the_cli_defaults() {
        let cli = get_mnv::cli::Args::parse_from([
            "get_mnv",
            "--vcf",
            "/tmp/v.vcf",
            "--fasta",
            "/tmp/ref.fasta",
            "--gff",
            "/tmp/ref.gff",
        ]);
        let gui = minimal_config("/tmp/v.vcf").into_args();

        assert_eq!(gui.frameshift_min_freq, cli.frameshift_min_freq);
        assert_eq!(gui.indel_anchor_depth, cli.indel_anchor_depth);
        assert_eq!(gui.phased_indel_min_reads, cli.phased_indel_min_reads);
        assert_eq!(gui.phased_indel_min_freq, cli.phased_indel_min_freq);
        assert_eq!(gui.count_mates_separately, cli.count_mates_separately);
        assert_eq!(gui.min_quality, cli.min_quality);
        assert_eq!(gui.min_mapq, cli.min_mapq);
        assert_eq!(gui.translation_table, cli.translation_table);
        assert_eq!(gui.min_snp_reads, cli.min_snp_reads);
        assert_eq!(gui.min_mnv_reads, cli.min_mnv_reads);
    }

    /// Reads `DEFAULT_CONFIG` out of the front end and compares it with the CLI.
    ///
    /// This is the check the fallback test above cannot make. The form ships
    /// with a value in nearly every field, so `into_args` never reaches its
    /// `unwrap_or`, and a test that only exercises the fallback path will pass
    /// while the app a user opens answers differently from `get_mnv` on the same
    /// files. The app is deliberately more conservative in a few places;
    /// `INTENDED_DIFFERENCES` is that list, and anything outside it is drift.
    ///
    /// The list is checked in both directions: an entry that no longer differs
    /// also fails, so aligning a default forces the note to be removed rather
    /// than left behind to mislead the next reader.
    const INTENDED_DIFFERENCES: &[(&str, &str)] = &[
        (
            "minMapq",
            "the form asks for MAPQ 20, so a run started by clicking Run \
                     ignores multi-mapping reads the CLI would count",
        ),
        ("minSnpReads", "the form asks for 2 supporting reads"),
        ("minMnvReads", "the form asks for 2 supporting reads"),
        (
            "normalizeAlleles",
            "the form trims shared REF/ALT context by default",
        ),
        (
            "splitMultiallelic",
            "the form splits multiallelic records rather than \
                               failing the run in front of the user",
        ),
    ];

    fn frontend_defaults() -> std::collections::HashMap<String, String> {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../frontend/src/types.ts");
        let source =
            std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
        let start = source
            .find("DEFAULT_CONFIG: AnalysisConfig = {")
            .expect("DEFAULT_CONFIG not found in types.ts");
        let body = &source[start..];
        let end = body.find("\n};").expect("unterminated DEFAULT_CONFIG");
        let mut out = std::collections::HashMap::new();
        for line in body[..end].lines().skip(1) {
            let line = line.trim().trim_end_matches(',');
            if let Some((key, value)) = line.split_once(':') {
                out.insert(key.trim().to_string(), value.trim().to_string());
            }
        }
        assert!(
            out.len() > 20,
            "parsed only {} fields out of DEFAULT_CONFIG; the literal's shape changed",
            out.len()
        );
        out
    }

    #[test]
    fn shipped_frontend_defaults_match_the_cli_except_where_intended() {
        let cli = get_mnv::cli::Args::parse_from([
            "get_mnv",
            "--vcf",
            "/tmp/v.vcf",
            "--fasta",
            "/tmp/ref.fasta",
            "--gff",
            "/tmp/ref.gff",
        ]);
        let gui = frontend_defaults();

        let numeric: Vec<(&str, f64)> = vec![
            ("minQuality", cli.min_quality as f64),
            ("minMapq", cli.min_mapq as f64),
            ("minSnpReads", cli.min_snp_reads as f64),
            ("minMnvReads", cli.min_mnv_reads as f64),
            ("minSnpFrequency", cli.min_snp_frequency),
            ("minMnvFrequency", cli.min_mnv_frequency),
            ("minSnpStrandReads", cli.min_snp_strand_reads as f64),
            ("minMnvStrandReads", cli.min_mnv_strand_reads as f64),
            ("minStrandBiasP", cli.min_strand_bias_p),
            ("frameshiftMinFreq", cli.frameshift_min_freq),
            ("phasedIndelMinReads", cli.phased_indel_min_reads as f64),
            ("phasedIndelMinFreq", cli.phased_indel_min_freq),
            ("translationTable", cli.translation_table as f64),
        ];
        let boolean: Vec<(&str, bool)> = vec![
            ("indelAnchorDepth", cli.indel_anchor_depth),
            ("normalizeAlleles", cli.normalize_alleles),
            ("splitMultiallelic", cli.split_multiallelic),
            ("strict", cli.strict),
            ("emitFiltered", cli.emit_filtered),
            ("strandBiasInfo", cli.strand_bias_info),
            ("keepOriginalInfo", cli.keep_original_info),
            ("excludeIntergenic", cli.exclude_intergenic),
        ];

        let mut differing = Vec::new();
        for (key, expected) in &numeric {
            let raw = gui
                .get(*key)
                .unwrap_or_else(|| panic!("{key} missing from DEFAULT_CONFIG"));
            let actual: f64 = raw
                .parse()
                .unwrap_or_else(|_| panic!("{key} is not a number in types.ts: {raw}"));
            if (actual - expected).abs() > f64::EPSILON {
                differing.push((*key, format!("form {actual}, CLI {expected}")));
            }
        }
        for (key, expected) in &boolean {
            let raw = gui
                .get(*key)
                .unwrap_or_else(|| panic!("{key} missing from DEFAULT_CONFIG"));
            let actual = raw == "true";
            if actual != *expected {
                differing.push((*key, format!("form {actual}, CLI {expected}")));
            }
        }

        let intended: Vec<&str> = INTENDED_DIFFERENCES.iter().map(|(k, _)| *k).collect();
        for (key, detail) in &differing {
            assert!(
                intended.contains(key),
                "{key} drifted away from the CLI ({detail}). Either restore the CLI \
                 default or add it to INTENDED_DIFFERENCES with the reason, and say \
                 so in docs/gui.md."
            );
        }
        for key in &intended {
            assert!(
                differing.iter().any(|(k, _)| k == key),
                "{key} is listed in INTENDED_DIFFERENCES but now matches the CLI; \
                 drop the entry and the note in docs/gui.md."
            );
        }
    }
}
