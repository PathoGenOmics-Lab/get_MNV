//! Tauri commands for the get_MNV desktop application.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::Header;
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
#[allow(dead_code)]
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
    pub exclude_intergenic: Option<bool>,
    pub gff_features: Option<String>,
    pub output_dir: Option<String>,
    pub convert: Option<bool>,
    pub both: Option<bool>,
    pub summary_json: Option<bool>,
    pub run_manifest: Option<bool>,
    pub error_json: Option<bool>,
    pub translation_table: Option<u8>,
}

impl AnalysisConfig {
    pub fn into_args(self) -> Vec<String> {
        let mut args = vec![
            "get_mnv".to_string(),
            "--vcf".to_string(),
            self.vcf_file,
            "--fasta".to_string(),
            self.fasta_file,
        ];

        // Determine annotation type by extension
        let genes_lower = self.genes_file.to_lowercase();
        if genes_lower.ends_with(".gff")
            || genes_lower.ends_with(".gff3")
            || genes_lower.ends_with(".gtf")
        {
            args.push("--gff".to_string());
        } else {
            args.push("--tsv-genes".to_string());
        }
        args.push(self.genes_file);

        if let Some(bam) = self.bam_file {
            args.push("--bam".to_string());
            args.push(bam);
        }
        if let Some(sample) = self.sample {
            if !sample.is_empty() {
                args.push("--sample".to_string());
                args.push(sample);
            }
        }
        if let Some(chrom) = self.chrom {
            if !chrom.is_empty() {
                args.push("--chrom".to_string());
                args.push(chrom);
            }
        }
        if let Some(true) = self.normalize_alleles {
            args.push("--normalize-alleles".to_string());
        }
        if let Some(q) = self.min_quality {
            args.push("--quality".to_string());
            args.push(q.to_string());
        }
        if let Some(m) = self.min_mapq {
            args.push("--mapq".to_string());
            args.push(m.to_string());
        }
        if let Some(t) = self.threads {
            args.push("--threads".to_string());
            args.push(t.to_string());
        }
        if let Some(r) = self.min_snp_reads {
            args.push("--min-snp-reads".to_string());
            args.push(r.to_string());
        }
        if let Some(r) = self.min_mnv_reads {
            args.push("--min-mnv-reads".to_string());
            args.push(r.to_string());
        }
        if let Some(s) = self.min_snp_strand_reads {
            args.push("--min-snp-strand".to_string());
            args.push(s.to_string());
        }
        if let Some(s) = self.min_mnv_strand_reads {
            args.push("--min-mnv-strand".to_string());
            args.push(s.to_string());
        }
        if let Some(p) = self.min_strand_bias_p {
            if p > 0.0 {
                args.push("--min-strand-bias-p".to_string());
                args.push(p.to_string());
            }
        }
        if let Some(true) = self.dry_run {
            args.push("--dry-run".to_string());
        }
        if let Some(true) = self.strict {
            args.push("--strict".to_string());
        }
        if let Some(true) = self.split_multiallelic {
            args.push("--split-multiallelic".to_string());
        }
        if let Some(true) = self.emit_filtered {
            args.push("--emit-filtered".to_string());
        }
        if let Some(true) = self.vcf_gz {
            args.push("--vcf-gz".to_string());
        }
        if let Some(true) = self.index_vcf_gz {
            args.push("--index-vcf-gz".to_string());
        }
        if let Some(true) = self.strand_bias_info {
            args.push("--strand-bias-info".to_string());
        }
        if let Some(true) = self.keep_original_info {
            args.push("--keep-original-info".to_string());
        }
        if let Some(true) = self.exclude_intergenic {
            args.push("--exclude-intergenic".to_string());
        }
        if let Some(ref features) = self.gff_features {
            if !features.is_empty() {
                args.push("--gff-features".to_string());
                args.push(features.clone());
            }
        }
        if let Some(ref dir) = self.output_dir {
            if !dir.is_empty() {
                args.push("--output-dir".to_string());
                args.push(dir.clone());
            }
        }

        // Output mode: --both or --convert
        let both = self.both.unwrap_or(false);
        let convert = self.convert.unwrap_or(false);
        if both {
            args.push("--both".to_string());
        } else if convert {
            args.push("--convert".to_string());
        }

        if let Some(true) = self.summary_json {
            args.push("--summary-json".to_string());
        }
        if let Some(true) = self.run_manifest {
            args.push("--run-manifest".to_string());
        }
        if let Some(true) = self.error_json {
            args.push("--error-json".to_string());
        }
        if let Some(t) = self.translation_table {
            if t != 11 {
                args.push("--translation-table".to_string());
                args.push(t.to_string());
            }
        }

        args
    }
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AnalysisResult {
    pub success: bool,
    pub message: String,
    pub output_tsv: Option<String>,
    pub output_vcf: Option<String>,
    pub output_bcf: Option<String>,
    pub summary_json: Option<String>,
    pub run_manifest: Option<String>,
    pub error_json: Option<String>,
}

/// Progress event payload.
#[derive(Debug, Clone, Serialize)]
pub struct ProgressEvent {
    pub stage: String,
    pub message: String,
    pub percent: Option<f64>,
}

#[tauri::command]
pub fn run_analysis(
    config: AnalysisConfig,
    window: tauri::Window,
) -> Result<AnalysisResult, String> {
    let args_vec = config.into_args();
    let args_str: Vec<&str> = args_vec.iter().map(|s| s.as_str()).collect();

    let _ = window.emit("analysis-progress", ProgressEvent {
        stage: "parsing".to_string(),
        message: "Parsing command line arguments...".to_string(),
        percent: Some(5.0),
    });

    let args = get_mnv::cli::Args::try_parse_from(&args_str)
        .map_err(|e| format!("Invalid arguments: {e}"))?;

    let _ = window.emit("analysis-progress", ProgressEvent {
        stage: "running".to_string(),
        message: "Running analysis...".to_string(),
        percent: Some(10.0),
    });

    match get_mnv::run_pipeline(args) {
        Ok(summary) => {
            let _ = window.emit("analysis-progress", ProgressEvent {
                stage: "complete".to_string(),
                message: "Analysis complete!".to_string(),
                percent: Some(100.0),
            });

            Ok(AnalysisResult {
                success: true,
                message: format!(
                    "Processed {} contigs: {} variants ({} MNVs)",
                    summary.contigs_processed,
                    summary.total_variants,
                    summary.total_mnvs
                ),
                output_tsv: summary.output_tsv,
                output_vcf: summary.output_vcf,
                output_bcf: summary.output_bcf,
                summary_json: summary.summary_json,
                run_manifest: summary.run_manifest,
                error_json: None,
            })
        }
        Err(e) => {
            let _ = window.emit("analysis-progress", ProgressEvent {
                stage: "error".to_string(),
                message: format!("Analysis failed: {e}"),
                percent: None,
            });

            Ok(AnalysisResult {
                success: false,
                message: format!("Error: {e}"),
                output_tsv: None,
                output_vcf: None,
                output_bcf: None,
                summary_json: None,
                run_manifest: None,
                error_json: None,
            })
        }
    }
}

// ---------------------------------------------------------------------------
// BAM Viewer
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BamViewRequest {
    pub bam_path: String,
    pub fasta_path: String,
    pub chrom: String,
    pub window_start: u64,
    pub window_end: u64,
    pub positions: Vec<u64>,
    pub ref_bases: Vec<String>,
    pub alt_bases: Vec<String>,
    pub min_mapq: Option<u8>,
    pub min_base_quality: Option<u8>,
    pub max_reads: Option<usize>,
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
fn base_at_ref_pos_noodles(
    record: &bam::Record,
    ref_pos: u64,
    min_bq: u8,
) -> Option<String> {
    let read_start = record.alignment_start()
        .and_then(|p| p.ok())
        .map(|p| { let v: usize = p.into(); v as u64 - 1 }) // 0-based
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
fn read_fasta_region(fasta_path: &str, chrom: &str, start: u64, end: u64) -> Result<String, String> {
    // Read .fai index to find the sequence
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
            seq_len = fields[1].parse().map_err(|_| "Bad .fai length".to_string())?;
            offset = fields[2].parse().map_err(|_| "Bad .fai offset".to_string())?;
            line_bases = fields[3].parse().map_err(|_| "Bad .fai line_bases".to_string())?;
            line_width = fields[4].parse().map_err(|_| "Bad .fai line_width".to_string())?;
            found = true;
            break;
        }
    }
    if !found {
        return Err(format!("Chromosome '{}' not found in FASTA index", chrom));
    }

    let actual_end = end.min(seq_len);
    let window_len = (actual_end - start) as usize;

    // Calculate byte offsets
    let start_line = start / line_bases;
    let start_byte = offset + start_line * line_width + (start % line_bases);

    let end_line = actual_end / line_bases;
    let end_byte = offset + end_line * line_width + (actual_end % line_bases);

    let bytes_to_read = (end_byte - start_byte) as usize;

    use std::io::{Read, Seek, SeekFrom};
    let mut fasta_file = std::fs::File::open(fasta_path)
        .map_err(|e| format!("Cannot open FASTA: {}", e))?;
    fasta_file.seek(SeekFrom::Start(start_byte))
        .map_err(|e| format!("FASTA seek error: {}", e))?;

    let mut buf = vec![0u8; bytes_to_read + 100]; // extra for newlines
    let n = fasta_file.read(&mut buf)
        .map_err(|e| format!("FASTA read error: {}", e))?;

    let seq: String = buf[..n].iter()
        .filter(|&&b| b != b'\n' && b != b'\r')
        .map(|&b| b.to_ascii_uppercase() as char)
        .take(window_len)
        .collect();

    // Pad if needed
    let mut result = seq;
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

    // Build variant site lookup
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

    // Read reference sequence
    let window_len = (request.window_end - request.window_start) as usize;
    let reference = read_fasta_region(
        &request.fasta_path,
        &request.chrom,
        request.window_start,
        request.window_end,
    )?;

    // Open BAM with noodles
    let mut bam_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&request.bam_path)
        .map_err(|e| format!("Cannot open BAM: {}", e))?;
    let header: Header = bam_reader.read_header()
        .map_err(|e| format!("Cannot read BAM header: {}", e))?;

    // Query region (1-based inclusive)
    let region_str = format!("{}:{}-{}", request.chrom, request.window_start + 1, request.window_end);
    let region: noodles::core::Region = region_str.parse()
        .map_err(|e| format!("Invalid region '{}': {}", region_str, e))?;
    let mut query = bam_reader.query(&header, &region)
        .map_err(|e| format!("BAM query error: {}", e))?;

    let mut all_reads: Vec<BamReadView> = Vec::new();
    let mut counts = BamSupportCounts::default();
    let mut total_reads: usize = 0;

    let mut record = bam::Record::default();
    while query.read_record(&mut record).map_err(|e| format!("BAM read error: {}", e))? != 0 {
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_qc_fail()
            || flags.is_duplicate()
        {
            continue;
        }

        let mapq = record.mapping_quality()
            .map(|q: noodles::sam::alignment::record::MappingQuality| q.get())
            .unwrap_or(255);
        if mapq < min_mapq {
            continue;
        }

        total_reads += 1;

        // Extract bases at each window position
        let mut bases: Vec<String> = Vec::with_capacity(window_len);
        for offset in 0..window_len {
            let rp = request.window_start + offset as u64;
            match base_at_ref_pos_noodles(&record, rp, min_bq) {
                Some(b) => bases.push(b),
                None => bases.push(" ".to_string()),
            }
        }

        // Classify read support
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

        match support.as_str() {
            "mnv" => counts.mnv += 1,
            "partial" => counts.partial += 1,
            "reference" => counts.reference += 1,
            _ => counts.other += 1,
        }

        let strand = if flags.is_reverse_complemented() { "-" } else { "+" }.to_string();

        let read_name = record.name()
            .map(|n| String::from_utf8_lossy(<_ as AsRef<[u8]>>::as_ref(&n)).to_string())
            .unwrap_or_else(|| "?".to_string());

        // Compute read end from alignment start + aligned length
        let read_start = record.alignment_start()
            .and_then(|p| p.ok())
            .map(|p| { let v: usize = p.into(); v as u64 - 1 })
            .unwrap_or(0);
        let aligned_len: u64 = record.cigar().iter()
            .filter_map(|r| r.ok())
            .map(|op| match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
                | Kind::Deletion | Kind::Skip => op.len() as u64,
                _ => 0,
            })
            .sum();
        let read_end = read_start + aligned_len;

        all_reads.push(BamReadView {
            name: read_name,
            strand,
            support,
            start: read_start,
            end: read_end,
            mapq,
            bases,
        });
    }

    counts.total = total_reads;

    // Per-position coverage from ALL reads
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
