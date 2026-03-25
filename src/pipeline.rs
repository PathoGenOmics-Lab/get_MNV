//! Pipeline orchestration: input validation, per-contig parallel processing,
//! variant emission, summary/manifest generation, and progress reporting.

use crate::cli::Args;
use crate::error::{AppError, AppResult, ErrorCode};
use crate::io::{self, AnnotationFormat, ReferenceMap, VcfPosition};
use crate::output;
use crate::read_count::{self, ReadCountSummary, RegionObservationCache};
use crate::variants::{self, Gene, VariantInfo, VariantType};
use log::info;
use rayon::prelude::*;
use rust_htslib::bam::IndexedReader;
use serde::Serialize;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Read as IoRead, Write};
use std::path::Path;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Debug, Clone, Default, Serialize)]
pub struct ContigSummary {
    pub contig: String,
    pub snp_records_in_vcf: usize,
    pub mapped_genes: usize,
    pub produced_variants: usize,
    pub snp_variants: usize,
    pub mnv_variants: usize,
    pub snp_mnv_variants: usize,
    pub indel_variants: usize,
    pub intergenic_variants: usize,
    pub region_cache_hits: usize,
    pub region_cache_misses: usize,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct GlobalSummary {
    pub contig_count: usize,
    pub snp_records_in_vcf: usize,
    pub mapped_genes: usize,
    pub produced_variants: usize,
    pub snp_variants: usize,
    pub mnv_variants: usize,
    pub snp_mnv_variants: usize,
    pub indel_variants: usize,
    pub intergenic_variants: usize,
    pub region_cache_hits: usize,
    pub region_cache_misses: usize,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct RunSummary {
    pub schema_version: String,
    pub sample: Option<String>,
    pub dry_run: bool,
    pub bam_provided: bool,
    pub inputs: RunInputs,
    pub output_tsv: Option<String>,
    pub output_vcf: Option<String>,
    pub output_bcf: Option<String>,
    pub contigs: Vec<ContigSummary>,
    pub timings: RunTimings,
    pub global: GlobalSummary,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct RunTimings {
    pub parse_inputs_ms: f64,
    pub process_ms: f64,
    pub emit_ms: f64,
    pub total_ms: f64,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct RunInputs {
    pub vcf: String,
    pub fasta: String,
    pub annotation: String,
    pub bam: Option<String>,
    pub checksums: InputChecksums,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct InputChecksums {
    pub vcf_sha256: String,
    pub fasta_sha256: String,
    pub annotation_sha256: String,
    pub bam_sha256: Option<String>,
}

/// Progress event emitted during pipeline execution (for desktop GUI).
#[derive(Debug, Clone, Serialize)]
pub struct ProgressEvent {
    pub phase: String,
    pub contig: Option<String>,
    pub current: usize,
    pub total: usize,
}

#[derive(Debug)]
struct ParsedInputs {
    base_name: String,
    references: ReferenceMap,
    snp_by_contig: HashMap<String, Vec<VcfPosition>>,
    contigs: Vec<String>,
    command_line: String,
    preloaded_gff: Option<HashMap<String, Vec<Gene>>>,
    original_info_headers: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct RegionCacheKey {
    contig: String,
    start: usize,
    end: usize,
    positions: Vec<usize>,
    min_mapq: u8,
}

/// Simple LRU cache backed by a `Vec` of entries sorted by access order.
/// Eviction is O(1) (pop front) instead of O(n) scan. Lookup is O(n) which
/// is acceptable for small capacities (e.g. 64). For larger sizes, consider
/// the `lru` crate.
struct SimpleLruCache<K, V> {
    capacity: usize,
    /// Entries ordered from least-recently-used (front) to most-recently-used (back).
    entries: Vec<(K, V)>,
}

impl<K, V> SimpleLruCache<K, V>
where
    K: Eq + std::hash::Hash + Clone,
    V: Clone,
{
    fn new(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(1),
            entries: Vec::with_capacity(capacity.max(1)),
        }
    }

    fn get_cloned(&mut self, key: &K) -> Option<V> {
        if let Some(idx) = self.entries.iter().position(|(k, _)| k == key) {
            // Move to back (most recently used)
            let entry = self.entries.remove(idx);
            let value = entry.1.clone();
            self.entries.push(entry);
            Some(value)
        } else {
            None
        }
    }

    fn insert(&mut self, key: K, value: V) {
        // If key exists, remove old entry first
        if let Some(idx) = self.entries.iter().position(|(k, _)| *k == key) {
            self.entries.remove(idx);
        }
        // Evict LRU (front) if at capacity
        if self.entries.len() >= self.capacity {
            self.entries.remove(0);
        }
        self.entries.push((key, value));
    }
}

struct WorkerState {
    bam: Option<IndexedReader>,
    region_cache: SimpleLruCache<RegionCacheKey, RegionObservationCache>,
}

#[derive(Default)]
struct WorkerResult {
    variants: Vec<VariantInfo>,
    region_cache_hits: usize,
    region_cache_misses: usize,
}

fn sort_variants(variants: &mut [VariantInfo]) {
    variants.sort_by(|a, b| {
        let a_min = a.positions.iter().copied().min().unwrap_or(usize::MAX);
        let b_min = b.positions.iter().copied().min().unwrap_or(usize::MAX);
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| a_min.cmp(&b_min))
            .then_with(|| a.positions.cmp(&b.positions))
            .then_with(|| a.variant_type.cmp(&b.variant_type))
            .then_with(|| a.gene.cmp(&b.gene))
    });
}

fn configure_threads(threads: Option<usize>) -> AppResult<()> {
    if let Some(threads) = threads {
        if threads == 0 {
            return Err(AppError::config("--threads must be >= 1"));
        }
        // build_global() succeeds only once per process; ignore
        // "already initialized" on subsequent runs in the desktop app.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }
    Ok(())
}

fn selected_contigs(
    args: &Args,
    snp_by_contig: &HashMap<String, Vec<VcfPosition>>,
) -> AppResult<Vec<String>> {
    let mut contigs = if let Some(chrom_arg) = args.chrom.as_ref() {
        chrom_arg
            .split(',')
            .map(str::trim)
            .filter(|c| !c.is_empty())
            .map(ToOwned::to_owned)
            .collect::<Vec<_>>()
    } else {
        snp_by_contig.keys().cloned().collect::<Vec<_>>()
    };

    contigs.sort();
    contigs.dedup();

    if contigs.is_empty() {
        return Err(AppError::validation("No contigs selected for processing"));
    }

    Ok(contigs)
}

fn validate_contig_inputs(
    contigs: &[String],
    references: &ReferenceMap,
    snp_by_contig: &HashMap<String, Vec<VcfPosition>>,
    annotation_format: AnnotationFormat,
) -> AppResult<()> {
    if matches!(annotation_format, AnnotationFormat::Tsv) && contigs.len() > 1 {
        return Err(AppError::validation(
            "TSV annotation does not include contig names; for multi-contig VCF use --gff or restrict with --chrom",
        ));
    }

    let missing_in_vcf = contigs
        .iter()
        .filter(|contig| !snp_by_contig.contains_key(*contig))
        .cloned()
        .collect::<Vec<_>>();
    let missing_in_fasta = contigs
        .iter()
        .filter(|contig| !references.contains_key(*contig))
        .cloned()
        .collect::<Vec<_>>();

    if !missing_in_vcf.is_empty() || !missing_in_fasta.is_empty() {
        return Err(AppError::validation(format!(
            "Contig validation failed. Missing in VCF: [{}]. Missing in FASTA: [{}].",
            if missing_in_vcf.is_empty() {
                "none".to_string()
            } else {
                missing_in_vcf.join(", ")
            },
            if missing_in_fasta.is_empty() {
                "none".to_string()
            } else {
                missing_in_fasta.join(", ")
            }
        )));
    }

    Ok(())
}

fn sanitized_command_line() -> String {
    let command_line_args = std::env::args()
        .skip(1)
        .map(|arg| {
            if arg.contains('/') {
                Path::new(&arg)
                    .file_name()
                    .and_then(|value| value.to_str())
                    .map(ToOwned::to_owned)
                    .unwrap_or(arg)
            } else {
                arg
            }
        })
        .collect::<Vec<_>>();
    format!("get_mnv {}", command_line_args.join(" "))
}

fn parse_inputs(args: &Args, sample_override: Option<&str>) -> AppResult<ParsedInputs> {
    let base_name = io::get_base_name(&args.vcf_file).map_err(reclassify_generic_as_validation)?;
    let references =
        io::load_references(&args.fasta_file).map_err(reclassify_generic_as_validation)?;
    let snp_by_contig = io::load_vcf_positions_by_contig(
        &args.vcf_file,
        sample_override,
        args.split_multiallelic,
        args.normalize_alleles,
        args.keep_original_info,
    )
    .map_err(reclassify_generic_as_validation)?;
    let annotation_format =
        io::detect_annotation_format(&args.genes_file).map_err(reclassify_generic_as_validation)?;
    let contigs = selected_contigs(args, &snp_by_contig)?;
    validate_strict_original_metrics(&contigs, &snp_by_contig, args.strict)?;
    validate_contig_inputs(&contigs, &references, &snp_by_contig, annotation_format)?;

    let preloaded_gff = if annotation_format == AnnotationFormat::Gff {
        Some(
            io::preload_gff_genes(&args.genes_file, &args.gff_features)
                .map_err(reclassify_generic_as_validation)?,
        )
    } else {
        None
    };

    let original_info_headers = if args.keep_original_info {
        io::extract_original_info_headers(&args.vcf_file)
            .map_err(reclassify_generic_as_validation)?
    } else {
        Vec::new()
    };

    Ok(ParsedInputs {
        base_name,
        references,
        snp_by_contig,
        contigs,
        command_line: sanitized_command_line(),
        preloaded_gff,
        original_info_headers,
    })
}

fn reclassify_generic_as_validation(error: AppError) -> AppError {
    if error.code == ErrorCode::Generic {
        AppError::validation(error.message)
    } else {
        error
    }
}

fn validate_strict_original_metrics(
    contigs: &[String],
    snp_by_contig: &HashMap<String, Vec<VcfPosition>>,
    strict: bool,
) -> AppResult<()> {
    if !strict {
        return Ok(());
    }

    let mut missing_dp = 0usize;
    let mut missing_freq = 0usize;
    let mut examples: Vec<String> = Vec::new();

    for contig in contigs {
        if let Some(positions) = snp_by_contig.get(contig) {
            for site in positions {
                let mut missing_fields: Vec<&str> = Vec::new();
                if site.original_dp.is_none() {
                    missing_dp += 1;
                    missing_fields.push("ODP");
                }
                if site.original_freq.is_none() {
                    missing_freq += 1;
                    missing_fields.push("OFREQ");
                }
                if !missing_fields.is_empty() && examples.len() < 5 {
                    examples.push(format!(
                        "{}:{}({})",
                        contig,
                        site.position,
                        missing_fields.join(",")
                    ));
                }
            }
        }
    }

    if missing_dp > 0 || missing_freq > 0 {
        return Err(AppError::validation(format!(
            "--strict enabled, but original VCF metrics are missing (ODP missing in {} records, OFREQ missing in {} records). First affected sites: {}",
            missing_dp,
            missing_freq,
            if examples.is_empty() {
                "none".to_string()
            } else {
                examples.join(", ")
            }
        )));
    }

    Ok(())
}

fn apply_read_summary(variant: &mut VariantInfo, summary: &ReadCountSummary) {
    variant.snp_reads = Some(summary.snp_counts.clone());
    variant.snp_forward_reads = Some(summary.snp_forward_counts.clone());
    variant.snp_reverse_reads = Some(summary.snp_reverse_counts.clone());
    variant.mnv_reads = Some(summary.mnv_count);
    variant.mnv_forward_reads = Some(summary.mnv_forward_count);
    variant.mnv_reverse_reads = Some(summary.mnv_reverse_count);
    variant.mnv_total_reads = Some(summary.mnv_total_reads);
    variant.total_reads = Some(summary.total_reads.clone());
    variant.total_forward_reads = Some(summary.total_forward_reads.clone());
    variant.total_reverse_reads = Some(summary.total_reverse_reads.clone());
    variant.mnv_total_forward_reads = Some(summary.mnv_total_forward_reads);
    variant.mnv_total_reverse_reads = Some(summary.mnv_total_reverse_reads);
}

fn annotate_variants_for_gene(
    gene: &Gene,
    snp_list: &[VcfPosition],
    reference: &io::Reference,
    contig: &str,
) -> Vec<VariantInfo> {
    variants::get_mnv_variants_for_gene(gene, snp_list, reference, contig)
}

fn count_gene_variant_reads(
    state: &mut WorkerState,
    args: &Args,
    contig: &str,
    gene: &Gene,
    variants: &mut [VariantInfo],
) -> AppResult<(usize, usize)> {
    if args.dry_run || state.bam.is_none() {
        return Ok((0, 0));
    }

    let mut target_positions = variants
        .iter()
        .filter(|variant| variant.variant_type != VariantType::Indel)
        .flat_map(|variant| variant.positions.iter().copied())
        .collect::<Vec<_>>();
    if target_positions.is_empty() {
        return Ok((0, 0));
    }

    target_positions.sort_unstable();
    target_positions.dedup();

    let cache_key = RegionCacheKey {
        contig: contig.to_string(),
        start: gene.start,
        end: gene.end,
        positions: target_positions.clone(),
        min_mapq: args.min_mapq,
    };

    let (cache, cache_hits, cache_misses) =
        if let Some(cached) = state.region_cache.get_cloned(&cache_key) {
            (cached, 1, 0)
        } else {
            let bam = state.bam.as_mut().expect("checked Some above; qed");
            let built = read_count::build_region_observation_cache(
                bam,
                contig,
                gene.start,
                gene.end,
                &target_positions,
                args.min_mapq,
            )
            .map_err(|e| {
                AppError::validation(format!(
                    "Failed building read cache for contig '{}' gene '{}' at interval {}-{}: {}",
                    contig, gene.name, gene.start, gene.end, e
                ))
            })?;

            state.region_cache.insert(cache_key, built.clone());
            (built, 0, 1)
        };

    for variant in variants {
        if variant.variant_type == VariantType::Indel {
            continue;
        }
        let summary = read_count::count_reads_from_cache(
            &cache,
            &variant.positions,
            &variant.base_changes,
            args.min_quality,
        )
        .map_err(|e| {
            AppError::validation(format!(
                "Failed counting reads from cache for contig '{}' gene '{}' at positions {:?}: {}",
                contig, gene.name, variant.positions, e
            ))
        })?;
        apply_read_summary(variant, &summary);
    }

    Ok((cache_hits, cache_misses))
}

fn summarize_contig_variants(
    contig: &str,
    snp_records_in_vcf: usize,
    mapped_genes: usize,
    variants: &[VariantInfo],
    region_cache_hits: usize,
    region_cache_misses: usize,
) -> ContigSummary {
    let mut summary = ContigSummary {
        contig: contig.to_string(),
        snp_records_in_vcf,
        mapped_genes,
        produced_variants: variants.len(),
        region_cache_hits,
        region_cache_misses,
        ..ContigSummary::default()
    };
    for variant in variants {
        if variant.gene == "intergenic" {
            summary.intergenic_variants += 1;
        } else {
            match variant.variant_type {
                VariantType::Snp => summary.snp_variants += 1,
                VariantType::Mnv => summary.mnv_variants += 1,
                VariantType::SnpMnv => summary.snp_mnv_variants += 1,
                VariantType::Indel => summary.indel_variants += 1,
            }
        }
    }
    summary
}

fn process_contig(
    args: &Args,
    contig: &str,
    references: &ReferenceMap,
    snp_by_contig: &HashMap<String, Vec<VcfPosition>>,
    preloaded_gff: Option<&HashMap<String, Vec<Gene>>>,
) -> AppResult<(Vec<VariantInfo>, ContigSummary)> {
    let reference =
        io::reference_for_chrom(references, contig).map_err(reclassify_generic_as_validation)?;
    let snp_list = snp_by_contig
        .get(contig)
        .ok_or_else(|| AppError::validation(format!("Missing VCF data for contig '{}'", contig)))?;
    io::validate_vcf_reference_alleles(contig, snp_list, &reference)
        .map_err(reclassify_generic_as_validation)?;

    let genes = if let Some(gff_genes) = preloaded_gff {
        let all_contig_genes = gff_genes.get(contig).cloned().unwrap_or_default();
        let filtered = io::filter_genes_with_snps(&all_contig_genes, snp_list);
        log::info!(
            "GFF/GFF3 contig '{}': {} gene entries, {} mapped to SNPs, {} without SNPs",
            contig,
            all_contig_genes.len(),
            filtered.len(),
            all_contig_genes.len() - filtered.len()
        );
        filtered
    } else {
        io::load_genes(&args.genes_file, snp_list, Some(contig), &args.gff_features)
            .map_err(reclassify_generic_as_validation)?
    };
    info!(
        "Contig '{}' -> {} SNP/variant records in VCF, {} mapped genes",
        contig,
        snp_list.len(),
        genes.len()
    );

    let bam_path = args.bam_file.clone();
    let should_count_reads = args.bam_file.is_some() && !args.dry_run;
    let worker_results: Result<Vec<WorkerResult>, AppError> = genes
        .par_iter()
        .map_init(
            || -> AppResult<WorkerState> {
                let bam = if should_count_reads {
                    if let Some(path) = bam_path.as_ref() {
                        Some(IndexedReader::from_path(path).map_err(|e| {
                            AppError::validation(format!(
                                "Failed to open BAM '{}' in worker thread: {}",
                                path, e
                            ))
                        })?)
                    } else {
                        None
                    }
                } else {
                    None
                };
                Ok(WorkerState {
                    bam,
                    region_cache: SimpleLruCache::new(64),
                })
            },
            |state_result, gene| -> AppResult<WorkerResult> {
                let state = state_result
                    .as_mut()
                    .map_err(|err| AppError::validation(err.to_string()))?;
                let mut variants = annotate_variants_for_gene(gene, snp_list, &reference, contig);
                if variants.is_empty() {
                    return Ok(WorkerResult::default());
                }
                let (cache_hits, cache_misses) =
                    count_gene_variant_reads(state, args, contig, gene, &mut variants)?;
                Ok(WorkerResult {
                    variants,
                    region_cache_hits: cache_hits,
                    region_cache_misses: cache_misses,
                })
            },
        )
        .collect();

    let mut all_variants = Vec::new();
    let mut cache_hits = 0usize;
    let mut cache_misses = 0usize;
    for result in worker_results.map_err(|e| {
        AppError::validation(format!("Error while processing contig '{}': {}", contig, e))
    })? {
        cache_hits += result.region_cache_hits;
        cache_misses += result.region_cache_misses;
        all_variants.extend(result.variants);
    }

    // Collect intergenic variants (positions not covered by any gene).
    if !args.exclude_intergenic {
        let mut covered: HashSet<usize> = HashSet::with_capacity(snp_list.len());
        for gene in &genes {
            for snp in snp_list {
                if snp.position >= gene.start && snp.position <= gene.end {
                    covered.insert(snp.position);
                }
            }
        }
        let mut intergenic_count = 0usize;
        for snp in snp_list {
            if !covered.contains(&snp.position) {
                all_variants.push(variants::build_intergenic_variant(contig, snp));
                intergenic_count += 1;
            }
        }
        if intergenic_count > 0 {
            info!(
                "Contig '{}' -> {} intergenic variant(s) added",
                contig, intergenic_count
            );
        }
    }

    sort_variants(&mut all_variants);
    let summary = summarize_contig_variants(
        contig,
        snp_list.len(),
        genes.len(),
        &all_variants,
        cache_hits,
        cache_misses,
    );
    Ok((all_variants, summary))
}

fn update_global_summary(global: &mut GlobalSummary, contig_summary: &ContigSummary) {
    global.contig_count += 1;
    global.snp_records_in_vcf += contig_summary.snp_records_in_vcf;
    global.mapped_genes += contig_summary.mapped_genes;
    global.produced_variants += contig_summary.produced_variants;
    global.snp_variants += contig_summary.snp_variants;
    global.mnv_variants += contig_summary.mnv_variants;
    global.snp_mnv_variants += contig_summary.snp_mnv_variants;
    global.indel_variants += contig_summary.indel_variants;
    global.intergenic_variants += contig_summary.intergenic_variants;
    global.region_cache_hits += contig_summary.region_cache_hits;
    global.region_cache_misses += contig_summary.region_cache_misses;
}

fn emit_contig_variants(
    tsv_writer: Option<&mut output::TsvWriter>,
    vcf_writer: Option<&mut output::VcfWriter>,
    variants: &[VariantInfo],
    references: &ReferenceMap,
) -> AppResult<()> {
    if let Some(writer) = tsv_writer {
        writer.write_variants(variants)?;
    }
    if let Some(writer) = vcf_writer {
        writer.write_variants(variants, references)?;
    }
    Ok(())
}

#[cfg(test)]
fn summary_to_json(summary: &RunSummary) -> String {
    serde_json::to_string_pretty(summary).unwrap_or_else(|e| {
        format!(
            "{{\"schema_version\":\"1.0.0\",\"error\":\"Failed to serialize summary: {}\"}}",
            e
        )
    }) + "\n"
}

pub fn write_summary_json(path: &str, summary: &RunSummary) -> AppResult<()> {
    let mut file = File::create(path)?;
    serde_json::to_writer_pretty(&mut file, summary)
        .map_err(|e| AppError::msg(format!("Failed to serialize summary JSON: {}", e)))?;
    file.write_all(b"\n")?;
    Ok(())
}

fn compute_sha256(path: &str) -> AppResult<String> {
    let file = File::open(path).map_err(|e| {
        AppError::validation(format!("Failed to open '{}' for checksum: {}", path, e))
    })?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 8192];
    loop {
        let read_bytes = reader.read(&mut buffer).map_err(|e| {
            AppError::validation(format!("Failed reading '{}' for checksum: {}", path, e))
        })?;
        if read_bytes == 0 {
            break;
        }
        hasher.update(&buffer[..read_bytes]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn build_input_metadata(args: &Args) -> AppResult<RunInputs> {
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

fn sanitize_sample_for_path(sample: &str) -> String {
    sample
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '_' || ch == '-' {
                ch
            } else {
                '_'
            }
        })
        .collect::<String>()
}

fn append_sample_suffix(base_name: &str, sample_suffix: Option<&str>) -> String {
    match sample_suffix {
        Some(sample) => format!("{}.sample_{}", base_name, sanitize_sample_for_path(sample)),
        None => base_name.to_string(),
    }
}

fn build_run_manifest_value(summary: &RunSummary, command_line: &str) -> AppResult<Value> {
    let timestamp_unix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| AppError::msg(format!("Failed to compute current time: {}", e)))?
        .as_secs();

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

fn write_json_value(path: &str, value: &Value) -> AppResult<()> {
    let mut file = File::create(path)?;
    serde_json::to_writer_pretty(&mut file, value)
        .map_err(|e| AppError::msg(format!("Failed to serialize JSON to '{}': {}", path, e)))?;
    file.write_all(b"\n")?;
    Ok(())
}

fn log_toggle(label: &str, enabled: bool) {
    info!(
        "{}: {}",
        label,
        if enabled { "enabled" } else { "disabled" }
    );
}

fn log_run_configuration(args: &Args, sample_override: Option<&str>) {
    info!("VCF file: {}", args.vcf_file);
    if let Some(ref bam) = args.bam_file {
        info!("BAM file: {}", bam);
    } else {
        info!("No BAM file provided. Output will be generated without read count fields.");
    }
    info!("FASTA file: {}", args.fasta_file);
    info!("Gene annotation file: {}", args.genes_file);
    info!("GFF feature types: {}", args.gff_features.join(", "));
    if let Ok(AnnotationFormat::Tsv) = io::detect_annotation_format(&args.genes_file) {
        let default_features = vec!["gene".to_string(), "pseudogene".to_string()];
        if args.gff_features != default_features {
            log::warn!("--gff-features is ignored when using TSV annotation format (--genes)");
        }
    }
    if let Some(sample) = sample_override {
        info!("Target sample for original FORMAT metrics: {}", sample);
    } else {
        info!("Target sample for original FORMAT metrics: first sample");
    }
    if let Some(chrom) = args.chrom.as_ref() {
        info!("Target chromosome(s): {}", chrom);
    } else {
        info!("Target chromosome(s): all contigs in input VCF");
    }
    if let Some(threads) = args.threads {
        info!("Threads: {}", threads);
    } else {
        info!("Threads: Rayon default");
    }
    info!("Minimum Phred quality: {}", args.min_quality);
    info!("Minimum mapping quality (MAPQ): {}", args.min_mapq);
    info!("Minimum SNP reads: {}", args.min_snp_reads);
    info!("Minimum MNV reads: {}", args.min_mnv_reads);
    info!("Minimum strand-bias p-value: {:.4}", args.min_strand_bias_p);
    info!(
        "Minimum SNP reads per strand: {}",
        args.min_snp_strand_reads
    );
    info!(
        "Minimum MNV reads per strand: {}",
        args.min_mnv_strand_reads
    );
    log_toggle("Strict original metrics", args.strict);
    log_toggle("Split multiallelic records", args.split_multiallelic);
    log_toggle("Compressed VCF output (--vcf-gz)", args.vcf_gz);
    log_toggle("Build index for .vcf.gz output", args.index_vcf_gz);
    log_toggle("Strand bias INFO fields", args.strand_bias_info);
    log_toggle("Emit filtered records", args.emit_filtered);
    log_toggle("Exclude intergenic variants", args.exclude_intergenic);
    log_toggle("BCF output", args.bcf);
    log_toggle("Dry run", args.dry_run);
    if args.dry_run && args.bam_file.is_some() {
        info!("Dry-run mode active: BAM read counting is skipped.");
    }
}

fn run_single(
    args: &Args,
    sample_override: Option<&str>,
    sample_suffix: Option<&str>,
    write_reports: bool,
    on_progress: &dyn Fn(ProgressEvent),
) -> AppResult<RunSummary> {
    let total_start = Instant::now();
    log_run_configuration(args, sample_override);

    on_progress(ProgressEvent {
        phase: "parsing".into(),
        contig: None,
        current: 0,
        total: 0,
    });

    let parse_start = Instant::now();
    let parsed = parse_inputs(args, sample_override)?;
    let inputs = build_input_metadata(args)?;
    let parse_inputs_ms = parse_start.elapsed().as_secs_f64() * 1000.0;
    let base_stem = append_sample_suffix(&parsed.base_name, sample_suffix);
    // Allow overriding the filename prefix (desktop app).
    let stem_name = match &args.output_prefix {
        Some(prefix) => prefix.clone(),
        None => base_stem,
    };
    // When output_dir is set (e.g. desktop app), write output files there
    // instead of the process CWD.
    let output_stem = match &args.output_dir {
        Some(dir) => {
            let dir_path = std::path::Path::new(dir);
            dir_path.join(&stem_name).to_string_lossy().into_owned()
        }
        None => stem_name,
    };

    let output_tsv = if args.dry_run {
        None
    } else if args.both || !args.convert {
        Some(format!("{}.MNV.tsv", output_stem))
    } else {
        None
    };
    let output_vcf = if args.dry_run {
        None
    } else if args.both || args.convert {
        Some(if args.vcf_gz {
            format!("{}.MNV.vcf.gz", output_stem)
        } else {
            format!("{}.MNV.vcf", output_stem)
        })
    } else {
        None
    };
    let output_bcf = if args.dry_run || !args.bcf {
        None
    } else {
        Some(format!("{}.MNV.bcf", output_stem))
    };

    let mut tsv_writer = if output_tsv.is_some() {
        Some(output::TsvWriter::new(
            &output_stem,
            args.bam_file.is_some(),
        )?)
    } else {
        None
    };
    let mut vcf_writer = if output_vcf.is_some() {
        Some(output::VcfWriter::new(
            &output_stem,
            args.bam_file.is_some(),
            args.min_snp_reads,
            args.min_mnv_reads,
            args.min_quality,
            args.min_mapq,
            &parsed.command_line,
            &parsed.contigs,
            args.vcf_gz,
            args.min_snp_strand_reads,
            args.min_mnv_strand_reads,
            args.min_strand_bias_p,
            args.emit_filtered,
            args.strand_bias_info,
            &parsed.original_info_headers,
        )?)
    } else {
        None
    };

    let mut process_ms = 0.0f64;
    let mut emit_ms = 0.0f64;
    let mut summary = RunSummary {
        schema_version: "1.0.0".to_string(),
        sample: sample_override.map(|sample| sample.to_string()),
        dry_run: args.dry_run,
        bam_provided: args.bam_file.is_some(),
        inputs,
        output_tsv,
        output_vcf,
        output_bcf,
        ..RunSummary::default()
    };

    let total_contigs = parsed.contigs.len();
    for (contig_idx, contig) in parsed.contigs.iter().enumerate() {
        on_progress(ProgressEvent {
            phase: "processing".into(),
            contig: Some(contig.clone()),
            current: contig_idx + 1,
            total: total_contigs,
        });

        let process_start = Instant::now();
        let (contig_variants, contig_summary) = process_contig(
            args,
            contig,
            &parsed.references,
            &parsed.snp_by_contig,
            parsed.preloaded_gff.as_ref(),
        )?;
        process_ms += process_start.elapsed().as_secs_f64() * 1000.0;

        let emit_start = Instant::now();
        emit_contig_variants(
            tsv_writer.as_mut(),
            vcf_writer.as_mut(),
            &contig_variants,
            &parsed.references,
        )?;
        emit_ms += emit_start.elapsed().as_secs_f64() * 1000.0;

        info!(
            "Summary '{}' -> variants={} (SNP={}, MNV={}, SNP/MNV={}, INDEL={}, intergenic={})",
            contig_summary.contig,
            contig_summary.produced_variants,
            contig_summary.snp_variants,
            contig_summary.mnv_variants,
            contig_summary.snp_mnv_variants,
            contig_summary.indel_variants,
            contig_summary.intergenic_variants
        );
        if !args.dry_run && args.bam_file.is_some() {
            info!(
                "Region cache '{}' -> hits={}, misses={}",
                contig_summary.contig,
                contig_summary.region_cache_hits,
                contig_summary.region_cache_misses
            );
        }

        update_global_summary(&mut summary.global, &contig_summary);
        summary.contigs.push(contig_summary);
    }

    on_progress(ProgressEvent {
        phase: "complete".into(),
        contig: None,
        current: total_contigs,
        total: total_contigs,
    });

    summary.timings.parse_inputs_ms = parse_inputs_ms;
    summary.timings.process_ms = process_ms;
    summary.timings.emit_ms = emit_ms;
    summary.timings.total_ms = total_start.elapsed().as_secs_f64() * 1000.0;

    drop(tsv_writer);
    drop(vcf_writer);

    if args.index_vcf_gz {
        if let Some(vcf_path) = summary.output_vcf.as_deref() {
            output::build_tabix_index(vcf_path)?;
            info!("Built Tabix index for {}", vcf_path);
        }
    }
    if let (Some(vcf_path), Some(bcf_path)) =
        (summary.output_vcf.as_deref(), summary.output_bcf.as_deref())
    {
        output::convert_vcf_to_bcf(vcf_path, bcf_path)?;
        info!("Converted {} to {}", vcf_path, bcf_path);
    }

    info!(
        "Global summary -> contigs={}, VCF records={}, mapped genes={}, emitted variants={} (SNP={}, MNV={}, SNP/MNV={}, INDEL={}, intergenic={})",
        summary.global.contig_count,
        summary.global.snp_records_in_vcf,
        summary.global.mapped_genes,
        summary.global.produced_variants,
        summary.global.snp_variants,
        summary.global.mnv_variants,
        summary.global.snp_mnv_variants,
        summary.global.indel_variants,
        summary.global.intergenic_variants
    );
    if !args.dry_run && args.bam_file.is_some() {
        info!(
            "Global region cache -> hits={}, misses={}",
            summary.global.region_cache_hits, summary.global.region_cache_misses
        );
    }

    if write_reports {
        if let Some(summary_json_path) = args.summary_json.as_deref() {
            write_summary_json(summary_json_path, &summary)?;
            info!("Summary JSON written to {}", summary_json_path);
        }
        if let Some(run_manifest_path) = args.run_manifest.as_deref() {
            let manifest = build_run_manifest_value(&summary, &parsed.command_line)?;
            write_json_value(run_manifest_path, &manifest)?;
            info!("Run manifest written to {}", run_manifest_path);
        }
    }

    if args.dry_run {
        info!("Dry-run completed successfully. No output files were written.");
    } else {
        info!("Processing complete. Output files generated successfully.");
    }

    Ok(summary)
}

/// Run the pipeline (CLI entry point, no progress reporting).
pub fn run(args: &Args) -> AppResult<RunSummary> {
    run_with_progress(args, &|_| {})
}

/// Run the pipeline with a progress callback (used by the desktop GUI).
pub fn run_with_progress(
    args: &Args,
    on_progress: &dyn Fn(ProgressEvent),
) -> AppResult<RunSummary> {
    if args.vcf_gz && !(args.convert || args.both) {
        return Err(AppError::config("--vcf-gz requires --convert or --both"));
    }
    if args.index_vcf_gz && !args.vcf_gz {
        return Err(AppError::config("--index-vcf-gz requires --vcf-gz"));
    }
    if args.bcf && !(args.convert || args.both) {
        return Err(AppError::config("--bcf requires --convert or --both"));
    }
    if args.keep_original_info && !(args.convert || args.both) {
        return Err(AppError::config(
            "--keep-original-info requires --convert or --both",
        ));
    }
    if !(0.0..=1.0).contains(&args.min_strand_bias_p) {
        return Err(AppError::config(
            "--min-strand-bias-p must be between 0 and 1",
        ));
    }

    configure_threads(args.threads)?;

    if args.sample.as_deref() != Some("all") {
        return run_single(args, args.sample.as_deref(), None, true, on_progress);
    }

    let sample_names =
        io::list_vcf_samples(&args.vcf_file).map_err(reclassify_generic_as_validation)?;
    if sample_names.is_empty() {
        return Err(AppError::validation(
            "Requested --sample all but input VCF has no sample columns",
        ));
    }

    let mut sample_summaries = Vec::new();
    for sample in &sample_names {
        info!("Processing sample '{}' in --sample all mode", sample);
        sample_summaries.push(run_single(
            args,
            Some(sample),
            Some(sample),
            false,
            on_progress,
        )?);
    }

    let mut aggregate = RunSummary {
        schema_version: "1.0.0".to_string(),
        sample: Some("all".to_string()),
        dry_run: args.dry_run,
        bam_provided: args.bam_file.is_some(),
        inputs: build_input_metadata(args)?,
        output_tsv: None,
        output_vcf: None,
        output_bcf: None,
        ..RunSummary::default()
    };

    for sample_summary in &sample_summaries {
        aggregate.global.contig_count += sample_summary.global.contig_count;
        aggregate.global.snp_records_in_vcf += sample_summary.global.snp_records_in_vcf;
        aggregate.global.mapped_genes += sample_summary.global.mapped_genes;
        aggregate.global.produced_variants += sample_summary.global.produced_variants;
        aggregate.global.snp_variants += sample_summary.global.snp_variants;
        aggregate.global.mnv_variants += sample_summary.global.mnv_variants;
        aggregate.global.snp_mnv_variants += sample_summary.global.snp_mnv_variants;
        aggregate.global.indel_variants += sample_summary.global.indel_variants;
        aggregate.global.intergenic_variants += sample_summary.global.intergenic_variants;
        aggregate.global.region_cache_hits += sample_summary.global.region_cache_hits;
        aggregate.global.region_cache_misses += sample_summary.global.region_cache_misses;
        aggregate.timings.parse_inputs_ms += sample_summary.timings.parse_inputs_ms;
        aggregate.timings.process_ms += sample_summary.timings.process_ms;
        aggregate.timings.emit_ms += sample_summary.timings.emit_ms;
        aggregate.timings.total_ms += sample_summary.timings.total_ms;
    }

    if let Some(summary_json_path) = args.summary_json.as_deref() {
        let payload = json!({
            "schema_version": "1.0.0",
            "mode": "sample_all",
            "sample_count": sample_summaries.len(),
            "sample_names": sample_names,
            "aggregate": aggregate,
            "samples": sample_summaries
        });
        write_json_value(summary_json_path, &payload)?;
        info!("Summary JSON written to {}", summary_json_path);
    }
    if let Some(run_manifest_path) = args.run_manifest.as_deref() {
        let timestamp_unix = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|e| AppError::msg(format!("Failed to compute current time: {}", e)))?
            .as_secs();
        let payload = json!({
            "schema_version": "1.0.0",
            "mode": "sample_all",
            "timestamp_unix": timestamp_unix,
            "tool_version": env!("CARGO_PKG_VERSION"),
            "command_line": sanitized_command_line(),
            "sample_names": sample_names,
            "aggregate": aggregate,
            "samples": sample_summaries
        });
        write_json_value(run_manifest_path, &payload)?;
        info!("Run manifest written to {}", run_manifest_path);
    }

    Ok(aggregate)
}

#[cfg(test)]
mod tests {
    use super::{
        summary_to_json, ContigSummary, GlobalSummary, InputChecksums, RunInputs, RunSummary,
        RunTimings,
    };

    #[test]
    fn test_summary_to_json_contains_expected_keys() {
        let summary = RunSummary {
            schema_version: "1.0.0".to_string(),
            sample: None,
            dry_run: true,
            bam_provided: false,
            inputs: RunInputs {
                vcf: "in.vcf".to_string(),
                fasta: "ref.fasta".to_string(),
                annotation: "genes.gff3".to_string(),
                bam: None,
                checksums: InputChecksums {
                    vcf_sha256: "vcf".to_string(),
                    fasta_sha256: "fasta".to_string(),
                    annotation_sha256: "ann".to_string(),
                    bam_sha256: None,
                },
            },
            output_tsv: None,
            output_vcf: Some("out.vcf".to_string()),
            output_bcf: None,
            contigs: vec![ContigSummary {
                contig: "chr1".to_string(),
                snp_records_in_vcf: 3,
                mapped_genes: 1,
                produced_variants: 3,
                snp_variants: 2,
                mnv_variants: 1,
                snp_mnv_variants: 0,
                indel_variants: 0,
                intergenic_variants: 0,
                region_cache_hits: 0,
                region_cache_misses: 1,
            }],
            timings: RunTimings {
                parse_inputs_ms: 1.0,
                process_ms: 2.0,
                emit_ms: 3.0,
                total_ms: 6.0,
            },
            global: GlobalSummary {
                contig_count: 1,
                snp_records_in_vcf: 3,
                mapped_genes: 1,
                produced_variants: 3,
                snp_variants: 2,
                mnv_variants: 1,
                snp_mnv_variants: 0,
                indel_variants: 0,
                intergenic_variants: 0,
                region_cache_hits: 0,
                region_cache_misses: 1,
            },
        };

        let json = summary_to_json(&summary);
        assert!(json.contains("\"dry_run\": true"));
        assert!(json.contains("\"schema_version\": \"1.0.0\""));
        assert!(json.contains("\"output_vcf\": \"out.vcf\""));
        assert!(json.contains("\"contig\": \"chr1\""));
        assert!(json.contains("\"timings\""));
        assert!(json.contains("\"parse_inputs_ms\": 1.0"));
        assert!(json.contains("\"inputs\""));
        assert!(json.contains("\"global\""));
    }
}
