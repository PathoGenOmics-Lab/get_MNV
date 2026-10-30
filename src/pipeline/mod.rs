//! Pipeline orchestration: input validation, per-contig parallel processing,
//! variant emission, summary/manifest generation, and progress reporting.

pub mod config;
pub mod manifest;
pub(crate) mod output_paths;
pub mod processing;
pub(crate) mod shared;
pub mod summary;

// Re-export public types.
pub use summary::{
    write_summary_json, ContigSummary, GlobalSummary, InputChecksums, ProgressEvent, RunInputs,
    RunSummary, RunTimings,
};

use config::{
    configure_threads, log_run_configuration, sanitize_sample_for_path, sanitized_command_line,
};
use manifest::{
    build_input_metadata, build_run_manifest_value, output_checksums_for, write_json_value,
};
use processing::{
    emit_contig_variants, parse_inputs, process_contig, reclassify_generic_as_validation,
    resolve_variant_input_format, SampleRun,
};
use shared::SharedInputs;
use summary::update_global_summary;

use crate::cli::{Args, VariantInputFormat};
use crate::error::{AppError, AppResult, IoResultExt};
use crate::io;
use crate::output;
use log::info;
use serde_json::json;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

fn run_single(
    args: &Args,
    sample_override: Option<&str>,
    sample_suffix: Option<&str>,
    write_reports: bool,
    on_progress: &dyn Fn(ProgressEvent),
    run: SampleRun<'_>,
) -> AppResult<RunSummary> {
    let total_start = Instant::now();
    let genetic_code =
        crate::genetic_code::GeneticCode::new(args.translation_table).ok_or_else(|| {
            AppError::config(format!(
                "Unsupported translation table: {}",
                args.translation_table
            ))
        })?;
    info!("Using genetic code: {genetic_code}");
    log_run_configuration(args, sample_override);

    on_progress(ProgressEvent {
        phase: "parsing".into(),
        contig: None,
        current: 0,
        total: 0,
    });

    let parse_start = Instant::now();
    let parsed = parse_inputs(args, sample_override, run.shared)?;
    let parse_inputs_ms = parse_start.elapsed().as_secs_f64() * 1000.0;
    // SHA-256 checksums are deferred until after parsing to keep parse_inputs_ms
    // measuring only I/O + validation, not hashing.
    let needs_checksums = args.summary_json.is_some() || args.run_manifest.is_some();
    let inputs = match run.shared {
        // The same files for every sample of a cohort, hashed once.
        Some(shared) if needs_checksums => shared.input_metadata(args)?,
        _ => build_input_metadata(args, needs_checksums)?,
    };
    let paths = output_paths::OutputPaths::resolve(args, &parsed.base_name, sample_suffix);

    // Fail fast on a missing/un-indexed BAM before any output file is created.
    if !args.dry_run {
        if let Some(bam_path) = args.bam_file.as_deref() {
            config::validate_bam(bam_path, &args.fasta_file)?;
        }
    }

    // The same fail-fast idea for --output-dir: create it here, before the run,
    // rather than discovering at the moment of writing that the directory does
    // not exist. On a cohort that moment arrives after every sample has been
    // annotated, and the work is lost for the sake of a missing folder.
    //
    // create_dir_all rather than create_dir, so `--output-dir results/run3`
    // works without the caller having made `results` first, and it succeeds
    // quietly when the directory is already there.
    if !args.dry_run {
        if let Some(dir) = args.output_dir.as_deref() {
            std::fs::create_dir_all(dir).with_path(dir)?;
        }
    }

    // Transactional output: if the run errors after the output files are created,
    // remove the partial files on unwind so downstream tooling never consumes a
    // truncated VCF/TSV. Disarmed once the outputs are fully written.
    struct OutputCleanup {
        paths: Vec<String>,
        armed: bool,
    }
    impl Drop for OutputCleanup {
        fn drop(&mut self) {
            if self.armed {
                for path in &self.paths {
                    let _ = std::fs::remove_file(path);
                }
            }
        }
    }
    let mut output_cleanup = OutputCleanup {
        paths: [paths.tsv.clone(), paths.vcf.clone(), paths.bcf.clone()]
            .into_iter()
            .flatten()
            .collect(),
        armed: true,
    };

    let mut tsv_writer = if paths.tsv.is_some() {
        Some(output::TsvWriter::new(output::TsvWriterConfig {
            filename: &paths.output_stem,
            bam_provided: args.bam_file.is_some(),
            min_snp_reads: args.min_snp_reads,
            min_snp_frequency: args.min_snp_frequency,
            min_mnv_reads: args.min_mnv_reads,
            min_mnv_frequency: args.min_mnv_frequency,
            min_snp_strand_reads: args.min_snp_strand_reads,
            min_mnv_strand_reads: args.min_mnv_strand_reads,
        })?)
    } else {
        None
    };
    let mut vcf_writer = if paths.vcf.is_some() {
        Some(output::VcfWriter::new(output::VcfWriterConfig {
            filename: &paths.output_stem,
            bam_provided: args.bam_file.is_some(),
            min_snp_reads: args.min_snp_reads,
            min_snp_frequency: args.min_snp_frequency,
            min_mnv_reads: args.min_mnv_reads,
            min_mnv_frequency: args.min_mnv_frequency,
            min_quality: args.min_quality,
            min_mapq: args.min_mapq,
            command_line: &parsed.command_line,
            contigs: &parsed.contigs,
            bgzf_output: args.vcf_gz,
            min_snp_strand_reads: args.min_snp_strand_reads,
            min_mnv_strand_reads: args.min_mnv_strand_reads,
            min_strand_bias_p: args.min_strand_bias_p,
            emit_filtered: args.emit_filtered,
            include_strand_bias_info: args.strand_bias_info,
            original_info_headers: &parsed.original_info_headers,
        })?)
    } else {
        None
    };

    let mut process_ms = 0.0f64;
    let mut emit_ms = 0.0f64;
    let mut summary = RunSummary {
        schema_version: "1.0.0".to_string(),
        sample: sample_override.map(std::string::ToString::to_string),
        dry_run: args.dry_run,
        bam_provided: args.bam_file.is_some(),
        translation_table: args.translation_table,
        inputs,
        output_tsv: paths.tsv,
        output_vcf: paths.vcf,
        output_bcf: paths.bcf,
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
        let (contig_variants, contig_summary) =
            process_contig(args, contig, &parsed, genetic_code, run)?;
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

    if let Some(writer) = vcf_writer.as_mut() {
        writer.flush()?;
    }
    drop(tsv_writer);
    drop(vcf_writer);
    // Outputs are fully written; keep them even if a post-write step (tabix/BCF
    // index, report JSON) fails.
    output_cleanup.armed = false;

    if args.index_vcf_gz {
        if let Some(vcf_path) = summary.output_vcf.as_deref() {
            if output::build_tabix_index(vcf_path)? == output::ExternalStep::Ran {
                info!("Built Tabix index for {vcf_path}");
            }
        }
    }
    if let (Some(vcf_path), Some(bcf_path)) =
        (summary.output_vcf.clone(), summary.output_bcf.clone())
    {
        if output::convert_vcf_to_bcf(&vcf_path, &bcf_path)? == output::ExternalStep::Ran {
            info!("Converted {vcf_path} to {bcf_path}");
        } else {
            // bcftools was not there, so the file was never written. Anything
            // downstream that trusts this path has to see the absence: the
            // summary JSON would otherwise name a file that does not exist, and
            // --run-manifest hashes every output the summary names.
            summary.output_bcf = None;
        }
    }

    info!(
        "Global summary -> contigs={}, VCF records={}, mapped genes={}, produced variants={} (SNP={}, MNV={}, SNP/MNV={}, INDEL={}, intergenic={})",
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
            info!("Summary JSON written to {summary_json_path}");
        }
        if let Some(run_manifest_path) = args.run_manifest.as_deref() {
            let manifest = build_run_manifest_value(&summary, &parsed.command_line)?;
            write_json_value(run_manifest_path, &manifest)?;
            info!("Run manifest written to {run_manifest_path}");
        }
        // Single-sample report. In `--sample all` mode `write_reports` is false
        // for each sample and one combined report is written by the caller.
        if args.report.is_some() || args.matrix.is_some() {
            let sample_label = summary.sample.clone().unwrap_or_else(|| {
                std::path::Path::new(&parsed.base_name)
                    .file_name()
                    .map(|name| name.to_string_lossy().into_owned())
                    .unwrap_or_else(|| "sample".to_string())
            });
            let tsv_path = summary.output_tsv.as_deref().ok_or_else(|| {
                AppError::config("--report and --matrix need the TSV output; use --both instead of --convert, or drop --convert")
            })?;
            let mut builder = output::ReportBuilder::new();
            builder.add_tsv(tsv_path, Some(&sample_label))?;
            write_report_outputs(args, &builder, &parsed.command_line)?;
        }
    }

    if args.dry_run {
        info!("Dry-run completed successfully. No output files were written.");
    } else {
        info!("Processing complete. Output files generated successfully.");
    }

    Ok(summary)
}

/// Civil date from a count of days since the Unix epoch (Howard Hinnant's
/// `civil_from_days`). Avoids pulling in a date dependency for the one
/// human-readable timestamp the report header shows.
fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let month = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if month <= 2 { year + 1 } else { year }, month, day)
}

/// Current UTC time as `YYYY-MM-DD HH:MM UTC`, for the report header.
fn report_timestamp() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    format!(
        "{year:04}-{month:02}-{day:02} {:02}:{:02} UTC",
        rem / 3600,
        (rem % 3600) / 60
    )
}

/// Write whichever of the HTML report and the variant matrix were asked for,
/// both from the same parsed TSVs.
fn write_report_outputs(
    args: &Args,
    builder: &output::ReportBuilder,
    command_line: &str,
) -> AppResult<()> {
    let samples = builder.sample_count();
    if let Some(report_path) = args.report.as_deref() {
        builder.write(report_path, command_line, &report_timestamp())?;
        info!(
            "HTML report written to {report_path} ({} variants across {samples} sample(s))",
            builder.len()
        );
    }
    if let Some(matrix_path) = args.matrix.as_deref() {
        let variants = builder.write_matrix(matrix_path)?;
        info!(
            "Variant matrix written to {matrix_path} ({variants} variants x {samples} sample(s))"
        );
    }
    Ok(())
}

/// Build the HTML report, the variant matrix or both from existing get_MNV TSV
/// files, without running the pipeline. Each input file becomes one sample.
fn run_report_from_tsvs(args: &Args) -> AppResult<RunSummary> {
    if args.report.is_none() && args.matrix.is_none() {
        return Err(AppError::config(
            "--report-from needs --report, --matrix or both to name an output",
        ));
    }

    // Label the whole input set at once: file names alone collide when a
    // per-sample pipeline writes `results/<sample>/variants.MNV.tsv`.
    let labels = output::sample_labels(&args.report_from);
    let mut builder = output::ReportBuilder::new();
    for (tsv_path, label) in args.report_from.iter().zip(labels.iter()) {
        builder.add_tsv(tsv_path, Some(label))?;
        info!("Report input: {tsv_path} (sample '{label}')");
    }
    if builder.is_empty() {
        return Err(AppError::validation(
            "The provided TSV files contain no variant rows, so there is nothing to report",
        ));
    }
    write_report_outputs(args, &builder, &sanitized_command_line())?;

    Ok(RunSummary {
        schema_version: "1.0.0".to_string(),
        ..RunSummary::default()
    })
}

/// What one sample of a cohort said and did, held until the samples before it
/// have been written.
struct SampleOutcome {
    index: usize,
    held: Vec<crate::logging::Held>,
    /// Each progress event with the number of held lines that came before it.
    progress: Vec<(usize, ProgressEvent)>,
    summary: AppResult<RunSummary>,
}

impl SampleOutcome {
    /// Write the sample's lines and report its progress, in the order the
    /// sample produced them.
    fn release(self, on_progress: &dyn Fn(ProgressEvent)) -> AppResult<RunSummary> {
        let mut written = 0;
        for (position, event) in self.progress {
            for line in &self.held[written..position] {
                line.write();
            }
            written = position;
            on_progress(event);
        }
        for line in &self.held[written..] {
            line.write();
        }
        self.summary
    }
}

/// Annotate every sample of a `--sample all` cohort, and return their
/// summaries in sample order, with the inputs they shared if they shared any.
///
/// The samples used to run one after another, each reading every input again.
/// They now share the inputs and run side by side, one sample per thread, and
/// each sample's log lines and progress events are held until the samples
/// before it have been written, so the log and the progress read exactly as
/// they did. Both lean on the logger [`crate::logging::install`] sets up;
/// without it the samples run one after another as before. A cohort with a BAM
/// and fewer samples than threads also goes one sample at a time, since its
/// cost is in counting reads gene by gene, which spreads over the threads
/// better than its few samples do.
fn run_samples(
    args: &Args,
    sample_names: &[String],
    on_progress: &dyn Fn(ProgressEvent),
) -> AppResult<(Vec<RunSummary>, Option<SharedInputs>)> {
    if !crate::logging::can_hold() {
        let summaries = samples_one_at_a_time(args, sample_names, on_progress, SampleRun::alone())?;
        return Ok((summaries, None));
    }
    let shared = SharedInputs::new(args);
    let threads = rayon::current_num_threads();
    let workers = threads.min(sample_names.len());
    let summaries = if workers < 2 || (args.bam_file.is_some() && sample_names.len() < threads) {
        let run = SampleRun {
            shared: Some(&shared),
            parallel_genes: true,
        };
        samples_one_at_a_time(args, sample_names, on_progress, run)?
    } else {
        samples_side_by_side(args, sample_names, on_progress, &shared, workers)?
    };
    Ok((summaries, Some(shared)))
}

fn samples_one_at_a_time(
    args: &Args,
    sample_names: &[String],
    on_progress: &dyn Fn(ProgressEvent),
    run: SampleRun<'_>,
) -> AppResult<Vec<RunSummary>> {
    let mut summaries = Vec::with_capacity(sample_names.len());
    for sample in sample_names {
        info!("Processing sample '{sample}' in --sample all mode");
        summaries.push(run_single(
            args,
            Some(sample),
            Some(sample),
            false,
            on_progress,
            run,
        )?);
    }
    Ok(summaries)
}

/// Run the samples on `workers` threads, a sample at a time on each, and
/// write what each said in sample order as soon as the samples before it have
/// been written.
///
/// Samples are started in order. When one fails, no sample after it is
/// started, the samples before it finish, and its error is the run's, as it
/// was when they ran one after another. Samples after it that were already
/// running finish and keep their output.
fn samples_side_by_side(
    args: &Args,
    sample_names: &[String],
    on_progress: &dyn Fn(ProgressEvent),
    shared: &SharedInputs,
    workers: usize,
) -> AppResult<Vec<RunSummary>> {
    let next = AtomicUsize::new(0);
    let first_failure = AtomicUsize::new(usize::MAX);
    std::thread::scope(|scope| -> AppResult<Vec<RunSummary>> {
        let (sender, outcomes) = std::sync::mpsc::channel::<SampleOutcome>();
        let mut started = 0;
        for _ in 0..workers {
            let sender = sender.clone();
            let (next, first_failure) = (&next, &first_failure);
            let worker = std::thread::Builder::new()
                .stack_size(config::WORKER_STACK_SIZE)
                .spawn_scoped(scope, move || loop {
                    let index = next.fetch_add(1, Ordering::SeqCst);
                    if index >= sample_names.len() || index > first_failure.load(Ordering::SeqCst) {
                        break;
                    }
                    let sample = sample_names[index].as_str();
                    let progress = std::cell::RefCell::new(Vec::new());
                    let report = |event: ProgressEvent| {
                        progress
                            .borrow_mut()
                            .push((crate::logging::held_so_far(), event));
                    };
                    let (summary, held) = crate::logging::holding(|| {
                        info!("Processing sample '{sample}' in --sample all mode");
                        run_single(
                            args,
                            Some(sample),
                            Some(sample),
                            false,
                            &report,
                            SampleRun {
                                shared: Some(shared),
                                parallel_genes: false,
                            },
                        )
                    });
                    if summary.is_err() {
                        first_failure.fetch_min(index, Ordering::SeqCst);
                    }
                    let outcome = SampleOutcome {
                        index,
                        held,
                        progress: progress.into_inner(),
                        summary,
                    };
                    if sender.send(outcome).is_err() {
                        break;
                    }
                });
            if worker.is_err() {
                break;
            }
            started += 1;
        }
        drop(sender);
        if started == 0 {
            let run = SampleRun {
                shared: Some(shared),
                parallel_genes: true,
            };
            return samples_one_at_a_time(args, sample_names, on_progress, run);
        }

        let mut waiting = std::collections::BTreeMap::new();
        let mut summaries = Vec::with_capacity(sample_names.len());
        for (index, sample) in sample_names.iter().enumerate() {
            let outcome = loop {
                if let Some(outcome) = waiting.remove(&index) {
                    break outcome;
                }
                match outcomes.recv() {
                    Ok(outcome) => {
                        waiting.insert(outcome.index, outcome);
                    }
                    // Every worker has stopped and this sample never came back,
                    // which only a worker that panicked leaves behind. Leaving
                    // the scope passes its panic on.
                    Err(_) => {
                        return Err(AppError::msg(format!(
                            "Sample '{sample}' was not processed"
                        )))
                    }
                }
            };
            summaries.push(outcome.release(on_progress)?);
        }
        Ok(summaries)
    })
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
    // Standalone report mode: aggregate existing TSVs, no pipeline run.
    if !args.report_from.is_empty() {
        return run_report_from_tsvs(args);
    }
    // Enforced here rather than by the clap ArgGroup, which cannot express
    // "unless --report-from" (that mode needs no variant input).
    if args.vcf_file.is_none() && args.tsv_file.is_none() {
        return Err(AppError::config(
            "one of --vcf or --tsv is required (or use --report-from to build a report from existing TSV files)",
        ));
    }
    if args.report.is_some() || args.matrix.is_some() {
        if args.dry_run {
            return Err(AppError::config(
                "--report and --matrix cannot be used with --dry-run",
            ));
        }
        if args.convert && !args.both {
            return Err(AppError::config(
                "--report and --matrix need the TSV output; use --both instead of --convert, or drop --convert",
            ));
        }
    }
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
    if !(0.0..=1.0).contains(&args.min_snp_frequency) {
        return Err(AppError::config(
            "--min-snp-frequency must be between 0 and 1",
        ));
    }
    if !(0.0..=1.0).contains(&args.min_mnv_frequency) {
        return Err(AppError::config(
            "--min-mnv-frequency must be between 0 and 1",
        ));
    }
    if !(0.0..=1.0).contains(&args.frameshift_min_freq) {
        return Err(AppError::config(
            "--frameshift-min-freq must be between 0 and 1 (it is an allele frequency, not a percentage)",
        ));
    }
    if !(0.0..=1.0).contains(&args.phased_indel_min_freq) {
        return Err(AppError::config(
            "--phased-indel-min-freq must be between 0 and 1 (it is an allele frequency, not a percentage)",
        ));
    }
    if (args.min_snp_frequency > 0.0 || args.min_mnv_frequency > 0.0) && args.bam_file.is_none() {
        return Err(AppError::config(
            "--min-snp-frequency and --min-mnv-frequency require --bam because frequencies are calculated from read support",
        ));
    }
    if crate::genetic_code::GeneticCode::new(args.translation_table).is_none() {
        return Err(AppError::config(format!(
            "--translation-table {} is not supported. Supported tables: {:?}",
            args.translation_table,
            crate::genetic_code::SUPPORTED_TABLES
        )));
    }

    configure_threads(args.threads)?;

    let input_format = resolve_variant_input_format(args)?;
    if input_format == VariantInputFormat::Tsv && args.sample.is_some() {
        return Err(AppError::config(
            "--sample is only supported for VCF input; TSV input is treated as a single sample",
        ));
    }

    if args.sample.as_deref() != Some("all") {
        return run_single(
            args,
            args.sample.as_deref(),
            None,
            true,
            on_progress,
            SampleRun::alone(),
        );
    }

    let variant_file = args.variant_file();
    let sample_names = if io::vcf_fast::use_fast_parser(variant_file) {
        io::vcf_fast::list_text_vcf_samples(variant_file)
            .map_err(reclassify_generic_as_validation)?
    } else {
        io::list_vcf_samples(variant_file).map_err(reclassify_generic_as_validation)?
    };
    if sample_names.is_empty() {
        return Err(AppError::validation(
            "Requested --sample all but input VCF has no sample columns",
        ));
    }

    // Every sample writes to a file named after it, with the characters a file
    // name cannot hold replaced. Two samples that come out the same, whether the
    // header repeats a name or `a/b` and `a_b` both sanitize to `a_b`, would
    // write to one file and the second would overwrite the first, leaving a
    // cohort short one sample with nothing said about it.
    let mut stems: std::collections::HashMap<String, &String> = std::collections::HashMap::new();
    for sample in &sample_names {
        let stem = sanitize_sample_for_path(sample);
        if let Some(previous) = stems.insert(stem.clone(), sample) {
            return Err(AppError::validation(format!(
                "Samples '{previous}' and '{sample}' both name their output 'sample_{stem}', so one \
                 would overwrite the other. Rename one of them in the VCF header, or process them \
                 one at a time with --sample."
            )));
        }
    }

    let (sample_summaries, shared) = run_samples(args, &sample_names, on_progress)?;

    let mut aggregate = RunSummary {
        schema_version: "1.0.0".to_string(),
        sample: Some("all".to_string()),
        dry_run: args.dry_run,
        bam_provided: args.bam_file.is_some(),
        translation_table: args.translation_table,
        inputs: match &shared {
            Some(shared) => shared.input_metadata(args)?,
            None => build_input_metadata(args, true)?,
        },
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

    // One combined report and matrix covering every sample of the multi-sample VCF.
    if args.report.is_some() || args.matrix.is_some() {
        let mut builder = output::ReportBuilder::new();
        for (sample, sample_summary) in sample_names.iter().zip(sample_summaries.iter()) {
            let tsv_path = sample_summary.output_tsv.as_deref().ok_or_else(|| {
                AppError::config("--report and --matrix need the TSV output; use --both instead of --convert, or drop --convert")
            })?;
            builder.add_tsv(tsv_path, Some(sample))?;
        }
        write_report_outputs(args, &builder, &sanitized_command_line())?;
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
        info!("Summary JSON written to {summary_json_path}");
    }
    if let Some(run_manifest_path) = args.run_manifest.as_deref() {
        let timestamp_unix = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs();
        // Each sample carries the checksums of the files it wrote. Listing the
        // paths without them, as this payload used to, left `--run-manifest`
        // promising checksums and delivering none in the one mode that writes
        // more than one output pair.
        let mut samples = Vec::with_capacity(sample_summaries.len());
        for sample_summary in &sample_summaries {
            let mut value = serde_json::to_value(sample_summary)?;
            if let Some(object) = value.as_object_mut() {
                object.insert(
                    "output_checksums".to_string(),
                    output_checksums_for(sample_summary)?,
                );
            }
            samples.push(value);
        }
        let payload = json!({
            "schema_version": "1.0.0",
            "mode": "sample_all",
            "timestamp_unix": timestamp_unix,
            "tool_version": env!("CARGO_PKG_VERSION"),
            "command_line": sanitized_command_line(),
            "sample_names": sample_names,
            "aggregate": aggregate,
            "samples": samples
        });
        write_json_value(run_manifest_path, &payload)?;
        info!("Run manifest written to {run_manifest_path}");
    }

    Ok(aggregate)
}

#[cfg(test)]
mod tests {
    use super::summary::summary_to_json;
    use super::{ContigSummary, GlobalSummary, InputChecksums, RunInputs, RunSummary, RunTimings};

    /// A cohort run side by side must say, report and write exactly what the
    /// same cohort run one sample at a time does: the same files byte for
    /// byte, the same log lines in the same order and the same progress
    /// events, whatever order its samples happen to finish in.
    #[test]
    fn a_cohort_side_by_side_matches_one_sample_at_a_time() {
        use super::{
            samples_one_at_a_time, samples_side_by_side, ProgressEvent, SampleRun, SharedInputs,
        };
        use clap::Parser;
        use std::cell::RefCell;
        use std::collections::BTreeMap;

        crate::logging::tests::install_for_tests();
        let root = std::env::temp_dir().join(format!("get_mnv_cohort_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        let example = concat!(env!("CARGO_MANIFEST_DIR"), "/example");

        // Eight samples over the example's first 300 records: genotypes that
        // differ from sample to sample, and a depth and frequency that some
        // samples give and others leave missing.
        let source = std::fs::read_to_string(format!("{example}/G35894.var.snp.vcf")).unwrap();
        let samples: Vec<String> = (0..8).map(|i| format!("S{i}")).collect();
        let mut vcf = String::new();
        for line in source.lines().filter(|line| line.starts_with("##")) {
            vcf.push_str(line);
            vcf.push('\n');
        }
        vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
        vcf.push_str(&samples.join("\t"));
        vcf.push('\n');
        for (row, line) in source
            .lines()
            .filter(|line| !line.starts_with('#'))
            .take(300)
            .enumerate()
        {
            let fixed: Vec<&str> = line.split('\t').take(8).collect();
            vcf.push_str(&fixed.join("\t"));
            vcf.push_str("\tGT:DP:FREQ");
            for sample in 0..samples.len() {
                let carried = (row * 7 + sample * 3) % 5 != 0;
                if sample % 3 == 2 {
                    vcf.push_str(&format!("\t{}:.:.", u8::from(carried)));
                } else {
                    vcf.push_str(&format!(
                        "\t{}:{}:{}%",
                        u8::from(carried),
                        10 + row % 40,
                        20 + sample * 9
                    ));
                }
            }
            vcf.push('\n');
        }
        let vcf_path = root.join("cohort.vcf");
        std::fs::write(&vcf_path, vcf).unwrap();

        let run = |label: &str, side_by_side: bool| {
            let out = root.join(label);
            let args = crate::cli::Args::parse_from([
                "get_MNV".to_string(),
                "--vcf".to_string(),
                vcf_path.to_string_lossy().into_owned(),
                "--fasta".to_string(),
                format!("{example}/MTB_ancestor.fas"),
                "--genes".to_string(),
                format!("{example}/anot_genes.txt"),
                "--sample".to_string(),
                "all".to_string(),
                "--both".to_string(),
                "--output-dir".to_string(),
                out.to_string_lossy().into_owned(),
            ]);
            std::fs::create_dir_all(&out).unwrap();
            let events = RefCell::new(Vec::new());
            // Each event is logged as well, so the lines also show where every
            // event fell among them.
            let on_progress = |event: ProgressEvent| {
                let said = format!(
                    "{} {:?} {}/{}",
                    event.phase, event.contig, event.current, event.total
                );
                log::info!(target: "progress", "{said}");
                events.borrow_mut().push(said);
            };
            let before = crate::logging::tests::written_here().len();
            let summaries = if side_by_side {
                let shared = SharedInputs::new(&args);
                samples_side_by_side(&args, &samples, &on_progress, &shared, 3)
            } else {
                samples_one_at_a_time(&args, &samples, &on_progress, SampleRun::alone())
            }
            .expect("the cohort runs");
            let lines: Vec<String> = crate::logging::tests::written_here()[before..]
                .iter()
                .map(|line| line.replace(&out.to_string_lossy().into_owned(), "<out>"))
                .collect();
            let files: BTreeMap<String, Vec<u8>> = std::fs::read_dir(&out)
                .unwrap()
                .map(|entry| {
                    let entry = entry.unwrap();
                    (
                        entry.file_name().to_string_lossy().into_owned(),
                        std::fs::read(entry.path()).unwrap(),
                    )
                })
                .collect();
            let counts: Vec<usize> = summaries
                .iter()
                .map(|summary| summary.global.produced_variants)
                .collect();
            (lines, events.into_inner(), files, counts)
        };

        let (lines, events, files, counts) = run("one_at_a_time", false);
        let (side_lines, side_events, side_files, side_counts) = run("side_by_side", true);
        assert_eq!(files.len(), 2 * samples.len());
        assert!(
            lines.len() > 10 * samples.len(),
            "only {} lines",
            lines.len()
        );
        assert_eq!(side_lines, lines);
        assert_eq!(side_events, events);
        assert_eq!(side_counts, counts);
        assert_eq!(
            side_files.keys().collect::<Vec<_>>(),
            files.keys().collect::<Vec<_>>()
        );
        for (name, bytes) in &files {
            assert!(side_files[name] == *bytes, "{name} differs");
        }
        std::fs::remove_dir_all(&root).ok();
    }

    #[test]
    fn test_summary_to_json_contains_expected_keys() {
        let summary = RunSummary {
            schema_version: "1.0.0".to_string(),
            sample: None,
            dry_run: true,
            bam_provided: false,
            translation_table: 11,
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
