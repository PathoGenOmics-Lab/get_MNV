//! Pipeline orchestration: input validation, per-contig parallel processing,
//! variant emission, summary/manifest generation, and progress reporting.

pub mod config;
pub mod manifest;
pub mod processing;
pub mod summary;

// Re-export public types.
pub use summary::{
    write_summary_json, ContigSummary, GlobalSummary, InputChecksums, ProgressEvent, RunInputs,
    RunSummary, RunTimings,
};

use config::{
    append_sample_suffix, configure_threads, log_run_configuration, sanitized_command_line,
};
use manifest::{build_input_metadata, build_run_manifest_value, write_json_value};
use processing::{
    emit_contig_variants, parse_inputs, process_contig, reclassify_generic_as_validation,
};
use summary::update_global_summary;

use crate::cli::Args;
use crate::error::{AppError, AppResult};
use crate::io;
use crate::output;
use log::info;
use serde_json::json;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

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
    let parse_inputs_ms = parse_start.elapsed().as_secs_f64() * 1000.0;
    // SHA-256 checksums are deferred until after parsing to keep parse_inputs_ms
    // measuring only I/O + validation, not hashing.
    let inputs = build_input_metadata(args)?;
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
        Some(format!("{output_stem}.MNV.tsv"))
    } else {
        None
    };
    let output_vcf = if args.dry_run {
        None
    } else if args.both || args.convert {
        Some(if args.vcf_gz {
            format!("{output_stem}.MNV.vcf.gz")
        } else {
            format!("{output_stem}.MNV.vcf")
        })
    } else {
        None
    };
    let output_bcf = if args.dry_run || !args.bcf {
        None
    } else {
        Some(format!("{output_stem}.MNV.bcf"))
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
        Some(output::VcfWriter::new(output::VcfWriterConfig {
            filename: &output_stem,
            bam_provided: args.bam_file.is_some(),
            min_snp_reads: args.min_snp_reads,
            min_mnv_reads: args.min_mnv_reads,
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
            info!("Built Tabix index for {vcf_path}");
        }
    }
    if let (Some(vcf_path), Some(bcf_path)) =
        (summary.output_vcf.as_deref(), summary.output_bcf.as_deref())
    {
        output::convert_vcf_to_bcf(vcf_path, bcf_path)?;
        info!("Converted {vcf_path} to {bcf_path}");
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
            info!("Summary JSON written to {summary_json_path}");
        }
        if let Some(run_manifest_path) = args.run_manifest.as_deref() {
            let manifest = build_run_manifest_value(&summary, &parsed.command_line)?;
            write_json_value(run_manifest_path, &manifest)?;
            info!("Run manifest written to {run_manifest_path}");
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

    let sample_names = if io::vcf_fast::use_fast_parser(&args.vcf_file) {
        io::vcf_fast::list_text_vcf_samples(&args.vcf_file)
            .map_err(reclassify_generic_as_validation)?
    } else {
        io::list_vcf_samples(&args.vcf_file).map_err(reclassify_generic_as_validation)?
    };
    if sample_names.is_empty() {
        return Err(AppError::validation(
            "Requested --sample all but input VCF has no sample columns",
        ));
    }

    let mut sample_summaries = Vec::new();
    for sample in &sample_names {
        info!("Processing sample '{sample}' in --sample all mode");
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
        info!("Summary JSON written to {summary_json_path}");
    }
    if let Some(run_manifest_path) = args.run_manifest.as_deref() {
        let timestamp_unix = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs();
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
        info!("Run manifest written to {run_manifest_path}");
    }

    Ok(aggregate)
}

#[cfg(test)]
mod tests {
    use super::summary::summary_to_json;
    use super::{ContigSummary, GlobalSummary, InputChecksums, RunInputs, RunSummary, RunTimings};

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
