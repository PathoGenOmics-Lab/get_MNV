//! Pipeline configuration: thread setup, contig selection, input validation,
//! and run logging.

use crate::cli::Args;
use crate::error::{AppError, AppResult};
use crate::io::{self, AnnotationFormat, ReferenceMap, VcfPosition};
use log::info;
use std::collections::HashMap;
use std::path::Path;

pub(crate) fn configure_threads(threads: Option<usize>) -> AppResult<()> {
    if let Some(threads) = threads {
        if threads == 0 {
            return Err(AppError::config("--threads must be >= 1"));
        }
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }
    Ok(())
}

pub(crate) fn selected_contigs(
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

pub(crate) fn validate_contig_inputs(
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

pub(crate) fn sanitized_command_line() -> String {
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

pub(crate) fn validate_strict_original_metrics(
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

pub(crate) fn log_toggle(label: &str, enabled: bool) {
    info!(
        "{}: {}",
        label,
        if enabled { "enabled" } else { "disabled" }
    );
}

pub(crate) fn log_run_configuration(args: &Args, sample_override: Option<&str>) {
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

pub(crate) fn sanitize_sample_for_path(sample: &str) -> String {
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

pub(crate) fn append_sample_suffix(base_name: &str, sample_suffix: Option<&str>) -> String {
    match sample_suffix {
        Some(sample) => format!("{}.sample_{}", base_name, sanitize_sample_for_path(sample)),
        None => base_name.to_string(),
    }
}
