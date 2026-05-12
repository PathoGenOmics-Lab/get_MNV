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
            let cleaned = if arg.contains('/') {
                Path::new(&arg)
                    .file_name()
                    .and_then(|value| value.to_str())
                    .map(ToOwned::to_owned)
                    .unwrap_or(arg)
            } else {
                arg
            };
            // Escape control characters that could corrupt VCF header lines
            cleaned
                .replace('\t', "\\t")
                .replace('\n', "\\n")
                .replace('\r', "\\r")
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
    info!("Variant input file: {}", args.vcf_file);
    info!("Variant input format: {:?}", args.input_format);
    if let Some(ref bam) = args.bam_file {
        info!("BAM file: {bam}");
    } else {
        info!("No BAM file provided. Output will be generated without read count fields.");
    }
    info!("FASTA file: {}", args.fasta_file);
    info!("Gene annotation file: {}", args.genes_file());
    info!("GFF feature types: {}", args.gff_features().join(", "));
    if let Ok(AnnotationFormat::Tsv) = io::detect_annotation_format(args.genes_file()) {
        let default_features = vec!["gene".to_string(), "pseudogene".to_string()];
        if args.gff_features() != default_features {
            log::warn!("--gff-features is ignored when using TSV annotation format (--genes)");
        }
    }
    if let Some(sample) = sample_override {
        info!("Target sample for original FORMAT metrics: {sample}");
    } else {
        info!("Target sample for original FORMAT metrics: first sample");
    }
    if let Some(chrom) = args.chrom.as_ref() {
        info!("Target chromosome(s): {chrom}");
    } else {
        info!("Target chromosome(s): all contigs in input VCF");
    }
    if let Some(threads) = args.threads {
        info!("Threads: {threads}");
    } else {
        info!("Threads: Rayon default");
    }
    info!("Minimum Phred quality: {}", args.min_quality);
    info!("Minimum mapping quality (MAPQ): {}", args.min_mapq);
    info!("Minimum SNP reads: {}", args.min_snp_reads);
    info!("Minimum SNP frequency: {:.4}", args.min_snp_frequency);
    info!("Minimum MNV reads: {}", args.min_mnv_reads);
    info!("Minimum MNV frequency: {:.4}", args.min_mnv_frequency);
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_configure_threads_zero_is_error() {
        let result = configure_threads(Some(0));
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("must be >= 1"));
    }

    #[test]
    fn test_configure_threads_none_is_ok() {
        assert!(configure_threads(None).is_ok());
    }

    #[test]
    fn test_selected_contigs_from_chrom_arg() {
        let mut args = default_test_args();
        args.chrom = Some("chrB, chrA, chrA".to_string());
        let snps: HashMap<String, Vec<VcfPosition>> = HashMap::new();
        let contigs = selected_contigs(&args, &snps).unwrap();
        assert_eq!(contigs, vec!["chrA", "chrB"]); // sorted, deduped
    }

    #[test]
    fn test_selected_contigs_empty_is_error() {
        let mut args = default_test_args();
        args.chrom = Some(String::new());
        let snps: HashMap<String, Vec<VcfPosition>> = HashMap::new();
        let result = selected_contigs(&args, &snps);
        assert!(result.is_err());
    }

    #[test]
    fn test_sanitize_sample_for_path_special_chars() {
        assert_eq!(sanitize_sample_for_path("sample/1:bad"), "sample_1_bad");
        assert_eq!(sanitize_sample_for_path("ok-name_123"), "ok-name_123");
    }

    #[test]
    fn test_append_sample_suffix_with_and_without() {
        assert_eq!(append_sample_suffix("output", None), "output");
        assert_eq!(
            append_sample_suffix("output", Some("sample1")),
            "output.sample_sample1"
        );
    }

    #[test]
    fn test_validate_contig_inputs_tsv_multi_contig_error() {
        let contigs = vec!["chr1".to_string(), "chr2".to_string()];
        let refs: ReferenceMap = [
            ("chr1".into(), "ACGT".into()),
            ("chr2".into(), "TTTT".into()),
        ]
        .into_iter()
        .collect();
        let snps: HashMap<String, Vec<VcfPosition>> =
            [("chr1".into(), vec![]), ("chr2".into(), vec![])]
                .into_iter()
                .collect();
        let result = validate_contig_inputs(&contigs, &refs, &snps, AnnotationFormat::Tsv);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("TSV annotation"));
    }

    #[test]
    fn test_validate_strict_original_metrics_missing_dp() {
        let contigs = vec!["chr1".to_string()];
        let snps: HashMap<String, Vec<VcfPosition>> = [(
            "chr1".into(),
            vec![VcfPosition {
                position: 10,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                original_dp: None,
                original_freq: Some(0.5),
                original_info: None,
            }],
        )]
        .into_iter()
        .collect();
        let result = validate_strict_original_metrics(&contigs, &snps, true);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("ODP missing"));
    }

    #[test]
    fn test_validate_strict_disabled_always_ok() {
        let contigs = vec!["chr1".to_string()];
        let snps: HashMap<String, Vec<VcfPosition>> = [(
            "chr1".into(),
            vec![VcfPosition {
                position: 10,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            }],
        )]
        .into_iter()
        .collect();
        assert!(validate_strict_original_metrics(&contigs, &snps, false).is_ok());
    }

    fn default_test_args() -> crate::cli::Args {
        crate::cli::Args {
            vcf_file: String::new(),
            input_format: crate::cli::VariantInputFormat::Auto,
            bam_file: None,
            fasta_file: String::new(),
            genes_file_tsv: Some(String::new()),
            gff_file: None,
            gff_features_raw: None,
            sample: None,
            chrom: None,
            normalize_alleles: false,
            min_quality: 20,
            min_mapq: 20,
            threads: None,
            min_snp_reads: 1,
            min_snp_frequency: 0.0,
            min_mnv_reads: 1,
            min_mnv_frequency: 0.0,
            min_snp_strand_reads: 0,
            min_mnv_strand_reads: 0,
            min_strand_bias_p: 0.0,
            dry_run: false,
            strict: false,
            split_multiallelic: false,
            emit_filtered: false,
            vcf_gz: false,
            index_vcf_gz: false,
            strand_bias_info: false,
            keep_original_info: false,
            exclude_intergenic: false,
            bcf: false,
            summary_json: None,
            error_json: None,
            run_manifest: None,
            convert: false,
            both: false,
            translation_table: 11,
            output_dir: None,
            output_prefix: None,
        }
    }
}
