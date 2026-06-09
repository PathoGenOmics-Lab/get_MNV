//! Input parsing: variant/reference/annotation loading and input-format
//! resolution into a `ParsedInputs` bundle.

use super::reclassify_generic_as_validation;
use super::ParsedInputs;
use crate::cli::{Args, VariantInputFormat};
use crate::error::{AppError, AppResult};
use crate::io::{self, AnnotationFormat};
use crate::pipeline::config::{
    sanitized_command_line, selected_contigs, validate_contig_inputs,
    validate_strict_original_metrics,
};

pub(crate) fn parse_inputs(args: &Args, sample_override: Option<&str>) -> AppResult<ParsedInputs> {
    let variant_file = args.variant_file();
    let base_name = io::get_base_name(variant_file).map_err(reclassify_generic_as_validation)?;
    let references = io::load_references(&args.fasta_file, args.chrom.as_deref())
        .map_err(reclassify_generic_as_validation)?;
    let input_format = resolve_variant_input_format(args)?;
    if input_format == VariantInputFormat::Tsv && sample_override.is_some() {
        return Err(AppError::config(
            "--sample is only supported for VCF input; TSV input is treated as a single sample",
        ));
    }

    let snp_by_contig = match input_format {
        VariantInputFormat::Tsv => io::ivar::load_ivar_tsv(variant_file, &references)
            .map_err(reclassify_generic_as_validation)?,
        VariantInputFormat::Vcf | VariantInputFormat::Auto => {
            // Use the fast text parser for plain .vcf files and the
            // BGZF-aware parser for .vcf.gz. BCF input is rejected with a
            // conversion hint.
            if io::vcf_fast::use_fast_parser(variant_file) {
                io::vcf_fast::load_vcf_text(
                    variant_file,
                    sample_override,
                    args.split_multiallelic,
                    args.normalize_alleles,
                    args.keep_original_info,
                )
                .map_err(reclassify_generic_as_validation)?
            } else {
                io::load_vcf_positions_by_contig(
                    variant_file,
                    sample_override,
                    args.split_multiallelic,
                    args.normalize_alleles,
                    args.keep_original_info,
                )
                .map_err(reclassify_generic_as_validation)?
            }
        }
    };
    let annotation_format = io::detect_annotation_format(args.genes_file())
        .map_err(reclassify_generic_as_validation)?;
    let contigs = selected_contigs(args, &snp_by_contig)?;
    validate_strict_original_metrics(&contigs, &snp_by_contig, args.strict)?;
    validate_contig_inputs(&contigs, &references, &snp_by_contig, annotation_format)?;

    // Auto-select the phase/splice-aware CDS model when the user did not set
    // --gff-features and the GFF actually contains CDS features; otherwise a
    // eukaryotic GFF would be analysed over whole-gene spans (introns included),
    // silently mis-numbering amino acids. Explicitly passing --gff-features
    // (e.g. `gene`) keeps the legacy behaviour and the existing phase warning.
    let gff_features = if args.gff_features_raw.is_none()
        && annotation_format == AnnotationFormat::Gff
        && io::annotation::gff_has_cds_features(args.genes_file())
            .map_err(reclassify_generic_as_validation)?
    {
        log::info!(
            "GFF contains CDS features and --gff-features was not set; analysing 'CDS' \
             (phase- and splice-aware codon model). Pass --gff-features gene to use whole-gene spans."
        );
        vec!["CDS".to_string()]
    } else {
        args.gff_features()
    };

    let preloaded_gff = if annotation_format == AnnotationFormat::Gff {
        Some(
            io::preload_gff_genes(args.genes_file(), &gff_features)
                .map_err(reclassify_generic_as_validation)?,
        )
    } else {
        None
    };

    if args.keep_original_info && input_format == VariantInputFormat::Tsv {
        log::warn!(
            "--keep-original-info has no effect with iVar TSV input: the TSV format has no \
             INFO column to preserve, so no original INFO fields will be carried over."
        );
    }
    let original_info_headers =
        if args.keep_original_info && input_format == VariantInputFormat::Vcf {
            if io::vcf_fast::use_fast_parser(variant_file) {
                io::vcf_fast::extract_text_info_headers(variant_file)
                    .map_err(reclassify_generic_as_validation)?
            } else {
                io::extract_original_info_headers(variant_file)
                    .map_err(reclassify_generic_as_validation)?
            }
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

pub(crate) fn resolve_variant_input_format(args: &Args) -> AppResult<VariantInputFormat> {
    if args.tsv_file.is_some() && args.input_format == VariantInputFormat::Vcf {
        return Err(AppError::config(
            "--tsv cannot be combined with --input-format vcf",
        ));
    }

    let variant_file = args.variant_file();
    match args.effective_input_format() {
        VariantInputFormat::Auto => {
            if io::ivar::looks_like_ivar_tsv(variant_file)
                .map_err(reclassify_generic_as_validation)?
            {
                Ok(VariantInputFormat::Tsv)
            } else {
                Ok(VariantInputFormat::Vcf)
            }
        }
        explicit => Ok(explicit),
    }
}
