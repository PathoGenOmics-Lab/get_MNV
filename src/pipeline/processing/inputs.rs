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
        VariantInputFormat::Tsv => {
            io::ivar::load_ivar_tsv(variant_file, &references, args.left_align)
                .map_err(reclassify_generic_as_validation)?
        }
        VariantInputFormat::Vcf | VariantInputFormat::Auto => io::load_vcf(
            variant_file,
            sample_override,
            args.split_multiallelic,
            args.allele_normalization(),
            args.keep_original_info,
            Some(&references),
        )
        .map_err(reclassify_generic_as_validation)?,
    };
    let annotation_format = io::detect_annotation_format(args.genes_file())
        .map_err(reclassify_generic_as_validation)?;
    let contigs = selected_contigs(args, &snp_by_contig)?;
    validate_strict_original_metrics(&contigs, &snp_by_contig, args.strict)?;
    let tsv_rows_without_contig = if annotation_format == AnnotationFormat::Tsv && contigs.len() > 1
    {
        io::tsv_rows_without_contig(args.genes_file()).map_err(reclassify_generic_as_validation)?
    } else {
        0
    };
    validate_contig_inputs(
        &contigs,
        &references,
        &snp_by_contig,
        annotation_format,
        tsv_rows_without_contig,
    )?;
    // The seventh column says which contig each gene is on. A table with
    // something else there (a locus tag, a description) used to be read fine;
    // read as contigs, its values match no contig and every gene is quietly
    // dropped. A column that names no contig of the reference at all is not a
    // contig column, and the run says so rather than annotate nothing.
    if annotation_format == AnnotationFormat::Tsv {
        let named =
            io::tsv_named_contigs(args.genes_file()).map_err(reclassify_generic_as_validation)?;
        if !named.is_empty() && !named.iter().any(|contig| references.contains_key(contig)) {
            let mut known: Vec<&String> = references.keys().collect();
            known.sort();
            return Err(crate::error::AppError::validation(format!(
                "The gene table's seventh column names contigs ({}) that are not in the reference ({}). \
                 That column is read as the contig each gene is on; remove it, or write the contig there",
                named.iter().take(3).cloned().collect::<Vec<_>>().join(", "),
                known.iter().take(3).map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
            )));
        }
    }

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

    let catalogue = args
        .catalogue
        .as_deref()
        .map(crate::variants::catalogue::Catalogue::load)
        .transpose()?;

    Ok(ParsedInputs {
        base_name,
        references,
        snp_by_contig,
        contigs,
        command_line: sanitized_command_line(),
        preloaded_gff,
        original_info_headers,
        catalogue,
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
