//! The inputs every sample of a `--sample all` cohort shares, loaded once.
//!
//! A cohort used to run the whole pipeline once per sample, and every run read
//! the reference, the annotation and the variant file again. Loading them once
//! must not change what any sample says. Each load here keeps the lines it
//! logged, and every sample that asks for it logs them again at the point
//! where its own run would have loaded it. A load that failed fails every
//! sample that asks, with the same error.

use super::manifest::build_input_metadata;
use super::summary::RunInputs;
use crate::cli::Args;
use crate::error::AppResult;
use crate::io::{self, AnnotationFormat, ReferenceMap, VcfPosition};
use crate::logging::{self, Held};
use crate::variants::catalogue::Catalogue;
use crate::variants::Gene;
use std::collections::HashMap;
use std::sync::{Arc, OnceLock};

/// One input, loaded by whichever sample asks for it first.
struct Loaded<T> {
    cell: OnceLock<(Vec<Held>, AppResult<Arc<T>>)>,
}

impl<T> Loaded<T> {
    fn new() -> Self {
        Loaded {
            cell: OnceLock::new(),
        }
    }

    fn get(&self, load: impl FnOnce() -> AppResult<T>) -> AppResult<Arc<T>> {
        let (said, loaded) = self.cell.get_or_init(|| {
            let (loaded, said) = logging::holding(load);
            (said, loaded.map(Arc::new))
        });
        for line in said {
            line.write();
        }
        match loaded {
            Ok(value) => Ok(Arc::clone(value)),
            Err(error) => Err(error.replica()),
        }
    }
}

/// What a cohort's samples load once between them.
pub(crate) struct SharedInputs {
    /// The variant file held in memory, or `None` when it could not be held
    /// and every sample reads it as before.
    variants: Option<io::vcf_text::VcfText>,
    references: Loaded<ReferenceMap>,
    gff_has_cds_features: Loaded<bool>,
    gff_genes: Loaded<HashMap<String, Vec<Gene>>>,
    tsv_rows: Loaded<Vec<io::annotation::TsvGeneRow>>,
    tsv_rows_without_contig: Loaded<usize>,
    tsv_named_contigs: Loaded<Vec<String>>,
    original_info_headers: Loaded<Vec<String>>,
    catalogue: Loaded<Catalogue>,
    input_metadata: Loaded<RunInputs>,
}

impl SharedInputs {
    /// Read the variant file into memory; everything else waits until a sample
    /// asks for it. A file that cannot be held, or that fails to open now,
    /// is left for each sample to read, so it fails the way it always did.
    pub(crate) fn new(args: &Args) -> Self {
        SharedInputs {
            variants: io::vcf_text::VcfText::read(args.variant_file())
                .ok()
                .flatten(),
            references: Loaded::new(),
            gff_has_cds_features: Loaded::new(),
            gff_genes: Loaded::new(),
            tsv_rows: Loaded::new(),
            tsv_rows_without_contig: Loaded::new(),
            tsv_named_contigs: Loaded::new(),
            original_info_headers: Loaded::new(),
            catalogue: Loaded::new(),
            input_metadata: Loaded::new(),
        }
    }

    pub(crate) fn references(&self, args: &Args) -> AppResult<Arc<ReferenceMap>> {
        self.references
            .get(|| io::load_references(&args.fasta_file, args.chrom.as_deref()))
    }

    /// One sample's variants, taken from the file held in memory.
    pub(crate) fn variants(
        &self,
        args: &Args,
        sample: Option<&str>,
        references: &ReferenceMap,
    ) -> AppResult<HashMap<String, Vec<VcfPosition>>> {
        let Some(text) = &self.variants else {
            return io::load_vcf(
                args.variant_file(),
                sample,
                args.split_multiallelic,
                args.normalize_alleles,
                args.keep_original_info,
                Some(references),
            );
        };
        match text.splitting() {
            io::vcf_text::Splitting::AsWritten => io::vcf_fast::load_held_vcf_text(
                text,
                sample,
                args.split_multiallelic,
                args.normalize_alleles,
                args.keep_original_info,
                Some(references),
            ),
            io::vcf_text::Splitting::Trimmed => io::vcf::load_held_vcf_positions(
                text,
                sample,
                args.split_multiallelic,
                args.normalize_alleles,
                args.keep_original_info,
                Some(references),
            ),
        }
    }

    /// The INFO header lines a run with `--keep-original-info` carries over.
    pub(crate) fn original_info_headers(&self, args: &Args) -> AppResult<Arc<Vec<String>>> {
        self.original_info_headers.get(|| {
            let variant_file = args.variant_file();
            match &self.variants {
                Some(text) => {
                    let lines = text.replay().map(|line| line.map(str::to_string));
                    match text.splitting() {
                        io::vcf_text::Splitting::AsWritten => {
                            io::vcf_fast::text_info_headers_from(lines)
                        }
                        io::vcf_text::Splitting::Trimmed => {
                            io::vcf::original_info_headers_from(lines)
                        }
                    }
                }
                None if io::vcf_fast::use_fast_parser(variant_file) => {
                    io::vcf_fast::extract_text_info_headers(variant_file)
                }
                None => io::extract_original_info_headers(variant_file),
            }
        })
    }

    pub(crate) fn gff_has_cds_features(&self, args: &Args) -> AppResult<bool> {
        self.gff_has_cds_features
            .get(|| io::annotation::gff_has_cds_features(args.genes_file()))
            .map(|has| *has)
    }

    /// The GFF annotation, for the feature types every sample of the run
    /// settles on alike.
    pub(crate) fn gff_genes(
        &self,
        args: &Args,
        gff_features: &[String],
    ) -> AppResult<Arc<HashMap<String, Vec<Gene>>>> {
        self.gff_genes
            .get(|| io::preload_gff_genes(args.genes_file(), gff_features))
    }

    pub(crate) fn tsv_rows_without_contig(&self, args: &Args) -> AppResult<usize> {
        self.tsv_rows_without_contig
            .get(|| io::tsv_rows_without_contig(args.genes_file()))
            .map(|count| *count)
    }

    pub(crate) fn tsv_named_contigs(&self, args: &Args) -> AppResult<Arc<Vec<String>>> {
        self.tsv_named_contigs
            .get(|| io::tsv_named_contigs(args.genes_file()))
    }

    /// The genes of `contig` a sample annotates, as [`io::load_genes`] gives
    /// them, or every gene of the contig when `snp_list` is `None`, as
    /// [`io::load_all_tsv_genes`] does. The TSV table is read once; a GFF here
    /// is read as before.
    pub(crate) fn contig_genes(
        &self,
        args: &Args,
        snp_list: Option<&[VcfPosition]>,
        contig: &str,
    ) -> AppResult<Vec<Gene>> {
        let genes_file = args.genes_file();
        if let Some(snp_list) = snp_list {
            if io::detect_annotation_format(genes_file)? == AnnotationFormat::Gff {
                return io::load_genes(genes_file, snp_list, Some(contig), &args.gff_features());
            }
        }
        io::annotation::load_tsv_genes_with(
            genes_file,
            || {
                self.tsv_rows
                    .get(|| io::annotation::read_tsv_gene_rows(genes_file))
            },
            snp_list,
            Some(contig),
        )
    }

    pub(crate) fn catalogue(&self, path: &str) -> AppResult<Arc<Catalogue>> {
        self.catalogue.get(|| Catalogue::load(path))
    }

    /// The input files and their checksums, which every summary repeats.
    pub(crate) fn input_metadata(&self, args: &Args) -> AppResult<RunInputs> {
        self.input_metadata
            .get(|| build_input_metadata(args, true))
            .map(|inputs| (*inputs).clone())
    }
}
