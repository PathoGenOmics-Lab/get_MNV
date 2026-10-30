//! Input parsing: VCF/BCF loading, BAM/CRAM opening, FASTA reference reading,
//! GFF/TSV gene annotation, allele validation, and original metrics extraction.

pub mod alignment;
pub mod annotation;
pub mod fasta;
pub mod ivar;
pub mod validation;
pub(crate) mod variant_index;
pub mod vcf;
pub mod vcf_fast;
pub(crate) mod vcf_text;

// Re-export public types and functions for backward compatibility.
pub use alignment::{is_cram_path, AlignmentQuery, AlignmentReader};
pub use annotation::{
    detect_annotation_format, filter_genes_with_snps, gene_overlaps_variant, load_all_tsv_genes,
    load_genes, preload_gff_genes, tsv_named_contigs, tsv_rows_without_contig, AnnotationFormat,
};
pub use fasta::{
    drop_reference_calls, left_anchor_insertions, load_references, reference_for_chrom,
    validate_vcf_reference_alleles, Reference, ReferenceMap,
};
pub use validation::get_base_name;
pub(crate) use variant_index::VariantIndex;
pub use vcf::{
    extract_original_info_headers, list_vcf_samples, load_vcf, load_vcf_positions_by_contig,
    parse_declared_phase, DeclaredPhase, VcfPosition,
};

#[cfg(test)]
mod tests;
