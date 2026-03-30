//! Input parsing: VCF loading, FASTA reference reading, GFF/TSV gene
//! annotation, allele validation, and original metrics extraction.

pub mod annotation;
pub mod fasta;
pub mod validation;
pub mod vcf;
pub mod vcf_fast;

// Re-export public types and functions for backward compatibility.
pub use annotation::{
    detect_annotation_format, filter_genes_with_snps, load_genes, preload_gff_genes,
    AnnotationFormat,
};
pub use fasta::{
    load_references, reference_for_chrom, validate_vcf_reference_alleles, Reference, ReferenceMap,
};
pub use validation::get_base_name;
pub use vcf::{
    extract_original_info_headers, list_vcf_samples, load_vcf_positions_by_contig, VcfPosition,
};

#[cfg(test)]
mod tests;
