//! Variant data types (SNP, MNV, Indel), codon processing, and amino acid
//! change classification.

mod codon;
mod event;
pub mod hgvs;
pub mod splice;
mod types;

// Re-export all public types and functions
pub use codon::{
    build_intergenic_variant, build_intron_variant, build_non_coding_variant,
    build_phased_indel_haplotype_variants, build_splice_variant, get_mnv_variants_for_gene,
    get_mnv_variants_for_gene_with_config, process_codon, FrameshiftPhasing, IndelAnnotationConfig,
};
pub use event::{
    decompose_allele, parse_component_label, substitution_components, AlleleComponent,
    AlleleComponentKind, AlleleEvent, AlleleEventClass,
};
pub use types::*;
