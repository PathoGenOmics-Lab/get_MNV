//! Variant data types (SNP, MNV, Indel), codon processing, and amino acid
//! change classification.

mod codon;
mod event;
mod types;

// Re-export all public types and functions
pub use codon::{
    build_intergenic_variant, build_phased_indel_haplotype_variants, get_mnv_variants_for_gene,
    process_codon,
};
pub use event::{
    decompose_allele, parse_component_label, substitution_components, AlleleComponent,
    AlleleComponentKind, AlleleEvent, AlleleEventClass,
};
pub use types::*;
