//! Variant data types (SNP, MNV, Indel), codon processing, and amino acid
//! change classification.

mod codon;
pub mod consequence;
pub mod declared_phase;
mod event;
pub mod hgvs;
pub mod linkage;
#[cfg(test)]
mod property_tests;
pub mod splice;
mod types;

// Re-export all public types and functions
pub use codon::{
    build_intergenic_variant, build_intron_variant, build_non_coding_variant, build_splice_variant,
    get_mnv_variants_for_gene, get_mnv_variants_for_gene_with_config, local_haplotype_components,
    phased_haplotype_variant, process_codon, FrameshiftPhasing, IndelAnnotationConfig, PairLinkage,
};
pub use event::{
    decompose_allele, observation_ref_len, parse_component_label, substitution_components,
    AlleleComponent, AlleleComponentKind, AlleleEvent, AlleleEventClass,
};
pub use types::*;
