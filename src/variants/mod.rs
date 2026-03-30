//! Variant data types (SNP, MNV, Indel), codon processing, and amino acid
//! change classification.

mod codon;
mod types;

// Re-export all public types and functions
pub use codon::{build_intergenic_variant, get_mnv_variants_for_gene, process_codon};
pub use types::*;
