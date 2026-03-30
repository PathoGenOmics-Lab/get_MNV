//! **get_MNV** — Multi-Nucleotide Variant detection in genomic sequences.
//!
//! This crate identifies codon-level MNVs from VCF + reference + gene
//! annotation, recalculates amino acid changes considering all SNVs together,
//! and optionally quantifies read support from BAM files.

pub mod cli;
pub mod error;
pub mod genetic_code;
pub mod io;
pub mod output;
pub mod pipeline;
pub mod read_count;
pub mod utils;
pub mod variants;

/// Core library version, sourced from Cargo.toml at compile time.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
