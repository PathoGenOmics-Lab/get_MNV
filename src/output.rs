//! Output writers for TSV and VCF/BCF formats.

mod common;
pub(crate) mod stats;
mod tsv;
mod vcf;

pub use tsv::TsvWriter;
pub use vcf::{build_tabix_index, convert_vcf_to_bcf, VcfWriter, VcfWriterConfig};
