//! Output writers for TSV and VCF/BCF formats.

mod common;
mod report;
pub(crate) mod stats;
mod tsv;
mod vcf;

pub use report::{sample_labels, ReportBuilder};
pub use tsv::{TsvWriter, TsvWriterConfig};
pub use vcf::{build_tabix_index, convert_vcf_to_bcf, VcfWriter, VcfWriterConfig};
