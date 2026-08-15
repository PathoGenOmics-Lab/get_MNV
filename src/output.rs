//! Output writers for TSV and VCF/BCF formats.

mod common;

// The consequence a row states, which is what tells two rows about the same base
// apart when they name the same gene. A GFF gives every transcript of a gene the
// same Name, so the name alone cannot say whether a row has already been written.
pub(crate) use common::so_consequence;
mod report;
pub(crate) mod stats;
mod tsv;
mod vcf;

pub use report::{sample_labels, ReportBuilder};
pub use tsv::{TsvWriter, TsvWriterConfig};
pub use vcf::{build_tabix_index, convert_vcf_to_bcf, ExternalStep, VcfWriter, VcfWriterConfig};
