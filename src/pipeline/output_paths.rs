//! Output file path resolution: stem name, TSV/VCF/BCF paths based on CLI flags.

use crate::cli::Args;

/// Resolved output file paths for a pipeline run.
pub(crate) struct OutputPaths {
    pub output_stem: String,
    pub tsv: Option<String>,
    pub vcf: Option<String>,
    pub bcf: Option<String>,
}

impl OutputPaths {
    /// Resolve output paths from CLI args, base name, and optional sample suffix.
    pub fn resolve(
        args: &Args,
        base_name: &str,
        sample_suffix: Option<&str>,
    ) -> Self {
        let base_stem = super::config::append_sample_suffix(base_name, sample_suffix);
        let stem_name = match &args.output_prefix {
            Some(prefix) => prefix.clone(),
            None => base_stem,
        };
        let output_stem = match &args.output_dir {
            Some(dir) => {
                let dir_path = std::path::Path::new(dir);
                dir_path.join(&stem_name).to_string_lossy().into_owned()
            }
            None => stem_name,
        };

        let tsv = if args.dry_run {
            None
        } else if args.both || !args.convert {
            Some(format!("{output_stem}.MNV.tsv"))
        } else {
            None
        };

        let vcf = if args.dry_run {
            None
        } else if args.both || args.convert {
            Some(if args.vcf_gz {
                format!("{output_stem}.MNV.vcf.gz")
            } else {
                format!("{output_stem}.MNV.vcf")
            })
        } else {
            None
        };

        let bcf = if args.dry_run || !args.bcf {
            None
        } else {
            Some(format!("{output_stem}.MNV.bcf"))
        };

        Self {
            output_stem,
            tsv,
            vcf,
            bcf,
        }
    }
}
