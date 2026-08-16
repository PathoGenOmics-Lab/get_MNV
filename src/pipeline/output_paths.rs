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
    pub fn resolve(args: &Args, base_name: &str, sample_suffix: Option<&str>) -> Self {
        // The sample suffix applies to a caller-supplied prefix too. Without it
        // every sample of a `--sample all` run wrote to the same file name and
        // overwrote the one before, so a cohort of ten produced one file holding
        // the last sample, with nothing to say the other nine had been lost.
        let stem_name = match &args.output_prefix {
            Some(prefix) => super::config::append_sample_suffix(prefix, sample_suffix),
            None => super::config::append_sample_suffix(base_name, sample_suffix),
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

#[cfg(test)]
mod tests {
    use super::OutputPaths;
    use clap::Parser;

    /// Both of these were fields the desktop app set and clap never exposed,
    /// marked `#[arg(skip)]`. Parsing them from an argument list is the whole
    /// point, and a test that only set the fields directly would have passed
    /// throughout the years they were unreachable from a terminal.
    #[test]
    fn the_command_line_can_say_where_output_goes() {
        let args = crate::cli::Args::parse_from([
            "get_mnv",
            "--vcf",
            "v.vcf",
            "--fasta",
            "r.fasta",
            "--genes",
            "g.txt",
            "--output-dir",
            "results/run3",
            "--output-prefix",
            "sample1",
        ]);
        assert_eq!(args.output_dir.as_deref(), Some("results/run3"));
        assert_eq!(args.output_prefix.as_deref(), Some("sample1"));

        let paths = OutputPaths::resolve(&args, "ignored", None);
        assert_eq!(paths.tsv.as_deref(), Some("results/run3/sample1.MNV.tsv"));
    }

    /// A cohort writes one file per sample, and the directory must not flatten
    /// that back into one.
    #[test]
    fn a_directory_and_a_cohort_still_give_one_file_per_sample() {
        let args = crate::cli::Args::parse_from([
            "get_mnv",
            "--vcf",
            "v.vcf",
            "--fasta",
            "r.fasta",
            "--genes",
            "g.txt",
            "--output-dir",
            "out",
        ]);
        let one = OutputPaths::resolve(&args, "cohort", Some("S1"));
        let two = OutputPaths::resolve(&args, "cohort", Some("S2"));
        assert_eq!(one.tsv.as_deref(), Some("out/cohort.sample_S1.MNV.tsv"));
        assert_ne!(one.tsv, two.tsv);
    }

    #[test]
    fn test_output_prefix_keeps_the_sample_suffix() {
        let mut args = crate::cli::Args::parse_from([
            "get_mnv", "--vcf", "v.vcf", "--fasta", "r.fasta", "--genes", "g.txt",
        ]);
        args.output_prefix = Some("custom".to_string());

        let one = OutputPaths::resolve(&args, "cohort", Some("S1"));
        let two = OutputPaths::resolve(&args, "cohort", Some("S2"));
        assert_eq!(one.tsv.as_deref(), Some("custom.sample_S1.MNV.tsv"));
        assert_eq!(two.tsv.as_deref(), Some("custom.sample_S2.MNV.tsv"));
        assert_ne!(
            one.tsv, two.tsv,
            "two samples must not write to the same file"
        );

        // Sin muestra el prefijo se usa tal cual.
        let plain = OutputPaths::resolve(&args, "cohort", None);
        assert_eq!(plain.tsv.as_deref(), Some("custom.MNV.tsv"));
    }
}
