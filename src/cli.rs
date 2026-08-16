//! Command-line argument parsing and validation.

use clap::{ArgGroup, Parser, ValueEnum};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum VariantInputFormat {
    Auto,
    Vcf,
    #[value(name = "tsv", alias = "ivar")]
    #[serde(rename = "tsv", alias = "ivar")]
    Tsv,
}

#[derive(Debug, Clone, Parser)]
#[command(
    name = "get_mnv",
    version,
    author = "Paula Ruiz Rodriguez",
    about = "Identifies codon-level SNPs/MNVs, annotates indels and complex alleles, calculates amino acid changes, and outputs results in TSV/VCF format."
)]
#[command(
    group(
        // Not `required(true)`: `--report-from` builds a report from existing
        // TSVs without running the pipeline, so it needs no variant input. The
        // "one of --vcf/--tsv is required" check then lives in `validate`.
        ArgGroup::new("variant_input")
            .required(false)
            .args(["vcf_file", "tsv_file"])
    )
)]
pub struct Args {
    /// Variant input file containing SNVs/MNVs and indels in plain or BGZF-compressed VCF format
    #[arg(
        short = 'v',
        long = "vcf",
        value_name = "VCF_FILE",
        conflicts_with = "tsv_file"
    )]
    pub vcf_file: Option<String>,

    /// iVar variants TSV input file
    #[arg(long = "tsv", value_name = "TSV_FILE", conflicts_with = "vcf_file")]
    pub tsv_file: Option<String>,

    /// Legacy variant input format selector
    #[arg(long = "input-format", value_enum, default_value_t = VariantInputFormat::Auto, hide = true)]
    pub input_format: VariantInputFormat,

    /// BAM file with aligned reads (optional)
    #[arg(short = 'b', long = "bam")]
    pub bam_file: Option<String>,

    /// Reference FASTA file
    #[arg(
        short = 'f',
        long = "fasta",
        required_unless_present = "report_from",
        default_value = ""
    )]
    pub fasta_file: String,

    /// Gene annotation file in TSV format (gene,start,end,strand)
    #[arg(
        short = 'g',
        long = "genes",
        required_unless_present_any = ["gff_file", "report_from"],
        conflicts_with = "gff_file"
    )]
    pub genes_file_tsv: Option<String>,

    /// Gene annotation file in GFF/GFF3 format
    #[arg(long = "gff", required_unless_present_any = ["genes_file_tsv", "report_from"])]
    pub gff_file: Option<String>,

    /// Comma-separated GFF feature types to analyze (default: gene,pseudogene)
    #[arg(long = "gff-features")]
    pub gff_features_raw: Option<String>,

    /// Chromosome/contig to process (optional; default: all contigs in the variant input)
    #[arg(long)]
    pub chrom: Option<String>,

    /// Sample name to use for original FORMAT metrics in multi-sample VCF (default: first sample, use 'all' for all samples)
    #[arg(long)]
    pub sample: Option<String>,

    /// Minimum base Phred quality for BAM read support (default: 20)
    #[arg(short = 'q', long = "quality", default_value_t = 20)]
    pub min_quality: u8,

    /// Minimum mapping quality MAPQ (default: 0)
    #[arg(long = "min-mapq", alias = "mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Normalize REF/ALT alleles (trim shared prefix/suffix) before processing
    #[arg(long = "normalize-alleles")]
    pub normalize_alleles: bool,

    /// Number of threads to use (default: Rayon auto)
    #[arg(long)]
    pub threads: Option<usize>,

    /// Minimum read count for SNP (default: 0)
    #[arg(short = 's', long = "snp", default_value_t = 0)]
    pub min_snp_reads: usize,

    /// Minimum BAM-derived SNP allele frequency, from 0.0 to 1.0 (default: 0.0)
    #[arg(long = "min-snp-frequency", default_value_t = 0.0)]
    pub min_snp_frequency: f64,

    /// Minimum read count for MNV (default: 0)
    #[arg(short = 'm', long = "mnv", default_value_t = 0)]
    pub min_mnv_reads: usize,

    /// Minimum BAM-derived MNV haplotype frequency, from 0.0 to 1.0 (default: 0.0)
    #[arg(long = "min-mnv-frequency", default_value_t = 0.0)]
    pub min_mnv_frequency: f64,

    /// Minimum supporting reads required on each strand for SNP calls (default: 0)
    #[arg(long = "min-snp-strand", default_value_t = 0)]
    pub min_snp_strand_reads: usize,

    /// Minimum supporting reads required on each strand for MNV calls (default: 0)
    #[arg(long = "min-mnv-strand", default_value_t = 0)]
    pub min_mnv_strand_reads: usize,

    /// Minimum Fisher exact p-value accepted for strand-bias metrics (default: 0.0)
    #[arg(long = "min-strand-bias-p", default_value_t = 0.0)]
    pub min_strand_bias_p: f64,

    /// Minimum allele frequency (0.0-1.0) an upstream indel must reach to mark
    /// downstream SNV/MNV codons as frameshifted. Default 0.5 propagates the
    /// frame shift only from a consensus (majority) upstream indel, so a
    /// high-frequency downstream substitution is not relabelled as frameshifted
    /// because of a low-frequency upstream indel that is almost certainly on a
    /// different molecule (intra-host data). Set 0.0 to propagate from every
    /// indel. Indels without a known frequency always propagate.
    #[arg(long = "frameshift-min-freq", default_value_t = 0.5)]
    pub frameshift_min_freq: f64,

    /// By default the indel locus depth (EDP/EFREQ denominator) is counted from
    /// reads observing the anchor base, which avoids under-counting depth and
    /// EFREQ bias for multi-base deletions. Pass --legacy-indel-depth to restrict
    /// the denominator to reads that fully span the REF allele instead.
    #[arg(long = "legacy-indel-depth", action = clap::ArgAction::SetFalse)]
    pub indel_anchor_depth: bool,

    /// Minimum BAM-supporting reads required to emit a phased indel/complex
    /// haplotype row. Default 2: haplotypes are read off the molecules, so a
    /// single read carrying a sequencing error at a called position mints a
    /// combination of its own, and one read is not evidence of a haplotype.
    /// Set 1 to emit every combination any read shows.
    #[arg(long = "phased-indel-min-reads", default_value_t = 2)]
    pub phased_indel_min_reads: usize,

    /// Minimum BAM-derived frequency (0.0-1.0) required to emit a phased
    /// indel/complex haplotype row (default: 0.0).
    #[arg(long = "phased-indel-min-freq", default_value_t = 0.0)]
    pub phased_indel_min_freq: f64,

    /// Count the two mates of a paired-end fragment as two observations
    /// instead of one molecule. By default they are one: a fragment is a
    /// single DNA molecule sequenced from both ends, so counting the mates
    /// separately double-counts wherever they overlap, and treats a variant on
    /// each mate as unrelated when it is in fact proof that the two are on the
    /// same molecule. Single-end data is unaffected either way.
    #[arg(long = "count-mates-separately")]
    pub count_mates_separately: bool,

    /// Parse and validate inputs, print per-contig summary, and skip writing TSV/VCF outputs
    #[arg(long = "dry-run")]
    pub dry_run: bool,

    /// Fail if original depth/frequency metrics (ODP/OFREQ) are missing in input variant calls
    #[arg(long)]
    pub strict: bool,

    /// Split multiallelic VCF records into independent ALT alleles instead of failing
    #[arg(long = "split-multiallelic")]
    pub split_multiallelic: bool,

    /// In VCF output, emit records that fail read/frequency/strand/strand-bias thresholds with FILTER tags instead of skipping them
    #[arg(long = "emit-filtered")]
    pub emit_filtered: bool,

    /// Write VCF output as BGZF-compressed .vcf.gz (requires VCF output mode)
    #[arg(long = "vcf-gz")]
    pub vcf_gz: bool,

    /// Build a Tabix .tbi index for .vcf.gz output (requires --vcf-gz)
    #[arg(long = "index-vcf-gz")]
    pub index_vcf_gz: bool,

    /// Add Fisher exact strand-bias p-values to VCF INFO fields (SBP/MSBP)
    #[arg(long = "strand-bias-info")]
    pub strand_bias_info: bool,

    /// Preserve all original INFO fields from the input VCF in the output VCF (requires --convert or --both)
    #[arg(long = "keep-original-info")]
    pub keep_original_info: bool,

    /// Also write a BCF output converted from generated VCF (requires --convert or --both)
    #[arg(long)]
    pub bcf: bool,

    /// Write run summary as JSON
    #[arg(long = "summary-json")]
    pub summary_json: Option<String>,

    /// Write structured error details as JSON when the command fails
    #[arg(long = "error-json")]
    pub error_json: Option<String>,

    /// Write a reproducibility manifest (inputs, outputs, checksums, runtime metadata)
    #[arg(long = "run-manifest")]
    pub run_manifest: Option<String>,

    /// Write a self-contained interactive HTML report of the called variants.
    /// Requires TSV output, which is the default. --convert writes the VCF
    /// *instead* of the TSV, so use --both rather than --convert when you want
    /// a report alongside VCF output. With --sample all the report covers
    /// every sample.
    #[arg(long = "report", value_name = "HTML_FILE")]
    pub report: Option<String>,

    /// Build the HTML report from existing get_MNV TSV files instead of running
    /// the pipeline, for cohorts processed one sample per run. Each file becomes
    /// one sample, labelled by its file name. Requires --report for the output.
    #[arg(
        long = "report-from",
        value_name = "TSV",
        num_args = 1..,
        requires = "report"
    )]
    pub report_from: Vec<String>,

    /// NCBI translation table number for codon-to-amino-acid mapping
    /// (default: 11 = Bacterial/Archaeal/Plant Plastid).
    /// Supported: 1 (Standard), 2 (Vertebrate Mito), 3 (Yeast Mito),
    /// 4 (Mold/Protozoan Mito), 5 (Invertebrate Mito), 6 (Ciliate),
    /// 11 (Bacterial), 12 (Alt Yeast Nuclear), 25 (SR1/Gracilibacteria).
    #[arg(long = "translation-table", default_value_t = 11)]
    pub translation_table: u8,

    /// Exclude intergenic SNPs (variants outside annotated genes) from the output
    #[arg(long = "exclude-intergenic")]
    pub exclude_intergenic: bool,

    /// Output in VCF format (output file will have a .MNV.vcf extension)
    #[arg(long, conflicts_with = "both")]
    pub convert: bool,

    /// Output both TSV and VCF in a single run
    #[arg(long, conflicts_with = "convert")]
    pub both: bool,

    // These two existed before as fields the desktop app could set, marked
    // `#[arg(skip)]` so the command line could not reach them. The app could put
    // its results wherever the user chose and the command line could not, and
    // every pipeline integration paid for that: the Snakemake wrapper written
    // for this tool runs it in a scratch directory and moves the files
    // afterwards, purely because it could not say where they should go.
    /// Directory to write the TSV, VCF and BCF into. Created if it does not
    /// exist. Without it they go to the current directory
    #[arg(long = "output-dir", value_name = "DIR")]
    pub output_dir: Option<String>,

    /// Base name for the output files, replacing the one taken from the input.
    /// `--output-prefix sample1` writes `sample1.MNV.tsv`. With --sample all the
    /// sample name is still appended, so a cohort does not write every sample
    /// over the same file
    #[arg(long = "output-prefix", value_name = "PREFIX")]
    pub output_prefix: Option<String>,
}

impl Args {
    /// Resolved variant input file path (from --vcf or --tsv).
    pub fn variant_file(&self) -> &str {
        self.tsv_file
            .as_deref()
            .or(self.vcf_file.as_deref())
            .expect("either --vcf or --tsv must be provided")
    }

    /// Variant parser selected by the public input option plus legacy overrides.
    pub fn effective_input_format(&self) -> VariantInputFormat {
        if self.tsv_file.is_some() {
            VariantInputFormat::Tsv
        } else {
            self.input_format
        }
    }

    /// Resolved gene annotation file path (from --genes or --gff).
    pub fn genes_file(&self) -> &str {
        self.gff_file
            .as_deref()
            .or(self.genes_file_tsv.as_deref())
            .expect("either --genes or --gff must be provided")
    }

    /// Parsed GFF feature types.
    pub fn gff_features(&self) -> Vec<String> {
        self.gff_features_raw
            .as_deref()
            .map(|s| {
                s.split(',')
                    .map(|ft| ft.trim().to_string())
                    .filter(|ft| !ft.is_empty())
                    .collect()
            })
            .unwrap_or_else(|| vec!["gene".to_string(), "pseudogene".to_string()])
    }
}

pub fn parse_args() -> Args {
    Args::parse()
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn test_default_gff_features() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
        ])
        .unwrap();
        assert_eq!(args.gff_features(), vec!["gene", "pseudogene"]);
    }

    #[test]
    fn test_custom_gff_features() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--fasta",
            "ref.fa",
            "--gff",
            "ann.gff",
            "--gff-features",
            "CDS,mRNA",
        ])
        .unwrap();
        assert_eq!(args.gff_features(), vec!["CDS", "mRNA"]);
    }

    #[test]
    fn test_translation_table_default() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
        ])
        .unwrap();
        assert_eq!(args.translation_table, 11);
    }

    #[test]
    fn test_translation_table_custom() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
            "--translation-table",
            "2",
        ])
        .unwrap();
        assert_eq!(args.translation_table, 2);
    }

    #[test]
    fn test_frequency_filters_parse() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
            "--min-snp-frequency",
            "0.15",
            "--min-mnv-frequency",
            "0.25",
        ])
        .unwrap();
        assert_eq!(args.min_snp_frequency, 0.15);
        assert_eq!(args.min_mnv_frequency, 0.25);
    }

    #[test]
    fn test_mapq_alias() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
            "--mapq",
            "20",
        ])
        .unwrap();
        assert_eq!(args.min_mapq, 20);
    }

    #[test]
    fn test_input_format_tsv() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "variants.tsv",
            "--input-format",
            "tsv",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
        ])
        .unwrap();
        assert_eq!(args.input_format, VariantInputFormat::Tsv);
        assert_eq!(args.effective_input_format(), VariantInputFormat::Tsv);
    }

    #[test]
    fn test_tsv_option_selects_tsv_input() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--tsv",
            "variants.tsv",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
        ])
        .unwrap();
        assert_eq!(args.tsv_file.as_deref(), Some("variants.tsv"));
        assert_eq!(args.variant_file(), "variants.tsv");
        assert_eq!(args.effective_input_format(), VariantInputFormat::Tsv);
    }

    #[test]
    fn test_vcf_and_tsv_conflict() {
        let result = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--tsv",
            "variants.tsv",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_input_format_ivar_alias() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "variants.tsv",
            "--input-format",
            "ivar",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
        ])
        .unwrap();
        assert_eq!(args.input_format, VariantInputFormat::Tsv);
    }

    #[test]
    fn test_genes_or_gff_required() {
        // Neither --genes nor --gff → should fail
        let result = Args::try_parse_from(["get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa"]);
        assert!(result.is_err());
    }

    #[test]
    fn test_convert_and_both_conflict() {
        let result = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--fasta",
            "ref.fa",
            "--genes",
            "g.txt",
            "--convert",
            "--both",
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_genes_file_accessor() {
        let args = Args::try_parse_from([
            "get_mnv",
            "--vcf",
            "in.vcf",
            "--fasta",
            "ref.fa",
            "--genes",
            "genes.txt",
        ])
        .unwrap();
        assert_eq!(args.genes_file(), "genes.txt");

        let args_gff = Args::try_parse_from([
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa", "--gff", "ann.gff",
        ])
        .unwrap();
        assert_eq!(args_gff.genes_file(), "ann.gff");
    }
}
