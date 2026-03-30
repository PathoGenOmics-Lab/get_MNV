//! Command-line argument parsing and validation.

use clap::Parser;

#[derive(Debug, Clone, Parser)]
#[command(
    name = "get_mnv",
    version,
    author = "Paula Ruiz Rodriguez",
    about = "Identifies SNPs within codons, reclassifies multi-nucleotide variants (MNVs), calculates amino acid changes, and outputs results in TSV/VCF format."
)]
pub struct Args {
    /// VCF file containing SNPs
    #[arg(short = 'v', long = "vcf")]
    pub vcf_file: String,

    /// BAM file with aligned reads (optional)
    #[arg(short = 'b', long = "bam")]
    pub bam_file: Option<String>,

    /// Reference FASTA file
    #[arg(short = 'f', long = "fasta")]
    pub fasta_file: String,

    /// Gene annotation file in TSV format (gene,start,end,strand)
    #[arg(short = 'g', long = "genes", required_unless_present = "gff_file", conflicts_with = "gff_file")]
    pub genes_file_tsv: Option<String>,

    /// Gene annotation file in GFF/GFF3 format
    #[arg(long = "gff", required_unless_present = "genes_file_tsv")]
    pub gff_file: Option<String>,

    /// Comma-separated GFF feature types to analyze (default: gene,pseudogene)
    #[arg(long = "gff-features")]
    pub gff_features_raw: Option<String>,

    /// Chromosome/contig to process (optional; default: all contigs in the VCF)
    #[arg(long)]
    pub chrom: Option<String>,

    /// Sample name to use for original FORMAT metrics in multi-sample VCF (default: first sample, use 'all' for all samples)
    #[arg(long)]
    pub sample: Option<String>,

    /// Minimum Phred quality (default: 20)
    #[arg(short = 'q', long = "quality", default_value_t = 20)]
    pub min_quality: u8,

    /// Minimum mapping quality MAPQ (default: 0)
    #[arg(long, default_value_t = 0)]
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

    /// Minimum read count for MNV (default: 0)
    #[arg(short = 'm', long = "mnv", default_value_t = 0)]
    pub min_mnv_reads: usize,

    /// Minimum supporting reads required on each strand for SNP calls (default: 0)
    #[arg(long = "min-snp-strand", default_value_t = 0)]
    pub min_snp_strand_reads: usize,

    /// Minimum supporting reads required on each strand for MNV calls (default: 0)
    #[arg(long = "min-mnv-strand", default_value_t = 0)]
    pub min_mnv_strand_reads: usize,

    /// Minimum Fisher exact p-value accepted for strand-bias metrics (default: 0.0)
    #[arg(long = "min-strand-bias-p", default_value_t = 0.0)]
    pub min_strand_bias_p: f64,

    /// Parse and validate inputs, print per-contig summary, and skip writing TSV/VCF outputs
    #[arg(long = "dry-run")]
    pub dry_run: bool,

    /// Fail if original depth/frequency metrics (ODP/OFREQ) are missing in input VCF
    #[arg(long)]
    pub strict: bool,

    /// Split multiallelic VCF records into independent ALT alleles instead of failing
    #[arg(long = "split-multiallelic")]
    pub split_multiallelic: bool,

    /// Emit records that fail read/strand/strand-bias thresholds with FILTER tags instead of skipping them
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

    /// Preserve all original INFO fields from the input VCF in the output VCF
    #[arg(long = "keep-original-info")]
    pub keep_original_info: bool,

    /// Also write a BCF output converted from generated VCF
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

    /// Optional output directory (used by desktop app)
    #[arg(skip)]
    pub output_dir: Option<String>,

    /// Optional output filename prefix (used by desktop app)
    #[arg(skip)]
    pub output_prefix: Option<String>,
}

impl Args {
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
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa", "--genes", "genes.txt",
        ]).unwrap();
        assert_eq!(args.gff_features(), vec!["gene", "pseudogene"]);
    }

    #[test]
    fn test_custom_gff_features() {
        let args = Args::try_parse_from([
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa", "--gff", "ann.gff",
            "--gff-features", "CDS,mRNA",
        ]).unwrap();
        assert_eq!(args.gff_features(), vec!["CDS", "mRNA"]);
    }

    #[test]
    fn test_translation_table_default() {
        let args = Args::try_parse_from([
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa", "--genes", "genes.txt",
        ]).unwrap();
        assert_eq!(args.translation_table, 11);
    }

    #[test]
    fn test_translation_table_custom() {
        let args = Args::try_parse_from([
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa", "--genes", "genes.txt",
            "--translation-table", "2",
        ]).unwrap();
        assert_eq!(args.translation_table, 2);
    }

    #[test]
    fn test_genes_or_gff_required() {
        // Neither --genes nor --gff → should fail
        let result = Args::try_parse_from([
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa",
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_convert_and_both_conflict() {
        let result = Args::try_parse_from([
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa", "--genes", "g.txt",
            "--convert", "--both",
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_genes_file_accessor() {
        let args = Args::try_parse_from([
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa", "--genes", "genes.txt",
        ]).unwrap();
        assert_eq!(args.genes_file(), "genes.txt");

        let args_gff = Args::try_parse_from([
            "get_mnv", "--vcf", "in.vcf", "--fasta", "ref.fa", "--gff", "ann.gff",
        ]).unwrap();
        assert_eq!(args_gff.genes_file(), "ann.gff");
    }
}
