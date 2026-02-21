use clap::{value_parser, Arg, ArgGroup, Command};

#[derive(Debug, Clone)]
pub struct Args {
    pub vcf_file: String,
    pub bam_file: Option<String>,
    pub fasta_file: String,
    pub genes_file: String,
    pub sample: Option<String>,
    pub chrom: Option<String>,
    pub normalize_alleles: bool,
    pub min_quality: u8,
    pub min_mapq: u8,
    pub threads: Option<usize>,
    pub min_snp_reads: usize,
    pub min_mnv_reads: usize,
    pub min_snp_strand_reads: usize,
    pub min_mnv_strand_reads: usize,
    pub min_strand_bias_p: f64,
    pub dry_run: bool,
    pub strict: bool,
    pub split_multiallelic: bool,
    pub emit_filtered: bool,
    pub vcf_gz: bool,
    pub index_vcf_gz: bool,
    pub strand_bias_info: bool,
    pub bcf: bool,
    pub summary_json: Option<String>,
    pub error_json: Option<String>,
    pub run_manifest: Option<String>,
    pub convert: bool,
    pub both: bool,
}

pub fn parse_args() -> Args {
    let matches = Command::new("get_mnv")
        .version("1.0.1")
        .author("Paula Ruiz Rodriguez")
        .about("Identifies SNPs within codons, reclassifies multi-nucleotide variants (MNVs), calculates amino acid changes, and outputs results in TSV/VCF format.")
        .arg(Arg::new("vcf")
            .short('v')
            .long("vcf")
            .value_name("VCF_FILE")
            .help("VCF file containing SNPs")
            .required(true))
        .arg(Arg::new("bam")
            .short('b')
            .long("bam")
            .value_name("BAM_FILE")
            .help("BAM file with aligned reads (optional)")
            .required(false))
        .arg(Arg::new("fasta")
            .short('f')
            .long("fasta")
            .value_name("FASTA_FILE")
            .help("Reference FASTA file")
            .required(true))
        .arg(Arg::new("genes")
            .short('g')
            .long("genes")
            .value_name("GENES_FILE")
            .help("Gene annotation file in TSV format (gene,start,end,strand)")
            .required(false))
        .arg(Arg::new("gff")
            .long("gff")
            .value_name("GFF_FILE")
            .help("Gene annotation file in GFF/GFF3 format")
            .required(false))
        .arg(Arg::new("chrom")
            .long("chrom")
            .value_name("CHROM")
            .help("Chromosome/contig to process (optional; default: all contigs in the VCF)")
            .required(false))
        .arg(Arg::new("sample")
            .long("sample")
            .value_name("SAMPLE")
            .help("Sample name to use for original FORMAT metrics in multi-sample VCF (default: first sample, use 'all' for all samples)")
            .required(false))
        .group(
            ArgGroup::new("annotation")
                .args(["genes", "gff"])
                .required(true)
                .multiple(false),
        )
        .arg(Arg::new("quality")
            .short('q')
            .long("quality")
            .value_name("QUALITY")
            .help("Minimum Phred quality (default: 20)")
            .num_args(1)
            .default_value("20")
            .value_parser(value_parser!(u8)))
        .arg(Arg::new("mapq")
            .long("mapq")
            .value_name("MAPQ")
            .help("Minimum mapping quality MAPQ (default: 0)")
            .num_args(1)
            .default_value("0")
            .value_parser(value_parser!(u8)))
        .arg(Arg::new("normalize_alleles")
            .long("normalize-alleles")
            .help("Normalize REF/ALT alleles (trim shared prefix/suffix) before processing")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("threads")
            .long("threads")
            .value_name("N")
            .help("Number of threads to use (default: Rayon auto)")
            .num_args(1)
            .required(false)
            .value_parser(value_parser!(usize)))
        .arg(Arg::new("min_snp_reads")
            .short('s')
            .long("snp")
            .value_name("SNP")
            .help("Minimum read count for SNP (default: 0)")
            .num_args(1)
            .default_value("0")
            .value_parser(value_parser!(usize)))
        .arg(Arg::new("min_mnv_reads")
            .short('m')
            .long("mnv")
            .value_name("MNV")
            .help("Minimum read count for MNV (default: 0)")
            .num_args(1)
            .default_value("0")
            .value_parser(value_parser!(usize)))
        .arg(Arg::new("min_snp_strand_reads")
            .long("min-snp-strand")
            .value_name("N")
            .help("Minimum supporting reads required on each strand for SNP calls (default: 0)")
            .num_args(1)
            .default_value("0")
            .value_parser(value_parser!(usize)))
        .arg(Arg::new("min_mnv_strand_reads")
            .long("min-mnv-strand")
            .value_name("N")
            .help("Minimum supporting reads required on each strand for MNV calls (default: 0)")
            .num_args(1)
            .default_value("0")
            .value_parser(value_parser!(usize)))
        .arg(Arg::new("min_strand_bias_p")
            .long("min-strand-bias-p")
            .value_name("P")
            .help("Minimum Fisher exact p-value accepted for strand-bias metrics (default: 0.0)")
            .num_args(1)
            .default_value("0.0")
            .value_parser(value_parser!(f64)))
        .arg(Arg::new("dry_run")
            .long("dry-run")
            .help("Parse and validate inputs, print per-contig summary, and skip writing TSV/VCF outputs")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("strict")
            .long("strict")
            .help("Fail if original depth/frequency metrics (ODP/OFREQ) are missing in input VCF")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("split_multiallelic")
            .long("split-multiallelic")
            .help("Split multiallelic VCF records into independent ALT alleles instead of failing")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("emit_filtered")
            .long("emit-filtered")
            .help("Emit records that fail read/strand/strand-bias thresholds with FILTER tags instead of skipping them")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("vcf_gz")
            .long("vcf-gz")
            .help("Write VCF output as BGZF-compressed .vcf.gz (requires VCF output mode)")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("index_vcf_gz")
            .long("index-vcf-gz")
            .help("Build a Tabix .tbi index for .vcf.gz output (requires --vcf-gz)")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("strand_bias_info")
            .long("strand-bias-info")
            .help("Add Fisher exact strand-bias p-values to VCF INFO fields (SBP/MSBP)")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("bcf")
            .long("bcf")
            .help("Also write a BCF output converted from generated VCF")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("summary_json")
            .long("summary-json")
            .value_name("JSON_FILE")
            .help("Write run summary as JSON")
            .required(false))
        .arg(Arg::new("error_json")
            .long("error-json")
            .value_name("JSON_FILE")
            .help("Write structured error details as JSON when the command fails")
            .required(false))
        .arg(Arg::new("run_manifest")
            .long("run-manifest")
            .value_name("JSON_FILE")
            .help("Write a reproducibility manifest (inputs, outputs, checksums, runtime metadata)")
            .required(false))
        .arg(Arg::new("convert")
            .long("convert")
            .help("Output in VCF format (output file will have a .MNV.vcf extension)")
            .conflicts_with("both")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("both")
            .long("both")
            .help("Output both TSV and VCF in a single run")
            .conflicts_with("convert")
            .action(clap::ArgAction::SetTrue))
        .get_matches();

    let genes_file = matches
        .get_one::<String>("gff")
        .or_else(|| matches.get_one::<String>("genes"))
        .expect("Critical error: either --genes or --gff must be provided.")
        .to_string();

    Args {
        vcf_file: matches.get_one::<String>("vcf").unwrap().to_string(),
        bam_file: matches.get_one::<String>("bam").cloned(),
        fasta_file: matches.get_one::<String>("fasta").unwrap().to_string(),
        genes_file,
        sample: matches.get_one::<String>("sample").cloned(),
        chrom: matches.get_one::<String>("chrom").cloned(),
        normalize_alleles: matches.get_flag("normalize_alleles"),
        min_quality: *matches.get_one::<u8>("quality").unwrap(),
        min_mapq: *matches.get_one::<u8>("mapq").unwrap(),
        threads: matches.get_one::<usize>("threads").copied(),
        min_snp_reads: *matches.get_one::<usize>("min_snp_reads").unwrap(),
        min_mnv_reads: *matches.get_one::<usize>("min_mnv_reads").unwrap(),
        min_snp_strand_reads: *matches.get_one::<usize>("min_snp_strand_reads").unwrap(),
        min_mnv_strand_reads: *matches.get_one::<usize>("min_mnv_strand_reads").unwrap(),
        min_strand_bias_p: *matches.get_one::<f64>("min_strand_bias_p").unwrap(),
        dry_run: matches.get_flag("dry_run"),
        strict: matches.get_flag("strict"),
        split_multiallelic: matches.get_flag("split_multiallelic"),
        emit_filtered: matches.get_flag("emit_filtered"),
        vcf_gz: matches.get_flag("vcf_gz"),
        index_vcf_gz: matches.get_flag("index_vcf_gz"),
        strand_bias_info: matches.get_flag("strand_bias_info"),
        bcf: matches.get_flag("bcf"),
        summary_json: matches.get_one::<String>("summary_json").cloned(),
        error_json: matches.get_one::<String>("error_json").cloned(),
        run_manifest: matches.get_one::<String>("run_manifest").cloned(),
        convert: matches.get_flag("convert"),
        both: matches.get_flag("both"),
    }
}
