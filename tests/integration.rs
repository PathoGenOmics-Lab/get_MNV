//! End-to-end integration tests using the example dataset.
//!
//! These tests run the full pipeline (`pipeline::run`) against the bundled
//! `example/` data and verify output correctness at a structural level.

use get_mnv::cli::{Args, VariantInputFormat};
use get_mnv::pipeline;
use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn example_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("example")
}

fn temp_dir(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("get_mnv_test_{}_{}", prefix, nanos));
    fs::create_dir_all(&dir).expect("create temp dir");
    dir
}

fn base_args() -> Args {
    let ex = example_dir();
    Args {
        vcf_file: Some(ex.join("G35894.var.snp.vcf").to_string_lossy().into()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: ex.join("MTB_ancestor.fas").to_string_lossy().into(),
        genes_file_tsv: Some(ex.join("anot_genes.txt").to_string_lossy().into()),
        gff_file: None,
        gff_features_raw: None,
        sample: None,
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: None,
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: false,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        output_dir: None,
        translation_table: 11,
        output_prefix: None,
    }
}

// ---------------------------------------------------------------------------
// T1: Full pipeline E2E — TSV output
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_tsv_output_matches_expected_variant_counts() {
    let tmp = temp_dir("e2e_tsv");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_out".to_string());

    let summary = pipeline::run(&args).expect("pipeline should succeed");

    // The example has 950 VCF records (1 contig: MTB_anc)
    assert_eq!(summary.global.contig_count, 1, "should process 1 contig");
    assert!(
        summary.global.snp_records_in_vcf >= 900,
        "expected ~950 VCF records, got {}",
        summary.global.snp_records_in_vcf
    );

    // Variant counts: 797 SNP + 9 MNV + 1 SNP/MNV = 807 genic + intergenic
    let total = summary.global.produced_variants;
    assert!(total >= 800, "expected ≥800 variants, got {total}");
    assert!(
        summary.global.snp_variants >= 700,
        "expected ≥700 SNPs, got {}",
        summary.global.snp_variants
    );
    // MNVs can be counted as mnv_variants or snp_mnv_variants depending
    // on whether individual SNP reads were also observed
    let mnv_total = summary.global.mnv_variants + summary.global.snp_mnv_variants;
    assert!(
        mnv_total >= 5,
        "expected ≥5 MNV+SNP/MNV, got {mnv_total} (mnv={}, snp_mnv={})",
        summary.global.mnv_variants,
        summary.global.snp_mnv_variants
    );

    // TSV file should exist
    let tsv_path = tmp.join("test_out.MNV.tsv");
    assert!(tsv_path.exists(), "TSV output file should exist");
    let tsv_content = fs::read_to_string(&tsv_path).expect("read TSV");
    let tsv_lines: Vec<&str> = tsv_content.lines().collect();
    // Header + data lines
    assert!(
        tsv_lines.len() >= 800,
        "TSV should have ≥800 lines (header + data), got {}",
        tsv_lines.len()
    );

    // Header should contain expected columns (no BAM)
    let header = tsv_lines[0];
    assert!(header.contains("Chromosome"), "header missing Chromosome");
    assert!(header.contains("Gene"), "header missing Gene");
    assert!(
        header.contains("Variant Type"),
        "header missing Variant Type"
    );
    assert!(header.contains("Change Type"), "header missing Change Type");

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_ivar_tsv_input_with_tsv_option() {
    let tmp = temp_dir("e2e_ivar");
    let ivar_path = tmp.join("sample_variants.tsv");
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");

    fs::write(
        &ivar_path,
        "REGION\tPOS\tREF\tALT\tREF_DP\tALT_DP\tALT_FREQ\tTOTAL_DP\tPASS\n\
chr1\t4\tG\tA\t1\t9\t0.9\t10\tTRUE\n\
chr1\t5\tG\tA\t1\t9\t0.9\t10\tTRUE\n\
chr1\t6\tA\tC\t1\t9\t0.9\t10\tTRUE\n\
chr1\t7\tC\tT\t1\t9\t0.9\t10\tFALSE\n\
chr1\t8\tC\t+A\t1\t9\t0.9\t10\tTRUE\n",
    )
    .unwrap();
    fs::write(&ref_path, ">chr1\nATGGGACCCTAA\n").unwrap();
    fs::write(&genes_path, "gene1\t1\t12\t+\n").unwrap();

    let mut args = base_args();
    args.vcf_file = None;
    args.tsv_file = Some(ivar_path.to_string_lossy().into());
    args.input_format = VariantInputFormat::Auto;
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("ivar_out".to_string());

    let summary = pipeline::run(&args).expect("iVar TSV pipeline should succeed");
    assert_eq!(summary.global.snp_records_in_vcf, 3);
    assert_eq!(summary.global.produced_variants, 1);
    assert_eq!(summary.global.snp_mnv_variants, 1);

    let tsv_content = fs::read_to_string(tmp.join("ivar_out.MNV.tsv")).expect("read TSV");
    assert!(tsv_content.contains("4, 5, 6"));
    assert!(tsv_content.contains("SNP/MNV"));

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1b: Full pipeline E2E — VCF output
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_vcf_output_is_valid() {
    let tmp = temp_dir("e2e_vcf");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_vcf".to_string());
    args.convert = true;

    let summary = pipeline::run(&args).expect("pipeline should succeed");

    let vcf_path = tmp.join("test_vcf.MNV.vcf");
    assert!(vcf_path.exists(), "VCF output file should exist");

    let vcf_content = fs::read_to_string(&vcf_path).expect("read VCF");
    let lines: Vec<&str> = vcf_content.lines().collect();

    // Check VCF header
    assert!(
        lines[0].starts_with("##fileformat=VCFv4"),
        "missing VCF header"
    );
    let header_line = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .expect("missing #CHROM line");
    assert!(header_line.contains("INFO"), "header missing INFO column");

    // Data lines (non-header, non-comment)
    let data_lines: Vec<&&str> = lines.iter().filter(|l| !l.starts_with('#')).collect();
    assert!(
        data_lines.len() >= 800,
        "VCF should have ≥800 data lines, got {}",
        data_lines.len()
    );

    // Every data line should have the contig
    for line in &data_lines {
        assert!(
            line.starts_with("MTB_anc\t"),
            "unexpected contig in VCF line"
        );
    }

    // Check INFO fields contain expected tags
    let sample_line = data_lines[0];
    assert!(sample_line.contains("GENE="), "missing GENE in INFO");
    assert!(sample_line.contains("TYPE="), "missing TYPE in INFO");
    assert!(sample_line.contains("CT="), "missing CT in INFO");

    assert!(
        summary.output_vcf.is_some(),
        "summary should report VCF output"
    );
    assert!(
        summary.output_tsv.is_none(),
        "should not produce TSV in convert mode"
    );

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1c: Both TSV + VCF output
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_both_outputs() {
    let tmp = temp_dir("e2e_both");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_both".to_string());
    args.both = true;

    let summary = pipeline::run(&args).expect("pipeline should succeed");

    assert!(tmp.join("test_both.MNV.tsv").exists(), "TSV should exist");
    assert!(tmp.join("test_both.MNV.vcf").exists(), "VCF should exist");
    assert!(summary.output_tsv.is_some());
    assert!(summary.output_vcf.is_some());

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1d: Dry run produces no output files
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_dry_run_no_files() {
    let tmp = temp_dir("e2e_dry");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_dry".to_string());
    args.dry_run = true;

    let summary = pipeline::run(&args).expect("pipeline should succeed");

    assert!(summary.dry_run, "summary should report dry-run");
    assert!(
        !tmp.join("test_dry.MNV.tsv").exists(),
        "TSV should NOT exist in dry-run"
    );
    assert!(
        summary.global.produced_variants >= 800,
        "dry-run should still count variants"
    );

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1e: Summary JSON output
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_summary_json() {
    let tmp = temp_dir("e2e_json");
    let summary_path = tmp.join("summary.json");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_json".to_string());
    args.summary_json = Some(summary_path.to_string_lossy().into());

    pipeline::run(&args).expect("pipeline should succeed");

    assert!(summary_path.exists(), "summary JSON should exist");
    let json_str = fs::read_to_string(&summary_path).expect("read JSON");
    let json: serde_json::Value = serde_json::from_str(&json_str).expect("parse JSON");
    assert_eq!(json["schema_version"], "1.0.0");
    assert!(json["global"]["produced_variants"].as_u64().unwrap() >= 800);
    assert!(json["timings"]["total_ms"].as_f64().unwrap() > 0.0);
    assert_eq!(json["contigs"].as_array().unwrap().len(), 1);

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1f: Exclude intergenic
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_exclude_intergenic() {
    let tmp = temp_dir("e2e_intergen");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_ig".to_string());

    // First run WITH intergenic
    let summary_with = pipeline::run(&args).expect("pipeline with intergenic");

    // Second run WITHOUT intergenic
    args.exclude_intergenic = true;
    args.output_prefix = Some("test_no_ig".to_string());
    let summary_without = pipeline::run(&args).expect("pipeline without intergenic");

    assert!(
        summary_with.global.intergenic_variants > 0,
        "should have some intergenic variants"
    );
    assert_eq!(
        summary_without.global.intergenic_variants, 0,
        "should have 0 intergenic when excluded"
    );
    assert!(
        summary_with.global.produced_variants > summary_without.global.produced_variants,
        "excluding intergenic should reduce total variants"
    );

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1g: Run manifest with checksums
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_run_manifest() {
    let tmp = temp_dir("e2e_manifest");
    let manifest_path = tmp.join("manifest.json");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_mf".to_string());
    args.run_manifest = Some(manifest_path.to_string_lossy().into());

    pipeline::run(&args).expect("pipeline should succeed");

    assert!(manifest_path.exists(), "manifest should exist");
    let json_str = fs::read_to_string(&manifest_path).expect("read manifest");
    let json: serde_json::Value = serde_json::from_str(&json_str).expect("parse manifest");
    assert_eq!(json["schema_version"], "1.0.0");
    assert!(json["tool_version"].as_str().unwrap().contains('.'));
    assert!(json["timestamp_unix"].as_u64().unwrap() > 0);
    // Output checksums
    let checksums = &json["output_checksums"];
    assert!(checksums["output_tsv_sha256"].is_string());

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1h: Normalize alleles flag
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_normalize_alleles() {
    let tmp = temp_dir("e2e_normalize");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_norm".to_string());
    args.normalize_alleles = true;

    let summary = pipeline::run(&args).expect("pipeline with normalize should succeed");
    assert!(summary.global.produced_variants >= 800);

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1i: Error JSON on invalid input
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_error_json_on_bad_input() {
    let tmp = temp_dir("e2e_error");
    let error_path = tmp.join("error.json");
    let mut args = base_args();
    args.vcf_file = Some("/nonexistent/file.vcf".to_string());
    args.tsv_file = None;
    args.error_json = Some(error_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());

    let result = pipeline::run(&args);
    assert!(result.is_err(), "should fail with nonexistent VCF");

    // Error JSON is written by main(), not pipeline::run(), so we just
    // verify the error is returned correctly
    let err = result.unwrap_err();
    assert!(!err.to_string().is_empty());

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1j: Keep original INFO fields
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_keep_original_info_in_vcf() {
    let tmp = temp_dir("e2e_keep_info");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_info".to_string());
    args.convert = true;
    args.keep_original_info = true;

    pipeline::run(&args).expect("pipeline should succeed");

    let vcf_path = tmp.join("test_info.MNV.vcf");
    let content = fs::read_to_string(&vcf_path).expect("read VCF");

    // The example VCF has ANN= from SnpEff — should be preserved
    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    let has_ann = data_lines.iter().any(|l| l.contains("ANN="));
    assert!(
        has_ann,
        "original ANN= INFO field should be preserved with --keep-original-info"
    );

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T1k: Split multiallelic
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_split_multiallelic_flag() {
    let tmp = temp_dir("e2e_split");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("test_split".to_string());
    args.split_multiallelic = true;

    let summary = pipeline::run(&args).expect("pipeline with split-multiallelic should succeed");
    assert!(summary.global.produced_variants >= 800);

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T3: --sample all with multi-sample VCF
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_sample_all_multisample_vcf() {
    let tmp = temp_dir("e2e_sample_all");

    // Create a minimal multi-sample VCF
    // Reference has A at pos 100, C at pos 200, G at pos 300 (1-based)
    let vcf_content = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B
chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:DP\t1/1:20\t0/0:15
chr1\t200\t.\tC\tG\t.\tPASS\t.\tGT:DP\t0/0:18\t1/1:22
chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT:DP\t1/1:25\t1/1:30
";
    let ref_content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGGACGTACGTACGTACGTACGT\n";
    let genes_content = "gene1\t98\t302\t+\n";

    let vcf_path = tmp.join("test.vcf");
    let ref_path = tmp.join("ref.fas");
    let genes_path = tmp.join("genes.txt");
    fs::write(&vcf_path, vcf_content).unwrap();
    fs::write(&ref_path, ref_content).unwrap();
    fs::write(&genes_path, genes_content).unwrap();

    let args = Args {
        vcf_file: Some(vcf_path.to_string_lossy().into()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: ref_path.to_string_lossy().into(),
        genes_file_tsv: Some(genes_path.to_string_lossy().into()),
        gff_file: None,
        gff_features_raw: None,
        sample: Some("all".to_string()),
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: None,
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: false,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        output_dir: Some(tmp.to_string_lossy().into()),
        translation_table: 11,
        output_prefix: None,
    };

    let summary = pipeline::run(&args).expect("--sample all should succeed");
    assert_eq!(summary.sample.as_deref(), Some("all"));

    // Should produce per-sample output files
    let sample_a_tsv = tmp.join("test.sample_SAMPLE_A.MNV.tsv");
    let sample_b_tsv = tmp.join("test.sample_SAMPLE_B.MNV.tsv");
    assert!(sample_a_tsv.exists(), "SAMPLE_A TSV should exist");
    assert!(sample_b_tsv.exists(), "SAMPLE_B TSV should exist");

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T3b: Single sample selection
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_sample_selection() {
    let tmp = temp_dir("e2e_sample_select");

    let vcf_content = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B
chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:DP\t1/1:20\t0/0:15
";
    let ref_content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGGACGTACGTACGTACGTACGT\n";
    let genes_content = "gene1\t98\t102\t+\n";

    let vcf_path = tmp.join("test.vcf");
    let ref_path = tmp.join("ref.fas");
    let genes_path = tmp.join("genes.txt");
    fs::write(&vcf_path, vcf_content).unwrap();
    fs::write(&ref_path, ref_content).unwrap();
    fs::write(&genes_path, genes_content).unwrap();

    let args = Args {
        vcf_file: Some(vcf_path.to_string_lossy().into()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: ref_path.to_string_lossy().into(),
        genes_file_tsv: Some(genes_path.to_string_lossy().into()),
        gff_file: None,
        gff_features_raw: None,
        sample: Some("SAMPLE_B".to_string()),
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: None,
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: false,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        output_dir: Some(tmp.to_string_lossy().into()),
        translation_table: 11,
        output_prefix: Some("selected".to_string()),
    };

    let summary = pipeline::run(&args).expect("--sample SAMPLE_B should succeed");
    assert_eq!(summary.sample.as_deref(), Some("SAMPLE_B"));
    assert!(tmp.join("selected.MNV.tsv").exists());

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T3c: Invalid sample name should error
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_invalid_sample_errors() {
    let tmp = temp_dir("e2e_bad_sample");

    let vcf_content = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A
chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t1/1
";
    let ref_content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGGACGTACGTACGTACGTACGT\n";
    let genes_content = "gene1\t98\t102\t+\n";

    fs::write(tmp.join("test.vcf"), vcf_content).unwrap();
    fs::write(tmp.join("ref.fas"), ref_content).unwrap();
    fs::write(tmp.join("genes.txt"), genes_content).unwrap();

    let args = Args {
        vcf_file: Some(tmp.join("test.vcf").to_string_lossy().into()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: tmp.join("ref.fas").to_string_lossy().into(),
        genes_file_tsv: Some(tmp.join("genes.txt").to_string_lossy().into()),
        gff_file: None,
        gff_features_raw: None,
        sample: Some("NONEXISTENT".to_string()),
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: None,
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: false,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        output_dir: Some(tmp.to_string_lossy().into()),
        translation_table: 11,
        output_prefix: None,
    };

    let result = pipeline::run(&args);
    assert!(result.is_err(), "nonexistent sample should fail");
    let err = result.unwrap_err().to_string();
    assert!(
        err.contains("NONEXISTENT") || err.contains("not found"),
        "error should mention sample name: {err}"
    );

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T4: Strict mode rejects missing ODP/OFREQ
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_strict_mode_with_minimal_vcf() {
    let tmp = temp_dir("e2e_strict");

    // Create a VCF without any depth/frequency fields — strict should reject it
    let vcf_content = "\
##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\tPASS\t.
";
    let ref_content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGGACGTACGTACGTACGTACGT\n";
    let genes_content = "gene1\t98\t102\t+\n";

    fs::write(tmp.join("test.vcf"), vcf_content).unwrap();
    fs::write(tmp.join("ref.fas"), ref_content).unwrap();
    fs::write(tmp.join("genes.txt"), genes_content).unwrap();

    let args = Args {
        vcf_file: Some(tmp.join("test.vcf").to_string_lossy().into()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: tmp.join("ref.fas").to_string_lossy().into(),
        genes_file_tsv: Some(tmp.join("genes.txt").to_string_lossy().into()),
        gff_file: None,
        gff_features_raw: None,
        sample: None,
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: None,
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: false,
        strict: true,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        output_dir: Some(tmp.to_string_lossy().into()),
        translation_table: 11,
        output_prefix: None,
    };

    let result = pipeline::run(&args);
    assert!(
        result.is_err(),
        "strict mode should fail when DP/FREQ missing"
    );
    let err = result.unwrap_err().to_string();
    assert!(
        err.contains("strict") || err.contains("ODP") || err.contains("missing"),
        "error should mention strict/metrics: {err}"
    );

    fs::remove_dir_all(&tmp).ok();
}

// ---------------------------------------------------------------------------
// T5: Pipeline config validation
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_vcf_gz_without_convert_errors() {
    let mut args = base_args();
    args.vcf_gz = true;
    // vcf_gz without convert or both should error
    let result = pipeline::run(&args);
    assert!(result.is_err(), "--vcf-gz without --convert should fail");
}

#[test]
fn test_e2e_index_without_vcf_gz_errors() {
    let mut args = base_args();
    args.index_vcf_gz = true;
    let result = pipeline::run(&args);
    assert!(
        result.is_err(),
        "--index-vcf-gz without --vcf-gz should fail"
    );
}

#[test]
fn test_e2e_bcf_without_convert_errors() {
    let mut args = base_args();
    args.bcf = true;
    let result = pipeline::run(&args);
    assert!(result.is_err(), "--bcf without --convert should fail");
}

#[test]
fn test_e2e_frequency_filters_require_bam() {
    let mut args = base_args();
    args.min_snp_frequency = 0.05;
    let err = pipeline::run(&args).expect_err("frequency filtering without BAM should fail");
    assert!(err.to_string().contains("require --bam"));
}

// ---------------------------------------------------------------------------
// T6: Fast text parser produces same output as htslib
// ---------------------------------------------------------------------------

#[test]
fn test_fast_parser_matches_htslib() {
    use get_mnv::io::{vcf, vcf_fast};

    let ex = example_dir();
    let vcf_file = ex.join("G35894.var.snp.vcf").to_string_lossy().to_string();

    // Fast text parser
    let fast_result = vcf_fast::load_vcf_text(&vcf_file, None, false, false, false)
        .expect("fast parser should succeed");

    // htslib parser
    let htslib_result = vcf::load_vcf_positions_by_contig(&vcf_file, None, false, false, false)
        .expect("htslib parser should succeed");

    // Same contigs
    assert_eq!(
        fast_result
            .keys()
            .collect::<std::collections::BTreeSet<_>>(),
        htslib_result
            .keys()
            .collect::<std::collections::BTreeSet<_>>(),
        "contigs should match"
    );

    for contig in fast_result.keys() {
        let fast_positions = &fast_result[contig];
        let htslib_positions = &htslib_result[contig];
        assert_eq!(
            fast_positions.len(),
            htslib_positions.len(),
            "contig '{}': position count mismatch (fast={}, htslib={})",
            contig,
            fast_positions.len(),
            htslib_positions.len()
        );
        for (i, (fp, hp)) in fast_positions
            .iter()
            .zip(htslib_positions.iter())
            .enumerate()
        {
            assert_eq!(
                fp.position, hp.position,
                "contig '{}' pos {}: position mismatch",
                contig, i
            );
            assert_eq!(
                fp.ref_allele, hp.ref_allele,
                "contig '{}' pos {}: ref mismatch",
                contig, i
            );
            assert_eq!(
                fp.alt_allele, hp.alt_allele,
                "contig '{}' pos {}: alt mismatch",
                contig, i
            );
            assert_eq!(
                fp.original_dp, hp.original_dp,
                "contig '{}' pos {} ({}): DP mismatch (fast={:?}, htslib={:?})",
                contig, i, fp.position, fp.original_dp, hp.original_dp
            );
            // Compare freq with tolerance for float rounding
            match (fp.original_freq, hp.original_freq) {
                (Some(f), Some(h)) => {
                    assert!(
                        (f - h).abs() < 1e-6,
                        "contig '{}' pos {} ({}): FREQ mismatch (fast={}, htslib={})",
                        contig,
                        i,
                        fp.position,
                        f,
                        h
                    );
                }
                (None, None) => {}
                _ => panic!(
                    "contig '{}' pos {} ({}): FREQ presence mismatch (fast={:?}, htslib={:?})",
                    contig, i, fp.position, fp.original_freq, hp.original_freq
                ),
            }
        }
    }
}

// ---------------------------------------------------------------------------
// T6b: Fast parser with --keep-original-info matches htslib
// ---------------------------------------------------------------------------

#[test]
fn test_fast_parser_keep_info_matches_htslib() {
    use get_mnv::io::{vcf, vcf_fast};

    let ex = example_dir();
    let vcf_file = ex.join("G35894.var.snp.vcf").to_string_lossy().to_string();

    let fast_result = vcf_fast::load_vcf_text(&vcf_file, None, false, false, true)
        .expect("fast parser should succeed");
    let htslib_result = vcf::load_vcf_positions_by_contig(&vcf_file, None, false, false, true)
        .expect("htslib parser should succeed");

    for contig in fast_result.keys() {
        let fast_positions = &fast_result[contig];
        let htslib_positions = &htslib_result[contig];
        for (i, (fp, hp)) in fast_positions
            .iter()
            .zip(htslib_positions.iter())
            .enumerate()
        {
            // Both should have or not have original_info
            assert_eq!(
                fp.original_info.is_some(),
                hp.original_info.is_some(),
                "contig '{}' pos {} ({}): original_info presence mismatch",
                contig,
                i,
                fp.position
            );
        }
    }
}

// ---------------------------------------------------------------------------
// T7: --translation-table support
// ---------------------------------------------------------------------------

#[test]
fn test_invalid_translation_table_fails() {
    let ex = example_dir();
    let out = std::env::temp_dir().join("get_mnv_test_invalid_tt");

    let args = Args {
        vcf_file: Some(ex.join("G35894.var.snp.vcf").to_string_lossy().to_string()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: ex.join("MTB_ancestor.fas").to_string_lossy().to_string(),
        genes_file_tsv: Some(ex.join("anot_genes.txt").to_string_lossy().to_string()),
        gff_file: None,
        gff_features_raw: None,
        sample: None,
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: Some(1),
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: true,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        translation_table: 99, // invalid
        output_dir: Some(out.to_string_lossy().to_string()),
        output_prefix: None,
    };

    let result = pipeline::run(&args);
    assert!(result.is_err(), "invalid translation table should fail");
    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("translation-table") || err_msg.contains("not supported"),
        "Error should mention translation table: {err_msg}"
    );
}

#[test]
fn test_translation_table_1_standard() {
    // Table 1 (Standard) should produce valid results — same as default for MTB
    let ex = example_dir();
    let out = std::env::temp_dir().join("get_mnv_test_table1");

    let args = Args {
        vcf_file: Some(ex.join("G35894.var.snp.vcf").to_string_lossy().to_string()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: ex.join("MTB_ancestor.fas").to_string_lossy().to_string(),
        genes_file_tsv: Some(ex.join("anot_genes.txt").to_string_lossy().to_string()),
        gff_file: None,
        gff_features_raw: None,
        sample: None,
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: Some(1),
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: true,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        translation_table: 1, // Standard
        output_dir: Some(out.to_string_lossy().to_string()),
        output_prefix: None,
    };

    let summary = pipeline::run(&args).expect("table 1 should succeed");
    // Tables 1 and 11 produce identical results for standard sense codons
    assert!(summary.global.produced_variants > 0);
}

// ---------------------------------------------------------------------------
// T8: Robustness — edge case inputs
// ---------------------------------------------------------------------------

#[test]
fn test_empty_vcf_no_records() {
    // A VCF with header but no data records should error (no contigs)
    use std::io::Write;
    let dir = std::env::temp_dir().join("get_mnv_robustness_empty_vcf");
    let _ = std::fs::create_dir_all(&dir);

    let vcf_path = dir.join("empty.vcf");
    {
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    }

    let fasta_path = dir.join("ref.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">chr1\nACGT").unwrap();
    }

    let ex = example_dir();
    let args = Args {
        vcf_file: Some(vcf_path.to_string_lossy().to_string()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: fasta_path.to_string_lossy().to_string(),
        genes_file_tsv: Some(ex.join("anot_genes.txt").to_string_lossy().to_string()),
        gff_file: None,
        gff_features_raw: None,
        sample: None,
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: Some(1),
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: true,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        translation_table: 11,
        output_dir: Some(dir.to_string_lossy().to_string()),
        output_prefix: None,
    };

    let result = pipeline::run(&args);
    assert!(result.is_err(), "Empty VCF should fail with no contigs");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn test_truncated_vcf_record() {
    // A VCF with a truncated record (fewer than 8 columns) should error
    use std::io::Write;
    let dir = std::env::temp_dir().join("get_mnv_robustness_truncated");
    let _ = std::fs::create_dir_all(&dir);

    let vcf_path = dir.join("truncated.vcf");
    {
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        // Valid record
        writeln!(f, "chr1\t100\t.\tA\tT\t.\tPASS\t.").unwrap();
        // Truncated record (only 3 columns)
        writeln!(f, "chr1\t200\t.").unwrap();
    }

    let fasta_path = dir.join("ref.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">chr1").unwrap();
        writeln!(f, "{}", "A".repeat(300)).unwrap();
    }

    let args = Args {
        vcf_file: Some(vcf_path.to_string_lossy().to_string()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: fasta_path.to_string_lossy().to_string(),
        genes_file_tsv: Some(
            example_dir()
                .join("anot_genes.txt")
                .to_string_lossy()
                .to_string(),
        ),
        gff_file: None,
        gff_features_raw: None,
        sample: None,
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: Some(1),
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: true,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        translation_table: 11,
        output_dir: Some(dir.to_string_lossy().to_string()),
        output_prefix: None,
    };

    let result = pipeline::run(&args);
    assert!(result.is_err(), "Truncated VCF record should error");
    let err = result.unwrap_err().to_string();
    assert!(
        err.contains("expected at least 8 columns") || err.contains("column"),
        "Error should mention columns: {err}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn test_vcf_no_header() {
    // A VCF file with no #CHROM header should error
    use std::io::Write;
    let dir = std::env::temp_dir().join("get_mnv_robustness_no_header");
    let _ = std::fs::create_dir_all(&dir);

    let vcf_path = dir.join("noheader.vcf");
    {
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        // No #CHROM line, directly data
        writeln!(f, "chr1\t100\t.\tA\tT\t.\tPASS\t.").unwrap();
    }

    let fasta_path = dir.join("ref.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">chr1\nACGT").unwrap();
    }

    let args = Args {
        vcf_file: Some(vcf_path.to_string_lossy().to_string()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: fasta_path.to_string_lossy().to_string(),
        genes_file_tsv: Some(
            example_dir()
                .join("anot_genes.txt")
                .to_string_lossy()
                .to_string(),
        ),
        gff_file: None,
        gff_features_raw: None,
        sample: None,
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: Some(1),
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: true,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: None,
        run_manifest: None,
        convert: false,
        both: false,
        translation_table: 11,
        output_dir: Some(dir.to_string_lossy().to_string()),
        output_prefix: None,
    };

    let result = pipeline::run(&args);
    assert!(result.is_err(), "VCF without header should fail");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn test_error_json_written_on_failure() {
    // Verify that --error-json produces a valid JSON file when the pipeline fails
    use std::io::Write;
    let dir = std::env::temp_dir().join("get_mnv_robustness_error_json");
    let _ = std::fs::create_dir_all(&dir);

    let vcf_path = dir.join("bad.vcf");
    {
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        // No #CHROM header → will error
    }

    let fasta_path = dir.join("ref.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">chr1\nACGT").unwrap();
    }

    let error_json_path = dir.join("error.json");

    let args = Args {
        vcf_file: Some(vcf_path.to_string_lossy().to_string()),
        tsv_file: None,
        input_format: VariantInputFormat::Auto,
        bam_file: None,
        fasta_file: fasta_path.to_string_lossy().to_string(),
        genes_file_tsv: Some(
            example_dir()
                .join("anot_genes.txt")
                .to_string_lossy()
                .to_string(),
        ),
        gff_file: None,
        gff_features_raw: None,
        sample: None,
        chrom: None,
        normalize_alleles: false,
        min_quality: 20,
        min_mapq: 0,
        threads: Some(1),
        min_snp_reads: 0,
        min_snp_frequency: 0.0,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        dry_run: true,
        strict: false,
        split_multiallelic: false,
        emit_filtered: false,
        vcf_gz: false,
        index_vcf_gz: false,
        strand_bias_info: false,
        keep_original_info: false,
        exclude_intergenic: false,
        bcf: false,
        summary_json: None,
        error_json: Some(error_json_path.to_string_lossy().to_string()),
        run_manifest: None,
        convert: false,
        both: false,
        translation_table: 11,
        output_dir: Some(dir.to_string_lossy().to_string()),
        output_prefix: None,
    };

    // Pipeline should fail
    let result = pipeline::run(&args);
    assert!(result.is_err());

    // But error-json should NOT be written by pipeline::run itself —
    // it's written by the main() wrapper. Verify the error propagates cleanly.
    let err = result.unwrap_err();
    let json_str = get_mnv::error::error_to_json(&err);
    let parsed: serde_json::Value = serde_json::from_str(&json_str).unwrap();
    assert!(
        parsed.get("code").is_some(),
        "Error JSON should have 'code' field"
    );
    assert!(
        parsed.get("message").is_some(),
        "Error JSON should have 'message' field"
    );
    assert!(
        parsed.get("schema_version").is_some(),
        "Error JSON should have 'schema_version'"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
