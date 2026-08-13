use super::{write_sorted_vcf_entries, VcfWriter, VcfWriterConfig};
use crate::variants::{ChangeType, VariantInfo, VariantType};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Read;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_temp_path(prefix: &str, suffix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before UNIX_EPOCH")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}_{suffix}"))
}

#[test]
fn test_write_sorted_vcf_entries_orders_by_position() {
    let path = unique_temp_path("get_mnv_vcf_sort", "vcf");
    let mut file = File::create(&path).expect("failed to create temp file");
    write_sorted_vcf_entries(
        &mut file,
        vec![
            (20, "chr\t20\t.\tA\tT\t.\tPASS\tTYPE=SNP".to_string()),
            (10, "chr\t10\t.\tA\tC\t.\tPASS\tTYPE=SNP".to_string()),
        ],
    )
    .expect("failed to write sorted entries");
    drop(file);

    let mut contents = String::new();
    File::open(&path)
        .expect("failed to open temp file")
        .read_to_string(&mut contents)
        .expect("failed to read temp file");
    let lines: Vec<&str> = contents.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].contains("\t10\t"));
    assert!(lines[1].contains("\t20\t"));
    let _ = fs::remove_file(path);
}

fn snp_variant_with_reads(support_reads: usize, depth: usize) -> VariantInfo {
    VariantInfo {
        chrom: "chr1".to_string(),
        gene: "geneA".to_string(),
        positions: vec![10],
        ref_bases: vec!["A".to_string()],
        base_changes: vec!["T".to_string()],
        aa_changes: vec!["Ala1Val".to_string()],
        snp_aa_changes: vec!["Ala1Val".to_string()],
        aa_changes_local: vec!["Ala1Val".to_string()],
        snp_aa_changes_local: vec!["Ala1Val".to_string()],
        variant_type: VariantType::Snp,
        change_type: ChangeType::NonSynonymous,
        snp_reads: Some(vec![support_reads]),
        snp_forward_reads: Some(vec![support_reads]),
        snp_reverse_reads: Some(vec![0]),
        mnv_reads: None,
        mnv_forward_reads: None,
        mnv_reverse_reads: None,
        mnv_total_reads: None,
        total_reads: Some(vec![depth]),
        total_forward_reads: Some(vec![depth]),
        total_reverse_reads: Some(vec![0]),
        mnv_total_forward_reads: None,
        mnv_total_reverse_reads: None,
        mnv_phasing_reads: None,
        ref_codon: Some("AAA".to_string()),
        snp_codon: Some("TAA".to_string()),
        mnv_codon: None,
        original_dp: None,
        original_freq: None,
        original_info: None,
        event_class: Some("snp".to_string()),
        event_components: vec!["SNV:10:A>T".to_string()],
        annotations: crate::variants::VariantAnnotations::default(),
    }
}

/// A codon-level call whose locus the BAM never reached: every count is zero,
/// which is the absence of evidence, not a failed threshold.
fn snp_mnv_variant_without_read_support() -> VariantInfo {
    VariantInfo {
        chrom: "chr1".to_string(),
        gene: "geneA".to_string(),
        positions: vec![10, 12],
        ref_bases: vec!["A".to_string(), "A".to_string()],
        base_changes: vec!["T".to_string(), "G".to_string()],
        aa_changes: vec!["Ala1Val".to_string()],
        snp_aa_changes: vec!["Ala1Val".to_string(), "Ala1Gly".to_string()],
        aa_changes_local: vec!["Ala1Val".to_string()],
        snp_aa_changes_local: vec!["Ala1Val".to_string(), "Ala1Gly".to_string()],
        variant_type: VariantType::SnpMnv,
        change_type: ChangeType::NonSynonymous,
        snp_reads: Some(vec![0, 0]),
        snp_forward_reads: Some(vec![0, 0]),
        snp_reverse_reads: Some(vec![0, 0]),
        mnv_reads: Some(0),
        mnv_forward_reads: Some(0),
        mnv_reverse_reads: Some(0),
        mnv_total_reads: Some(0),
        total_reads: Some(vec![0, 0]),
        total_forward_reads: Some(vec![0, 0]),
        total_reverse_reads: Some(vec![0, 0]),
        mnv_total_forward_reads: Some(0),
        mnv_total_reverse_reads: Some(0),
        mnv_phasing_reads: None,
        ref_codon: Some("AAA".to_string()),
        snp_codon: Some("TAA".to_string()),
        mnv_codon: Some("TAG".to_string()),
        original_dp: None,
        original_freq: None,
        original_info: None,
        event_class: Some("mnv".to_string()),
        event_components: vec!["SNV:10:A>T".to_string(), "SNV:12:A>G".to_string()],
        annotations: crate::variants::VariantAnnotations::default(),
    }
}

fn write_one_variant(variant: &VariantInfo, min_reads: usize, emit_filtered: bool) -> String {
    let stem = unique_temp_path("get_mnv_vcf_zero_support", "out");
    let stem_str = stem.to_string_lossy().into_owned();
    let out_path = format!("{stem_str}.MNV.vcf");
    let mut writer = VcfWriter::new(VcfWriterConfig {
        filename: &stem_str,
        bam_provided: true,
        min_snp_reads: min_reads,
        min_snp_frequency: 0.0,
        min_mnv_reads: min_reads,
        min_mnv_frequency: 0.0,
        min_quality: 20,
        min_mapq: 0,
        command_line: "get_mnv test",
        contigs: &["chr1".to_string()],
        bgzf_output: false,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        emit_filtered,
        include_strand_bias_info: false,
        original_info_headers: &[],
    })
    .expect("writer should build");
    let mut references = HashMap::new();
    references.insert("chr1".to_string(), "AAAAAAAAAAAAAAAAAAAA".to_string());
    writer
        .write_variants(std::slice::from_ref(variant), &references)
        .expect("variant should write");
    drop(writer);

    let mut contents = String::new();
    File::open(&out_path)
        .expect("failed to open VCF")
        .read_to_string(&mut contents)
        .expect("failed to read VCF");
    let _ = fs::remove_file(out_path);
    contents
}

/// A `--bam` run used to drop any allele with zero supporting reads before the
/// filter gate could see it, so the VCF silently lost calls the TSV still wrote
/// and `--emit-filtered` could not recover them. Zero support with no threshold
/// set is not a failed filter: the record belongs in the file.
#[test]
fn zero_support_records_are_written_when_no_threshold_asks_otherwise() {
    let contents = write_one_variant(&snp_mnv_variant_without_read_support(), 0, false);
    let records: Vec<&str> = contents
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(
        records.len(),
        3,
        "expected both SNP records and the MNV record, got:\n{contents}"
    );
    assert!(records.iter().all(|line| line.contains("\tPASS\t")));
}

/// The same row under a threshold it cannot meet: skipped by default, and
/// present with a `FILTER` tag when the user asks to see what was dropped.
/// Neither was reachable while the zero-support early exit stood.
#[test]
fn zero_support_records_obey_the_threshold_and_emit_filtered() {
    let variant = snp_mnv_variant_without_read_support();

    let skipped = write_one_variant(&variant, 2, false);
    assert!(
        !skipped.lines().any(|line| !line.starts_with('#')),
        "a threshold of 2 should skip every zero-support record, got:\n{skipped}"
    );

    let tagged = write_one_variant(&variant, 2, true);
    let records: Vec<&str> = tagged
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(records.len(), 3, "got:\n{tagged}");
    assert!(records.iter().all(|line| line.contains("LowSupport")));
}

#[test]
fn test_low_frequency_filter_marks_snp_when_emit_filtered() {
    let stem = unique_temp_path("get_mnv_vcf_frequency", "out");
    let stem_str = stem.to_string_lossy().into_owned();
    let out_path = format!("{stem_str}.MNV.vcf");
    let mut writer = VcfWriter::new(VcfWriterConfig {
        filename: &stem_str,
        bam_provided: true,
        min_snp_reads: 0,
        min_snp_frequency: 0.5,
        min_mnv_reads: 0,
        min_mnv_frequency: 0.0,
        min_quality: 20,
        min_mapq: 0,
        command_line: "get_mnv test",
        contigs: &["chr1".to_string()],
        bgzf_output: false,
        min_snp_strand_reads: 0,
        min_mnv_strand_reads: 0,
        min_strand_bias_p: 0.0,
        emit_filtered: true,
        include_strand_bias_info: false,
        original_info_headers: &[],
    })
    .expect("writer should build");
    writer
        .write_variants(&[snp_variant_with_reads(2, 10)], &HashMap::new())
        .expect("variant should write");
    drop(writer);

    let mut contents = String::new();
    File::open(&out_path)
        .expect("failed to open VCF")
        .read_to_string(&mut contents)
        .expect("failed to read VCF");
    assert!(contents.contains("##FILTER=<ID=LowFrequency"));
    assert!(contents.contains("\tLowFrequency\t"));
    assert!(contents.contains("FREQ=0.2000"));
    let _ = fs::remove_file(out_path);
}
