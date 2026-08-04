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
