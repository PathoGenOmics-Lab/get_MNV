//! End-to-end integration tests using the example dataset.
//!
//! These tests run the full pipeline (`pipeline::run`) against the bundled
//! `example/` data and verify output correctness at a structural level.

use get_mnv::cli::{Args, VariantInputFormat};
use get_mnv::pipeline;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
        convert: false,
        both: false,
        output_dir: None,
        translation_table: 11,
        output_prefix: None,
    }
}

struct SyntheticRead {
    name: &'static str,
    start: usize,
    cigar: &'static str,
    sequence: &'static str,
}

fn parse_test_cigar(cigar: &str) -> noodles::sam::alignment::record_buf::Cigar {
    use noodles::sam::alignment::record::cigar::op::{Kind, Op};

    let mut ops = Vec::new();
    let mut len = String::new();
    for symbol in cigar.chars() {
        if symbol.is_ascii_digit() {
            len.push(symbol);
            continue;
        }
        let count = len
            .parse::<usize>()
            .expect("synthetic CIGAR operation length");
        len.clear();
        let kind = match symbol {
            'M' => Kind::Match,
            'I' => Kind::Insertion,
            'D' => Kind::Deletion,
            'N' => Kind::Skip,
            'S' => Kind::SoftClip,
            'H' => Kind::HardClip,
            'P' => Kind::Pad,
            '=' => Kind::SequenceMatch,
            'X' => Kind::SequenceMismatch,
            _ => panic!("unsupported synthetic CIGAR operator: {symbol}"),
        };
        ops.push(Op::new(kind, count));
    }
    assert!(len.is_empty(), "dangling synthetic CIGAR length");
    ops.into_iter().collect()
}

fn write_synthetic_bam(path: &Path, reference_len: usize, reads: &[SyntheticRead]) {
    write_synthetic_bam_q(path, reference_len, reads, Some(40));
}

fn write_synthetic_bam_q(
    path: &Path,
    reference_len: usize,
    reads: &[SyntheticRead],
    quality: Option<u8>,
) {
    use noodles::core::Position;
    use noodles::sam::{
        self,
        alignment::{
            io::Write as _,
            record::{Flags, MappingQuality},
            record_buf::{QualityScores, Sequence},
            RecordBuf,
        },
        header::record::value::{
            map::{
                self,
                header::{sort_order::COORDINATE, tag::SORT_ORDER},
                ReferenceSequence,
            },
            Map,
        },
    };
    use std::num::NonZeroUsize;

    let header = sam::Header::builder()
        .set_header(
            Map::<map::Header>::builder()
                .insert(SORT_ORDER, COORDINATE)
                .build()
                .expect("synthetic BAM header"),
        )
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(
                NonZeroUsize::new(reference_len).expect("non-empty reference"),
            ),
        )
        .build();

    let file = fs::File::create(path).expect("create synthetic BAM");
    let mut writer = noodles::bam::io::Writer::new(file);
    writer.write_header(&header).expect("write BAM header");

    for read in reads {
        let cigar = parse_test_cigar(read.cigar);
        assert_eq!(
            cigar.read_length(),
            read.sequence.len(),
            "synthetic read {} has a sequence/CIGAR length mismatch",
            read.name
        );
        let qualities = match quality {
            Some(q) => QualityScores::from(vec![q; read.sequence.len()]),
            None => QualityScores::default(), // SAM QUAL='*' (no per-base qualities)
        };
        let record = RecordBuf::builder()
            .set_name(read.name)
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(read.start).expect("alignment start"))
            .set_mapping_quality(MappingQuality::new(60).expect("mapping quality"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(read.sequence.as_bytes()))
            .set_quality_scores(qualities)
            .build();
        writer
            .write_alignment_record(&header, &record)
            .expect("write BAM record");
    }
    writer.try_finish().expect("finish synthetic BAM");

    let index = noodles::bam::fs::index(path).expect("index synthetic BAM");
    noodles::bam::bai::fs::write(path.with_extension("bam.bai"), &index)
        .expect("write synthetic BAI");
}

fn write_phase_reference_files(tmp: &Path) -> (PathBuf, PathBuf) {
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");
    fs::write(&ref_path, ">chr1\nATGAAATTTCCC\n").unwrap();
    fs::write(&genes_path, "gene1\t1\t12\t+\n").unwrap();
    (ref_path, genes_path)
}

fn read_tsv_rows(path: &Path) -> Vec<HashMap<String, String>> {
    let content = fs::read_to_string(path).expect("read TSV output");
    let mut lines = content.lines();
    let header = lines
        .next()
        .expect("TSV header")
        .split('\t')
        .map(str::to_string)
        .collect::<Vec<_>>();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(key, value)| (key.clone(), value.to_string()))
                .collect()
        })
        .collect()
}

fn find_row<'a>(
    rows: &'a [HashMap<String, String>],
    event_class: &str,
    reference_bases: &str,
    base_changes: &str,
) -> Option<&'a HashMap<String, String>> {
    rows.iter().find(|row| {
        row.get("Event Class").map(String::as_str) == Some(event_class)
            && row.get("Reference Bases").map(String::as_str) == Some(reference_bases)
            && row.get("Base Changes").map(String::as_str) == Some(base_changes)
    })
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
    assert_eq!(summary.global.snp_records_in_vcf, 4);
    assert_eq!(summary.global.produced_variants, 2);
    assert_eq!(summary.global.snp_mnv_variants, 1);
    assert_eq!(summary.global.indel_variants, 1);

    let tsv_content = fs::read_to_string(tmp.join("ivar_out.MNV.tsv")).expect("read TSV");
    assert!(tsv_content.contains("4, 5, 6"));
    assert!(tsv_content.contains("SNP/MNV"));
    assert!(tsv_content.contains("insertion"));
    assert!(tsv_content.contains("INS:8:+A"));

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_vcf_phased_insertion_haplotype_is_supported_by_same_read() {
    let tmp = temp_dir("e2e_phased_ins_vcf");
    let (ref_path, genes_path) = write_phase_reference_files(&tmp);
    let vcf_path = tmp.join("phase.vcf");
    let bam_path = tmp.join("phase.bam");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=12>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tC\t.\tPASS\t.\n\
chr1\t6\t.\tA\tAT\t.\tPASS\t.\n",
    )
    .unwrap();
    write_synthetic_bam(
        &bam_path,
        12,
        &[
            SyntheticRead {
                name: "phased_ins",
                start: 1,
                cigar: "6M1I6M",
                sequence: "ATGCAATTTTCCC",
            },
            SyntheticRead {
                name: "ins_only",
                start: 1,
                cigar: "6M1I6M",
                sequence: "ATGAAATTTTCCC",
            },
            SyntheticRead {
                name: "snv_only",
                start: 1,
                cigar: "12M",
                sequence: "ATGCAATTTCCC",
            },
            SyntheticRead {
                name: "ref",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
        ],
    );

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.tsv_file = None;
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.threads = Some(1);
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("phase".to_string());

    let summary = pipeline::run(&args).expect("phased insertion pipeline should succeed");
    assert_eq!(summary.global.snp_records_in_vcf, 2);
    assert!(
        summary.global.indel_variants >= 2,
        "expected original insertion plus phased haplotype"
    );

    let rows = read_tsv_rows(&tmp.join("phase.MNV.tsv"));
    let compound = find_row(&rows, "complex_indel", "AAA", "CAAT")
        .expect("compound insertion+SNV haplotype row");
    assert_eq!(compound["Positions"], "4");
    assert_eq!(compound["Event Reads"], "1");
    assert_eq!(compound["Event Depth"], "4");
    assert!(compound["Event Components"].contains("SNV:4:A>C"));
    assert!(compound["Event Components"].contains("INS:6:+T"));

    let insertion = find_row(&rows, "insertion", "A", "AT").expect("original insertion row");
    assert_eq!(insertion["Event Reads"], "2");
    assert_eq!(insertion["Event Depth"], "4");

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_vcf_phased_mnv_plus_insertion_emits_full_complex_haplotype() {
    let tmp = temp_dir("e2e_phased_mnv_ins_vcf");
    let (ref_path, genes_path) = write_phase_reference_files(&tmp);
    let vcf_path = tmp.join("phase_mnv.vcf");
    let bam_path = tmp.join("phase_mnv.bam");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=12>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tC\t.\tPASS\t.\n\
chr1\t5\t.\tA\tG\t.\tPASS\t.\n\
chr1\t6\t.\tA\tAT\t.\tPASS\t.\n",
    )
    .unwrap();
    write_synthetic_bam(
        &bam_path,
        12,
        &[
            SyntheticRead {
                name: "full_complex",
                start: 1,
                cigar: "6M1I6M",
                sequence: "ATGCGATTTTCCC",
            },
            SyntheticRead {
                name: "mnv_only",
                start: 1,
                cigar: "12M",
                sequence: "ATGCGATTTCCC",
            },
            SyntheticRead {
                name: "ins_only",
                start: 1,
                cigar: "6M1I6M",
                sequence: "ATGAAATTTTCCC",
            },
            SyntheticRead {
                name: "ref",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
        ],
    );

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.tsv_file = None;
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.threads = Some(1);
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("phase_mnv".to_string());

    let summary = pipeline::run(&args).expect("phased MNV+insertion pipeline should succeed");
    assert_eq!(summary.global.snp_records_in_vcf, 3);
    assert!(
        summary.global.indel_variants >= 2,
        "expected original insertion plus phased complex haplotype"
    );

    let rows = read_tsv_rows(&tmp.join("phase_mnv.MNV.tsv"));
    let mnv_row = find_row(&rows, "mnv", "A, A", "C, G").expect("codon MNV row");
    assert_eq!(mnv_row["Positions"], "4, 5");
    assert_eq!(mnv_row["Change Type"], "Non-synonymous");
    assert_ne!(mnv_row["AA Changes"], "Unknown");

    let compound =
        find_row(&rows, "complex_indel", "AAA", "CGAT").expect("full MNV+insertion haplotype row");
    assert_eq!(compound["Positions"], "4");
    assert_eq!(compound["Event Reads"], "1");
    assert_eq!(compound["Event Depth"], "4");
    assert!(compound["Event Components"].contains("SNV:4:A>C"));
    assert!(compound["Event Components"].contains("SNV:5:A>G"));
    assert!(compound["Event Components"].contains("INS:6:+T"));

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_vcf_phased_two_indels_require_cigar_components() {
    let tmp = temp_dir("e2e_phased_two_indels_vcf");
    let (ref_path, genes_path) = write_phase_reference_files(&tmp);
    let vcf_path = tmp.join("phase_two_indels.vcf");
    let bam_path = tmp.join("phase_two_indels.bam");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=12>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tAG\t.\tPASS\t.\n\
chr1\t5\t.\tAA\tA\t.\tPASS\t.\n",
    )
    .unwrap();
    write_synthetic_bam(
        &bam_path,
        12,
        &[
            SyntheticRead {
                name: "full_complex",
                start: 1,
                cigar: "4M1I1M1D6M",
                sequence: "ATGAGATTTCCC",
            },
            SyntheticRead {
                name: "sequence_mimic_snp",
                start: 1,
                cigar: "12M",
                sequence: "ATGAGATTTCCC",
            },
            SyntheticRead {
                name: "ins_only",
                start: 1,
                cigar: "4M1I8M",
                sequence: "ATGAGAATTTCCC",
            },
            SyntheticRead {
                name: "del_only",
                start: 1,
                cigar: "5M1D6M",
                sequence: "ATGAATTTCCC",
            },
            SyntheticRead {
                name: "ref",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
        ],
    );

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.tsv_file = None;
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.threads = Some(1);
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("phase_two_indels".to_string());

    let summary = pipeline::run(&args).expect("phased two-indel pipeline should succeed");
    assert_eq!(summary.global.snp_records_in_vcf, 2);
    assert!(
        summary.global.indel_variants >= 3,
        "expected original indels plus phased complex haplotype"
    );

    let rows = read_tsv_rows(&tmp.join("phase_two_indels.MNV.tsv"));
    let compound = find_row(&rows, "complex_indel", "AAA", "AGA").expect("two-indel haplotype row");
    assert_eq!(compound["Event Reads"], "1");
    assert_eq!(compound["Event Depth"], "5");
    assert!(compound["Event Components"].contains("INS:4:+G"));
    assert!(compound["Event Components"].contains("DEL:6:A"));

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_vcf_phased_deletion_and_insertion_in_adjacent_codons_combine() {
    // A deletion in one codon and an insertion in the next codon, carried on the
    // same reads, are emitted as a single exact `complex_indel` haplotype, while
    // the individual indel rows are preserved. This is the cross-codon case:
    // codon 2 (pos 4-6) carries the deletion of pos 6, codon 3 (pos 7-9) carries
    // the insertion after pos 7. Only the read that carries BOTH supports the
    // combined event.
    let tmp = temp_dir("e2e_phased_crosscodon_delins");
    let (ref_path, genes_path) = write_phase_reference_files(&tmp);
    let vcf_path = tmp.join("crosscodon.vcf");
    let bam_path = tmp.join("crosscodon.bam");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=12>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t5\t.\tAA\tA\t.\tPASS\t.\n\
chr1\t7\t.\tT\tTG\t.\tPASS\t.\n",
    )
    .unwrap();
    write_synthetic_bam(
        &bam_path,
        12,
        &[
            // Carries the codon-2 deletion (pos 6) AND the codon-3 insertion
            // (G after pos 7) on the same molecule.
            SyntheticRead {
                name: "both",
                start: 1,
                cigar: "5M1D1M1I5M",
                sequence: "ATGAATGTTCCC",
            },
            // Deletion only.
            SyntheticRead {
                name: "del_only",
                start: 1,
                cigar: "5M1D6M",
                sequence: "ATGAATTTCCC",
            },
            // Insertion only.
            SyntheticRead {
                name: "ins_only",
                start: 1,
                cigar: "7M1I5M",
                sequence: "ATGAAATGTTCCC",
            },
            // Plain reference.
            SyntheticRead {
                name: "ref",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
        ],
    );

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.tsv_file = None;
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.threads = Some(1);
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("crosscodon".to_string());

    let summary = pipeline::run(&args).expect("cross-codon del+ins pipeline should succeed");
    assert_eq!(summary.global.snp_records_in_vcf, 2);
    assert!(
        summary.global.indel_variants >= 3,
        "expected the two original indels plus the phased cross-codon haplotype"
    );

    let rows = read_tsv_rows(&tmp.join("crosscodon.MNV.tsv"));

    // The combined complex_indel row carries both indel components and is
    // supported only by the single read that carries both indels.
    let compound = rows
        .iter()
        .find(|row| {
            row.get("Event Class").map(String::as_str) == Some("complex_indel")
                && row
                    .get("Event Components")
                    .is_some_and(|c| c.contains("DEL:6:A") && c.contains("INS:7:+G"))
        })
        .expect("a complex_indel row combining the codon-2 deletion and codon-3 insertion");
    assert_eq!(
        compound["Event Reads"], "1",
        "only the read carrying both the deletion and the insertion supports the combined event"
    );

    // The individual indels are preserved as their own rows.
    assert!(
        find_row(&rows, "deletion", "AA", "A").is_some(),
        "the standalone codon-2 deletion row should still be present"
    );
    assert!(
        find_row(&rows, "insertion", "T", "TG").is_some(),
        "the standalone codon-3 insertion row should still be present"
    );

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_vcf_close_indel_and_snv_do_not_phase_without_shared_read() {
    let tmp = temp_dir("e2e_unphased_ins_vcf");
    let (ref_path, genes_path) = write_phase_reference_files(&tmp);
    let vcf_path = tmp.join("unphased.vcf");
    let bam_path = tmp.join("unphased.bam");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=12>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tC\t.\tPASS\t.\n\
chr1\t6\t.\tA\tAT\t.\tPASS\t.\n",
    )
    .unwrap();
    write_synthetic_bam(
        &bam_path,
        12,
        &[
            SyntheticRead {
                name: "ins_only",
                start: 1,
                cigar: "6M1I6M",
                sequence: "ATGAAATTTTCCC",
            },
            SyntheticRead {
                name: "snv_only",
                start: 1,
                cigar: "12M",
                sequence: "ATGCAATTTCCC",
            },
            SyntheticRead {
                name: "ref",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
        ],
    );

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.tsv_file = None;
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.threads = Some(1);
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("unphased".to_string());

    let summary = pipeline::run(&args).expect("unphased insertion pipeline should succeed");
    assert_eq!(summary.global.indel_variants, 1);

    let rows = read_tsv_rows(&tmp.join("unphased.MNV.tsv"));
    assert!(
        find_row(&rows, "complex_indel", "AAA", "CAAT").is_none(),
        "compound haplotype must require a read carrying both events"
    );
    let insertion = find_row(&rows, "insertion", "A", "AT").expect("original insertion row");
    assert_eq!(insertion["Event Reads"], "1");
    assert_eq!(insertion["Event Depth"], "3");

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_vcf_phased_deletion_haplotype_is_supported_by_same_read() {
    let tmp = temp_dir("e2e_phased_del_vcf");
    let (ref_path, genes_path) = write_phase_reference_files(&tmp);
    let vcf_path = tmp.join("phase_del.vcf");
    let bam_path = tmp.join("phase_del.bam");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=12>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tC\t.\tPASS\t.\n\
chr1\t5\t.\tAA\tA\t.\tPASS\t.\n",
    )
    .unwrap();
    write_synthetic_bam(
        &bam_path,
        12,
        &[
            SyntheticRead {
                name: "phased_del",
                start: 1,
                cigar: "5M1D6M",
                sequence: "ATGCATTTCCC",
            },
            SyntheticRead {
                name: "del_only",
                start: 1,
                cigar: "5M1D6M",
                sequence: "ATGAATTTCCC",
            },
            SyntheticRead {
                name: "snv_only",
                start: 1,
                cigar: "12M",
                sequence: "ATGCAATTTCCC",
            },
            SyntheticRead {
                name: "ref",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
        ],
    );

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.tsv_file = None;
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.threads = Some(1);
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("phase_del".to_string());

    let summary = pipeline::run(&args).expect("phased deletion pipeline should succeed");
    assert_eq!(summary.global.snp_records_in_vcf, 2);
    assert!(
        summary.global.indel_variants >= 2,
        "expected original deletion plus phased haplotype"
    );

    let rows = read_tsv_rows(&tmp.join("phase_del.MNV.tsv"));
    let compound =
        find_row(&rows, "complex_indel", "AAA", "CA").expect("compound deletion+SNV row");
    assert_eq!(compound["Positions"], "4");
    assert_eq!(compound["Event Reads"], "1");
    assert_eq!(compound["Event Depth"], "4");
    assert!(compound["Event Components"].contains("SNV:4:A>C"));
    assert!(compound["Event Components"].contains("DEL:6:A"));

    let deletion = find_row(&rows, "deletion", "AA", "A").expect("original deletion row");
    assert_eq!(deletion["Event Reads"], "2");
    assert_eq!(deletion["Event Depth"], "4");

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_ivar_phased_insertion_haplotype_matches_vcf_semantics() {
    let tmp = temp_dir("e2e_phased_ins_ivar");
    let (ref_path, genes_path) = write_phase_reference_files(&tmp);
    let ivar_path = tmp.join("phase.tsv");
    let bam_path = tmp.join("phase.bam");

    fs::write(
        &ivar_path,
        "REGION\tPOS\tREF\tALT\tREF_DP\tALT_DP\tALT_FREQ\tTOTAL_DP\tPASS\n\
chr1\t4\tA\tC\t1\t9\t0.9\t10\tTRUE\n\
chr1\t6\tA\t+T\t1\t9\t0.9\t10\tTRUE\n",
    )
    .unwrap();
    write_synthetic_bam(
        &bam_path,
        12,
        &[
            SyntheticRead {
                name: "phased_ins",
                start: 1,
                cigar: "6M1I6M",
                sequence: "ATGCAATTTTCCC",
            },
            SyntheticRead {
                name: "ins_only",
                start: 1,
                cigar: "6M1I6M",
                sequence: "ATGAAATTTTCCC",
            },
            SyntheticRead {
                name: "snv_only",
                start: 1,
                cigar: "12M",
                sequence: "ATGCAATTTCCC",
            },
            SyntheticRead {
                name: "ref",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
        ],
    );

    let mut args = base_args();
    args.vcf_file = None;
    args.tsv_file = Some(ivar_path.to_string_lossy().into());
    args.input_format = VariantInputFormat::Auto;
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.threads = Some(1);
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("phase_ivar".to_string());

    let summary = pipeline::run(&args).expect("phased iVar pipeline should succeed");
    assert_eq!(summary.global.snp_records_in_vcf, 2);
    assert!(
        summary.global.indel_variants >= 2,
        "expected original iVar insertion plus phased haplotype"
    );

    let rows = read_tsv_rows(&tmp.join("phase_ivar.MNV.tsv"));
    let compound = find_row(&rows, "complex_indel", "AAA", "CAAT")
        .expect("compound iVar insertion+SNV haplotype row");
    assert_eq!(compound["Event Reads"], "1");
    assert_eq!(compound["Event Depth"], "4");
    assert!(compound["Event Components"].contains("SNV:4:A>C"));
    assert!(compound["Event Components"].contains("INS:6:+T"));

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

#[test]
fn test_e2e_deletion_anchored_before_gene_not_intergenic_duplicate() {
    let tmp = temp_dir("e2e_boundary_del");
    let vcf_path = tmp.join("boundary.vcf");
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tCA\tC\t60\tPASS\tDP=12;AF=0.5\n",
    )
    .unwrap();
    fs::write(&ref_path, ">chr1\nCATGAAATTT\n").unwrap();
    fs::write(&genes_path, "gene1\t2\t10\t+\n").unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("boundary_del".to_string());

    let summary = pipeline::run(&args).expect("pipeline should annotate boundary deletion");
    assert_eq!(summary.global.indel_variants, 1);
    assert_eq!(summary.global.intergenic_variants, 0);

    let rows = read_tsv_rows(&tmp.join("boundary_del.MNV.tsv"));
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0]["Gene"], "gene1");
    assert_eq!(rows[0]["Event Class"], "deletion");
    assert_eq!(rows[0]["Event Components"], "DEL:2:A");
    assert_ne!(rows[0]["AA Changes"], "Unknown");

    fs::remove_dir_all(&tmp).ok();
}

#[test]
fn test_e2e_bam_without_base_qualities_still_counts_reads() {
    // Regression: a BAM with QUAL='*' (no per-base qualities) must not lose all
    // read support. Missing base quality is treated as passing (mirroring the
    // missing-MAPQ default of 255); otherwise every base scored 0, failed the
    // default --quality 20 gate, and EDP/EFREQ collapsed to 0 for every read.
    let tmp = temp_dir("e2e_qual_star");
    let (ref_path, genes_path) = write_phase_reference_files(&tmp); // ATGAAATTTCCC, gene1 1-12 +
    let vcf_path = tmp.join("q.vcf");
    let bam_path = tmp.join("q.bam");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=12>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tC\t.\tPASS\t.\n",
    )
    .unwrap();
    // Four reads span position 4; two carry the C alt. All have QUAL='*'.
    write_synthetic_bam_q(
        &bam_path,
        12,
        &[
            SyntheticRead {
                name: "alt1",
                start: 1,
                cigar: "12M",
                sequence: "ATGCAATTTCCC",
            },
            SyntheticRead {
                name: "alt2",
                start: 1,
                cigar: "12M",
                sequence: "ATGCAATTTCCC",
            },
            SyntheticRead {
                name: "ref1",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
            SyntheticRead {
                name: "ref2",
                start: 1,
                cigar: "12M",
                sequence: "ATGAAATTTCCC",
            },
        ],
        None, // QUAL='*'
    );

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.tsv_file = None;
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.gff_file = None;
    args.threads = Some(1);
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("q".to_string());

    pipeline::run(&args).expect("pipeline should succeed on a QUAL='*' BAM");

    let rows = read_tsv_rows(&tmp.join("q.MNV.tsv"));
    let snp = rows
        .iter()
        .find(|r| r["Base Changes"] == "C" && r["Positions"] == "4")
        .expect("SNP row at position 4");
    assert_eq!(
        snp["Total Reads"], "4",
        "all four spanning reads must count toward depth even without base qualities"
    );
    assert_eq!(
        snp["SNP Reads"], "2",
        "the two reads carrying the C allele must be counted"
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
fn test_e2e_bcf_input_errors_with_conversion_hint() {
    let mut args = base_args();
    args.vcf_file = Some("input.bcf".to_string());

    let err = pipeline::run(&args).expect_err("BCF input should fail");
    let message = err.to_string();
    assert!(message.contains("BCF input is not supported"));
    assert!(message.contains("bcftools view"));
}

#[test]
fn test_e2e_frequency_filters_require_bam() {
    let mut args = base_args();
    args.min_snp_frequency = 0.05;
    let err = pipeline::run(&args).expect_err("frequency filtering without BAM should fail");
    assert!(err.to_string().contains("require --bam"));
}

#[test]
fn test_e2e_frameshift_min_freq_out_of_range_errors() {
    // A common mistake is passing a percentage (50) instead of a fraction (0.5),
    // which previously silently disabled all frameshift propagation. It must now
    // be rejected like the other frequency flags.
    let mut args = base_args();
    args.frameshift_min_freq = 50.0;
    let err =
        pipeline::run(&args).expect_err("out-of-range --frameshift-min-freq should be rejected");
    assert!(err.to_string().contains("frameshift-min-freq"));
}

#[test]
fn test_e2e_phased_indel_min_freq_out_of_range_errors() {
    let mut args = base_args();
    args.phased_indel_min_freq = 1.5;
    let err =
        pipeline::run(&args).expect_err("out-of-range --phased-indel-min-freq should be rejected");
    assert!(err.to_string().contains("phased-indel-min-freq"));
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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
        frameshift_min_freq: 0.0,
        indel_anchor_depth: false,
        phased_indel_min_reads: 1,
        phased_indel_min_freq: 0.0,
        count_mates_separately: false,
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
        report: None,
        report_from: Vec::new(),
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

// The HTML report must be produced from a real run, embed its data, and stay
// self-contained (no external requests) so it opens offline.
#[test]
fn test_e2e_html_report_is_self_contained() {
    let tmp = temp_dir("e2e_report");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("rep".to_string());
    let report_path = tmp.join("report.html");
    args.report = Some(report_path.to_string_lossy().into());

    pipeline::run(&args).expect("pipeline with --report should succeed");

    let html = fs::read_to_string(&report_path).expect("report should exist");
    assert!(
        !html.contains("__GET_MNV_REPORT_DATA__"),
        "the data placeholder should have been replaced"
    );
    // Assert on the document title, which is stable, rather than on the exact
    // masthead markup, which carries styling spans.
    assert!(
        html.contains("<title>get_MNV report</title>"),
        "missing document title"
    );
    assert!(
        html.to_lowercase().contains("variant report"),
        "missing masthead heading"
    );
    assert!(html.contains("alt=\"get_MNV\""), "missing the project logo");
    assert!(
        html.contains("\"rows\":["),
        "the embedded payload should carry variant rows"
    );
    // Self-contained: nothing is fetched when the file is opened. Hyperlinks are
    // exempt on purpose, since an <a href> loads nothing until it is clicked;
    // what must not appear is anything the browser resolves on load.
    for pattern in ["src=\"http", "@import", "cdn.", "url(http", "<link "] {
        assert!(
            !html.contains(pattern),
            "report must not load external resources, found {pattern}"
        );
    }
    // Outbound links are allowed, but only as links, and they must not hand the
    // opener over to the target page.
    for (href, label) in [
        ("https://github.com/PathoGenOmics-Lab/get_MNV", "GitHub"),
        (
            "https://pathogenomics-lab.github.io/get_MNV/",
            "documentation",
        ),
    ] {
        assert!(
            html.contains(&format!("href=\"{href}\"")),
            "missing the {label} link"
        );
    }
    assert_eq!(
        html.matches("target=\"_blank\"").count(),
        html.matches("rel=\"noopener noreferrer\"").count(),
        "every new-tab link must carry rel=noopener noreferrer"
    );
    let _ = std::fs::remove_dir_all(&tmp);
}

// --report-from aggregates TSVs from separate runs into one multi-sample report
// without running the pipeline again.
#[test]
fn test_e2e_report_from_aggregates_existing_tsvs() {
    let tmp = temp_dir("e2e_report_from");
    let mut args = base_args();
    args.output_dir = Some(tmp.to_string_lossy().into());

    // Two per-sample TSVs, as a one-sample-per-run pipeline would leave behind.
    let mut tsvs = Vec::new();
    for sample in ["CASE-01", "CASE-02"] {
        args.output_prefix = Some(sample.to_string());
        let summary = pipeline::run(&args).expect("pipeline should succeed");
        tsvs.push(summary.output_tsv.expect("TSV output path"));
    }

    let report_path = tmp.join("cohort.html");
    let mut report_args = base_args();
    report_args.report_from = tsvs;
    report_args.report = Some(report_path.to_string_lossy().into());
    pipeline::run(&report_args).expect("--report-from should succeed");

    let html = fs::read_to_string(&report_path).expect("cohort report should exist");
    assert!(
        html.contains("CASE-01") && html.contains("CASE-02"),
        "both samples should be labelled in the report"
    );
    assert!(!html.contains("__GET_MNV_REPORT_DATA__"));
    let _ = std::fs::remove_dir_all(&tmp);
}

// ---------------------------------------------------------------------------
// Non-coding features declared in the TSV annotation are not translated
// ---------------------------------------------------------------------------

/// A gene table entry marked as a non-coding RNA must not be translated. Without
/// the biotype column get_MNV assumed every feature was protein-coding and
/// invented an amino-acid change for a gene that is never translated.
#[test]
fn test_e2e_tsv_biotype_keeps_non_coding_rna_untranslated() {
    let tmp = temp_dir("e2e_biotype");
    let ref_path = tmp.join("ref.fasta");
    // Two features over the same 12 bases so the only difference between the two
    // runs below is the declared biotype.
    fs::write(&ref_path, ">chr1\nATGAAATTTCCC\n").unwrap();
    let vcf_path = tmp.join("in.vcf");
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tG\t.\tPASS\t.\n",
    )
    .unwrap();

    let run = |genes: &str, stem: &str| -> Vec<HashMap<String, String>> {
        let genes_path = tmp.join(format!("{stem}.txt"));
        fs::write(&genes_path, genes).unwrap();
        let mut args = base_args();
        args.vcf_file = Some(vcf_path.to_string_lossy().into());
        args.fasta_file = ref_path.to_string_lossy().into();
        args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
        args.output_prefix = Some(tmp.join(stem).to_string_lossy().into());
        pipeline::run(&args).expect("pipeline should succeed");
        read_tsv_rows(&tmp.join(format!("{stem}.MNV.tsv")))
    };

    // Four columns: the historical format, still translated.
    let coding = run("ncR\t1\t12\t+\n", "legacy");
    assert_eq!(coding.len(), 1);
    assert_eq!(coding[0]["SO Term"], "missense_variant");
    assert_eq!(coding[0]["Impact"], "MODERATE");
    assert_ne!(coding[0]["AA Changes"], "-");

    // Same feature declared as a non-coding RNA: reported against its gene, with
    // no amino-acid change invented for it.
    let non_coding = run("ncR\t1\t12\t+\t0\tncRNA\n", "ncrna");
    assert_eq!(non_coding.len(), 1);
    assert_eq!(non_coding[0]["Gene"], "ncR");
    assert_eq!(
        non_coding[0]["SO Term"],
        "non_coding_transcript_exon_variant"
    );
    assert_eq!(non_coding[0]["Impact"], "MODIFIER");
    assert_eq!(non_coding[0]["AA Changes"], "-");

    let _ = fs::remove_dir_all(&tmp);
}

/// An unrecognised biotype must fail loudly: silently choosing either way would
/// invent a protein or hide a real one.
#[test]
fn test_e2e_tsv_unknown_biotype_is_rejected() {
    let tmp = temp_dir("e2e_biotype_bad");
    let ref_path = tmp.join("ref.fasta");
    fs::write(&ref_path, ">chr1\nATGAAATTTCCC\n").unwrap();
    let vcf_path = tmp.join("in.vcf");
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tG\t.\tPASS\t.\n",
    )
    .unwrap();
    let genes_path = tmp.join("genes.txt");
    fs::write(&genes_path, "g1\t1\t12\t+\t0\tprotein-coding\n").unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.output_prefix = Some(tmp.join("out").to_string_lossy().into());

    let err = pipeline::run(&args).expect_err("an unknown biotype must be rejected");
    let message = err.to_string();
    assert!(message.contains("unknown biotype"), "{message}");
    assert!(
        message.contains("protein_coding"),
        "the error should list what is accepted: {message}"
    );

    let _ = fs::remove_dir_all(&tmp);
}

/// A read that skips an intron (CIGAR `N`) must not be counted as carrying a
/// deletion of those bases. Both operations advance the reference without
/// consuming query bases, so reconstructing the observed allele from the CIGAR
/// alone makes an intron skip look exactly like the deletion allele.
#[test]
fn test_intron_skip_is_not_counted_as_deletion_support() {
    let tmp = temp_dir("e2e_skip_vs_del");
    let ref_path = tmp.join("ref.fasta");
    // 60 bases: a 3 bp deletion is called at 10 (REF=ATG ALT=A deletes 11-12).
    // Positions 1-9 = C, 10-12 = ATG, 13-60 = C.
    let sequence: &'static str = "CCCCCCCCCATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    fs::write(&ref_path, format!(">chr1\n{sequence}\n")).unwrap();
    let genes_path = tmp.join("genes.txt");
    fs::write(&genes_path, "gene1\t1\t60\t+\n").unwrap();
    let vcf_path = tmp.join("in.vcf");
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t10\t.\tATG\tA\t.\tPASS\t.\n",
    )
    .unwrap();

    // The read spans position 10 and then SKIPS 11-12 as an intron. It does not
    // carry the deletion.
    // 10 aligned bases (1-10), the intron 11-12 skipped, then 18 more (13-30).
    let spliced: &'static str = "CCCCCCCCCACCCCCCCCCCCCCCCCCC";
    let reads = vec![SyntheticRead {
        name: "skip1",
        start: 1,
        cigar: "10M2N18M",
        sequence: spliced,
    }];

    let bam_path = tmp.join("reads.bam");
    write_synthetic_bam(&bam_path, sequence.len(), &reads);

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.bam_file = Some(bam_path.to_string_lossy().into());
    args.output_prefix = Some(tmp.join("out").to_string_lossy().into());
    pipeline::run(&args).expect("pipeline should succeed");

    let rows = read_tsv_rows(&tmp.join("out.MNV.tsv"));
    let deletion = rows
        .iter()
        .find(|row| row.get("Event Class").map(String::as_str) == Some("deletion"))
        .expect("the deletion must be reported");
    assert_eq!(
        deletion.get("Event Reads").map(String::as_str),
        Some("0"),
        "an intron skip is not a deletion: no read carries this allele. Row: {deletion:?}"
    );
    // The user-visible number: this deletion used to be reported at frequency
    // 1.0000 on the strength of a read that does not carry it.
    assert_eq!(
        deletion.get("Event Frequency").map(String::as_str),
        Some("0.0000"),
        "frequency must not be inflated by a spliced read. Row: {deletion:?}"
    );

    let _ = fs::remove_dir_all(&tmp);
}

/// `--run-manifest` promises "inputs, outputs, checksums". With `--sample all`
/// it listed the per-sample output paths and no checksum beside them, because
/// that mode builds its own payload and skipped the block that computes them.
/// A manifest without them cannot serve the purpose it exists for.
#[test]
fn test_sample_all_manifest_records_a_checksum_for_every_sample() {
    use std::io::Write;

    let tmp = temp_dir("manifest_sample_all");
    let ex = example_dir();

    // A two-sample VCF built from the first few example records.
    let source = std::fs::read_to_string(ex.join("G35894.var.snp.vcf")).expect("example VCF");
    let multi = tmp.join("multi.vcf");
    let mut out = std::fs::File::create(&multi).expect("create multi-sample VCF");
    let mut written = 0;
    for line in source.lines() {
        if line.starts_with("##") {
            writeln!(out, "{line}").unwrap();
        } else if line.starts_with("#CHROM") {
            let fields: Vec<&str> = line.split('\t').take(8).collect();
            writeln!(out, "{}\tFORMAT\tSAMPA\tSAMPB", fields.join("\t")).unwrap();
        } else if written < 20 {
            let fields: Vec<&str> = line.split('\t').take(8).collect();
            writeln!(out, "{}\tGT\t1/1\t1/1", fields.join("\t")).unwrap();
            written += 1;
        }
    }
    drop(out);

    let manifest_path = tmp.join("run.manifest.json");
    let mut args = base_args();
    args.vcf_file = Some(multi.to_string_lossy().into());
    args.sample = Some("all".to_string());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.run_manifest = Some(manifest_path.to_string_lossy().into());

    pipeline::run(&args).expect("pipeline should succeed");

    let manifest: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&manifest_path).expect("manifest written"))
            .expect("manifest is JSON");
    let samples = manifest["samples"].as_array().expect("samples array");
    assert_eq!(samples.len(), 2, "one entry per sample");

    for sample in samples {
        let checksums = sample
            .get("output_checksums")
            .unwrap_or_else(|| panic!("no output_checksums for {sample:?}"));
        let recorded = checksums["output_tsv_sha256"]
            .as_str()
            .expect("a TSV checksum");
        let path = sample["output_tsv"].as_str().expect("a TSV path");
        let bytes = std::fs::read(path).expect("the TSV the manifest names");
        let actual = format!("{:x}", <sha2::Sha256 as sha2::Digest>::digest(&bytes));
        assert_eq!(recorded, actual, "checksum must match the file written");
    }

    let _ = std::fs::remove_dir_all(&tmp);
}

/// An intergenic substitution spanning two bases must reach the VCF whole.
///
/// The TSV row for one of these carries an entry per changed base. The VCF
/// writer read only the first entry, so the second base appeared in the TSV and
/// nowhere in the VCF: two outputs of a single run disagreeing about what had
/// been called. Reference: a row with two positions writes two records.
#[test]
fn test_e2e_intergenic_multibase_writes_every_base_to_vcf() {
    let tmp = temp_dir("e2e_intergenic_multibase_vcf");
    let vcf_path = tmp.join("interg.vcf");
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tAT\tCG\t60\tPASS\tDP=12;AF=0.5\n",
    )
    .unwrap();
    fs::write(&ref_path, ">chr1\nCATGAAATTTGGGCCC\n").unwrap();
    // El gen queda lejos: las dos bases son intergenicas.
    fs::write(&genes_path, "gene1\t10\t16\t+\n").unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("interg".to_string());
    args.both = true;

    let summary = pipeline::run(&args).expect("pipeline should annotate intergenic substitution");
    assert_eq!(summary.global.intergenic_variants, 1);

    let rows = read_tsv_rows(&tmp.join("interg.MNV.tsv"));
    assert_eq!(rows.len(), 1, "one TSV row for the intergenic substitution");
    let tsv_positions = rows[0]
        .get("Positions")
        .expect("Positions column")
        .split(", ")
        .map(str::to_string)
        .collect::<Vec<_>>();
    assert_eq!(tsv_positions, vec!["2".to_string(), "3".to_string()]);

    let vcf_body = fs::read_to_string(tmp.join("interg.MNV.vcf")).expect("VCF written");
    let records: Vec<&str> = vcf_body
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .collect();
    let vcf_positions: Vec<String> = records
        .iter()
        .map(|line| line.split('\t').nth(1).unwrap_or_default().to_string())
        .collect();
    assert_eq!(
        vcf_positions, tsv_positions,
        "every base of the TSV row must have a VCF record: {records:?}"
    );
    let alts: Vec<&str> = records
        .iter()
        .map(|line| line.split('\t').nth(4).unwrap_or_default())
        .collect();
    assert_eq!(alts, vec!["C", "G"], "each record carries its own base");
}

/// Each sample of a cohort is annotated for the alleles its genotype carries.
///
/// A VCF record lists every ALT seen at that site across the cohort. Annotating
/// all of them for every sample made each sample carry every variant, so the
/// cohort matrix of `--sample all` claimed a call for a sample whose genotype
/// reads 0/0, and two adjacent sites merged into an MNV for a sample that had
/// only one of them. A sample carrying nothing must also leave the rest of the
/// cohort alone instead of aborting the run.
#[test]
fn test_e2e_sample_all_respects_each_genotype() {
    let tmp = temp_dir("e2e_sample_all_genotypes");
    let vcf_path = tmp.join("cohort.vcf");
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");

    // S1 lleva solo la 11, S2 lleva la 10 y la 11 (un MNV), S3 no lleva nada.
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t10\t.\tA\tC\t60\tPASS\tDP=30\tGT\t0/0\t1/1\t0/0\n\
chr1\t11\t.\tA\tC\t60\tPASS\tDP=30\tGT\t1/1\t1/1\t0/0\n",
    )
    .unwrap();
    fs::write(&ref_path, ">chr1\nATGCCTAAAAAGGGTTTCCCATGCCTAAAG\n").unwrap();
    fs::write(&genes_path, "gene1\t1\t30\t+\n").unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("cohort".to_string());
    args.sample = Some("all".to_string());

    pipeline::run(&args).expect("a sample carrying nothing must not abort the cohort");

    let positions = |sample: &str| -> Vec<String> {
        let path = tmp.join(format!("cohort.sample_{sample}.MNV.tsv"));
        assert!(path.exists(), "{sample} should get an output file");
        read_tsv_rows(&path)
            .iter()
            .map(|row| row.get("Positions").cloned().unwrap_or_default())
            .collect()
    };

    assert_eq!(
        positions("S1"),
        vec!["11".to_string()],
        "S1 carries only 11"
    );
    assert_eq!(
        positions("S2"),
        vec!["10, 11".to_string()],
        "S2 carries both, so they merge into one MNV row"
    );
    assert!(
        positions("S3").is_empty(),
        "S3 carries nothing and gets an empty file, not a row"
    );

    fs::remove_dir_all(&tmp).ok();
}

/// Two samples whose output files would share a name stop the run.
///
/// The per-sample file is named after the sample with the characters a file name
/// cannot hold replaced, so `a/b` and `a_b` both become `sample_a_b`. Writing
/// them in turn left one file holding the second sample and nothing to say the
/// first had been overwritten, which is a cohort silently short one member.
#[test]
fn test_e2e_sample_all_rejects_colliding_output_names() {
    let tmp = temp_dir("e2e_sample_name_collision");
    let vcf_path = tmp.join("collide.vcf");
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");

    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta/b\ta_b\n\
chr1\t10\t.\tA\tC\t60\tPASS\tDP=30\tGT\t1/1\t1/1\n",
    )
    .unwrap();
    fs::write(&ref_path, ">chr1\nATGCCTAAAAAGGGTTTCCCATGCCTAAAG\n").unwrap();
    fs::write(&genes_path, "gene1\t1\t30\t+\n").unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.sample = Some("all".to_string());

    let error = pipeline::run(&args).expect_err("colliding sample names must stop the run");
    let message = error.to_string();
    assert!(
        message.contains("sample_a_b"),
        "the error should name the file both samples want: {message}"
    );

    let outputs = fs::read_dir(&tmp)
        .unwrap()
        .filter_map(|entry| entry.ok())
        .filter(|entry| entry.file_name().to_string_lossy().contains(".MNV."))
        .count();
    assert_eq!(outputs, 0, "nothing should be written before the run stops");

    fs::remove_dir_all(&tmp).ok();
}

/// A frameshift's premature stop is compared in the numbering of the rows it judges.
///
/// The stop is found by translating the alternate CDS, so it arrives as a residue
/// of the mutant protein, while the rows it is compared against carry reference
/// codon numbers. A 31 nt deletion puts the stop at mutant residue 9 and reference
/// codon 19, and every substitution between those two numbers was reported as
/// untranslated although the mutant protein still translates it. Translating the
/// deletion + substitution haplotype by hand shows residue 16 changing, so the row
/// must keep an amino-acid change; a row past the stop must keep the label.
#[test]
fn test_e2e_frameshift_stop_compared_in_reference_numbering() {
    let tmp = temp_dir("e2e_frameshift_ptc_numbering");
    let vcf_path = tmp.join("fs.vcf");
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");

    fs::write(&ref_path, ">chr1\nGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGATGAAACTGGTACACGCCAGGGGGGATGTAAGTGGAACTTTCGATATTGCTGGGCTAAGTTTGGGTTGTTCGCAAATCGAACTACAGACCTCAACGCTATTTGTCCTTCGCCCTTGCGGCGAAGCTCAGTGTGTGGCATCTTATAAAAAACGGTATTATTTTGAATTTCTTGCTCTTTGTGAATGTCCCCCGTGGACCGAGAGAAAGCTTTGCGGGAGGTGGGGTCAGCCATTCCAACGCTGTGGATACCTGTCTCATCCAGCGAACAATTTTAACAGAGTTCATTCTCCACACTAATGTGTTACCCAGTTCGAGCGCAT\n").unwrap();
    fs::write(&genes_path, "geneA\t101\t400\t+\n").unwrap();
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t112\t.\tACACGCCAGGGGGGATGTAAGTGGAACTTTCG\tA\t100\tPASS\tDP=30\n\
chr1\t146\t.\tA\tC\t100\tPASS\tDP=30\n\
chr1\t170\t.\tT\tA\t100\tPASS\tDP=30\n",
    )
    .unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("fs".to_string());

    pipeline::run(&args).expect("pipeline should annotate the frameshift");

    let rows = read_tsv_rows(&tmp.join("fs.MNV.tsv"));
    let by_position = |position: &str| -> (String, String) {
        let row = rows
            .iter()
            .find(|row| row.get("Positions").map(String::as_str) == Some(position))
            .unwrap_or_else(|| panic!("no row at {position}: {rows:?}"));
        (
            row.get("AA Changes").cloned().unwrap_or_default(),
            row.get("Impact").cloned().unwrap_or_default(),
        )
    };

    const UNTRANSLATED: &str = "downstream of premature stop";

    let (deletion_aa, _) = by_position("112");
    assert!(
        deletion_aa.ends_with("fs"),
        "the deletion itself is the frameshift: {deletion_aa}"
    );

    // Reference codon 16: after the stop in mutant numbering (9), before it in
    // reference numbering (19), and the mutant protein still translates it.
    let (inside_aa, inside_impact) = by_position("146");
    assert_ne!(inside_aa, UNTRANSLATED, "the mutant protein translates it");
    assert!(
        inside_aa.ends_with("(fs)"),
        "it is a frameshifted codon: {inside_aa}"
    );
    assert_eq!(inside_impact, "HIGH");

    // Reference codon 24: past the stop in both numberings.
    let (past_aa, past_impact) = by_position("170");
    assert_eq!(past_aa, UNTRANSLATED);
    assert_eq!(past_impact, "MODIFIER");

    fs::remove_dir_all(&tmp).ok();
}

/// An NMD verdict needs a stop that is actually premature.
///
/// The alternate protein's first stop was compared with the reference's by index,
/// in two proteins of different length. Deleting one codon moves the same natural
/// stop one index down, so every in-frame deletion was called premature and got an
/// NMD verdict; adding codons moves it up, so a stop an insertion genuinely
/// introduced got none. Both rows come from one run of the same two-exon
/// transcript.
#[test]
fn test_e2e_nmd_verdict_needs_a_premature_stop() {
    let tmp = temp_dir("e2e_nmd_premature_stop");
    let vcf_path = tmp.join("nmd.vcf");
    let ref_path = tmp.join("ref.fasta");
    let gff_path = tmp.join("genes.gff3");

    // Exon 1 = 1..90, intron 91..190, exon 2 = 191..280. Spliced CDS is
    // ATG + 58 Ala + stop, so the last exon-exon junction sits at CDS offset 90.
    let exon1 = format!("ATG{}", "GCA".repeat(29));
    let intron = format!("GT{}AG", "T".repeat(96));
    let exon2 = format!("{}TAA", "GCA".repeat(29));
    fs::write(&ref_path, format!(">chr1\n{exon1}{intron}{exon2}\n")).unwrap();
    fs::write(
        &gff_path,
        "##gff-version 3\n\
chr1\tsyn\tgene\t1\t280\t.\t+\t.\tID=gene-g1;Name=g1\n\
chr1\tsyn\tmRNA\t1\t280\t.\t+\t.\tID=mrna-g1;Parent=gene-g1;Name=g1\n\
chr1\tsyn\tCDS\t1\t90\t.\t+\t0\tID=c1;Parent=mrna-g1;Name=g1\n\
chr1\tsyn\tCDS\t191\t280\t.\t+\t0\tID=c2;Parent=mrna-g1;Name=g1\n",
    )
    .unwrap();

    // An in-frame insertion of 59 codons plus a stop, and a plain 3 nt deletion.
    let insertion = format!("A{}TAA", "GCA".repeat(59));
    fs::write(
        &vcf_path,
        format!(
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=280>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t21\t.\tA\t{insertion}\t100\tPASS\t.\n\
chr1\t30\t.\tAGCA\tA\t100\tPASS\t.\n"
        ),
    )
    .unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = None;
    args.gff_file = Some(gff_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("nmd".to_string());

    pipeline::run(&args).expect("pipeline should annotate both indels");

    let rows = read_tsv_rows(&tmp.join("nmd.MNV.tsv"));
    let field = |position: &str, column: &str| -> String {
        rows.iter()
            .find(|row| row.get("Positions").map(String::as_str) == Some(position))
            .unwrap_or_else(|| panic!("no row at {position}: {rows:?}"))
            .get(column)
            .cloned()
            .unwrap_or_default()
    };

    // The insertion carries its own stop, 59 codons before the natural one.
    assert_eq!(field("21", "SO Term"), "stop_gained");
    assert_eq!(
        field("21", "NMD Prediction"),
        "NMD-triggering",
        "a stop the insertion introduced is premature"
    );

    // The deletion removes one codon and gains no stop at all.
    assert_eq!(field("30", "SO Term"), "inframe_deletion");
    assert_eq!(
        field("30", "NMD Prediction"),
        "-",
        "the natural stop simply moved one index down"
    );

    fs::remove_dir_all(&tmp).ok();
}

/// The bases a declared phase skips still count towards an indel's length.
///
/// A phase of 2 means the first two bases of the CDS row belong to a codon that
/// began in a neighbouring exon: they code, they just cannot be translated from
/// this row alone. Measuring the indel's length change over the phase-trimmed
/// window counted a deletion of three coding bases as two, so an in-frame
/// deletion was reported as a frameshift and dragged "(fs)" onto every downstream
/// codon of the feature, where bcftools csq reports plain missense.
#[test]
fn test_e2e_phase_skipped_bases_count_towards_indel_length() {
    let tmp = temp_dir("e2e_phase_skipped_delta");
    let vcf_path = tmp.join("phase.vcf");
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");

    fs::write(
        &ref_path,
        ">chr1\nTTTTTTTTTAAGCTGCAGCTGCAGCTGCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
    )
    .unwrap();
    // CDS 10-29 on the plus strand with phase 2: genomic 10 and 11 are skipped.
    fs::write(&genes_path, "gA\t10\t29\t+\t2\n").unwrap();
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=60>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tAAGC\tA\t100\tPASS\tDP=30\n\
chr1\t18\t.\tG\tC\t100\tPASS\tDP=30\n",
    )
    .unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("phase".to_string());

    pipeline::run(&args).expect("pipeline should annotate the deletion and the substitution");

    let rows = read_tsv_rows(&tmp.join("phase.MNV.tsv"));
    let field = |position: &str, column: &str| -> String {
        rows.iter()
            .find(|row| row.get("Positions").map(String::as_str) == Some(position))
            .unwrap_or_else(|| panic!("no row at {position}: {rows:?}"))
            .get(column)
            .cloned()
            .unwrap_or_default()
    };

    // Three coding bases removed: the frame downstream is unchanged.
    assert_eq!(field("10", "SO Term"), "inframe_deletion");
    assert_eq!(field("10", "Change Type"), "In-frame Indel");
    // The affected codon began in a neighbouring exon, so the residues that
    // changed cannot be named from this row.
    assert_eq!(field("10", "AA Changes"), "Unknown");

    // The substitution downstream is a plain missense, not a frameshifted codon.
    assert_eq!(field("18", "SO Term"), "missense_variant");
    assert_eq!(field("18", "Change Type"), "Non-synonymous");
    let downstream_aa = field("18", "AA Changes");
    assert!(
        !downstream_aa.contains("fs"),
        "no frameshift to propagate: {downstream_aa}"
    );

    fs::remove_dir_all(&tmp).ok();
}

/// A splice call follows the bases a record changes, not the base it starts on.
///
/// The same change, `SNV:46:G>A` on a splice donor, can be written three ways: as
/// itself, padded with reference bases on its left, or as a deletion that removes
/// the whole GT donor. Keying the splice and intron tests on the record's POS
/// asked about a base the record leaves alone, so the padded encoding was called
/// intergenic and `--exclude-intergenic` deleted a HIGH-impact call, while the
/// deletion of both essential bases was reported as a low-impact splice region.
#[test]
fn test_e2e_splice_call_follows_the_changed_bases() {
    let tmp = temp_dir("e2e_splice_changed_bases");
    let ref_path = tmp.join("ref.fasta");
    let gff_path = tmp.join("genes.gff3");

    // Exon 1 = 1..45, intron 46..145 opening with the GT donor, exon 2 = 146..190.
    let exon1 = format!("ATG{}", "GCT".repeat(14));
    let intron = format!("GT{}AG", "T".repeat(96));
    let exon2 = format!("{}TAA", "GCT".repeat(14));
    fs::write(&ref_path, format!(">c1\n{exon1}{intron}{exon2}\n")).unwrap();
    fs::write(
        &gff_path,
        "##gff-version 3\n\
c1\tsyn\tgene\t1\t190\t.\t+\t.\tID=gene-g1;Name=g1\n\
c1\tsyn\tmRNA\t1\t190\t.\t+\t.\tID=m1;Parent=gene-g1;Name=g1\n\
c1\tsyn\tCDS\t1\t45\t.\t+\t0\tID=cds-g1a;Parent=m1;Name=g1\n\
c1\tsyn\tCDS\t146\t190\t.\t+\t0\tID=cds-g1b;Parent=m1;Name=g1\n",
    )
    .unwrap();

    let encodings = [
        ("plain", "c1\t46\t.\tG\tA\t100\tPASS\tDP=30\n"),
        ("padded", "c1\t40\t.\tGCTGCTG\tGCTGCTA\t100\tPASS\tDP=30\n"),
        ("donor_deleted", "c1\t45\t.\tTGT\tT\t100\tPASS\tDP=30\n"),
    ];

    for (name, record) in encodings {
        let vcf_path = tmp.join(format!("{name}.vcf"));
        fs::write(
            &vcf_path,
            format!(
                "##fileformat=VCFv4.2\n\
##contig=<ID=c1,length=190>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}"
            ),
        )
        .unwrap();

        let mut args = base_args();
        args.vcf_file = Some(vcf_path.to_string_lossy().into());
        args.fasta_file = ref_path.to_string_lossy().into();
        args.genes_file_tsv = None;
        args.gff_file = Some(gff_path.to_string_lossy().into());
        args.gff_features_raw = Some("CDS".to_string());
        args.output_dir = Some(tmp.to_string_lossy().into());
        args.output_prefix = Some(name.to_string());

        pipeline::run(&args).unwrap_or_else(|e| panic!("{name} should annotate: {e}"));

        let rows = read_tsv_rows(&tmp.join(format!("{name}.MNV.tsv")));
        assert_eq!(rows.len(), 1, "{name} should produce one row");
        assert_eq!(
            rows[0].get("SO Term").map(String::as_str),
            Some("splice_donor_variant"),
            "{name} destroys the donor"
        );
        assert_eq!(rows[0].get("Impact").map(String::as_str), Some("HIGH"));
        assert_eq!(
            rows[0].get("Gene").map(String::as_str),
            Some("g1"),
            "{name} belongs to the gene whose donor it destroys"
        );

        // The call must survive the flag that only removes variants outside
        // every feature.
        args.exclude_intergenic = true;
        args.output_prefix = Some(format!("{name}_excluded"));
        pipeline::run(&args).unwrap_or_else(|e| panic!("{name} with the flag: {e}"));
        let kept = read_tsv_rows(&tmp.join(format!("{name}_excluded.MNV.tsv")));
        assert_eq!(
            kept.len(),
            1,
            "--exclude-intergenic must not remove {name}: it is not intergenic"
        );
    }

    fs::remove_dir_all(&tmp).ok();
}

/// A record straddling a gene's edge is annotated once, and only where it is.
///
/// The fallback ran for any record no gene fully covered and rebuilt the whole
/// record, so `74 AAC>GGT` on a gene ending at 75 produced a second row repeating
/// 74 and 75 with a contradictory "Unknown" call, and claimed base 76 for a gene
/// that ends before it. Only the bases nobody annotated belong in that row.
#[test]
fn test_e2e_record_straddling_a_gene_edge_is_not_duplicated() {
    let tmp = temp_dir("e2e_straddling_record");
    let vcf_path = tmp.join("straddle.vcf");
    let ref_path = tmp.join("ref.fasta");
    let genes_path = tmp.join("genes.txt");

    let coding = format!("ATG{}{}TAA", "GCT".repeat(10), "AAA".repeat(10));
    fs::write(
        &ref_path,
        format!(">c1\n{}{coding}CCCGGGTTT\n", "T".repeat(9)),
    )
    .unwrap();
    fs::write(&genes_path, "geneA\t10\t75\t+\n").unwrap();
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=c1,length=84>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
c1\t74\t.\tAAC\tGGT\t100\tPASS\tDP=30\n",
    )
    .unwrap();

    let mut args = base_args();
    args.vcf_file = Some(vcf_path.to_string_lossy().into());
    args.fasta_file = ref_path.to_string_lossy().into();
    args.genes_file_tsv = Some(genes_path.to_string_lossy().into());
    args.output_dir = Some(tmp.to_string_lossy().into());
    args.output_prefix = Some("straddle".to_string());

    pipeline::run(&args).expect("pipeline should annotate the straddling record");

    let rows = read_tsv_rows(&tmp.join("straddle.MNV.tsv"));
    let described: Vec<(String, String)> = rows
        .iter()
        .map(|row| {
            (
                row.get("Positions").cloned().unwrap_or_default(),
                row.get("Gene").cloned().unwrap_or_default(),
            )
        })
        .collect();
    assert_eq!(
        described,
        vec![
            ("74, 75".to_string(), "geneA".to_string()),
            ("76".to_string(), "intergenic".to_string()),
        ],
        "the coding bases belong to the gene, the one past its end does not"
    );
    assert_eq!(
        rows[0].get("SO Term").map(String::as_str),
        Some("stop_lost")
    );

    // The gene's own row must survive the flag that removes only what lies
    // outside every feature.
    args.exclude_intergenic = true;
    args.output_prefix = Some("straddle_excluded".to_string());
    pipeline::run(&args).expect("pipeline should run with the flag");
    let kept = read_tsv_rows(&tmp.join("straddle_excluded.MNV.tsv"));
    assert_eq!(kept.len(), 1, "only the intergenic base is removed");
    assert_eq!(kept[0].get("Gene").map(String::as_str), Some("geneA"));

    fs::remove_dir_all(&tmp).ok();
}
