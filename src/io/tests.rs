use super::annotation::{detect_annotation_format, preload_gff_genes, AnnotationFormat};
use super::fasta::{load_references, validate_vcf_reference_alleles, Reference};
use super::validation::validate_iupac_sequence;
use super::vcf::{load_vcf_positions_by_contig, VcfPosition};
use std::collections::HashMap;
use std::fs;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_temp_path(prefix: &str, suffix: &str) -> std::path::PathBuf {
    // The clock alone does not make a name unique: macOS reads it in
    // microseconds, and two tests running in parallel took the same file, one
    // deleting it while the other still read it.
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before UNIX_EPOCH")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "{prefix}_{}_{nanos}_{n}_{suffix}",
        std::process::id()
    ))
}

#[test]
fn test_parse_annotation_format_by_extension() {
    let gff_format = detect_annotation_format("/tmp/example_annotations.gff3")
        .expect("should detect GFF by extension");
    let tsv_format =
        detect_annotation_format("/tmp/example_annotations.tsv").expect("should detect TSV");
    assert_eq!(gff_format, AnnotationFormat::Gff);
    assert_eq!(tsv_format, AnnotationFormat::Tsv);
}

#[test]
fn test_parse_annotation_format_by_content() {
    let path = unique_temp_path("get_mnv_annotations", "noext");
    fs::write(&path, "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=gene-abc\n")
        .expect("failed to write temp file");
    let format = detect_annotation_format(path.to_string_lossy().as_ref())
        .expect("should detect GFF by content");
    assert_eq!(format, AnnotationFormat::Gff);
    let _ = fs::remove_file(path);
}

#[test]
fn test_gene_name_from_gff_attributes() {
    use super::annotation::gene_name_from_gff;
    let mut attrs = HashMap::new();
    attrs.insert("gene".to_string(), "dnaN".to_string());
    attrs.insert("locus_tag".to_string(), "Rv0002".to_string());
    assert_eq!(gene_name_from_gff(&attrs), "dnaN_Rv0002");
}

#[test]
fn test_load_genes_from_gff_filters_contigs() {
    use super::annotation::load_genes_from_gff;
    let gff = "\
##gff-version 3
chrA\tsrc\tgene\t10\t30\t.\t+\t.\tID=gene-one;gene=one;locus_tag=L1
chrB\tsrc\tgene\t10\t30\t.\t+\t.\tID=gene-two;gene=two;locus_tag=L2
";
    let path = unique_temp_path("get_mnv_gff", "gff3");
    fs::write(&path, gff).expect("failed to write temp gff");

    let snp_list = vec![
        VcfPosition {
            record_start: 15,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
            declared_phase: None,
        },
        VcfPosition {
            record_start: 20,
            ref_allele: "C".to_string(),
            alt_allele: "G".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
            declared_phase: None,
        },
    ];

    let default_features = vec!["gene".to_string(), "pseudogene".to_string()];
    let genes = load_genes_from_gff(
        path.to_string_lossy().as_ref(),
        &snp_list,
        Some("chrA"),
        &default_features,
    )
    .expect("should load genes");
    assert_eq!(genes.len(), 1);
    assert_eq!(genes[0].name, "one_L1");
    assert_eq!(genes[0].start, 10);
    assert_eq!(genes[0].end, 30);
    let _ = fs::remove_file(path);
}

#[test]
fn test_parse_gff_attributes_handles_encoded_and_gtf_styles() {
    use super::annotation::parse_gff_attributes;
    let attrs = parse_gff_attributes(
        "ID=gene-abc;Name=dnaN%2Fbeta;locus_tag=Rv0002;gene_name \"dnaN\";note=\"ATPase;essential\"",
    );
    assert_eq!(attrs.get("ID"), Some(&"gene-abc".to_string()));
    assert_eq!(attrs.get("Name"), Some(&"dnaN/beta".to_string()));
    assert_eq!(attrs.get("locus_tag"), Some(&"Rv0002".to_string()));
    assert_eq!(attrs.get("gene_name"), Some(&"dnaN".to_string()));
    assert_eq!(attrs.get("note"), Some(&"ATPase;essential".to_string()));
}

#[test]
fn test_load_references_multiple_contigs() {
    let fasta_content = ">chr1\nACTG\n>chr2\nTTAA\n";
    let path = unique_temp_path("get_mnv_fasta_multi", "fas");
    fs::write(&path, fasta_content).expect("failed to write temp fasta");
    let refs =
        load_references(path.to_string_lossy().as_ref(), None).expect("should load references");
    assert_eq!(refs.get("chr1"), Some(&"ACTG".to_string()));
    assert_eq!(refs.get("chr2"), Some(&"TTAA".to_string()));
    let _ = fs::remove_file(path);
}

#[test]
fn test_validate_vcf_reference_alleles_detects_mismatch() {
    let reference = Reference { sequence: "ACTG" };
    let snp_list = vec![VcfPosition {
        record_start: 2,
        ref_allele: "A".to_string(),
        alt_allele: "T".to_string(),
        original_dp: None,
        original_freq: None,
        original_info: None,
        declared_phase: None,
    }];
    let error = validate_vcf_reference_alleles("chr1", &snp_list, &reference)
        .expect_err("expected mismatch");
    assert!(error.to_string().contains("VCF REF/FASTA mismatch"));
}

#[test]
fn test_validate_iupac_sequence_rejects_invalid_base() {
    let error = validate_iupac_sequence("ACXT", "test context")
        .expect_err("should reject invalid IUPAC base");
    assert!(error.to_string().contains("Invalid base"));
}

#[test]
fn test_parse_optional_depth_handles_integer_and_float_strings() {
    use super::vcf::parse_optional_depth;
    assert_eq!(parse_optional_depth("22"), Some(22));
    assert_eq!(parse_optional_depth("22.0"), Some(22));
    assert_eq!(parse_optional_depth("."), None);
    assert_eq!(parse_optional_depth(""), None);
}

// ---- decode_percent_encoded tests ----

#[test]
fn test_decode_percent_basic_ascii() {
    use super::annotation::decode_percent_encoded;
    assert_eq!(decode_percent_encoded("hello%20world"), "hello world");
    assert_eq!(decode_percent_encoded("a%2Cb%3Bc"), "a,b;c");
    assert_eq!(decode_percent_encoded("100%25"), "100%");
}

#[test]
fn test_decode_percent_passthrough() {
    use super::annotation::decode_percent_encoded;
    assert_eq!(decode_percent_encoded("plain text"), "plain text");
    assert_eq!(decode_percent_encoded(""), "");
}

#[test]
fn test_decode_percent_truncated_sequence() {
    use super::annotation::decode_percent_encoded;
    assert_eq!(decode_percent_encoded("abc%"), "abc%");
    assert_eq!(decode_percent_encoded("abc%2"), "abc%2");
}

#[test]
fn test_decode_percent_invalid_hex() {
    use super::annotation::decode_percent_encoded;
    assert_eq!(decode_percent_encoded("%ZZ"), "%ZZ");
    assert_eq!(decode_percent_encoded("%GG"), "%GG");
}

#[test]
fn test_decode_percent_multibyte_utf8() {
    use super::annotation::decode_percent_encoded;
    assert_eq!(decode_percent_encoded("%C3%A9"), "\u{00E9}");
}

// ---- split_attribute_fields tests ----

#[test]
fn test_split_attribute_fields_basic() {
    use super::annotation::split_attribute_fields;
    let fields = split_attribute_fields("ID=gene-abc;Name=foo;locus_tag=L1");
    assert_eq!(fields, vec!["ID=gene-abc", "Name=foo", "locus_tag=L1"]);
}

#[test]
fn test_split_attribute_fields_quoted_semicolon() {
    use super::annotation::split_attribute_fields;
    let fields = split_attribute_fields(r#"ID=gene-abc;note="ATPase;essential";gene=foo"#);
    assert_eq!(
        fields,
        vec!["ID=gene-abc", r#"note="ATPase;essential""#, "gene=foo"]
    );
}

#[test]
fn test_split_attribute_fields_empty_and_trailing() {
    use super::annotation::split_attribute_fields;
    assert!(split_attribute_fields("").is_empty());
    let fields = split_attribute_fields("ID=abc;Name=foo;");
    assert_eq!(fields, vec!["ID=abc", "Name=foo"]);
}

// ---- gene_name_from_gff fallback chain tests ----

#[test]
fn test_gene_name_fallback_chain() {
    use super::annotation::gene_name_from_gff;
    let mut attrs = HashMap::new();
    attrs.insert("ID".to_string(), "gene-Rv0001".to_string());
    // No locus_tag and no human name → just the stripped ID, not duplicated.
    assert_eq!(gene_name_from_gff(&attrs), "Rv0001");

    let mut attrs = HashMap::new();
    attrs.insert("Name".to_string(), "dnaA".to_string());
    attrs.insert("locus_tag".to_string(), "Rv0001".to_string());
    attrs.insert("ID".to_string(), "gene-dnaA".to_string());
    assert_eq!(gene_name_from_gff(&attrs), "dnaA_Rv0001");

    let attrs = HashMap::new();
    assert_eq!(gene_name_from_gff(&attrs), "unknown_gene");
}

// ---- preload_gff_genes test ----

#[test]
fn test_preload_gff_genes_groups_by_contig() {
    let gff = "\
##gff-version 3
chrA\tsrc\tgene\t10\t30\t.\t+\t.\tID=gene-one;gene=one;locus_tag=L1
chrA\tsrc\tgene\t40\t60\t.\t-\t.\tID=gene-two;gene=two;locus_tag=L2
chrB\tsrc\tgene\t100\t200\t.\t+\t.\tID=gene-three;gene=three;locus_tag=L3
chrB\tsrc\tCDS\t100\t200\t.\t+\t.\tID=cds-one
";
    let path = unique_temp_path("get_mnv_preload", "gff3");
    fs::write(&path, gff).expect("failed to write temp gff");

    let default_features = vec!["gene".to_string(), "pseudogene".to_string()];
    let result = preload_gff_genes(path.to_string_lossy().as_ref(), &default_features)
        .expect("should preload genes");
    assert_eq!(result.len(), 2);
    assert_eq!(result["chrA"].len(), 2);
    assert_eq!(result["chrB"].len(), 1);
    assert_eq!(result["chrB"][0].name, "three_L3");
    let _ = fs::remove_file(path);
}

// ---- gff feature type filtering tests ----

#[test]
fn test_parse_gff_gene_records_filters_by_feature_type() {
    use super::annotation::parse_gff_gene_records;
    let gff = "\
##gff-version 3
chrA\tsrc\tgene\t10\t30\t.\t+\t.\tID=gene-one;gene=one;locus_tag=L1
chrA\tsrc\tCDS\t10\t30\t.\t+\t0\tID=cds-one;gene=one;locus_tag=L1
chrA\tsrc\ttRNA\t50\t80\t.\t-\t.\tID=trna-one;gene=trn1;locus_tag=L2
chrA\tsrc\tpseudogene\t100\t200\t.\t+\t.\tID=gene-pg;gene=pg1;locus_tag=L3
";
    let path = unique_temp_path("get_mnv_ft_filter", "gff3");
    fs::write(&path, gff).expect("failed to write temp gff");
    let path_str = path.to_string_lossy();

    let default_ft = vec!["gene".to_string(), "pseudogene".to_string()];
    let records = parse_gff_gene_records(&path_str, &default_ft).unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].gene.name, "one_L1");
    assert_eq!(records[1].gene.name, "pg1_L3");

    let cds_ft = vec!["CDS".to_string()];
    let records = parse_gff_gene_records(&path_str, &cds_ft).unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].gene.start, 10);

    let multi_ft = vec!["CDS".to_string(), "tRNA".to_string()];
    let records = parse_gff_gene_records(&path_str, &multi_ft).unwrap();
    assert_eq!(records.len(), 2);

    let _ = fs::remove_file(path);
}

fn next_u64(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    *seed
}

fn random_char(seed: &mut u64) -> char {
    const ALPHABET: &[u8] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_=;:%\" ";
    let idx = (next_u64(seed) as usize) % ALPHABET.len();
    ALPHABET[idx] as char
}

#[test]
fn test_parse_gff_attributes_property_randomized_input() {
    use super::annotation::parse_gff_attributes;
    let mut seed = 0xC0FFEE_u64;
    for _ in 0..500 {
        let len = (next_u64(&mut seed) as usize % 120) + 1;
        let mut raw = String::with_capacity(len);
        for _ in 0..len {
            raw.push(random_char(&mut seed));
        }
        let parsed = parse_gff_attributes(&raw);
        for (key, value) in parsed {
            assert!(!key.trim().is_empty());
            assert!(!value.trim().is_empty());
        }
    }
}

#[test]
fn test_load_vcf_positions_multiallelic_split_mode() {
    let path = unique_temp_path("get_mnv_vcf_multi_split", "vcf");
    fs::write(
        &path,
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t2\t.\tT\tC,G\t.\tPASS\tDP=10;AF=0.5\n",
    )
    .expect("failed to write temp vcf");

    let err = load_vcf_positions_by_contig(
        path.to_string_lossy().as_ref(),
        None,
        false,
        false,
        false,
        None,
    )
    .expect_err("multiallelic should fail without split mode");
    assert!(err.to_string().contains("Multiallelic"));

    let parsed = load_vcf_positions_by_contig(
        path.to_string_lossy().as_ref(),
        None,
        true,
        false,
        false,
        None,
    )
    .expect("split mode should parse");
    let positions = parsed.get("chr1").expect("missing chr1");
    assert_eq!(positions.len(), 2);
    assert_eq!(positions[0].record_start, 2);
    assert_eq!(positions[0].alt_allele, "C");
    assert_eq!(positions[1].alt_allele, "G");

    let _ = fs::remove_file(path);
}

// ---- BCF input ----

const BCF_FIXTURE_VCF: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t2\t.\tT\tC,G\t50\tPASS\tDP=10;AF=0.5,0.3\tGT:DP\t1/2:10\n\
chr1\t5\t.\tA\tAT\t50\tPASS\tDP=12;AF=0.25\tGT:DP\t0/1:12\n\
chr1\t9\t.\tGC\tTA\t50\tPASS\tDP=8;AF=1\tGT:DP\t1/1:8\n";

/// Writes `vcf_text` as BCF through noodles, the way a caller's own BCF would
/// arrive, so the reader is tested against a real binary encoding.
fn write_bcf(vcf_text: &str, path: &std::path::Path) {
    use noodles::vcf::variant::io::Write as _;
    let mut reader = noodles::vcf::io::Reader::new(std::io::Cursor::new(vcf_text.as_bytes()));
    let header = reader.read_header().expect("fixture VCF header");
    let mut writer = noodles::bcf::io::Writer::new(fs::File::create(path).expect("create BCF"));
    writer.write_header(&header).expect("write BCF header");
    // Close the header's BGZF block, so the records start a block of their own
    // and a test can damage them without touching the header.
    std::io::Write::flush(writer.get_mut()).expect("flush BCF header block");
    for record in reader.records() {
        let record = record.expect("fixture VCF record");
        writer
            .write_variant_record(&header, &record)
            .expect("write BCF record");
    }
    writer.try_finish().expect("finish BCF");
}

fn describe(positions: &HashMap<String, Vec<VcfPosition>>) -> Vec<String> {
    let mut rows: Vec<String> = positions
        .iter()
        .flat_map(|(contig, values)| {
            values.iter().map(move |v| {
                format!(
                    "{contig}:{} {}>{} dp={:?} freq={:?}",
                    v.record_start, v.ref_allele, v.alt_allele, v.original_dp, v.original_freq
                )
            })
        })
        .collect();
    rows.sort();
    rows
}

#[test]
fn test_bcf_input_reads_the_same_as_the_vcf_it_encodes() {
    let base = unique_temp_path("get_mnv_bcf_same", "in");
    let vcf = base.with_extension("vcf");
    let bcf = base.with_extension("bcf");
    fs::write(&vcf, BCF_FIXTURE_VCF).unwrap();
    write_bcf(BCF_FIXTURE_VCF, &bcf);

    let from_vcf =
        load_vcf_positions_by_contig(vcf.to_str().unwrap(), None, true, false, false, None)
            .expect("VCF loads");
    let from_bcf =
        load_vcf_positions_by_contig(bcf.to_str().unwrap(), None, true, false, false, None)
            .expect("BCF loads");
    assert_eq!(describe(&from_vcf).len(), 4, "{:?}", describe(&from_vcf));
    assert_eq!(describe(&from_bcf), describe(&from_vcf));

    assert_eq!(
        super::vcf::list_vcf_samples(bcf.to_str().unwrap()).unwrap(),
        vec!["S1".to_string()]
    );
    let headers = super::vcf::extract_original_info_headers(bcf.to_str().unwrap()).unwrap();
    assert!(
        headers.iter().any(|line| line.contains("ID=AF")),
        "{headers:?}"
    );

    fs::remove_file(&vcf).ok();
    fs::remove_file(&bcf).ok();
}

#[test]
fn test_damaged_bcf_fails_instead_of_reading_short() {
    let bcf = unique_temp_path("get_mnv_bcf_damaged", "in").with_extension("bcf");
    write_bcf(BCF_FIXTURE_VCF, &bcf);
    let bytes = fs::read(&bcf).unwrap();
    // Walk the BGZF blocks (BSIZE, bytes 16-17, is the block length minus one)
    // and cut the records block in half, after a header that decodes cleanly.
    let mut blocks = Vec::new();
    let mut offset = 0usize;
    while offset + 18 <= bytes.len() {
        let length = u16::from_le_bytes([bytes[offset + 16], bytes[offset + 17]]) as usize + 1;
        blocks.push((offset, length));
        offset += length;
    }
    let (records_start, records_length) = blocks[1];
    assert!(
        records_length > 28,
        "second block should hold the records: {blocks:?}"
    );
    fs::write(&bcf, &bytes[..records_start + records_length / 2]).unwrap();

    let samples =
        super::vcf::list_vcf_samples(bcf.to_str().unwrap()).expect("the header alone still reads");
    assert_eq!(samples, vec!["S1".to_string()]);
    let result =
        load_vcf_positions_by_contig(bcf.to_str().unwrap(), None, true, false, false, None);
    assert!(
        result.is_err(),
        "a truncated BCF must not read as a shorter VCF"
    );

    fs::remove_file(&bcf).ok();
}
