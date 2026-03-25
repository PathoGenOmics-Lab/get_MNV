use super::annotation::{detect_annotation_format, preload_gff_genes, AnnotationFormat};
use super::fasta::{load_references, validate_vcf_reference_alleles, Reference};
use super::validation::validate_iupac_sequence;
use super::vcf::{load_vcf_positions_by_contig, VcfPosition};
use std::collections::HashMap;
use std::fs;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_temp_path(prefix: &str, suffix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before UNIX_EPOCH")
        .as_nanos();
    std::env::temp_dir().join(format!("{}_{}_{}", prefix, nanos, suffix))
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
            position: 15,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        VcfPosition {
            position: 20,
            ref_allele: "C".to_string(),
            alt_allele: "G".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
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
        load_references(path.to_string_lossy().as_ref()).expect("should load references");
    assert_eq!(refs.get("chr1"), Some(&"ACTG".to_string()));
    assert_eq!(refs.get("chr2"), Some(&"TTAA".to_string()));
    let _ = fs::remove_file(path);
}

#[test]
fn test_validate_vcf_reference_alleles_detects_mismatch() {
    let reference = Reference {
        sequence: "ACTG".to_string(),
    };
    let snp_list = vec![VcfPosition {
        position: 2,
        ref_allele: "A".to_string(),
        alt_allele: "T".to_string(),
        original_dp: None,
        original_freq: None,
        original_info: None,
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
    assert_eq!(gene_name_from_gff(&attrs), "Rv0001_Rv0001");

    let mut attrs = HashMap::new();
    attrs.insert("Name".to_string(), "dnaA".to_string());
    attrs.insert("locus_tag".to_string(), "Rv0001".to_string());
    attrs.insert("ID".to_string(), "gene-dnaA".to_string());
    assert_eq!(gene_name_from_gff(&attrs), "dnaA_Rv0001");

    let attrs = HashMap::new();
    assert_eq!(gene_name_from_gff(&attrs), "unknown_gene_unknown_gene");
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
    )
    .expect_err("multiallelic should fail without split mode");
    assert!(err.to_string().contains("Multiallelic"));

    let parsed =
        load_vcf_positions_by_contig(path.to_string_lossy().as_ref(), None, true, false, false)
            .expect("split mode should parse");
    let positions = parsed.get("chr1").expect("missing chr1");
    assert_eq!(positions.len(), 2);
    assert_eq!(positions[0].position, 2);
    assert_eq!(positions[0].alt_allele, "C");
    assert_eq!(positions[1].alt_allele, "G");

    let _ = fs::remove_file(path);
}
