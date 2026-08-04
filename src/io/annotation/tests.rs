use super::*;
use crate::variants::{Biotype, CdsSegment};

// ---- decode_percent_encoded ----

#[test]
fn test_decode_percent_plain() {
    assert_eq!(decode_percent_encoded("hello"), "hello");
}

#[test]
fn test_decode_percent_space() {
    assert_eq!(decode_percent_encoded("hello%20world"), "hello world");
}

#[test]
fn test_decode_percent_multiple() {
    assert_eq!(decode_percent_encoded("%2C%3B"), ",;");
}

#[test]
fn test_decode_percent_trailing_incomplete() {
    // Incomplete percent sequence at end preserved literally
    assert_eq!(decode_percent_encoded("a%2"), "a%2");
    assert_eq!(decode_percent_encoded("a%"), "a%");
}

#[test]
fn test_decode_percent_invalid_hex() {
    // Invalid hex digits preserved literally
    assert_eq!(decode_percent_encoded("%ZZ"), "%ZZ");
}

// ---- split_attribute_fields ----

#[test]
fn test_split_simple() {
    let fields = split_attribute_fields("ID=gene1;Name=rpoB");
    assert_eq!(fields, vec!["ID=gene1", "Name=rpoB"]);
}

#[test]
fn test_split_quoted_semicolon() {
    // Semicolons inside quotes should not split
    let fields = split_attribute_fields("Note=\"a;b\";ID=gene1");
    assert_eq!(fields, vec!["Note=\"a;b\"", "ID=gene1"]);
}

#[test]
fn test_split_trailing_semicolon() {
    let fields = split_attribute_fields("ID=gene1;");
    assert_eq!(fields, vec!["ID=gene1"]);
}

#[test]
fn test_split_empty() {
    let fields = split_attribute_fields("");
    assert!(fields.is_empty());
}

// ---- parse_gff_attributes ----

#[test]
fn test_parse_gff_key_value() {
    let attrs = parse_gff_attributes("ID=gene1;Name=rpoB;locus_tag=Rv0667");
    assert_eq!(attrs.get("ID").unwrap(), "gene1");
    assert_eq!(attrs.get("Name").unwrap(), "rpoB");
    assert_eq!(attrs.get("locus_tag").unwrap(), "Rv0667");
}

#[test]
fn test_parse_gff_space_separated() {
    // GTF-style: key "value"
    let attrs = parse_gff_attributes("gene_id \"ENSG001\"; gene_name \"TP53\"");
    assert_eq!(attrs.get("gene_id").unwrap(), "ENSG001");
    assert_eq!(attrs.get("gene_name").unwrap(), "TP53");
}

#[test]
fn test_parse_gff_percent_encoded() {
    let attrs = parse_gff_attributes("Name=rpoB%20gene");
    assert_eq!(attrs.get("Name").unwrap(), "rpoB gene");
}

// ---- has_snp_in_interval ----

#[test]
fn test_has_snp_in_interval_found() {
    let snps = vec![
        VcfPosition {
            position: 100,
            ref_allele: "A".into(),
            alt_allele: "T".into(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        VcfPosition {
            position: 200,
            ref_allele: "G".into(),
            alt_allele: "C".into(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];
    assert!(has_snp_in_interval(&snps, 50, 150));
    assert!(has_snp_in_interval(&snps, 100, 100)); // exact match
    assert!(has_snp_in_interval(&snps, 150, 250));
}

#[test]
fn test_has_snp_in_interval_not_found() {
    let snps = vec![VcfPosition {
        position: 100,
        ref_allele: "A".into(),
        alt_allele: "T".into(),
        original_dp: None,
        original_freq: None,
        original_info: None,
    }];
    assert!(!has_snp_in_interval(&snps, 200, 300));
    assert!(!has_snp_in_interval(&snps, 50, 99));
}

#[test]
fn test_has_snp_in_interval_empty() {
    assert!(!has_snp_in_interval(&[], 1, 1000));
}

// ---- gene_name_from_gff ----

#[test]
fn test_gene_name_from_gff_full() {
    let mut attrs = HashMap::new();
    attrs.insert("gene".to_string(), "rpoB".to_string());
    attrs.insert("locus_tag".to_string(), "Rv0667".to_string());
    assert_eq!(gene_name_from_gff(&attrs), "rpoB_Rv0667");
}

#[test]
fn test_gene_name_from_gff_no_gene_uses_name() {
    let mut attrs = HashMap::new();
    attrs.insert("Name".to_string(), "gene-katG".to_string());
    // "gene-" prefix should be stripped; no locus_tag → no duplicated suffix
    assert_eq!(gene_name_from_gff(&attrs), "katG");
}

#[test]
fn test_gene_name_from_gff_empty() {
    let attrs = HashMap::new();
    assert_eq!(gene_name_from_gff(&attrs), "unknown_gene");
}

#[test]
fn test_gene_name_from_gff_prefers_gene_name_over_id() {
    // Regression: eukaryotic GTF/GFF3 rows often carry gene_name (human
    // readable) plus an ID coming from annotation tools like agat. We
    // should prefer gene_name over ID so the TSV shows GNAQ instead of
    // agat-cds-37838.
    let mut attrs = HashMap::new();
    attrs.insert("ID".to_string(), "agat-cds-37838".to_string());
    attrs.insert("gene_id".to_string(), "ENSG00000156052.11".to_string());
    attrs.insert("gene_name".to_string(), "GNAQ".to_string());
    assert_eq!(gene_name_from_gff(&attrs), "GNAQ");
}

#[test]
fn test_gene_name_from_gff_uses_gene_id_when_no_name() {
    let mut attrs = HashMap::new();
    attrs.insert("ID".to_string(), "agat-cds-1".to_string());
    attrs.insert("gene_id".to_string(), "ENSG00000156052.11".to_string());
    assert_eq!(gene_name_from_gff(&attrs), "ENSG00000156052.11");
}

#[test]
fn test_gene_name_from_gff_locus_tag_equal_to_primary_not_duplicated() {
    // When locus_tag is the ONLY name available, primary and locus_tag
    // are equal so we must not emit `primary_primary`.
    let mut attrs = HashMap::new();
    attrs.insert("locus_tag".to_string(), "Rv0007".to_string());
    assert_eq!(gene_name_from_gff(&attrs), "Rv0007");
}

#[test]
fn test_gene_name_from_gff_empty_locus_tag_no_trailing_underscore() {
    // Regression: an empty `locus_tag=""` (some annotation tools emit it)
    // used to produce "BRCA1_" — a trailing underscore that broke
    // downstream filtering by gene name.
    let mut attrs = HashMap::new();
    attrs.insert("gene_name".to_string(), "BRCA1".to_string());
    attrs.insert("locus_tag".to_string(), String::new());
    assert_eq!(gene_name_from_gff(&attrs), "BRCA1");
}

// ---- filter_genes_with_snps ----

// ---- assign_cds_protein_offsets / multi-exon aggregation ----

fn record(
    contig: &str,
    start: usize,
    end: usize,
    strand: crate::variants::Strand,
    phase: u8,
    ft: &str,
    tid: Option<&str>,
) -> GffGeneRecord {
    GffGeneRecord {
        contig: contig.to_string(),
        gene: Gene {
            name: "x".into(),
            start,
            end,
            strand,
            phase,
            protein_offset: 0,
            transcript_id: tid.map(|s| s.to_string()),
            cds_segments: Vec::new(),
            biotype: Biotype::ProteinCoding,
        },
        feature_type: ft.to_string(),
        transcript_id: tid.map(|s| s.to_string()),
    }
}

#[test]
fn test_assign_protein_offsets_gnaq_minus_strand() {
    use crate::variants::Strand;
    // Seven CDS rows of ENST00000286548.9 (GNAQ, minus strand) copied from
    // the minimal example. On minus strand, the first exon of the
    // transcript is the one with the HIGHEST genomic coordinate.
    let mut recs = vec![
        // (genomic order, not transcript order)
        record(
            "chr9",
            77_721_323,
            77_721_513,
            Strand::Minus,
            2,
            "CDS",
            Some("ENST00000286548.9"),
        ), // exon 7
        record(
            "chr9",
            77_728_514,
            77_728_667,
            Strand::Minus,
            0,
            "CDS",
            Some("ENST00000286548.9"),
        ), // exon 6
        record(
            "chr9",
            77_794_463,
            77_794_592,
            Strand::Minus,
            1,
            "CDS",
            Some("ENST00000286548.9"),
        ), // exon 5
        record(
            "chr9",
            77_797_520,
            77_797_648,
            Strand::Minus,
            1,
            "CDS",
            Some("ENST00000286548.9"),
        ), // exon 4
        record(
            "chr9",
            77_815_616,
            77_815_770,
            Strand::Minus,
            0,
            "CDS",
            Some("ENST00000286548.9"),
        ), // exon 3
        record(
            "chr9",
            77_922_161,
            77_922_345,
            Strand::Minus,
            2,
            "CDS",
            Some("ENST00000286548.9"),
        ), // exon 2
        record(
            "chr9",
            78_031_100,
            78_031_235,
            Strand::Minus,
            0,
            "CDS",
            Some("ENST00000286548.9"),
        ), // exon 1
    ];
    assign_cds_protein_offsets(&mut recs);
    // Build a map (start -> offset) so we don't depend on ordering.
    let by_start: std::collections::HashMap<usize, usize> = recs
        .iter()
        .map(|r| (r.gene.start, r.gene.protein_offset))
        .collect();
    // Expected cumulative codon offset before each exon, derived from the
    // GFF spec formula `(sum_prior_lengths + phase_i) / 3`. Hand-checked
    // against the canonical Ensembl GNAQ-201 protein numbering: the
    // famous oncogenic mutation Q209L (chr9:77794572) must end up at
    // codon 209, which only happens with these offsets.
    //   exon 1: S=0,    phase=0  → offset = 0      (codons 1..45)
    //   exon 2: S=136,  phase=2  → offset = 46     (codons 47..107; 46 is split 1↔2)
    //   exon 3: S=321,  phase=0  → offset = 107    (codons 108..158)
    //   exon 4: S=476,  phase=1  → offset = 159    (codons 160..201; 159 is split 3↔4)
    //   exon 5: S=605,  phase=1  → offset = 202    (codons 203..245; 202 is split 4↔5)
    //   exon 6: S=735,  phase=0  → offset = 245    (codons 246..296)
    //   exon 7: S=889,  phase=2  → offset = 297    (codons 298..360; 297 is split 6↔7)
    assert_eq!(by_start[&78_031_100], 0, "exon 1");
    assert_eq!(by_start[&77_922_161], 46, "exon 2");
    assert_eq!(by_start[&77_815_616], 107, "exon 3");
    assert_eq!(by_start[&77_797_520], 159, "exon 4");
    assert_eq!(by_start[&77_794_463], 202, "exon 5");
    assert_eq!(by_start[&77_728_514], 245, "exon 6");
    assert_eq!(by_start[&77_721_323], 297, "exon 7");
}

#[test]
fn test_assign_protein_offsets_gnaq_q209l_hotspot() {
    use crate::variants::Strand;
    // Hard-coded sanity check: GNAQ Q209L canonical mutation lives at
    // chr9:77794572 (middle base of codon 209). Verify the chain
    //   exon-5 effective end → local codon index → +protein_offset
    // produces codon 209 exactly.
    let mut recs = vec![
        record(
            "chr9",
            78_031_100,
            78_031_235,
            Strand::Minus,
            0,
            "CDS",
            Some("ENST00000286548.9"),
        ),
        record(
            "chr9",
            77_922_161,
            77_922_345,
            Strand::Minus,
            2,
            "CDS",
            Some("ENST00000286548.9"),
        ),
        record(
            "chr9",
            77_815_616,
            77_815_770,
            Strand::Minus,
            0,
            "CDS",
            Some("ENST00000286548.9"),
        ),
        record(
            "chr9",
            77_797_520,
            77_797_648,
            Strand::Minus,
            1,
            "CDS",
            Some("ENST00000286548.9"),
        ),
        record(
            "chr9",
            77_794_463,
            77_794_592,
            Strand::Minus,
            1,
            "CDS",
            Some("ENST00000286548.9"),
        ),
        record(
            "chr9",
            77_728_514,
            77_728_667,
            Strand::Minus,
            0,
            "CDS",
            Some("ENST00000286548.9"),
        ),
        record(
            "chr9",
            77_721_323,
            77_721_513,
            Strand::Minus,
            2,
            "CDS",
            Some("ENST00000286548.9"),
        ),
    ];
    assign_cds_protein_offsets(&mut recs);
    let exon5 = recs.iter().find(|r| r.gene.start == 77_794_463).unwrap();
    // The exon-5 effective end (after subtracting phase=1) is 77794591.
    // For position 77794572: offset = 77794591 - 77794572 = 19,
    // local codon index = 19/3 = 6, local_aa_pos = 7,
    // total = protein_offset(202) + 7 = 209.
    let pos = 77_794_572usize;
    let eff_end = exon5.gene.end - exon5.gene.phase as usize;
    let local_aa = (eff_end - pos) / 3 + 1;
    let aa = exon5.gene.protein_offset + local_aa;
    assert_eq!(
        aa, 209,
        "GNAQ Q209L canonical position must map to codon 209"
    );
}

#[test]
fn test_assign_protein_offsets_plus_strand_two_exons() {
    use crate::variants::Strand;
    // Plus strand, two exons. exon1 length 12, phase 0 → 4 complete codons
    // contributed; exon2 starts at protein_offset=4.
    let mut recs = vec![
        record("chrX", 100, 111, Strand::Plus, 0, "CDS", Some("T1")),
        record("chrX", 200, 217, Strand::Plus, 0, "CDS", Some("T1")),
    ];
    assign_cds_protein_offsets(&mut recs);
    let by_start: std::collections::HashMap<usize, usize> = recs
        .iter()
        .map(|r| (r.gene.start, r.gene.protein_offset))
        .collect();
    assert_eq!(by_start[&100], 0);
    assert_eq!(by_start[&200], 4);
}

#[test]
fn test_count_multi_exon_cds_transcripts() {
    use crate::variants::Strand;
    let recs = vec![
        record("chr1", 10, 30, Strand::Plus, 0, "CDS", Some("T1")),
        record("chr1", 50, 80, Strand::Plus, 0, "CDS", Some("T1")),
        record("chr1", 100, 130, Strand::Plus, 0, "CDS", Some("T2")),
        record("chr1", 200, 230, Strand::Plus, 0, "exon", Some("T3")),
        record("chr2", 10, 30, Strand::Plus, 0, "CDS", None),
    ];

    assert_eq!(count_multi_exon_cds_transcripts(&recs), 1);
}

#[test]
fn test_assign_protein_offsets_non_cds_untouched() {
    use crate::variants::Strand;
    let mut recs = vec![
        record("chr1", 10, 30, Strand::Plus, 0, "gene", None),
        record("chr1", 10, 30, Strand::Plus, 0, "exon", Some("T1")),
    ];
    assign_cds_protein_offsets(&mut recs);
    for rec in &recs {
        assert_eq!(rec.gene.protein_offset, 0);
    }
}

#[test]
fn test_build_transcript_cds_records_collapses_exons() {
    use crate::variants::Strand;
    let mut recs = vec![
        record("chr1", 1, 4, Strand::Plus, 0, "CDS", Some("tx1")),
        record("chr1", 10, 14, Strand::Plus, 2, "CDS", Some("tx1")),
    ];
    assign_cds_protein_offsets(&mut recs);
    let collapsed = build_transcript_cds_records(recs);

    assert_eq!(collapsed.len(), 1);
    assert_eq!(collapsed[0].gene.start, 1);
    assert_eq!(collapsed[0].gene.end, 14);
    assert_eq!(collapsed[0].gene.transcript_id.as_deref(), Some("tx1"));
    assert_eq!(collapsed[0].gene.cds_segments.len(), 2);
    assert_eq!(collapsed[0].gene.cds_segments[0].start, 1);
    assert_eq!(collapsed[0].gene.cds_segments[1].start, 10);
}

// ---- filter_genes_with_snps ----

#[test]
fn test_filter_genes_with_snps() {
    let genes = vec![
        Gene {
            name: "gene1".into(),
            start: 100,
            end: 200,
            strand: crate::variants::Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
            biotype: Biotype::ProteinCoding,
        },
        Gene {
            name: "gene2".into(),
            start: 500,
            end: 600,
            strand: crate::variants::Strand::Minus,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
            biotype: Biotype::ProteinCoding,
        },
    ];
    let snps = vec![VcfPosition {
        position: 150,
        ref_allele: "A".into(),
        alt_allele: "T".into(),
        original_dp: None,
        original_freq: None,
        original_info: None,
    }];
    let filtered = filter_genes_with_snps(&genes, &snps);
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].name, "gene1");
}

#[test]
fn test_transcript_gene_filter_keeps_intronic_variants() {
    // Two-exon transcript with a large intron (genomic 101..149).
    let genes = vec![Gene {
        name: "tx_gene".into(),
        start: 1,
        end: 200,
        strand: crate::variants::Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx1".to_string()),
        cds_segments: vec![
            CdsSegment { start: 1, end: 100 },
            CdsSegment {
                start: 150,
                end: 200,
            },
        ],
        biotype: Biotype::ProteinCoding,
    }];
    let snp = |position: usize| VcfPosition {
        position,
        ref_allele: "A".to_string(),
        alt_allele: "G".to_string(),
        original_dp: None,
        original_freq: None,
        original_info: None,
    };

    // A deep intronic variant (>8 nt from either junction) keeps the gene so it
    // can be reported as an intron_variant against the transcript it sits in.
    // It still does not overlap the CDS, so it is not coding.
    let deep = vec![snp(125)];
    assert_eq!(filter_genes_with_snps(&genes, &deep).len(), 1);
    assert!(!gene_overlaps_variant(&genes[0], &deep[0]));

    // An intronic variant in the essential donor site (first base of the intron)
    // keeps the gene so it can be annotated as a splice variant, even though it
    // still does not overlap the CDS.
    let splice = vec![snp(101)];
    assert_eq!(filter_genes_with_snps(&genes, &splice).len(), 1);
    assert!(!gene_overlaps_variant(&genes[0], &splice[0]));

    // A variant outside the transcript entirely is neither, and the gene goes.
    let outside = vec![snp(500)];
    assert!(filter_genes_with_snps(&genes, &outside).is_empty());
}
