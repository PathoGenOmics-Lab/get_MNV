use super::{build_phased_indel_haplotype_variants, codon_bounds_for_position, process_codon};
use crate::utils::reverse_complement;
use crate::variants::{CdsSegment, ChangeType, CodonInfo, Gene, Snp, Strand, VariantType};

fn next_u64(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    *seed
}

fn random_base(seed: &mut u64) -> char {
    match next_u64(seed) % 4 {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        _ => 'T',
    }
}

#[test]
fn test_property_reverse_complement_involution_for_random_sequences() {
    let mut seed = 123456789u64;
    for _ in 0..200 {
        let mut seq = String::new();
        for _ in 0..30 {
            seq.push(random_base(&mut seed));
        }
        let twice = reverse_complement(&reverse_complement(&seq));
        assert_eq!(seq, twice);
    }
}

#[test]
fn test_property_codon_bounds_plus_strand_cover_position() {
    let gene = Gene {
        name: "gene_plus".to_string(),
        start: 100,
        end: 399,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    for pos in gene.start..=gene.end {
        let bounds = codon_bounds_for_position(&gene, pos).expect("expected codon bounds");
        assert_eq!(bounds.1 - bounds.0 + 1, 3);
        assert!(bounds.0 <= pos && pos <= bounds.1);
        assert!(bounds.0 >= gene.start && bounds.1 <= gene.end);
    }
}

#[test]
fn test_property_codon_bounds_minus_strand_cover_position() {
    let gene = Gene {
        name: "gene_minus".to_string(),
        start: 100,
        end: 399,
        strand: Strand::Minus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    for pos in gene.start..=gene.end {
        let bounds = codon_bounds_for_position(&gene, pos).expect("expected codon bounds");
        assert_eq!(bounds.1 - bounds.0 + 1, 3);
        assert!(bounds.0 <= pos && pos <= bounds.1);
        assert!(bounds.0 >= gene.start && bounds.1 <= gene.end);
    }
}

#[test]
fn test_change_type_from_label_roundtrip() {
    let labels = [
        "Synonymous",
        "Non-synonymous",
        "Stop gained",
        "Stop lost",
        "Unknown",
        "Indel overlap",
        "Synonymous (frameshift)",
        "Non-synonymous (frameshift)",
        "Stop gained (frameshift)",
        "Stop lost (frameshift)",
        "Unknown (frameshift)",
        "Frameshift Indel",
        "In-frame Indel",
    ];
    for label in labels {
        let ct = ChangeType::from_label(label);
        assert_eq!(ct.to_string(), label);
    }
}

#[test]
fn test_change_type_from_label_unknown_fallback() {
    assert_eq!(ChangeType::from_label("garbage"), ChangeType::Unknown);
    assert_eq!(ChangeType::from_label(""), ChangeType::Unknown);
}

#[test]
fn test_with_frameshift_idempotent() {
    let fs = ChangeType::FrameshiftSynonymous;
    assert_eq!(fs.with_frameshift(), ChangeType::FrameshiftSynonymous);
    let base = ChangeType::StopGained;
    assert_eq!(base.with_frameshift(), ChangeType::FrameshiftStopGained);
}

#[test]
fn test_variant_type_display() {
    assert_eq!(VariantType::Snp.to_string(), "SNP");
    assert_eq!(VariantType::Mnv.to_string(), "MNV");
    assert_eq!(VariantType::SnpMnv.to_string(), "SNP/MNV");
    assert_eq!(VariantType::Indel.to_string(), "INDEL");
}

#[test]
fn test_strand_from_str() {
    assert_eq!("+".parse::<Strand>(), Ok(Strand::Plus));
    assert_eq!("-".parse::<Strand>(), Ok(Strand::Minus));
    assert!("?".parse::<Strand>().is_err());
}

#[test]
fn test_vcf_mnv_record_is_decomposed_into_codon_haplotype() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 4,
            ref_allele: "AA".to_string(),
            alt_allele: "CC".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].positions, vec![4, 5]);
    assert_eq!(variants[0].event_class.as_deref(), Some("mnv"));
    assert_eq!(
        variants[0].event_components,
        vec!["SNV:4:A>C".to_string(), "SNV:5:A>C".to_string()]
    );
}

#[test]
fn test_transcript_model_groups_mnv_across_exon_junction() {
    let gene = Gene {
        name: "tx_cds".to_string(),
        start: 1,
        end: 14,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx1".to_string()),
        cds_segments: vec![
            CdsSegment { start: 1, end: 4 },
            CdsSegment { start: 10, end: 14 },
        ],
    };
    let reference = crate::io::Reference {
        sequence: "ATGACCCCCAATTT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[
            crate::io::VcfPosition {
                position: 4,
                ref_allele: "A".to_string(),
                alt_allele: "C".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            crate::io::VcfPosition {
                position: 10,
                ref_allele: "A".to_string(),
                alt_allele: "G".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].variant_type, VariantType::SnpMnv);
    assert_eq!(variants[0].positions, vec![4, 10]);
    assert_eq!(variants[0].ref_codon.as_deref(), Some("AAA"));
    assert_eq!(variants[0].mnv_codon.as_deref(), Some("CGA"));
    assert_eq!(variants[0].aa_changes, vec!["Lys2Arg".to_string()]);
    assert_eq!(
        variants[0].event_components,
        vec!["SNV:4:A>C".to_string(), "SNV:10:A>G".to_string()]
    );
}

#[test]
fn test_transcript_synonymous_snp_is_synonymous() {
    // Regression: the transcript/CDS path classified the change type from the
    // three-letter IUPAC string ("Leu2Leu"), so `determine_change_type` compared
    // 'L' vs 'u' and reported every synonymous codon as Non-synonymous. A
    // synonymous Lys (AAA -> AAG) must classify as Synonymous.
    let gene = Gene {
        name: "tx_cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx1".to_string()),
        cds_segments: vec![CdsSegment { start: 1, end: 9 }],
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTT", // Met-Lys-Phe
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 6, // 3rd base of codon 2: AAA -> AAG, both Lys
            ref_allele: "A".to_string(),
            alt_allele: "G".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].change_type, ChangeType::Synonymous);
    assert_eq!(variants[0].aa_changes, vec!["Lys2Lys".to_string()]);
}

#[test]
fn test_transcript_start_lost_snp_is_start_lost() {
    // Regression: the same three-letter classification bug also hid Start-lost
    // calls on the transcript path (Met1 -> Thr was reported Non-synonymous).
    let gene = Gene {
        name: "tx_cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx1".to_string()),
        cds_segments: vec![CdsSegment { start: 1, end: 9 }],
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 2, // ATG -> ACG: Met1 -> Thr (start lost)
            ref_allele: "T".to_string(),
            alt_allele: "C".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].change_type, ChangeType::StartLost);
    assert_eq!(variants[0].aa_changes, vec!["Met1Thr".to_string()]);
}

#[test]
fn test_transcript_model_restored_frame_does_not_mark_downstream_snp_frameshift() {
    let gene = Gene {
        name: "tx_cds".to_string(),
        start: 1,
        end: 25,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx1".to_string()),
        cds_segments: vec![
            CdsSegment { start: 1, end: 6 },
            CdsSegment { start: 20, end: 25 },
        ],
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAAGGGGGGGGGGGGGTTTCCC",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[
            crate::io::VcfPosition {
                position: 3,
                ref_allele: "G".to_string(),
                alt_allele: "GA".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            crate::io::VcfPosition {
                position: 20,
                ref_allele: "TT".to_string(),
                alt_allele: "T".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            crate::io::VcfPosition {
                position: 23,
                ref_allele: "C".to_string(),
                alt_allele: "A".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    let snp_row = variants
        .iter()
        .find(|variant| variant.positions == vec![23])
        .expect("downstream SNP row");
    assert!(!snp_row.aa_changes[0].contains("(fs)"));
    assert!(!matches!(
        snp_row.change_type,
        ChangeType::FrameshiftSynonymous
            | ChangeType::FrameshiftNonSynonymous
            | ChangeType::FrameshiftStopGained
            | ChangeType::FrameshiftStopLost
            | ChangeType::FrameshiftUnknown
    ));
}

#[test]
fn test_indel_reports_event_components_and_frameshift_protein_effect() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 6,
            ref_allele: "A".to_string(),
            alt_allele: "AT".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].variant_type, VariantType::Indel);
    assert_eq!(variants[0].event_class.as_deref(), Some("insertion"));
    assert_eq!(variants[0].event_components, vec!["INS:6:+T".to_string()]);
    assert!(variants[0].aa_changes[0].contains("fs"));
}

// H8: an in-frame insertion on a minus-strand gene must be recognised as
// in-frame (length change is strand-invariant) and must not be mislabelled
// as a frameshift. Genomic AAATTTCAT reverse-complements to the CDS
// ATGAAATTT (Met-Lys-Phe).
#[test]
fn test_minus_strand_inframe_insertion_is_inframe() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Minus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "AAATTTCAT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 4,
            ref_allele: "T".to_string(),
            alt_allele: "TGGG".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].variant_type, VariantType::Indel);
    assert_eq!(variants[0].event_class.as_deref(), Some("insertion"));
    assert_eq!(variants[0].event_components, vec!["INS:4:+GGG".to_string()]);
    assert_eq!(variants[0].change_type, ChangeType::InFrameIndel);
    assert!(!variants[0].aa_changes[0].contains("fs"));
}

#[test]
fn test_minus_strand_legacy_deletion_reports_consistent_codons() {
    // Regression (BUG 2): on the legacy (no transcript model) minus strand, the
    // indel codon display sliced the length-changed alt CDS with a
    // reference-anchored offset, reporting the wrong reference codon and
    // mnv_codon = None. CDS = revcomp("CCCAAATTTCAT") = "ATGAAATTTGGG"
    // (Met-Lys-Phe-Gly); deleting genomic 4-6 ("AAA") removes Phe (codon 3).
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 12,
        strand: Strand::Minus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "CCCAAATTTCAT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 3,
            ref_allele: "CAAA".to_string(),
            alt_allele: "C".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].variant_type, VariantType::Indel);
    // Reference codon is the deleted Phe (codon 3); the alt codon at the same CDS
    // offset is the residue that shifts into its place. Both must be present.
    assert_eq!(variants[0].ref_codon.as_deref(), Some("TTT"));
    assert_eq!(variants[0].mnv_codon.as_deref(), Some("GGG"));
}

#[test]
fn test_plus_strand_insertion_at_internal_exon_junction_is_coding() {
    // Regression: an insertion anchored at the last base of an internal exon
    // lands at the exon-exon junction inside the spliced CDS, so it must be
    // annotated as a coding indel, not dropped (and later reported intergenic).
    // exon1 = [1,3] "ATG", exon2 = [10,15] "AAATTT"; CDS = "ATGAAATTT"
    // (Met-Lys-Phe). Inserting GGG after pos 3 -> spliced CDS "ATGGGGAAATTT".
    let gene = Gene {
        name: "tx_cds".to_string(),
        start: 1,
        end: 15,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx1".to_string()),
        cds_segments: vec![
            CdsSegment { start: 1, end: 3 },
            CdsSegment { start: 10, end: 15 },
        ],
    };
    let reference = crate::io::Reference {
        sequence: "ATGCCCCCCAAATTT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 3,
            ref_allele: "G".to_string(),
            alt_allele: "GGGG".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(variants.len(), 1);
    assert_ne!(variants[0].gene, "intergenic");
    assert_eq!(variants[0].variant_type, VariantType::Indel);
    assert_eq!(variants[0].change_type, ChangeType::InFrameIndel);
    assert_eq!(variants[0].event_components, vec!["INS:3:+GGG".to_string()]);
}

#[test]
fn test_minus_strand_insertion_at_internal_exon_junction_is_coding() {
    // Same junction case on the minus strand. cds_segments are in transcript
    // order (descending genomic): exon0 = [10,15], exon1 = [1,3]. CDS =
    // RC("TTTCAT") + RC("AAA") = "ATGAAATTT". An insertion of GGG after genomic
    // pos 3 (exon1's genomic end = transcript-5' boundary, an internal junction)
    // becomes RC("GGG")="CCC" at the exon0|exon1 junction -> "ATGAAACCCTTT".
    let gene = Gene {
        name: "tx_cds".to_string(),
        start: 1,
        end: 15,
        strand: Strand::Minus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx1".to_string()),
        cds_segments: vec![
            CdsSegment { start: 10, end: 15 },
            CdsSegment { start: 1, end: 3 },
        ],
    };
    let reference = crate::io::Reference {
        sequence: "AAACCCCCCTTTCAT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 3,
            ref_allele: "A".to_string(),
            alt_allele: "AGGG".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(variants.len(), 1);
    assert_ne!(variants[0].gene, "intergenic");
    assert_eq!(variants[0].variant_type, VariantType::Indel);
    assert_eq!(variants[0].change_type, ChangeType::InFrameIndel);
    assert_eq!(variants[0].event_components, vec!["INS:3:+GGG".to_string()]);
}

// H5: an in-frame insertion that introduces a stop codon must be classified
// as Stop gained rather than the generic In-frame Indel. Inserting TAA after
// the start codon of ATGAAATTT yields ATG-TAA-AAA-TTT (Met-*-Lys-Phe).
#[test]
fn test_inframe_insertion_creating_stop_is_stop_gained() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 3,
            ref_allele: "G".to_string(),
            alt_allele: "GTAA".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].variant_type, VariantType::Indel);
    assert_eq!(variants[0].change_type, ChangeType::StopGained);
}

// H1: a low-frequency upstream frameshift indel must not relabel a
// high-frequency downstream SNV as frameshifted once the frequency gate is
// raised above the indel frequency.
#[test]
fn test_frameshift_frequency_gate_skips_low_freq_upstream_indel() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 12,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTTGGG",
    };
    // Upstream 1bp deletion (pos 4, low freq) + downstream SNV in codon 4
    // (pos 11, high freq).
    let variants_in = [
        crate::io::VcfPosition {
            position: 3,
            ref_allele: "GA".to_string(),
            alt_allele: "G".to_string(),
            original_dp: None,
            original_freq: Some(0.05),
            original_info: None,
        },
        crate::io::VcfPosition {
            position: 11,
            ref_allele: "G".to_string(),
            alt_allele: "A".to_string(),
            original_dp: None,
            original_freq: Some(0.95),
            original_info: None,
        },
    ];

    let find_snp = |variants: &[crate::variants::VariantInfo]| -> crate::variants::VariantInfo {
        variants
            .iter()
            .find(|v| v.variant_type == VariantType::Snp)
            .expect("downstream SNP row")
            .clone()
    };

    // Default config (no gate): the downstream SNV is marked frameshifted.
    let default_variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &variants_in,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    let default_snp = find_snp(&default_variants);
    assert!(
        default_snp.aa_changes[0].contains("fs"),
        "without gate the downstream SNP should be frameshifted, got {:?}",
        default_snp.aa_changes
    );

    // Gate above the indel frequency: the downstream SNV is NOT frameshifted.
    let gated_variants = crate::variants::get_mnv_variants_for_gene_with_config(
        &gene,
        &variants_in,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
        &crate::variants::IndelAnnotationConfig {
            frameshift_min_freq: 0.5,
        },
        &crate::variants::FrameshiftPhasing::default(),
    );
    let gated_snp = find_snp(&gated_variants);
    assert!(
        !gated_snp.aa_changes[0].contains("fs"),
        "with gate the downstream SNP should not be frameshifted, got {:?}",
        gated_snp.aa_changes
    );

    // Phasing evidence that the indel (pos 3) and the downstream SNV (pos 11) are
    // in trans suppresses propagation even with the default (0.0) frequency gate.
    let trans = crate::variants::FrameshiftPhasing::from_trans_pairs(
        std::collections::HashSet::from([(3usize, 11usize)]),
    );
    let phased_variants = crate::variants::get_mnv_variants_for_gene_with_config(
        &gene,
        &variants_in,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
        &crate::variants::IndelAnnotationConfig::default(),
        &trans,
    );
    let phased_snp = find_snp(&phased_variants);
    assert!(
        !phased_snp.aa_changes[0].contains("fs"),
        "an indel in trans with the downstream SNV should not frameshift it, got {:?}",
        phased_snp.aa_changes
    );
}

#[test]
fn test_deletion_anchored_before_cds_keeps_protein_effect() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 2,
        end: 10,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "CATGAAATTT",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 1,
            ref_allele: "CA".to_string(),
            alt_allele: "C".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].variant_type, VariantType::Indel);
    assert_eq!(variants[0].event_class.as_deref(), Some("deletion"));
    assert_eq!(variants[0].event_components, vec!["DEL:2:A".to_string()]);
    assert_eq!(variants[0].change_type, ChangeType::FrameshiftIndel);
    assert_ne!(variants[0].aa_changes, vec!["Unknown".to_string()]);
    assert!(variants[0].aa_changes[0].contains("fs"));
}

#[test]
fn test_insertion_after_cds_end_is_not_coding() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTTC",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 9,
            ref_allele: "T".to_string(),
            alt_allele: "TA".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert!(variants.is_empty());
}

#[test]
fn test_insertion_after_codon_end_does_not_mask_that_codon_snp() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 12,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTTCCC",
    };
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[
            crate::io::VcfPosition {
                position: 7,
                ref_allele: "T".to_string(),
                alt_allele: "C".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            crate::io::VcfPosition {
                position: 9,
                ref_allele: "T".to_string(),
                alt_allele: "TA".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    let snp_row = variants
        .iter()
        .find(|variant| variant.variant_type == VariantType::Snp)
        .expect("SNP row");

    assert_ne!(snp_row.change_type, ChangeType::IndelOverlap);
    assert_ne!(snp_row.aa_changes, vec!["Unknown".to_string()]);
}

#[test]
fn test_build_phased_indel_haplotype_combines_nearby_snv() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTT",
    };
    let variants = vec![
        crate::io::VcfPosition {
            position: 4,
            ref_allele: "A".to_string(),
            alt_allele: "C".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        crate::io::VcfPosition {
            position: 6,
            ref_allele: "A".to_string(),
            alt_allele: "AT".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];

    let phased = build_phased_indel_haplotype_variants(
        &gene,
        &variants,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(phased.len(), 1);
    assert_eq!(phased[0].variant_type, VariantType::Indel);
    assert_eq!(phased[0].event_class.as_deref(), Some("complex_indel"));
    assert_eq!(phased[0].positions, vec![4]);
    assert_eq!(phased[0].ref_bases, vec!["AAA".to_string()]);
    assert_eq!(phased[0].base_changes, vec!["CAAT".to_string()]);
    assert_eq!(
        phased[0].event_components,
        vec!["SNV:4:A>C".to_string(), "INS:6:+T".to_string()]
    );
}

#[test]
fn test_build_phased_indel_haplotype_preserves_deletion_component_coordinate() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTT",
    };
    let variants = vec![
        crate::io::VcfPosition {
            position: 4,
            ref_allele: "A".to_string(),
            alt_allele: "C".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        crate::io::VcfPosition {
            position: 5,
            ref_allele: "AA".to_string(),
            alt_allele: "A".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];

    let phased = build_phased_indel_haplotype_variants(
        &gene,
        &variants,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(phased.len(), 1);
    assert_eq!(phased[0].event_class.as_deref(), Some("complex_indel"));
    assert_eq!(phased[0].positions, vec![4]);
    assert_eq!(phased[0].ref_bases, vec!["AAA".to_string()]);
    assert_eq!(phased[0].base_changes, vec!["CA".to_string()]);
    assert_eq!(
        phased[0].event_components,
        vec!["SNV:4:A>C".to_string(), "DEL:6:A".to_string()]
    );
}

#[test]
fn test_build_phased_indel_haplotype_combines_two_indels() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTT",
    };
    let variants = vec![
        crate::io::VcfPosition {
            position: 4,
            ref_allele: "A".to_string(),
            alt_allele: "AG".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        crate::io::VcfPosition {
            position: 5,
            ref_allele: "AA".to_string(),
            alt_allele: "A".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];

    let phased = build_phased_indel_haplotype_variants(
        &gene,
        &variants,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert_eq!(phased.len(), 1);
    assert_eq!(phased[0].event_class.as_deref(), Some("complex_indel"));
    assert_eq!(phased[0].positions, vec![4]);
    assert_eq!(phased[0].ref_bases, vec!["AAA".to_string()]);
    assert_eq!(phased[0].base_changes, vec!["AGA".to_string()]);
    assert_eq!(
        phased[0].event_components,
        vec!["INS:4:+G".to_string(), "DEL:6:A".to_string()]
    );
}

#[test]
fn test_build_phased_indel_haplotype_ignores_distant_snv() {
    let gene = Gene {
        name: "cds".to_string(),
        start: 1,
        end: 12,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = crate::io::Reference {
        sequence: "ATGAAATTTCCC",
    };
    let variants = vec![
        crate::io::VcfPosition {
            position: 10,
            ref_allele: "C".to_string(),
            alt_allele: "G".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        crate::io::VcfPosition {
            position: 6,
            ref_allele: "A".to_string(),
            alt_allele: "AT".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];

    let phased = build_phased_indel_haplotype_variants(
        &gene,
        &variants,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    assert!(phased.is_empty());
}

#[test]
fn test_codon_bounds_phase_skipped_position_returns_none() {
    // Plus strand, phase=2: positions 100 and 101 are in the
    // phase-skipped region and must NOT be reported as belonging to a
    // codon (the codon they sit in spans into the previous exon).
    let gene = Gene {
        name: "cds_plus_phase2".to_string(),
        start: 100,
        end: 120,
        strand: Strand::Plus,
        phase: 2,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    assert!(codon_bounds_for_position(&gene, 100).is_none());
    assert!(codon_bounds_for_position(&gene, 101).is_none());
    // 102 is the first base of the first complete codon.
    assert_eq!(codon_bounds_for_position(&gene, 102), Some((102, 104)));

    // Minus strand symmetric: phase=2 means the LAST 2 bases of the
    // exon are skipped (they belong to a codon ending in the next exon
    // in transcript order). 120 and 119 → None, 118 → first complete
    // codon.
    let minus = Gene {
        name: "cds_minus_phase2".to_string(),
        start: 100,
        end: 120,
        strand: Strand::Minus,
        phase: 2,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    assert!(codon_bounds_for_position(&minus, 120).is_none());
    assert!(codon_bounds_for_position(&minus, 119).is_none());
    assert_eq!(codon_bounds_for_position(&minus, 118), Some((116, 118)));
}

#[test]
fn test_codon_bounds_plus_strand_with_phase_1() {
    // Phase=1 means skip 1 base from gene.start: first codon is [start+1, start+3]
    let gene = Gene {
        name: "cds_plus_phase1".to_string(),
        start: 100,
        end: 120,
        strand: Strand::Plus,
        phase: 1,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    // Position 100 is in the skipped region, no codon
    assert!(codon_bounds_for_position(&gene, 100).is_none());
    // 101,102,103 → first codon
    assert_eq!(codon_bounds_for_position(&gene, 101), Some((101, 103)));
    assert_eq!(codon_bounds_for_position(&gene, 103), Some((101, 103)));
    // 104 → next codon
    assert_eq!(codon_bounds_for_position(&gene, 104), Some((104, 106)));
}

#[test]
fn test_codon_bounds_minus_strand_phase_1_gnaq_regression() {
    // Regression for issue #12: GNAQ CDS chr9:77794463-77794592, minus strand, phase=1.
    // Without phase fix, 77794516 and 77794517 were grouped together and 77794518 left alone.
    // With phase=1 applied, 77794517 and 77794518 must share a codon, and 77794516 must be on its own.
    let gene = Gene {
        name: "GNAQ_cds".to_string(),
        start: 77_794_463,
        end: 77_794_592,
        strand: Strand::Minus,
        phase: 1,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let b516 = codon_bounds_for_position(&gene, 77_794_516).expect("bounds for 516");
    let b517 = codon_bounds_for_position(&gene, 77_794_517).expect("bounds for 517");
    let b518 = codon_bounds_for_position(&gene, 77_794_518).expect("bounds for 518");
    assert_eq!(b517, b518, "517 and 518 must share a codon under phase=1");
    assert_ne!(
        b516, b517,
        "516 must NOT share a codon with 517/518 under phase=1"
    );
}

#[test]
fn test_codon_bounds_outside_gene_returns_none() {
    let gene = Gene {
        name: "g".to_string(),
        start: 100,
        end: 199,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    assert!(codon_bounds_for_position(&gene, 50).is_none());
    assert!(codon_bounds_for_position(&gene, 300).is_none());
}

#[test]
fn test_get_mnv_variants_for_gene_mixed_snps_and_indels() {
    use super::get_mnv_variants_for_gene;
    use crate::io::{Reference, VcfPosition};

    let gene = Gene {
        name: "testGene".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let reference = Reference {
        sequence: "ATGATGATG",
    };

    let snps = vec![
        VcfPosition {
            position: 2,
            ref_allele: "T".to_string(),
            alt_allele: "C".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        VcfPosition {
            position: 5,
            ref_allele: "TG".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];

    let variants = get_mnv_variants_for_gene(
        &gene,
        &snps,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert!(
        variants.len() >= 2,
        "expected SNP + indel, got {}",
        variants.len()
    );
    let has_snp = variants.iter().any(|v| v.variant_type == VariantType::Snp);
    let has_indel = variants
        .iter()
        .any(|v| v.variant_type == VariantType::Indel);
    assert!(has_snp, "missing SNP variant");
    assert!(has_indel, "missing indel variant");
}

#[test]
fn test_build_intergenic_variant_snp() {
    use super::build_intergenic_variant;
    use crate::io::VcfPosition;

    let pos = VcfPosition {
        position: 42,
        ref_allele: "A".to_string(),
        alt_allele: "G".to_string(),
        original_dp: Some(30),
        original_freq: Some(0.8),
        original_info: None,
    };
    let v = build_intergenic_variant("chrX", &pos);
    assert_eq!(v.gene, "intergenic");
    assert_eq!(v.variant_type, VariantType::Snp);
    assert_eq!(v.positions, vec![42]);
    assert_eq!(v.original_dp, Some(vec![30]));
}

#[test]
fn test_build_intergenic_variant_indel() {
    use super::build_intergenic_variant;
    use crate::io::VcfPosition;

    let pos = VcfPosition {
        position: 10,
        ref_allele: "AT".to_string(),
        alt_allele: "A".to_string(),
        original_dp: None,
        original_freq: None,
        original_info: None,
    };
    let v = build_intergenic_variant("chr1", &pos);
    assert_eq!(v.variant_type, VariantType::Indel);
}

#[test]
fn test_process_codon_mnv_two_snps_in_codon() {
    let codon_info = CodonInfo {
        codon_list: vec![
            Snp {
                index: 100,
                position: 100,
                ref_base: "A".to_string(),
                base: "T".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
            Snp {
                index: 101,
                position: 101,
                ref_base: "T".to_string(),
                base: "C".to_string(),
                original_dp: None,
                original_freq: None,
                original_info: None,
            },
        ],
        original_codon: "ATG".to_string(),
        gene_name: "gene1".to_string(),
        gene_start: 100,
        gene_end: 399,
        codon_start: 100,
        codon_end: 102,
        protein_offset: 0,
    };
    let result = process_codon(
        codon_info,
        Strand::Plus,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(result.variant_type, VariantType::SnpMnv);
    assert_eq!(result.positions.len(), 2);
    assert!(result.mnv_codon.is_some());
}

#[test]
fn test_process_codon_emits_expected_variant_type_for_single_snp() {
    let codon_info = CodonInfo {
        codon_list: vec![Snp {
            index: 101,
            position: 101,
            ref_base: "T".to_string(),
            base: "C".to_string(),
            original_dp: Some(20),
            original_freq: Some(0.5),
            original_info: None,
        }],
        original_codon: "ATG".to_string(),
        gene_name: "gene1".to_string(),
        gene_start: 100,
        gene_end: 399,
        codon_start: 100,
        codon_end: 102,
        protein_offset: 0,
    };

    let result = process_codon(
        codon_info,
        Strand::Plus,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(result.variant_type, VariantType::Snp);
    assert!(matches!(
        result.change_type,
        ChangeType::Synonymous
            | ChangeType::NonSynonymous
            | ChangeType::StartLost
            | ChangeType::StopGained
            | ChangeType::StopLost
            | ChangeType::Unknown
    ));
}

#[test]
fn test_construct_codon_uppercase_alt() {
    // ALT alleles should be uppercased to match the codon lookup table.
    let codon_info = CodonInfo {
        codon_list: vec![Snp {
            index: 101,
            position: 101,
            ref_base: "A".to_string(),
            base: "t".to_string(), // lowercase ALT
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        original_codon: "ATG".to_string(),
        gene_name: "test".to_string(),
        gene_start: 100,
        gene_end: 120,
        codon_start: 100,
        codon_end: 102,
        protein_offset: 0,
    };
    let result = process_codon(
        codon_info,
        Strand::Plus,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    // Position 101 is the 2nd base (index 1). ATG → ATG with pos 101 T→t
    // The codon becomes "ATt" → uppercased to "ATT" → Ile (I), not X.
    assert_ne!(
        result.aa_changes[0], "X",
        "Lowercase ALT should not produce unknown amino acid"
    );
}

#[test]
fn test_multiallelic_position_emits_one_row_per_alt() {
    // After --split-multiallelic, two ALTs at the same position must
    // each produce an independent VariantInfo row (one per alt). True
    // duplicates (same position + same alt) remain deduplicated.
    use crate::io::{Reference, VcfPosition};
    use crate::variants::get_mnv_variants_for_gene;

    let gene = Gene {
        name: "geneA".to_string(),
        start: 100,
        end: 111,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    // Two VCF records at position 101 with different ALTs (T and G)
    let snps = vec![
        VcfPosition {
            position: 101,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        VcfPosition {
            position: 101,
            ref_allele: "A".to_string(),
            alt_allele: "G".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];
    let seq = "N".repeat(99) + "ATGATGATGATG";
    let reference = Reference { sequence: &seq };
    let variants = get_mnv_variants_for_gene(
        &gene,
        &snps,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(
        variants.len(),
        2,
        "Two multi-allelic alts at the same position must produce two rows"
    );
    let mut base_changes: Vec<String> =
        variants.iter().map(|v| v.base_changes[0].clone()).collect();
    base_changes.sort();
    assert_eq!(base_changes, vec!["G".to_string(), "T".to_string()]);
}

#[test]
fn test_true_duplicate_position_still_dedup() {
    // Two records with the same position AND same alt remain a single row.
    use crate::io::{Reference, VcfPosition};
    use crate::variants::get_mnv_variants_for_gene;

    let gene = Gene {
        name: "geneA".to_string(),
        start: 100,
        end: 111,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let snps = vec![
        VcfPosition {
            position: 101,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        VcfPosition {
            position: 101,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];
    let seq = "N".repeat(99) + "ATGATGATGATG";
    let reference = Reference { sequence: &seq };
    let variants = get_mnv_variants_for_gene(
        &gene,
        &snps,
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(variants.len(), 1, "True duplicate must be deduplicated");
    assert_eq!(variants[0].base_changes[0], "T");
}

// A nonsense SNV in the first exon of a two-exon transcript, far upstream of the
// last exon-exon junction, must be flagged NMD-triggering end-to-end (through the
// spliced-transcript codon path), not just by the unit-level classifier.
#[test]
fn test_nonsense_snv_far_upstream_of_last_junction_is_nmd_triggering() {
    // exon1 = 120 nt (offsets 0..120), intron = 20 nt, exon2 = 60 nt.
    // Spliced CDS = 180 nt, last exon-exon junction at CDS offset 120.
    let mut exon1 = String::from("ATG"); // Met, codon 0
    exon1.push_str(&"GCT".repeat(9)); // Ala x9 (offsets 3..30)
    exon1.push_str("CAA"); // Gln, codon 10 at CDS offset 30
    exon1.push_str(&"GCT".repeat(29)); // Ala x29 (offsets 33..120)
    let intron = "T".repeat(20);
    let mut exon2 = "GCT".repeat(19); // Ala x19
    exon2.push_str("TAA"); // natural stop at the protein C-terminus
    let seq = format!("{exon1}{intron}{exon2}");

    let gene = Gene {
        name: "tx".to_string(),
        start: 1,
        end: seq.len(),
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx".to_string()),
        cds_segments: vec![
            CdsSegment { start: 1, end: 120 },
            CdsSegment {
                start: 141,
                end: 200,
            },
        ],
    };
    let reference = crate::io::Reference { sequence: &seq };

    // Genomic 31 C>T turns codon CAA (CDS offset 30) into TAA: a premature stop
    // 90 nt (> 50) upstream of the last junction.
    let variants = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 31,
            ref_allele: "C".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );

    let stop = variants
        .iter()
        .find(|v| v.change_type == ChangeType::StopGained)
        .expect("nonsense SNV should be classified Stop gained");
    assert_eq!(
        stop.annotations.nmd,
        Some(crate::variants::NmdPrediction::Triggering),
        "a PTC 90 nt upstream of the last junction should trigger NMD"
    );
}

// The HGVS coding descriptor must number from the CDS start and use coding-strand
// bases, so a minus-strand variant reads opposite to its genomic base. Verifies
// the wiring end-to-end for both strands (unit tests cover the string formatting).
#[test]
fn test_hgvs_coding_descriptor_plus_and_minus_strand() {
    let gc = crate::genetic_code::GeneticCode::default();

    // Plus strand, single-feature gene 1..9. An MNV at genomic 1 and 2 (C>T)
    // maps directly to CDS positions 1 and 2.
    let plus = Gene {
        name: "gp".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let plus_ref = crate::io::Reference {
        sequence: "CCATTTGGG",
    };
    let mnv = [
        crate::io::VcfPosition {
            position: 1,
            ref_allele: "C".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
        crate::io::VcfPosition {
            position: 2,
            ref_allele: "C".to_string(),
            alt_allele: "T".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        },
    ];
    let plus_variants =
        crate::variants::get_mnv_variants_for_gene(&plus, &mnv, &plus_ref, "chr1", gc);
    assert_eq!(
        plus_variants[0].annotations.hgvs_c.as_deref(),
        Some("c.[1C>T;2C>T]")
    );

    // Minus strand, gene 1..9 whose CDS is the reverse complement (ATGTTTAAA).
    // Genomic position 9 is CDS position 1; a genomic T>A there is a coding A>T.
    let minus = Gene {
        name: "gm".to_string(),
        start: 1,
        end: 9,
        strand: Strand::Minus,
        phase: 0,
        protein_offset: 0,
        transcript_id: None,
        cds_segments: Vec::new(),
    };
    let minus_ref = crate::io::Reference {
        sequence: "TTTAAACAT",
    };
    let snv = [crate::io::VcfPosition {
        position: 9,
        ref_allele: "T".to_string(),
        alt_allele: "A".to_string(),
        original_dp: None,
        original_freq: None,
        original_info: None,
    }];
    let minus_variants =
        crate::variants::get_mnv_variants_for_gene(&minus, &snv, &minus_ref, "chr1", gc);
    assert_eq!(
        minus_variants[0].annotations.hgvs_c.as_deref(),
        Some("c.1A>T"),
        "minus-strand coding descriptor should use CDS numbering and coding-strand bases"
    );
}

// An exonic coding variant within 3 nt of an internal exon-exon junction must be
// flagged as a splice-region variant (folded into the SO term downstream). Uses
// the real transcript path via get_mnv_variants_for_gene.
#[test]
fn test_exonic_variant_near_junction_is_splice_region() {
    // Two-exon plus-strand CDS: exon1 = 120 nt (genomic 1..120), exon2 141..200.
    let mut exon1 = String::from("ATG");
    exon1.push_str(&"GCT".repeat(39)); // 3 + 117 = 120 nt
    let intron = "T".repeat(20);
    let exon2 = "GCT".repeat(20); // 60 nt
    let seq = format!("{exon1}{intron}{exon2}");
    let gene = Gene {
        name: "tx".to_string(),
        start: 1,
        end: seq.len(),
        strand: Strand::Plus,
        phase: 0,
        protein_offset: 0,
        transcript_id: Some("tx".to_string()),
        cds_segments: vec![
            CdsSegment { start: 1, end: 120 },
            CdsSegment {
                start: 141,
                end: 200,
            },
        ],
    };
    let reference = crate::io::Reference { sequence: &seq };
    // Genomic 120 is exon1's last base: exonic, at the donor boundary.
    let near = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 120,
            ref_allele: "T".to_string(),
            alt_allele: "A".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(
        near[0].annotations.splice,
        Some(crate::variants::SpliceConsequence::Region),
        "a coding variant at the exon boundary should be in the splice region"
    );

    // Genomic 60 is deep inside exon1: coding only, no splice region.
    let deep = crate::variants::get_mnv_variants_for_gene(
        &gene,
        &[crate::io::VcfPosition {
            position: 60,
            ref_allele: "T".to_string(),
            alt_allele: "A".to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
        }],
        &reference,
        "chr1",
        crate::genetic_code::GeneticCode::default(),
    );
    assert_eq!(deep[0].annotations.splice, None);
}
