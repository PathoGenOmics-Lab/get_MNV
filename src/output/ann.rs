//! The `ANN` INFO field: each record's annotation in the format snpEff writes
//! (the VCF annotation format, version 1.0), which `bcftools +split-vep`,
//! SnpSift and most VCF tooling already read.
//!
//! One entry per record, about the record's own ALT. On the record of a single
//! base of a codon with two or more substitutions that is the base on its own,
//! as the record's `AA` already reads; the codon's combined change belongs to
//! the codon's record, which says `INFO_COMPOUND_ANNOTATION`, snpEff's mark for
//! an effect that comes from combining variants.
//!
//! An `ANN` the input already carried (`--keep-original-info`) annotated the
//! input's records, often with more than get_MNV says: other transcripts, a
//! lab's own categories. Its entries for the record's ALT follow get_MNV's.

use super::common::{so_consequence, so_consequence_for};
use crate::variants::{NonCodingReason, VariantInfo};

/// snpEff's own header line, character for character: tools read the field
/// names out of its description.
pub(crate) const ANN_HEADER: &str = "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' \">";

/// Percent-encode whatever would end a field, an entry or the INFO column, and
/// the whitespace the format does not allow.
fn field(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for c in value.chars() {
        match c {
            '%' => out.push_str("%25"),
            '|' => out.push_str("%7C"),
            ',' => out.push_str("%2C"),
            ';' => out.push_str("%3B"),
            '=' => out.push_str("%3D"),
            '&' => out.push_str("%26"),
            ' ' => out.push_str("%20"),
            '\t' => out.push_str("%09"),
            '\n' => out.push_str("%0A"),
            '\r' => out.push_str("%0D"),
            other => out.push(other),
        }
    }
    out
}

/// The protein change as HGVS.p, when `AA` holds one: `Ser450Leu`, `Met1?`,
/// `Gln5*`, `Lys2_Phe3insGly`, `Ser450fs`. A label (`Synonymous`, `Unknown`,
/// `downstream of premature stop`), a residue that could not be read (`X`) or a
/// change read in a shifted frame (`Gly92Glu (fs)`) is not one.
fn hgvs_p(aa: &str) -> Option<String> {
    if aa.contains(|c: char| c.is_whitespace() || c == ',') {
        return None;
    }
    let b = aa.as_bytes();
    let residue = b.len() > 3
        && b[0].is_ascii_uppercase()
        && b[1].is_ascii_lowercase()
        && b[2].is_ascii_lowercase()
        && b[3].is_ascii_digit();
    if residue {
        return Some(format!("p.{aa}"));
    }
    // A lost stop (`*23Arg`) is an extension in HGVS, to a new stop that a
    // substitution's row does not place: `p.Ter23Argext*?`, as snpEff writes it.
    let lost_stop = aa.strip_prefix('*').and_then(|rest| {
        let digits = rest.chars().take_while(char::is_ascii_digit).count();
        let residue = &rest[digits..];
        (digits > 0 && residue.len() == 3 && residue.chars().all(|c| c.is_ascii_alphabetic()))
            .then(|| format!("p.Ter{}{residue}ext*?", &rest[..digits]))
    });
    lost_stop
}

/// The `ANN` value of one record. `alt` is the record's ALT, `base` the index of
/// the row's base it carries (`None` for the record of the whole row) and `aa`
/// the protein change its `AA` shows.
pub(crate) fn ann_value(
    variant: &VariantInfo,
    alt: &str,
    base: Option<usize>,
    aa: Option<&str>,
) -> String {
    let one_of_several = base.filter(|_| variant.positions.len() > 1);
    let ((term, impact), hgvs_c) = match one_of_several {
        Some(index) => match variant.annotations.single_bases.get(index) {
            Some(effect) => (
                so_consequence_for(
                    variant,
                    effect.in_context(variant.change_type),
                    effect.splice,
                ),
                effect.hgvs_c.as_deref(),
            ),
            // Not a codon: every base of the row is in the same kind of place,
            // but the row's descriptor names all of them, not this one.
            None => (so_consequence(variant), None),
        },
        None => (
            so_consequence(variant),
            variant.annotations.hgvs_c.as_deref(),
        ),
    };

    let intergenic = variant.gene == "intergenic";
    let gene = if intergenic {
        String::new()
    } else {
        field(&variant.gene)
    };
    let (feature_type, biotype) = if intergenic {
        ("intergenic_region", "")
    } else if variant.annotations.non_coding == Some(NonCodingReason::NonCodingTranscript) {
        ("transcript", "Noncoding")
    } else {
        ("transcript", "protein_coding")
    };
    // A codon's record joins its bases, and a local haplotype joins the input
    // records it was read from: either effect comes from combining variants.
    let compound =
        (base.is_none() && variant.positions.len() > 1) || variant.annotations.combines_records;

    [
        field(alt),
        term,
        impact.to_string(),
        gene.clone(),
        gene.clone(),
        feature_type.to_string(),
        gene,
        biotype.to_string(),
        // Rank, then the cDNA, CDS and protein positions and the distance.
        String::new(),
        hgvs_c.map(field).unwrap_or_default(),
        aa.and_then(hgvs_p).map(|p| field(&p)).unwrap_or_default(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        if compound {
            "INFO_COMPOUND_ANNOTATION".to_string()
        } else {
            String::new()
        },
    ]
    .join("|")
}

/// The whole `ANN` value of one record: get_MNV's entry, then those of the
/// input's `ANN` (in `original_info`) that are about the same allele.
pub(crate) fn ann_field(
    variant: &VariantInfo,
    alt: &str,
    base: Option<usize>,
    aa: Option<&str>,
    original_info: Option<&str>,
) -> String {
    let own = ann_value(variant, alt, base, aa);
    let kept: Vec<&str> = original_info
        .into_iter()
        .flat_map(|info| info.split(';'))
        .filter_map(|field| field.strip_prefix("ANN="))
        .flat_map(|value| value.split(','))
        .filter(|entry| {
            entry
                .split('|')
                .next()
                .is_some_and(|allele| allele.eq_ignore_ascii_case(alt))
                && *entry != own
        })
        .collect();
    std::iter::once(own.as_str())
        .chain(kept)
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::{ChangeType, SingleBaseEffect, VariantAnnotations, VariantType};

    fn codon(change_type: ChangeType, singles: &[ChangeType]) -> VariantInfo {
        VariantInfo {
            chrom: "chr1".to_string(),
            gene: "rpoB".to_string(),
            positions: vec![10, 11],
            ref_bases: vec!["T".to_string(), "C".to_string()],
            base_changes: vec!["G".to_string(), "A".to_string()],
            aa_changes: vec!["Ser4Glu".to_string()],
            snp_aa_changes: vec!["Ser4Ala".to_string(), "Ser4Ser".to_string()],
            aa_changes_local: vec!["Ser4Glu".to_string()],
            snp_aa_changes_local: vec!["Ser4Ala".to_string(), "Ser4Ser".to_string()],
            variant_type: VariantType::SnpMnv,
            change_type,
            snp_reads: None,
            snp_forward_reads: None,
            snp_reverse_reads: None,
            mnv_reads: None,
            mnv_forward_reads: None,
            mnv_reverse_reads: None,
            mnv_total_reads: None,
            total_reads: None,
            total_forward_reads: None,
            total_reverse_reads: None,
            mnv_total_forward_reads: None,
            mnv_total_reverse_reads: None,
            mnv_phasing_reads: None,
            ref_codon: None,
            snp_codon: None,
            mnv_codon: None,
            original_dp: None,
            original_freq: None,
            original_info: None,
            event_class: None,
            event_components: Vec::new(),
            annotations: VariantAnnotations {
                hgvs_c: Some("c.[10T>G;11C>A]".to_string()),
                single_bases: singles
                    .iter()
                    .zip(["c.10T>G", "c.11C>A"])
                    .map(|(&change_type, c)| SingleBaseEffect {
                        change_type,
                        hgvs_c: Some(c.to_string()),
                        splice: None,
                        original_info: None,
                    })
                    .collect(),
                ..Default::default()
            },
        }
    }

    fn fields(value: &str) -> Vec<&str> {
        let parts: Vec<&str> = value.split('|').collect();
        assert_eq!(parts.len(), 16, "{value}");
        parts
    }

    #[test]
    fn each_base_of_a_codon_is_annotated_on_its_own() {
        let row = codon(
            ChangeType::NonSynonymous,
            &[ChangeType::NonSynonymous, ChangeType::Synonymous],
        );
        let first = ann_value(&row, "G", Some(0), Some("Ser4Ala"));
        assert_eq!(
            fields(&first),
            [
                "G",
                "missense_variant",
                "MODERATE",
                "rpoB",
                "rpoB",
                "transcript",
                "rpoB",
                "protein_coding",
                "",
                "c.10T>G",
                "p.Ser4Ala",
                "",
                "",
                "",
                "",
                ""
            ]
        );
        // The second base alone changes nothing, whatever the codon does.
        let second = ann_value(&row, "A", Some(1), Some("Ser4Ser"));
        assert_eq!(fields(&second)[1..3], ["synonymous_variant", "LOW"]);
        assert_eq!(fields(&second)[9..11], ["c.11C>A", "p.Ser4Ser"]);
    }

    #[test]
    fn the_codon_record_carries_the_combined_change_marked_compound() {
        let row = codon(
            ChangeType::NonSynonymous,
            &[ChangeType::NonSynonymous, ChangeType::Synonymous],
        );
        let value = ann_value(&row, "GA", None, Some("Ser4Glu"));
        let parts = fields(&value);
        assert_eq!(parts[0], "GA");
        assert_eq!(parts[1], "missense_variant");
        // The allele list's `;` cannot stand in an INFO value.
        assert_eq!(parts[9], "c.[10T>G%3B11C>A]");
        assert_eq!(parts[10], "p.Ser4Glu");
        assert_eq!(parts[15], "INFO_COMPOUND_ANNOTATION");
    }

    #[test]
    fn a_codon_read_in_a_shifted_frame_keeps_that_label_on_every_base() {
        let row = codon(
            ChangeType::FrameshiftNonSynonymous,
            &[ChangeType::NonSynonymous, ChangeType::Synonymous],
        );
        let second = ann_value(&row, "A", Some(1), Some("Ser4Ser (fs)"));
        assert_eq!(fields(&second)[1], "frameshift_variant");
        assert_eq!(fields(&second)[10], "");
        let overlapped = codon(
            ChangeType::IndelOverlap,
            &[ChangeType::NonSynonymous, ChangeType::Synonymous],
        );
        let first = ann_value(&overlapped, "G", Some(0), Some("Unknown"));
        assert_eq!(
            fields(&first)[1..3],
            ["coding_sequence_variant", "MODIFIER"]
        );
    }

    #[test]
    fn an_intergenic_base_names_no_gene_or_transcript() {
        let mut row = codon(ChangeType::Unknown, &[]);
        row.gene = "intergenic".to_string();
        row.annotations.hgvs_c = None;
        let value = ann_value(&row, "G", Some(0), None);
        assert_eq!(
            fields(&value),
            [
                "G",
                "intergenic_variant",
                "MODIFIER",
                "",
                "",
                "intergenic_region",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                ""
            ]
        );
    }

    #[test]
    fn only_protein_changes_become_hgvs_p() {
        for (aa, expected) in [
            ("Ser450Leu", Some("p.Ser450Leu")),
            ("Met1?", Some("p.Met1?")),
            ("Gln5*", Some("p.Gln5*")),
            ("*100Gln", Some("p.Ter100Glnext*?")),
            ("*100*", None),
            ("Lys2_Phe3insGly", Some("p.Lys2_Phe3insGly")),
            ("Ser450fs", Some("p.Ser450fs")),
            ("Gly92Glu (fs)", None),
            ("Synonymous", None),
            ("Unknown", None),
            ("downstream of premature stop", None),
            ("X12Ala", None),
            ("-", None),
            ("Ser1Leu,Ala2Gly", None),
        ] {
            assert_eq!(hgvs_p(aa).as_deref(), expected, "{aa}");
        }
    }

    #[test]
    fn an_input_ann_follows_for_this_allele_only() {
        let row = codon(
            ChangeType::NonSynonymous,
            &[ChangeType::NonSynonymous, ChangeType::Synonymous],
        );
        let input = "DP=9;ANN=G|missense_variant|MODERATE|rpoB|g1|transcript|t1|Coding|1/1|c.10T>G|p.Ser4Ala||||||,G|custom|MODIFIER|||CUSTOM|essential||||||||||,T|missense_variant|MODERATE|rpoB|g1|transcript|t1|Coding|1/1|c.10T>A|p.Ser4Thr||||||;XX";
        let value = ann_field(&row, "G", Some(0), Some("Ser4Ala"), Some(input));
        let entries: Vec<&str> = value.split(',').collect();
        assert_eq!(entries.len(), 3, "{value}");
        assert!(entries[0].starts_with("G|missense_variant|MODERATE|rpoB|rpoB|"));
        assert!(entries[1].ends_with("|t1|Coding|1/1|c.10T>G|p.Ser4Ala||||||"));
        assert!(entries[2].starts_with("G|custom|"));
        // The codon's record has another allele: none of them describe it.
        let codon_record = ann_field(&row, "GA", None, Some("Ser4Glu"), Some(input));
        assert_eq!(codon_record.split(',').count(), 1, "{codon_record}");
        // A run on get_MNV's own output does not repeat its entry.
        let again = ann_field(
            &row,
            "G",
            Some(0),
            Some("Ser4Ala"),
            Some(&format!("ANN={}", entries[0])),
        );
        assert_eq!(again, entries[0]);
    }

    #[test]
    fn a_field_cannot_break_the_entry_or_the_info_column() {
        assert_eq!(field("a|b,c;d=e&f g%"), "a%7Cb%2Cc%3Bd%3De%26f%20g%25");
    }
}
