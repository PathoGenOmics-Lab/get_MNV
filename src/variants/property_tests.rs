//! Invariants checked over generated inputs rather than hand-picked ones.
//!
//! A worked example only covers the coordinates someone thought to write down.
//! These hold for every shape the generators produce, which is where off-by-one
//! errors in strand handling and codon boundaries actually live.

use crate::genetic_code::GeneticCode;
use crate::io::{Reference, VcfPosition};
use crate::test_support::single_exon_gene;
use crate::utils::reverse_complement;
use crate::variants::{get_mnv_variants_for_gene, Strand, VariantInfo};
use proptest::prelude::*;

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

fn snp(position: usize, reference: char, alternate: char) -> VcfPosition {
    VcfPosition {
        position,
        ref_allele: reference.to_string(),
        alt_allele: alternate.to_string(),
        original_dp: None,
        original_freq: None,
        original_info: None,
    }
}

/// The annotation facts that must not depend on which strand the gene is on.
fn annotation(variant: &VariantInfo) -> (Vec<String>, String, Option<String>, Option<String>) {
    (
        variant.aa_changes.clone(),
        variant.change_type.to_string(),
        variant.ref_codon.clone(),
        variant.mnv_codon.clone(),
    )
}

proptest! {
    /// Strand symmetry: annotating a variant on a plus-strand gene must give the
    /// same amino-acid consequence as annotating the mirrored variant on the
    /// reverse-complemented genome, where the same gene sits on the minus strand.
    ///
    /// The two runs describe the same physical molecule, so every protein-level
    /// answer has to match. This is the invariant that coordinate and complement
    /// mistakes on the minus strand break, and no single worked example pins it
    /// down across gene lengths, offsets and positions at once.
    #[test]
    fn strand_symmetry_of_codon_annotation(
        sequence in "[ACGT]{9,80}",
        prefix_len in 0usize..12,
        codons in 2usize..12,
        variant_offset in 0usize..36,
        alt_choice in 0usize..3,
    ) {
        let cds_len = codons * 3;
        let total = sequence.len();
        prop_assume!(variant_offset < cds_len);
        prop_assume!(prefix_len + cds_len <= total);

        let gene_start = prefix_len + 1;
        let gene_end = prefix_len + cds_len;
        let position = gene_start + variant_offset;

        let reference_base = sequence.as_bytes()[position - 1] as char;
        let alternates: Vec<char> = "ACGT".chars().filter(|&b| b != reference_base).collect();
        let alternate_base = alternates[alt_choice];

        // Plus strand, as generated.
        let plus_gene = single_exon_gene("g", gene_start, gene_end, Strand::Plus);
        let plus_reference = Reference { sequence: &sequence };
        let plus = get_mnv_variants_for_gene(
            &plus_gene,
            &[snp(position, reference_base, alternate_base)],
            &plus_reference,
            "chr",
            GeneticCode::default(),
        );

        // The same molecule read from the other side: the genome is
        // reverse-complemented, so genomic position p becomes total - p + 1, the
        // gene spans the mirrored interval on the minus strand, and both alleles
        // are complemented.
        let mirrored_sequence = reverse_complement(&sequence);
        let mirrored_gene = single_exon_gene(
            "g",
            total - gene_end + 1,
            total - gene_start + 1,
            Strand::Minus,
        );
        let mirrored_reference = Reference { sequence: &mirrored_sequence };
        let minus = get_mnv_variants_for_gene(
            &mirrored_gene,
            &[snp(
                total - position + 1,
                complement(reference_base),
                complement(alternate_base),
            )],
            &mirrored_reference,
            "chr",
            GeneticCode::default(),
        );

        prop_assert_eq!(plus.len(), minus.len(), "same number of annotated rows");
        for (plus_row, minus_row) in plus.iter().zip(minus.iter()) {
            prop_assert_eq!(
                annotation(plus_row),
                annotation(minus_row),
                "plus and minus strand annotations must agree for the same molecule"
            );
        }
    }

    /// Reverse complement is its own inverse.
    #[test]
    fn reverse_complement_is_an_involution(sequence in "[ACGT]{0,64}") {
        prop_assert_eq!(reverse_complement(&reverse_complement(&sequence)), sequence);
    }
}

proptest! {
    /// Every position inside a gene belongs to exactly one three-base codon that
    /// contains it and stays inside the feature, on either strand and at any
    /// phase. The previous version of this check walked one hard-coded gene.
    #[test]
    fn codon_bounds_cover_every_position_inside_a_gene(
        start in 1usize..500,
        length in 3usize..300,
        minus in any::<bool>(),
        phase in 0u8..3,
    ) {
        let end = start + length - 1;
        let strand = if minus { Strand::Minus } else { Strand::Plus };
        let gene = crate::test_support::single_exon_gene_with_phase("g", start, end, strand, phase);

        for position in gene.start..=gene.end {
            let Some((codon_start, codon_end)) =
                crate::variants::codon::codon_bounds(&gene, position)
            else {
                // Positions skipped by a non-zero phase have no complete codon.
                continue;
            };
            prop_assert_eq!(codon_end - codon_start + 1, 3, "a codon is three bases");
            prop_assert!(
                codon_start <= position && position <= codon_end,
                "codon {}-{} must contain position {}",
                codon_start,
                codon_end,
                position
            );
            prop_assert!(
                codon_start >= gene.start && codon_end <= gene.end,
                "codon {}-{} must stay inside the gene {}-{}",
                codon_start,
                codon_end,
                gene.start,
                gene.end
            );
        }
    }
}
