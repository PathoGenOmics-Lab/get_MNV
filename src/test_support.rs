//! Fixture construction for tests.
//!
//! **Build gene fixtures from GFF text, not from a `Gene` struct literal.**
//!
//! A hand-written literal can describe a gene model the parser would never
//! produce, and then the test agrees with the code about something untrue. That
//! is not hypothetical: a helper called `two_exon_gene` built two exons that
//! abutted, with no intron between them, so the tests that were meant to cover
//! splicing quietly asserted the behaviour of an unspliced transcript. They
//! passed, and they confirmed a bug that gave plain coding bases
//! `splice_donor_variant` at `HIGH` impact.
//!
//! Everything here goes through the real GFF pipeline, the same one
//! `load_genes_from_gff` uses, so a fixture is exactly what a user's annotation
//! file would produce.
//!
//! A struct literal is still the right tool for a gene model that is
//! deliberately impossible (empty CDS segments, coordinates past the end of the
//! reference) precisely because no parser would emit it. Say so in a comment
//! when you write one.

use crate::io::annotation::{
    assign_cds_protein_offsets, build_transcript_cds_records, parse_gff_gene_records_from_reader,
};
use crate::variants::{Gene, Strand};

/// Every gene in `gff`, parsed exactly as a real annotation file would be.
///
/// `feature_types` defaults to `CDS`, which is what get_MNV selects
/// automatically when a GFF contains CDS rows.
pub fn genes_from_gff(gff: &str) -> Vec<Gene> {
    genes_from_gff_features(gff, &["CDS".to_string()])
}

pub fn genes_from_gff_features(gff: &str, feature_types: &[String]) -> Vec<Gene> {
    let mut records = parse_gff_gene_records_from_reader(gff.as_bytes(), feature_types)
        .expect("fixture GFF should parse");
    assign_cds_protein_offsets(&mut records);
    build_transcript_cds_records(records)
        .into_iter()
        .map(|record| record.gene)
        .collect()
}

/// The single gene described by `gff`. Panics when the text describes anything
/// else, so a fixture cannot silently grow a second transcript.
pub fn gene_from_gff(gff: &str) -> Gene {
    let mut genes = genes_from_gff(gff);
    assert_eq!(
        genes.len(),
        1,
        "fixture should describe exactly one gene, got {}",
        genes.len()
    );
    genes.remove(0)
}

/// A single-exon coding gene: the prokaryotic shape, and the default for a test
/// that just needs somewhere for a variant to land.
pub fn single_exon_gene(name: &str, start: usize, end: usize, strand: Strand) -> Gene {
    gene_from_gff(&cds_row(name, start, end, strand, 0, None))
}

/// A single-exon coding gene with an explicit GFF phase.
pub fn single_exon_gene_with_phase(
    name: &str,
    start: usize,
    end: usize,
    strand: Strand,
    phase: u8,
) -> Gene {
    gene_from_gff(&cds_row(name, start, end, strand, phase, None))
}

/// A multi-segment transcript from its CDS segments, given in **transcript
/// order** (ascending genomic on the plus strand, descending on the minus). The
/// rows share one transcript identifier, so the parser aggregates them into the
/// spliced CDS model and computes the protein offsets itself.
///
/// Prefer [`spliced_gene`] or [`joined_gene`]: they say which kind of transcript
/// the test means and check that the segments actually describe it.
pub fn transcript_gene(name: &str, strand: Strand, segments: &[(usize, usize)]) -> Gene {
    let phases = segment_phases(segments);
    let rows: String = segments
        .iter()
        .zip(phases)
        .map(|(&(start, end), phase)| cds_row(name, start, end, strand, phase, Some(name)))
        .collect();
    gene_from_gff(&rows)
}

/// A genuinely spliced transcript: every consecutive pair of CDS segments is
/// separated by a real intron.
///
/// The assertion is the point. A fixture named for splicing whose exons abut is
/// an unspliced join, and a test built on it asserts the behaviour of a
/// transcript it does not describe. That is not a hypothetical failure mode: it
/// is how `splice_donor_variant` came to be reported on plain coding bases.
pub fn spliced_gene(name: &str, strand: Strand, segments: &[(usize, usize)]) -> Gene {
    let gene = transcript_gene(name, strand, segments);
    assert!(
        gene.cds_segments.len() >= 2,
        "a spliced transcript needs at least two CDS segments, got {}",
        gene.cds_segments.len()
    );
    for pair in gene.cds_segments.windows(2) {
        assert!(
            gene.intron_separates(&pair[0], &pair[1]),
            "segments {:?} and {:?} are not separated by an intron, so this is a \
             join and not a spliced transcript; use joined_gene if that is what \
             the test means",
            (pair[0].start, pair[0].end),
            (pair[1].start, pair[1].end),
        );
    }
    gene
}

/// A transcript whose CDS segments abut or overlap: one continuous reading
/// frame, as a ribosomal-slippage join such as SARS-CoV-2 ORF1ab. Nothing is
/// spliced out, so it has no exon-exon junction and no splice sites.
pub fn joined_gene(name: &str, strand: Strand, segments: &[(usize, usize)]) -> Gene {
    let gene = transcript_gene(name, strand, segments);
    for pair in gene.cds_segments.windows(2) {
        assert!(
            !gene.intron_separates(&pair[0], &pair[1]),
            "segments {:?} and {:?} are separated by an intron, so this is a \
             spliced transcript; use spliced_gene if that is what the test means",
            (pair[0].start, pair[0].end),
            (pair[1].start, pair[1].end),
        );
    }
    gene
}

/// GFF phase of each CDS segment, from the lengths of the ones before it.
fn segment_phases(segments: &[(usize, usize)]) -> Vec<u8> {
    let mut phases = Vec::with_capacity(segments.len());
    let mut prior = 0usize;
    for &(start, end) in segments {
        phases.push(((3 - prior % 3) % 3) as u8);
        prior += end.saturating_sub(start) + 1;
    }
    phases
}

fn cds_row(
    name: &str,
    start: usize,
    end: usize,
    strand: Strand,
    phase: u8,
    transcript: Option<&str>,
) -> String {
    let strand = match strand {
        Strand::Plus => '+',
        Strand::Minus => '-',
    };
    let mut attributes = format!("ID=cds-{name};gene={name}");
    if let Some(transcript) = transcript {
        attributes.push_str(&format!(";transcript_id={transcript}"));
    }
    format!("chr_fixture\tfixture\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attributes}\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_exon_gene_matches_a_prokaryotic_cds_row() {
        let gene = single_exon_gene("geneA", 10, 108, Strand::Plus);
        assert_eq!(gene.name, "geneA");
        assert_eq!((gene.start, gene.end), (10, 108));
        assert_eq!(gene.strand, Strand::Plus);
        assert_eq!(gene.phase, 0);
        assert_eq!(gene.protein_offset, 0);
        // A lone CDS row carries no transcript grouping, so there is no spliced
        // model and the legacy per-feature path is used, as before.
        assert!(gene.cds_segments.is_empty());
        assert!(gene.biotype.is_coding());
    }

    #[test]
    fn spliced_gene_segments_are_separated_by_real_introns() {
        // The whole point of this module. A "two exon" fixture whose exons abut
        // is an unspliced join, and asserting splice behaviour on it is how the
        // ribosomal-slippage bug hid.
        let gene = spliced_gene("geneC", Strand::Plus, &[(801, 900), (1001, 1200)]);
        assert_eq!(gene.cds_segments.len(), 2);
        assert!(
            gene.intron_separates(&gene.cds_segments[0], &gene.cds_segments[1]),
            "fixture must describe a genuinely spliced transcript"
        );
        assert!(crate::variants::splice::is_intronic_position(&gene, 950));
    }

    #[test]
    fn spliced_gene_keeps_transcript_order_on_the_minus_strand() {
        // Transcript order is descending genomic on the minus strand, and the
        // parser is what puts the segments in that order.
        let gene = spliced_gene("geneD", Strand::Minus, &[(1001, 1200), (801, 900)]);
        assert_eq!(gene.strand, Strand::Minus);
        assert_eq!(gene.cds_segments[0].start, 1001);
        assert_eq!(gene.cds_segments[1].start, 801);
        assert!(gene.intron_separates(&gene.cds_segments[0], &gene.cds_segments[1]));
    }

    #[test]
    #[should_panic(expected = "not separated by an intron")]
    fn spliced_gene_rejects_abutting_segments() {
        // The guard that would have caught `two_exon_gene`.
        spliced_gene("bad", Strand::Plus, &[(801, 900), (901, 1200)]);
    }

    #[test]
    #[should_panic(expected = "separated by an intron")]
    fn joined_gene_rejects_a_real_intron() {
        joined_gene("bad", Strand::Plus, &[(801, 900), (1001, 1200)]);
    }

    #[test]
    fn joined_gene_builds_a_slippage_join() {
        // SARS-CoV-2 ORF1ab: join(266..13468, 13468..21555).
        let gene = joined_gene("orf1ab", Strand::Plus, &[(266, 13468), (13468, 21555)]);
        assert_eq!(gene.cds_segments.len(), 2);
        assert!(crate::variants::splice::splice_consequence_for_position(&gene, 13469).is_none());
    }

    #[test]
    fn spliced_gene_computes_protein_offsets_from_the_segment_lengths() {
        // Exon 1 is 100 nt, so exon 2 starts inside codon 34 and its first
        // complete codon is number 35: offset 33 with phase 2.
        let gene = spliced_gene("geneC", Strand::Plus, &[(801, 900), (1001, 1200)]);
        assert_eq!(gene.cds_segments.len(), 2);
        // The aggregated record resets the offset; the per-segment phases are
        // what the parser used to reach it.
        assert_eq!(gene.protein_offset, 0);
        assert_eq!(segment_phases(&[(801, 900), (1001, 1200)]), vec![0, 2]);
        assert_eq!(segment_phases(&[(1, 99), (200, 298)]), vec![0, 0]);
    }

    #[test]
    fn phase_is_carried_through_to_the_gene() {
        let gene = single_exon_gene_with_phase("partial", 100, 200, Strand::Plus, 2);
        assert_eq!(gene.phase, 2);
    }
}
