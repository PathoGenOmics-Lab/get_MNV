//! Splice-site classification near internal exon-exon junctions.
//!
//! Uses the multi-exon transcript CDS model (`cds_segments`, in transcript
//! order): the two intronic bases flanking each internal junction are the
//! essential donor (5' of the intron) / acceptor (3' of the intron) sites; the
//! exon's first/last 3 bases and the intron's 3rd-8th bases are the wider splice
//! region. Single-exon transcripts (prokaryotes, most viruses) have no junction
//! and are never annotated.

use crate::variants::{Gene, SpliceConsequence, Strand};

/// Which half of a junction a boundary describes, so the two essential intronic
/// bases map to the correct donor / acceptor term.
#[derive(Clone, Copy)]
enum HalfKind {
    Donor,
    Acceptor,
}

/// Classify a genomic `position` against one junction boundary: `boundary` is
/// the exon's terminal base, `intron_dir` is +1 when the intron lies at higher
/// genomic coordinates than the exon and -1 when lower. Distance into the intron
/// of 1-2 is the essential site, 3-8 the intronic splice region, and 0 to -2 the
/// exonic splice region (the exon's last three bases).
fn classify_boundary(
    position: usize,
    boundary: usize,
    intron_dir: isize,
    kind: HalfKind,
) -> Option<SpliceConsequence> {
    let into_intron = (position as isize - boundary as isize) * intron_dir;
    match into_intron {
        1..=2 => Some(match kind {
            HalfKind::Donor => SpliceConsequence::Donor,
            HalfKind::Acceptor => SpliceConsequence::Acceptor,
        }),
        // Intronic splice region (3-8 nt in) and exonic splice region (the exon's
        // last three bases, distances 0 to -2).
        -2..=0 | 3..=8 => Some(SpliceConsequence::Region),
        _ => None,
    }
}

/// Splice consequence of a genomic `position` for `gene`, or `None` when it is
/// not within a splice region of any internal exon-exon junction (or the gene is
/// single-exon). The most severe match across junctions is returned.
pub fn splice_consequence_for_position(gene: &Gene, position: usize) -> Option<SpliceConsequence> {
    if gene.cds_segments.len() < 2 {
        return None;
    }
    let mut best: Option<SpliceConsequence> = None;
    // Only real introns have splice sites, and `Gene::introns` is the single
    // place that decides where they are: abutting or overlapping segments are a
    // continuous reading frame (a ribosomal-slippage join such as SARS-CoV-2
    // ORF1ab) and yield no intron at all.
    for intron in gene.introns() {
        // The donor is at the 5' end of the intron in transcript orientation and
        // the acceptor at its 3' end, which on the minus strand means the higher
        // and lower genomic coordinates respectively.
        let (donor_dir, acceptor_dir) = match gene.strand {
            Strand::Plus => (1, -1),
            Strand::Minus => (-1, 1),
        };
        let candidates = [
            classify_boundary(position, intron.exon_before_end, donor_dir, HalfKind::Donor),
            classify_boundary(
                position,
                intron.exon_after_start,
                acceptor_dir,
                HalfKind::Acceptor,
            ),
        ];
        for candidate in candidates.into_iter().flatten() {
            best = Some(match best {
                Some(current) if current.severity() >= candidate.severity() => current,
                _ => candidate,
            });
        }
    }
    best
}

/// Whether a genomic `position` falls inside a real intron of `gene`.
///
/// Only pairs of CDS segments that a genuine intron separates are considered,
/// so an unspliced transcript (a single segment, or a ribosomal-slippage join)
/// has no intronic positions at all. Splice sites are intronic too, but they are
/// classified by [`splice_consequence_for_position`], which is more specific;
/// this answers the wider question of "inside the transcript but not coding".
pub fn is_intronic_position(gene: &Gene, position: usize) -> bool {
    gene.introns()
        .any(|intron| position >= intron.start && position <= intron.end)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{joined_gene, spliced_gene, transcript_gene};
    use crate::variants::Biotype;

    #[test]
    fn plus_strand_donor_acceptor_and_regions() {
        // exon1 801-900, intron 901-1000, exon2 1001-1200.
        let g = spliced_gene("t", Strand::Plus, &[(801, 900), (1001, 1200)]);
        // Donor: first two intronic bases after exon1.
        assert_eq!(
            splice_consequence_for_position(&g, 901),
            Some(SpliceConsequence::Donor)
        );
        assert_eq!(
            splice_consequence_for_position(&g, 902),
            Some(SpliceConsequence::Donor)
        );
        // Acceptor: last two intronic bases before exon2.
        assert_eq!(
            splice_consequence_for_position(&g, 1000),
            Some(SpliceConsequence::Acceptor)
        );
        assert_eq!(
            splice_consequence_for_position(&g, 999),
            Some(SpliceConsequence::Acceptor)
        );
        // Intronic splice region (3-8 bp in) and exonic splice region (last 3 bp).
        assert_eq!(
            splice_consequence_for_position(&g, 905),
            Some(SpliceConsequence::Region)
        );
        assert_eq!(
            splice_consequence_for_position(&g, 900),
            Some(SpliceConsequence::Region)
        );
        assert_eq!(
            splice_consequence_for_position(&g, 1003),
            Some(SpliceConsequence::Region)
        );
        // Deep intron and deep exon: nothing.
        assert_eq!(splice_consequence_for_position(&g, 950), None);
        assert_eq!(splice_consequence_for_position(&g, 850), None);
    }

    #[test]
    fn minus_strand_flips_donor_and_acceptor() {
        // Transcript order (descending genomic): exon 1001-1200 then 801-900.
        // Intron 901-1000. In transcript orientation the donor is the low end of
        // the intron (just below 1001) and the acceptor the high end (just above
        // 900).
        let g = spliced_gene("t", Strand::Minus, &[(1001, 1200), (801, 900)]);
        assert_eq!(
            splice_consequence_for_position(&g, 1000),
            Some(SpliceConsequence::Donor)
        );
        assert_eq!(
            splice_consequence_for_position(&g, 999),
            Some(SpliceConsequence::Donor)
        );
        assert_eq!(
            splice_consequence_for_position(&g, 901),
            Some(SpliceConsequence::Acceptor)
        );
        assert_eq!(
            splice_consequence_for_position(&g, 902),
            Some(SpliceConsequence::Acceptor)
        );
        assert_eq!(splice_consequence_for_position(&g, 950), None);
    }

    #[test]
    fn single_exon_has_no_splice_sites() {
        let g = transcript_gene("t", Strand::Plus, &[(801, 900)]);
        assert_eq!(splice_consequence_for_position(&g, 901), None);
        // Written by hand on purpose: no annotation file produces a CDS record
        // with zero segments, so this one cannot come from the parser.
        let empty = Gene {
            name: "t".to_string(),
            start: 801,
            end: 900,
            strand: Strand::Plus,
            phase: 0,
            protein_offset: 0,
            transcript_id: Some("tx".to_string()),
            cds_segments: Vec::new(),
            biotype: Biotype::ProteinCoding,
        };
        assert_eq!(splice_consequence_for_position(&empty, 901), None);
    }

    #[test]
    fn abutting_segments_are_not_a_junction() {
        // Two CDS rows describing one continuous reading frame: no intron, so
        // the coding bases around the boundary are not splice sites.
        let g = joined_gene("t", Strand::Plus, &[(801, 900), (901, 1200)]);
        for pos in [899, 900, 901, 902, 903, 905] {
            assert_eq!(
                splice_consequence_for_position(&g, pos),
                None,
                "position {pos} is plain coding sequence"
            );
        }
    }

    #[test]
    fn ribosomal_slippage_join_is_not_a_junction() {
        // SARS-CoV-2 ORF1ab is annotated `join(266..13468,13468..21555)`: the
        // segments overlap by the re-read base, and nothing is spliced out.
        let g = joined_gene("t", Strand::Plus, &[(266, 13468), (13468, 21555)]);
        for pos in [13466, 13467, 13469, 13470, 13473] {
            assert_eq!(
                splice_consequence_for_position(&g, pos),
                None,
                "position {pos} is plain coding sequence"
            );
        }
    }

    #[test]
    fn a_real_intron_beside_a_slippage_join_still_yields_splice_sites() {
        // exon 1 (801-900), a real intron, exon 2 (1001-1200) continued by a
        // joined segment (1201-1400). Only the real intron has splice sites.
        let g = transcript_gene("t", Strand::Plus, &[(801, 900), (1001, 1200), (1201, 1400)]);
        assert_eq!(
            splice_consequence_for_position(&g, 901),
            Some(SpliceConsequence::Donor)
        );
        assert_eq!(
            splice_consequence_for_position(&g, 1000),
            Some(SpliceConsequence::Acceptor)
        );
        assert_eq!(splice_consequence_for_position(&g, 1201), None);
        assert_eq!(splice_consequence_for_position(&g, 1202), None);
    }

    #[test]
    fn intronic_positions_are_inside_a_real_intron() {
        // exon1 801-900, intron 901-1000, exon2 1001-1200.
        let g = spliced_gene("t", Strand::Plus, &[(801, 900), (1001, 1200)]);
        for pos in [901, 950, 1000] {
            assert!(is_intronic_position(&g, pos), "{pos} is in the intron");
        }
        for pos in [900, 1001, 850, 1100] {
            assert!(!is_intronic_position(&g, pos), "{pos} is exonic");
        }

        // Minus strand: transcript order is descending, the intron is the same
        // genomic interval.
        let m = spliced_gene("t", Strand::Minus, &[(1001, 1200), (801, 900)]);
        for pos in [901, 950, 1000] {
            assert!(is_intronic_position(&m, pos), "{pos} is in the intron");
        }
        assert!(!is_intronic_position(&m, 900));
        assert!(!is_intronic_position(&m, 1001));
    }

    #[test]
    fn an_unspliced_transcript_has_no_intronic_positions() {
        // Single segment, abutting segments and a slippage join all lack introns,
        // so no position between or around them is intronic.
        let single = transcript_gene("t", Strand::Plus, &[(801, 900)]);
        assert!(!is_intronic_position(&single, 950));
        let abutting = joined_gene("t", Strand::Plus, &[(801, 900), (901, 1200)]);
        assert!(!is_intronic_position(&abutting, 900));
        assert!(!is_intronic_position(&abutting, 901));
        let slippage = joined_gene("t", Strand::Plus, &[(266, 13468), (13468, 21555)]);
        assert!(!is_intronic_position(&slippage, 13468));
        assert!(!is_intronic_position(&slippage, 13469));
    }

    #[test]
    fn minus_strand_slippage_join_is_not_a_junction() {
        // Transcript order descending: 1001-1200 then 800-1000 (abutting).
        let g = joined_gene("t", Strand::Minus, &[(1001, 1200), (800, 1000)]);
        for pos in [999, 1000, 1001, 1002] {
            assert_eq!(
                splice_consequence_for_position(&g, pos),
                None,
                "position {pos} is plain coding sequence"
            );
        }
    }
}
