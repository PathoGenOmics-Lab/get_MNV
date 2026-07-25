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
    for pair in gene.cds_segments.windows(2) {
        let (earlier, later) = (&pair[0], &pair[1]);
        // `cds_segments` is in transcript order, so `earlier` is the 5' exon of
        // the junction. On the plus strand its 3' end is the higher genomic base
        // and the intron lies above it; on the minus strand it is the lower base
        // and the intron lies below.
        let (donor_boundary, donor_dir, acceptor_boundary, acceptor_dir) = match gene.strand {
            Strand::Plus => (earlier.end, 1, later.start, -1),
            Strand::Minus => (earlier.start, -1, later.end, 1),
        };
        let candidates = [
            classify_boundary(position, donor_boundary, donor_dir, HalfKind::Donor),
            classify_boundary(
                position,
                acceptor_boundary,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::CdsSegment;

    fn gene(strand: Strand, segments: &[(usize, usize)]) -> Gene {
        Gene {
            name: "t".to_string(),
            start: segments.iter().map(|s| s.0).min().unwrap_or(0),
            end: segments.iter().map(|s| s.1).max().unwrap_or(0),
            strand,
            phase: 0,
            protein_offset: 0,
            transcript_id: Some("tx".to_string()),
            cds_segments: segments
                .iter()
                .map(|&(start, end)| CdsSegment { start, end })
                .collect(),
        }
    }

    #[test]
    fn plus_strand_donor_acceptor_and_regions() {
        // exon1 801-900, intron 901-1000, exon2 1001-1200.
        let g = gene(Strand::Plus, &[(801, 900), (1001, 1200)]);
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
        let g = gene(Strand::Minus, &[(1001, 1200), (801, 900)]);
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
        let g = gene(Strand::Plus, &[(801, 900)]);
        assert_eq!(splice_consequence_for_position(&g, 901), None);
        let empty = gene(Strand::Plus, &[]);
        assert_eq!(splice_consequence_for_position(&empty, 901), None);
    }
}
