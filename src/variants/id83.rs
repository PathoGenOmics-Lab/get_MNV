//! COSMIC ID83 classes of small insertions and deletions.
//!
//! The scheme SigProfiler uses, labelled the way it labels them
//! (`size:Type:Subtype:count`):
//!
//! * one-base indels by the base, written from the pyrimidine of the pair (`C`
//!   or `T`), and the length of the homopolymer they sit in: `1:Del:C:0` is a C
//!   deleted where no other C is next to it, `1:Ins:T:5` a T inserted next to
//!   five or more;
//! * longer indels by size (2, 3, 4, 5 or more) and by how many further copies
//!   of the indel's own sequence sit next to it (`R`, 0 to 5 or more);
//! * longer deletions with no such copy, by the microhomology at their
//!   breakpoint (`M`): how many bases at one end of the deleted sequence are
//!   repeated on the other side of it.
//!
//! The copies are counted on both sides of the event, so the class is the same
//! wherever in a repeat the caller placed it, as it would be after
//! left-alignment.

fn same(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

/// How many whole copies of `unit` sit immediately before 1-based position
/// `start` (exclusive) and immediately after `end` (exclusive), in `seq`.
fn flanking_copies(seq: &[u8], unit: &[u8], start: usize, end: usize) -> usize {
    let length = unit.len();
    let mut copies = 0;
    // To the right: positions end+1 ... end+length, then the next copy.
    let mut at = end;
    while at + length <= seq.len()
        && seq[at..at + length]
            .iter()
            .zip(unit)
            .all(|(&a, &b)| same(a, b))
    {
        copies += 1;
        at += length;
    }
    // To the left: positions start-length ... start-1.
    let mut at = start.saturating_sub(1);
    while at >= length
        && seq[at - length..at]
            .iter()
            .zip(unit)
            .all(|(&a, &b)| same(a, b))
    {
        copies += 1;
        at -= length;
    }
    copies
}

/// The longest stretch of `deleted` repeated across one of its breakpoints:
/// its start against the bases after it, or its end against the bases before.
fn microhomology(seq: &[u8], deleted: &[u8], start: usize, end: usize) -> usize {
    let length = deleted.len();
    let mut best = 0;
    for k in (1..length).rev() {
        let after = end + k <= seq.len()
            && seq[end..end + k]
                .iter()
                .zip(&deleted[..k])
                .all(|(&a, &b)| same(a, b));
        let before = start > k
            && seq[start - 1 - k..start - 1]
                .iter()
                .zip(&deleted[length - k..])
                .all(|(&a, &b)| same(a, b));
        if after || before {
            best = k;
            break;
        }
    }
    best
}

fn pyrimidine(base: u8) -> Option<char> {
    match base.to_ascii_uppercase() {
        b'C' | b'G' => Some('C'),
        b'T' | b'A' => Some('T'),
        _ => None,
    }
}

/// The ID83 class of a pure deletion or insertion given as an anchored REF/ALT
/// at 1-based `pos` in `sequence`, or `None` for anything else (a substitution,
/// a deletion-insertion, a symbolic allele, an ambiguous base).
pub fn id83_class(
    sequence: &str,
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<String> {
    let seq = sequence.as_bytes();
    let (start, removed, added) = crate::variants::hgvs::indel_core(pos, ref_allele, alt_allele);
    let (kind, unit, first, last) = match (removed.is_empty(), added.is_empty()) {
        (false, true) => {
            let end = start + removed.len() - 1;
            if end > seq.len()
                || !seq[start - 1..end]
                    .iter()
                    .zip(removed.bytes())
                    .all(|(&a, b)| same(a, b))
            {
                return None;
            }
            ("Del", removed, start, end)
        }
        // An insertion falls between `start - 1` and `start`: no base of the
        // reference is inside it.
        (true, false) => ("Ins", added, start, start - 1),
        _ => return None,
    };
    if !unit.bytes().all(|b| b"ACGT".contains(&b)) || first == 0 {
        return None;
    }
    let size = unit.len();
    let copies = flanking_copies(seq, unit.as_bytes(), first, last);
    if size == 1 {
        let base = pyrimidine(unit.as_bytes()[0])?;
        return Some(format!("1:{kind}:{base}:{}", copies.min(5)));
    }
    let size_label = size.min(5);
    if kind == "Del" && copies == 0 {
        let homology = microhomology(seq, unit.as_bytes(), first, last);
        if homology > 0 {
            return Some(format!("{size_label}:Del:M:{}", homology.min(5)));
        }
    }
    Some(format!("{size_label}:{kind}:R:{}", copies.min(5)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_base_indels_by_base_and_homopolymer() {
        //               123456789
        let sequence = "AGCCCATTA";
        // One C of three deleted: two more Cs beside it.
        assert_eq!(
            id83_class(sequence, 2, "GC", "G").as_deref(),
            Some("1:Del:C:2")
        );
        // Wherever in the run the caller put it.
        assert_eq!(
            id83_class(sequence, 4, "CC", "C").as_deref(),
            Some("1:Del:C:2")
        );
        // A G is written from the pyrimidine: C. No other G next to it.
        assert_eq!(
            id83_class(sequence, 1, "AG", "A").as_deref(),
            Some("1:Del:C:0")
        );
        // A T inserted between the two Ts: two copies already there.
        assert_eq!(
            id83_class(sequence, 7, "T", "TT").as_deref(),
            Some("1:Ins:T:2")
        );
        // An A inserted where no A is next to it reads as a T, zero copies.
        assert_eq!(
            id83_class(sequence, 3, "C", "CA").as_deref(),
            Some("1:Ins:T:0")
        );
    }

    #[test]
    fn longer_indels_by_size_repeats_and_microhomology() {
        //               1234567890123
        let sequence = "GACACACAGTTTG";
        // AC three times at 2-7: deleting one leaves two more beside it.
        assert_eq!(
            id83_class(sequence, 1, "GAC", "G").as_deref(),
            Some("2:Del:R:2")
        );
        // Inserting AC next to the run: three copies already there.
        assert_eq!(
            id83_class(sequence, 1, "G", "GAC").as_deref(),
            Some("2:Ins:R:3")
        );
        // Deleting GCTA (3-6) from TTGCTAGCAA: no copy of GCTA beside it, but
        // its first two bases, GC, are repeated right after it.
        assert_eq!(
            id83_class("TTGCTAGCAA", 2, "TGCTA", "T").as_deref(),
            Some("4:Del:M:2")
        );
        // Five or more bases: the 5+ size class, here one more copy beside it.
        assert_eq!(
            id83_class("GACGTAACGTAC", 1, "GACGTA", "G").as_deref(),
            Some("5:Del:R:1")
        );
        // Neither a repeat nor microhomology.
        assert_eq!(
            id83_class("AAGCTTTT", 2, "AGC", "A").as_deref(),
            Some("2:Del:R:0")
        );
    }

    #[test]
    fn substitutions_and_delins_have_no_indel_class() {
        assert_eq!(id83_class("ACGT", 2, "C", "T"), None);
        assert_eq!(id83_class("ACGT", 2, "CG", "A"), None);
        assert_eq!(id83_class("ACGT", 2, "C", "<DEL>"), None);
    }
}
