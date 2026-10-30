//! Left-alignment of indels against the reference, the normalization that
//! `bcftools norm -f ref.fa` performs.
//!
//! An insertion or deletion inside a homopolymer or a tandem repeat can be
//! written at any of several places that all describe the same change.
//! Aligners put the gap at the leftmost one, and get_MNV counts an indel's read
//! support against the exact gap in each read's CIGAR, so an allele written
//! further right than the aligner put it found no supporting reads at all.
//!
//! The algorithm is the one bcftools runs (Tan, Abecasis and Kang 2015): while
//! both alleles end in the same base, drop that base, and when an allele runs
//! out, extend both by the reference base to their left; then drop the bases
//! both alleles start with, keeping at least one in each. What comes out:
//!
//! - an insertion or deletion at its leftmost equivalent place, anchored on the
//!   base before it, or on the base after it when it starts at the contig's
//!   first base and there is none before;
//! - a substitution (SNV or MNV) or a complex allele (a delins, or an indel and
//!   a substitution in one record) trimmed of the bases both alleles share, but
//!   not moved: once trimmed, its two alleles end in different bases, so there
//!   is no other place that describes it;
//! - and left as written: a symbolic or breakend allele, `*`, an allele with a
//!   base other than A, C, G, T or N (bcftools refuses those), a REF that
//!   disagrees with the reference (the run stops on that mismatch further on,
//!   as bcftools does), and a record whose ALT repeats its REF.
//!
//! Bases compare without regard to case and are written in upper case. One
//! deliberate difference from bcftools: an ambiguity code in the reference,
//! such as `R`, matches only itself here. bcftools reads it as its first base
//! (`R` as `A`), moves an indel through it and writes that base into REF, which
//! then disagrees with the FASTA the run checks every REF against.

use std::collections::VecDeque;

/// Whether every base of an allele is one left-alignment can compare: A, C, G,
/// T or N, in either case.
fn plain_bases(allele: &str) -> bool {
    !allele.is_empty()
        && allele
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

/// The record at 1-based `pos` with alleles `ref_allele`/`alt_allele`,
/// normalized as `bcftools norm -f` writes it against `sequence`, the record's
/// whole contig. `None` when it stays exactly as written (see the module notes
/// for what is never touched).
pub(crate) fn left_align(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    sequence: &str,
) -> Option<(usize, String, String)> {
    if pos == 0
        || !plain_bases(ref_allele)
        || !plain_bases(alt_allele)
        || ref_allele.eq_ignore_ascii_case(alt_allele)
    {
        return None;
    }
    let seq = sequence.as_bytes();
    // REF as a stretch of the reference, 0-based and half-open. It grows to
    // the left as the record moves, so it never has to be copied.
    let mut start = pos - 1;
    let mut end = start + ref_allele.len();
    if end > seq.len() || !seq[start..end].eq_ignore_ascii_case(ref_allele.as_bytes()) {
        return None;
    }
    let mut alt: VecDeque<u8> = alt_allele.bytes().map(|b| b.to_ascii_uppercase()).collect();

    // Both alleles end in the same base: it is context, not change, so drop it.
    // An allele left empty takes the reference base before the record, and so
    // does the other, which walks an indel leftwards through its repeat.
    while let Some(&last) = alt.back() {
        if !seq[end - 1].eq_ignore_ascii_case(&last) {
            break;
        }
        // At the contig's first base there is nothing to the left to take.
        if (end - start).min(alt.len()) <= 1 && start == 0 {
            break;
        }
        end -= 1;
        alt.pop_back();
        if end == start || alt.is_empty() {
            start -= 1;
            alt.push_front(seq[start].to_ascii_uppercase());
        }
    }
    // Then the bases both start with, keeping one in each allele.
    while end - start > 1 && alt.len() > 1 && seq[start].eq_ignore_ascii_case(&alt[0]) {
        start += 1;
        alt.pop_front();
    }

    let aligned_ref = sequence[start..end].to_ascii_uppercase();
    let aligned_alt: String = alt.into_iter().map(char::from).collect();
    let unchanged = start + 1 == pos
        && aligned_ref.eq_ignore_ascii_case(ref_allele)
        && aligned_alt.eq_ignore_ascii_case(alt_allele);
    (!unchanged).then_some((start + 1, aligned_ref, aligned_alt))
}

#[cfg(test)]
mod tests {
    use super::left_align;

    //               1        10        20        30
    //               |        |         |         |
    const SEQ: &str = "AAAACGTCACACACGTTTTTNNNNNCGRAAACGTAGGGG";

    fn aligned(pos: usize, ref_allele: &str, alt_allele: &str) -> (usize, String, String) {
        left_align(pos, ref_allele, alt_allele, SEQ)
            .unwrap_or_else(|| (pos, ref_allele.to_string(), alt_allele.to_string()))
    }

    fn record(pos: usize, ref_allele: &str, alt_allele: &str) -> (usize, String, String) {
        (pos, ref_allele.to_string(), alt_allele.to_string())
    }

    // Every expectation here is what bcftools 1.24 `norm -f` writes for the
    // same record against the same sequence.

    #[test]
    fn a_deletion_in_a_homopolymer_moves_to_its_left_end() {
        assert_eq!(aligned(19, "TT", "T"), record(15, "GT", "G"));
    }

    #[test]
    fn an_insertion_in_a_homopolymer_moves_to_its_left_end() {
        assert_eq!(aligned(17, "T", "TT"), record(15, "G", "GT"));
        assert_eq!(aligned(39, "G", "GG"), record(35, "A", "AG"));
    }

    #[test]
    fn a_repeat_unit_moves_whole_and_comes_out_rotated() {
        // CACACAC at 8-14: deleting one CA anywhere in it removes 8-9.
        assert_eq!(aligned(11, "ACA", "A"), record(7, "TCA", "T"));
        // The same deletion with its context written out on both sides.
        assert_eq!(aligned(11, "ACAC", "AC"), record(7, "TCA", "T"));
    }

    #[test]
    fn an_insertion_anchored_on_its_right_moves_onto_the_base_before_it() {
        assert_eq!(aligned(18, "T", "AT"), record(17, "T", "TA"));
    }

    #[test]
    fn an_insertion_that_is_not_a_copy_of_its_neighbour_stays() {
        assert_eq!(left_align(13, "A", "AAC", SEQ), None);
        assert_eq!(left_align(20, "T", "TN", SEQ), None);
    }

    #[test]
    fn at_the_contig_start_the_anchor_is_the_first_base() {
        assert_eq!(aligned(2, "AA", "A"), record(1, "AA", "A"));
        assert_eq!(aligned(2, "AAA", "A"), record(1, "AAA", "A"));
        assert_eq!(aligned(3, "A", "AA"), record(1, "A", "AA"));
        assert_eq!(left_align(1, "A", "CA", SEQ), None);
    }

    #[test]
    fn substitutions_and_complex_alleles_are_trimmed_but_never_moved() {
        assert_eq!(aligned(36, "GGG", "GTG"), record(37, "G", "T"));
        assert_eq!(aligned(11, "ACAC", "TGC"), record(11, "ACA", "TG"));
        assert_eq!(left_align(37, "G", "T", SEQ), None);
    }

    #[test]
    fn n_is_a_base_like_any_other() {
        assert_eq!(aligned(24, "NN", "N"), record(20, "TN", "T"));
    }

    #[test]
    fn case_is_ignored_and_the_result_is_upper_case() {
        assert_eq!(aligned(11, "aca", "a"), record(7, "TCA", "T"));
    }

    #[test]
    fn what_bcftools_refuses_or_cannot_move_is_left_as_written() {
        // A REF the reference does not have, past its end, or at position 0.
        assert_eq!(left_align(11, "ACG", "A", SEQ), None);
        assert_eq!(left_align(39, "GG", "G", SEQ), None);
        assert_eq!(left_align(0, "A", "AA", SEQ), None);
        // Symbolic, spanning-deletion and ambiguity-code alleles.
        assert_eq!(left_align(15, "G", "<DEL>", SEQ), None);
        assert_eq!(left_align(15, "G", "*", SEQ), None);
        assert_eq!(left_align(31, "A", "AR", SEQ), None);
        assert_eq!(left_align(28, "RA", "R", SEQ), None);
        // A record that describes no change.
        assert_eq!(left_align(19, "TT", "TT", SEQ), None);
    }

    #[test]
    fn an_ambiguity_code_in_the_reference_stops_the_walk() {
        // bcftools reads R as A and writes `27 GA>G`, a REF the FASTA does
        // not have; the record stops on the R instead.
        assert_eq!(aligned(30, "AA", "A"), record(28, "RA", "R"));
    }

    /// A contig made of short units repeated a few times each, so that most
    /// places sit in a homopolymer or a tandem repeat.
    fn repetitive_sequence() -> impl proptest::strategy::Strategy<Value = String> {
        use proptest::prelude::*;
        let unit = prop::sample::select(vec!["A", "C", "G", "T", "AC", "AG", "CAG", "ACGT"]);
        prop::collection::vec((unit, 1usize..6), 2..12).prop_map(|runs| {
            runs.iter()
                .map(|(unit, times)| unit.repeat(*times))
                .collect::<String>()
        })
    }

    /// What a record turns the reference into.
    fn haplotype(sequence: &str, pos: usize, ref_allele: &str, alt_allele: &str) -> String {
        format!(
            "{}{alt_allele}{}",
            &sequence[..pos - 1],
            &sequence[pos - 1 + ref_allele.len()..]
        )
    }

    /// Where a pure insertion or deletion belongs, found the slow way: every
    /// place a change of that length could produce the same haplotype, tried
    /// from the left, and the first written as VCF anchors it.
    fn leftmost_by_search(sequence: &str, haplotype: &str) -> (usize, String, String) {
        let (s, h) = (sequence.as_bytes(), haplotype.as_bytes());
        if h.len() < s.len() {
            let length = s.len() - h.len();
            let at = (0..=h.len())
                .find(|&p| s[..p] == h[..p] && s[p + length..] == h[p..])
                .expect("a deletion produced this haplotype");
            if at == 0 {
                record(1, &sequence[..=length], &sequence[length..=length])
            } else {
                record(at, &sequence[at - 1..at + length], &sequence[at - 1..at])
            }
        } else {
            let length = h.len() - s.len();
            let at = (0..=s.len())
                .find(|&p| h[..p] == s[..p] && h[p + length..] == s[p..])
                .expect("an insertion produced this haplotype");
            let inserted = &haplotype[at..at + length];
            if at == 0 {
                record(1, &sequence[..1], &format!("{inserted}{}", &sequence[..1]))
            } else {
                let anchor = &sequence[at - 1..at];
                record(at, anchor, &format!("{anchor}{inserted}"))
            }
        }
    }

    proptest::proptest! {
        /// Any spelling of a pure indel, anchored left or right and padded with
        /// any amount of context, comes out at the leftmost place a search over
        /// every placement finds.
        #[test]
        fn every_spelling_of_an_indel_reaches_the_leftmost_place(
            sequence in repetitive_sequence(),
            deletion in proptest::bool::ANY,
            at in 0usize..1000,
            length in 1usize..7,
            copy_from in 0usize..1000,
            pad_left in 0usize..4,
            pad_right in 0usize..4,
        ) {
            let n = sequence.len();
            // The gap: bases at..at+removed come out, `inserted` goes in at `at`,
            // which for an insertion may be past the last base.
            let at = if deletion { at % n } else { at % (n + 1) };
            let (removed, inserted) = if deletion {
                (length.min(n - at), String::new())
            } else {
                // Mostly a copy of nearby bases, which is what repeats grow by.
                let from = copy_from % n;
                let copied = &sequence[from..(from + length).min(n)];
                (0, copied.to_string())
            };
            proptest::prop_assume!(removed < n);
            // A window of the reference around the gap, wide enough that the
            // shorter allele keeps a base.
            let mut first = at.saturating_sub(pad_left);
            let mut last = (at + removed + pad_right).min(n);
            if last - first <= removed {
                if first > 0 { first -= 1 } else { last += 1 }
            }
            proptest::prop_assume!(last <= n && last > first);
            let ref_allele = &sequence[first..last];
            let alt_allele = format!(
                "{}{inserted}{}",
                &sequence[first..at],
                &sequence[at + removed..last]
            );
            let pos = first + 1;
            let wanted = leftmost_by_search(
                &sequence,
                &haplotype(&sequence, pos, ref_allele, &alt_allele),
            );
            let got = aligned_against(&sequence, pos, ref_allele, &alt_allele);
            proptest::prop_assert_eq!(
                &got, &wanted,
                "{}:{}>{} on {}", pos, ref_allele, alt_allele, sequence
            );
        }

        /// Whatever the alleles, the result describes the same change, cannot
        /// be trimmed or moved any further, and a second pass leaves it alone.
        #[test]
        fn any_record_keeps_its_change_and_comes_out_normalized(
            sequence in repetitive_sequence(),
            at in 0usize..1000,
            ref_length in 1usize..6,
            alt_allele in "[ACGT]{1,6}",
        ) {
            let n = sequence.len();
            let first = at % n;
            let ref_allele = &sequence[first..(first + ref_length).min(n)];
            proptest::prop_assume!(!ref_allele.eq_ignore_ascii_case(&alt_allele));
            let pos = first + 1;
            let (new_pos, new_ref, new_alt) =
                aligned_against(&sequence, pos, ref_allele, &alt_allele);

            proptest::prop_assert_eq!(
                haplotype(&sequence, new_pos, &new_ref, &new_alt),
                haplotype(&sequence, pos, ref_allele, &alt_allele)
            );
            let shared_last = new_ref.as_bytes().last() == new_alt.as_bytes().last();
            proptest::prop_assert!(
                !shared_last || (new_pos == 1 && new_ref.len().min(new_alt.len()) == 1),
                "{}:{}>{} still ends in a shared base", new_pos, new_ref, new_alt
            );
            let shared_first = new_ref.as_bytes()[0] == new_alt.as_bytes()[0];
            proptest::prop_assert!(
                !shared_first || new_ref.len().min(new_alt.len()) == 1,
                "{}:{}>{} still starts with a shared base", new_pos, new_ref, new_alt
            );
            if new_ref.len() == new_alt.len() {
                // A substitution only ever loses context.
                proptest::prop_assert!(
                    new_pos >= pos && new_pos + new_ref.len() <= pos + ref_allele.len()
                );
            }
            proptest::prop_assert_eq!(left_align(new_pos, &new_ref, &new_alt, &sequence), None);
        }
    }

    fn aligned_against(
        sequence: &str,
        pos: usize,
        ref_allele: &str,
        alt_allele: &str,
    ) -> (usize, String, String) {
        left_align(pos, ref_allele, alt_allele, sequence)
            .unwrap_or_else(|| record(pos, ref_allele, alt_allele))
    }

    #[test]
    fn a_long_repeat_is_walked_in_one_pass() {
        // The last AT of a 50,000-unit run, deleted, belongs after the G.
        let sequence = format!("G{}C", "AT".repeat(50_000));
        let anchor = sequence.len() - 3;
        assert_eq!(
            left_align(anchor, "TAT", "T", &sequence),
            Some(record(1, "GAT", "G"))
        );
    }
}
