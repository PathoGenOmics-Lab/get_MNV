//! COSMIC single-base-substitution (SBS96) classes.
//!
//! Each single-base substitution is named by its trinucleotide context with the
//! mutated base in the middle, written from the pyrimidine of the base pair:
//! `A[C>T]G`. A purine reference (`G` or `A`) is read from the other strand, so
//! the context is reverse-complemented with it. Six substitution types times 16
//! flanking pairs make the 96 classes of a mutational spectrum.
//!
//! The class is attached to every row as one entry per substitution component,
//! in the order of the row's Event Components, so a report can pair them and
//! count each substituted base once however many rows mention it.

use super::{parse_component_label, AlleleComponentKind, VariantInfo};

/// Placeholder for a component whose context cannot be read: at a contig end,
/// beside an ambiguous base, or with a reference that is not a single base.
pub const NO_CONTEXT: &str = "-";

fn complement(base: u8) -> Option<u8> {
    match base {
        b'A' => Some(b'T'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        b'T' => Some(b'A'),
        _ => None,
    }
}

/// The SBS96 class of a substitution of `ref_base` by `alt_base` at 1-based
/// `position` of `sequence`, or `None` when either flank is missing or is not
/// one of A, C, G and T, or the reference there is not `ref_base`.
pub fn sbs96_class(sequence: &str, position: usize, ref_base: u8, alt_base: u8) -> Option<String> {
    let bytes = sequence.as_bytes();
    if position < 2 || position >= bytes.len() {
        return None;
    }
    let upper = |base: u8| base.to_ascii_uppercase();
    let (before, at, after) = (
        upper(bytes[position - 2]),
        upper(bytes[position - 1]),
        upper(bytes[position]),
    );
    let (ref_base, alt_base) = (upper(ref_base), upper(alt_base));
    if at != ref_base || ref_base == alt_base {
        return None;
    }
    // Every base must be one of the four for the class to exist.
    for base in [before, at, after, alt_base] {
        complement(base)?;
    }
    let (five, from, to, three) = if matches!(ref_base, b'C' | b'T') {
        (before, ref_base, alt_base, after)
    } else {
        (
            complement(after)?,
            complement(ref_base)?,
            complement(alt_base)?,
            complement(before)?,
        )
    };
    Some(format!(
        "{}[{}>{}]{}",
        five as char, from as char, to as char, three as char
    ))
}

/// One SBS96 class (or [`NO_CONTEXT`]) per substitution component of the row,
/// in Event Components order. Empty when the row has no substitution.
pub fn sbs96_for_variant(variant: &VariantInfo, sequence: &str) -> Vec<String> {
    variant
        .event_components
        .iter()
        .filter_map(|label| parse_component_label(label))
        .filter(|component| component.kind == AlleleComponentKind::Snp)
        .map(|component| {
            match (
                component.ref_allele.as_bytes(),
                component.alt_allele.as_bytes(),
            ) {
                ([from], [to]) => sbs96_class(sequence, component.position, *from, *to),
                _ => None,
            }
            .unwrap_or_else(|| NO_CONTEXT.to_string())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_pyrimidine_reference_keeps_its_own_strand() {
        //                1234567
        let sequence = "AACGTTA";
        assert_eq!(
            sbs96_class(sequence, 3, b'C', b'T').as_deref(),
            Some("A[C>T]G")
        );
        assert_eq!(
            sbs96_class(sequence, 5, b'T', b'G').as_deref(),
            Some("G[T>G]T")
        );
    }

    #[test]
    fn a_purine_reference_is_read_from_the_other_strand() {
        // G>A in A-G-T is C>T in A-C-T on the opposite strand.
        let sequence = "AAGTTA";
        assert_eq!(
            sbs96_class(sequence, 3, b'G', b'A').as_deref(),
            Some("A[C>T]T")
        );
        // Lowercase (soft-masked) reference reads the same.
        assert_eq!(
            sbs96_class("aagtta", 3, b'g', b'a').as_deref(),
            Some("A[C>T]T")
        );
    }

    #[test]
    fn no_class_without_two_clean_flanks_and_a_matching_reference() {
        assert_eq!(sbs96_class("ACG", 1, b'A', b'T'), None);
        assert_eq!(sbs96_class("ACG", 3, b'G', b'T'), None);
        assert_eq!(sbs96_class("ANCG", 3, b'C', b'T'), None);
        assert_eq!(sbs96_class("ACG", 2, b'A', b'T'), None);
        assert_eq!(sbs96_class("ACG", 2, b'C', b'C'), None);
    }
}
