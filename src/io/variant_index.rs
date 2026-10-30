//! Where each variant of a contig acts, worked out once, so that a gene looks
//! only at the variants near it.
//!
//! A bacterial genome has thousands of genes and a sample hundreds of variants.
//! Asking every gene about every variant, and decomposing the allele again for
//! each question, was nearly all of a run's time, and a cohort paid it once per
//! sample.

use super::vcf::VcfPosition;
use crate::variants::splice::SPLICE_REACH;
use crate::variants::Gene;

/// Spans of bases kept in the order of the list they came from, which can say
/// what stretch of that list touches a given span.
///
/// The list does not have to be sorted. A query returns a stretch of the list,
/// so whoever filters it still sees the items in the order the list holds them;
/// a sorted list just makes that stretch short.
struct Stretches {
    /// For each item, the highest base reached by it or any item before it.
    /// Never decreases, so it can be searched.
    furthest_end: Vec<usize>,
    /// For each item, the lowest base reached by it or any item after it.
    /// Never decreases either.
    nearest_start: Vec<usize>,
}

impl Stretches {
    fn new(spans: &[(usize, usize)]) -> Self {
        let mut end = 0;
        let furthest_end = spans
            .iter()
            .map(|&(_, span_end)| {
                end = end.max(span_end);
                end
            })
            .collect();
        let mut nearest_start = vec![0; spans.len()];
        let mut start = usize::MAX;
        for (slot, &(span_start, _)) in nearest_start.iter_mut().zip(spans).rev() {
            start = start.min(span_start);
            *slot = start;
        }
        Stretches {
            furthest_end,
            nearest_start,
        }
    }

    /// The stretch holding every item whose span touches a base from `start`
    /// to `end` inclusive. Everything before it ends below `start`, and
    /// everything after it begins above `end`.
    fn stretch(&self, start: usize, end: usize) -> std::ops::Range<usize> {
        let first = self.furthest_end.partition_point(|&reach| reach < start);
        let last = self.nearest_start.partition_point(|&reach| reach <= end);
        first..last.max(first)
    }
}

/// The variants of one contig in their list order, with the bases each acts on.
pub struct VariantIndex<'a> {
    variants: &'a [VcfPosition],
    stretches: Stretches,
}

impl<'a> VariantIndex<'a> {
    pub fn new(variants: &'a [VcfPosition]) -> Self {
        let spans: Vec<(usize, usize)> = variants.iter().map(VcfPosition::acting_span).collect();
        Self {
            variants,
            stretches: Stretches::new(&spans),
        }
    }

    /// Every variant that acts on a base from `start` to `end` inclusive, as the
    /// stretch of the list that holds them. The stretch can also hold variants
    /// that act elsewhere, so callers still ask their own question of each; it
    /// never leaves one out.
    pub fn near(&self, start: usize, end: usize) -> &'a [VcfPosition] {
        &self.variants[self.stretches.stretch(start, end)]
    }

    /// Every variant a gene could claim: one on its span or its CDS segments,
    /// an insertion at one of its junctions, or one in the splice region or the
    /// intron of a junction.
    pub fn near_gene(&self, gene: &Gene) -> &'a [VcfPosition] {
        let (start, end) = gene_reach(gene);
        self.near(start, end)
    }
}

/// The genes of one contig in their list order, with the bases on which a
/// variant can still belong to each, so a variant asks only the genes near it.
pub(crate) struct GeneIndex<'a> {
    genes: &'a [Gene],
    stretches: Stretches,
}

impl<'a> GeneIndex<'a> {
    pub(crate) fn new(genes: &'a [Gene]) -> Self {
        let spans: Vec<(usize, usize)> = genes.iter().map(gene_reach).collect();
        Self {
            genes,
            stretches: Stretches::new(&spans),
        }
    }

    /// Every gene that could claim a base from `start` to `end` inclusive, as
    /// the stretch of the list that holds them, in list order. A gene outside
    /// it has no coding base, splice site or intron there.
    pub(crate) fn near(&self, start: usize, end: usize) -> &'a [Gene] {
        &self.genes[self.stretches.stretch(start, end)]
    }
}

/// The bases on which a variant can still belong to a gene: the gene's span
/// and its CDS segments, widened by the splice region, which runs up to
/// `SPLICE_REACH` bases from the exon base bounding a junction. Intronic
/// positions lie between segments, so they are inside already.
pub(crate) fn gene_reach(gene: &Gene) -> (usize, usize) {
    let mut low = gene.start.min(gene.end);
    let mut high = gene.start.max(gene.end);
    for segment in &gene.cds_segments {
        low = low.min(segment.start).min(segment.end);
        high = high.max(segment.start).max(segment.end);
    }
    (
        low.saturating_sub(SPLICE_REACH),
        high.saturating_add(SPLICE_REACH),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::spliced_gene;
    use crate::variants::Strand;

    fn variant(position: usize, ref_allele: &str, alt_allele: &str) -> VcfPosition {
        VcfPosition {
            record_start: position,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            original_dp: None,
            original_freq: None,
            original_info: None,
            declared_phase: None,
        }
    }

    /// The span is where the components act, which is not always the record's
    /// own span: a right-anchored insertion acts on the base before its POS, and
    /// a padded record acts only on the bases it changes.
    #[test]
    fn a_record_acts_where_its_components_do() {
        assert_eq!(variant(10, "A", "G").acting_span(), (10, 10));
        assert_eq!(variant(10, "ATG", "A").acting_span(), (11, 12));
        assert_eq!(variant(30, "T", "AT").acting_span(), (29, 29));
        assert_eq!(variant(40, "GCTGCTG", "GCTGCTA").acting_span(), (46, 46));
        assert_eq!(variant(50, "AC", "GT").acting_span(), (50, 51));
        assert_eq!(variant(60, "A", "A").acting_span(), (60, 60));
    }

    /// Whatever the list holds, sorted or not, with long deletions or not, a
    /// query must hand back every variant whose span touches the interval, in
    /// list order, so that filtering the stretch gives what filtering the whole
    /// list gave.
    #[test]
    fn a_query_never_leaves_out_a_variant_that_acts_there() {
        let variants = vec![
            variant(5, "A", "G"),
            variant(20, "ACGTACGTACGTACGTACGTA", "A"),
            variant(18, "C", "T"),
            variant(30, "T", "AT"),
            variant(31, "G", "C"),
            variant(45, "A", "A"),
            variant(44, "AC", "GT"),
            variant(90, "G", "GTT"),
            variant(70, "C", "A"),
        ];
        let index = VariantIndex::new(&variants);
        for start in 0..100 {
            for end in start..100 {
                let stretch = index.stretches.stretch(start, end);
                for (i, v) in variants.iter().enumerate() {
                    let (lo, hi) = v.acting_span();
                    if lo <= end && hi >= start {
                        assert!(
                            stretch.contains(&i),
                            "{start}-{end} left out {}:{}>{}",
                            v.record_start,
                            v.ref_allele,
                            v.alt_allele
                        );
                    }
                }
            }
        }
        assert!(index.near(95, 99).is_empty());
        // The list is not sorted: the insertion at 90 is listed before the
        // substitution at 70, so the stretch for 70 holds it too.
        assert_eq!(index.near(70, 70).len(), 2);
    }

    /// A gene's reach covers its CDS segments and the splice region past them,
    /// even where that region runs beyond the gene's own coordinates.
    #[test]
    fn a_gene_reaches_its_splice_region() {
        let gene = spliced_gene("t", Strand::Plus, &[(801, 900), (902, 903)]);
        assert_eq!(gene_reach(&gene), (801 - SPLICE_REACH, 903 + SPLICE_REACH));
        let variants = vec![
            variant(700, "A", "G"),
            variant(908, "C", "T"),
            variant(990, "G", "A"),
        ];
        let index = VariantIndex::new(&variants);
        let near = index.near_gene(&gene);
        assert_eq!(near.len(), 1);
        assert_eq!(near[0].record_start, 908);
    }
}
