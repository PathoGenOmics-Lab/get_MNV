//! Local haplotypes read straight off the molecules.
//!
//! Given a window and the variants called inside it, this reports which
//! combinations the reads actually carry, and how many reads carry each. It
//! replaces asking the BAM "does this combination exist?" once per candidate
//! combination, which needs `2^n` questions and, worse, answers yes to
//! combinations no molecule is: a read carrying A, B and C matches the compound
//! allele of A and B as well, because that allele only spans as far as B.
//!
//! Reading the molecules instead has no combinatorial limit, so a window with
//! more variants than the old cap no longer has to fall back to reporting only
//! pairs and the full set.

use crate::error::{AppError, AppResult};
use crate::variants::AlleleComponent;
use noodles::bam;
use noodles::sam::Header;
use std::collections::{BTreeMap, HashMap};

use super::observation::{
    build_molecule_key, build_region, observed_allele_for_ref_span, observed_supports_components,
    MoleculeKey,
};

/// One candidate variant of the window, and the reference span a read has to
/// cover before it can say anything about it.
pub struct LocalVariant {
    pub start: usize,
    pub ref_len: usize,
    pub components: Vec<AlleleComponent>,
}

pub struct LocalHaplotypeRequest<'a> {
    pub chrom: &'a str,
    /// The window to query, spanning every candidate variant.
    pub start: usize,
    pub end: usize,
    /// The candidate variants, in the caller's order. The `carried` vector of
    /// every returned haplotype uses that same order.
    pub variants: &'a [LocalVariant],
    pub min_phred_quality: u8,
    pub min_mapq: u8,
    /// Whether the two mates of a fragment are one molecule. Must match how the
    /// resulting haplotypes are then counted, or discovery and counting
    /// disagree about what a molecule is.
    pub pair_aware: bool,
}

#[derive(Debug, Clone)]
pub struct ObservedHaplotype {
    /// `carried[i]` is true when the reads of this haplotype carry variant `i`.
    pub carried: Vec<bool>,
    pub reads: usize,
}

#[derive(Debug, Clone, Default)]
pub struct LocalHaplotypeObservations {
    /// One entry per distinct combination seen, including the all-reference
    /// one, ordered deterministically.
    pub haplotypes: Vec<ObservedHaplotype>,
    /// Molecules that observed every candidate variant and were assigned a
    /// combination.
    pub spanning_reads: usize,
    /// Molecules that reached the window but left at least one candidate
    /// unobserved, or fell below the base-quality floor there. A molecule seen
    /// in part is not a haplotype: it carries what it carries, but placing it
    /// in a combination would assert it lacks the variants nobody looked at.
    pub partial_reads: usize,
}

pub fn observe_local_haplotypes(
    bam_reader: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: &Header,
    request: LocalHaplotypeRequest<'_>,
) -> AppResult<LocalHaplotypeObservations> {
    let LocalHaplotypeRequest {
        chrom,
        start,
        end,
        variants,
        min_phred_quality,
        min_mapq,
        pair_aware,
    } = request;

    if start == 0 || end < start {
        return Err(AppError::validation(format!(
            "Invalid local haplotype window {chrom}:{start}-{end}"
        )));
    }
    if variants.is_empty() {
        return Ok(LocalHaplotypeObservations::default());
    }

    let region = build_region(chrom, start, end)?;
    let mut query = bam_reader.query(header, &region).map_err(|e| {
        AppError::validation(format!("BAM query failed for {chrom}:{start}-{end}: {e}"))
    })?;

    // Each variant is judged on its own reference span, so a molecule read from
    // both ends can answer for a variant its first mate covers and another its
    // second does. Demanding one contiguous observation across the whole window
    // would throw that away, and would leave discovery disagreeing with the
    // pair-aware counting that follows about what a molecule is.
    struct MoleculeView {
        observed: Vec<bool>,
        carried: Vec<bool>,
        /// Positions where the molecule's two mates read different things. One
        /// of them is wrong and there is no telling which, so the molecule has
        /// observed nothing there.
        conflicted: Vec<bool>,
    }
    let mut views: HashMap<MoleculeKey, MoleculeView> = HashMap::new();
    let mut record_index = 0usize;
    let mut record = bam::Record::default();
    while query
        .read_record(&mut record)
        .map_err(|e| AppError::validation(format!("BAM read error: {e}")))?
        != 0
    {
        record_index += 1;
        let flags = record.flags();
        if flags.is_duplicate() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let mapq = record
            .mapping_quality()
            .map(|q: noodles::sam::alignment::record::MappingQuality| q.get())
            .unwrap_or(255);
        if mapq < min_mapq {
            continue;
        }

        let per_variant = variants
            .iter()
            .map(|variant| {
                observed_allele_for_ref_span(&record, variant.start, variant.ref_len)
                    .filter(|observed| observed.min_quality >= min_phred_quality)
                    .map(|observed| observed_supports_components(&observed, &variant.components))
            })
            .collect::<Vec<_>>();
        if per_variant.iter().all(Option::is_none) {
            continue;
        }

        let view = views
            .entry(build_molecule_key(&record, pair_aware, record_index))
            .or_insert_with(|| MoleculeView {
                observed: vec![false; variants.len()],
                carried: vec![false; variants.len()],
                conflicted: vec![false; variants.len()],
            });
        for (index, support) in per_variant.into_iter().enumerate() {
            if let Some(carries) = support {
                if view.observed[index] && view.carried[index] != carries {
                    view.conflicted[index] = true;
                }
                view.observed[index] = true;
                view.carried[index] |= carries;
            }
        }
    }

    let mut by_combination: BTreeMap<Vec<bool>, usize> = BTreeMap::new();
    let mut partial_reads = 0usize;
    for view in views.into_values() {
        if view.conflicted.iter().any(|clash| *clash) {
            partial_reads += 1;
        } else if view.observed.iter().all(|seen| *seen) {
            *by_combination.entry(view.carried).or_default() += 1;
        } else {
            partial_reads += 1;
        }
    }

    let spanning_reads = by_combination.values().sum();
    let haplotypes = by_combination
        .into_iter()
        .map(|(carried, reads)| ObservedHaplotype { carried, reads })
        .collect();

    Ok(LocalHaplotypeObservations {
        haplotypes,
        spanning_reads,
        partial_reads,
    })
}
