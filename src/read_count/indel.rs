//! Direct indel / complex-allele read counting from a BAM region query.

use crate::error::{AppError, AppResult};
use noodles::bam;
use noodles::sam::Header;
use std::collections::HashSet;
use std::rc::Rc;

use super::observation::{
    anchor_base_quality, build_molecule_key, build_region, observed_allele_for_ref_span,
    observed_supports_components, MoleculeKey,
};
use super::{IndelReadCountRequest, ReadCountSummary};

pub fn count_indel_reads(
    bam_reader: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: &Header,
    request: IndelReadCountRequest<'_>,
) -> AppResult<ReadCountSummary> {
    let IndelReadCountRequest {
        chrom,
        position,
        ref_allele,
        alt_allele,
        required_components,
        min_phred_quality,
        min_mapq,
        anchor_depth,
        pair_aware,
    } = request;

    if position == 0 || ref_allele.is_empty() {
        return Err(AppError::validation(format!(
            "Invalid indel allele for read counting at {chrom}:{position} REF='{ref_allele}' ALT='{alt_allele}'"
        )));
    }

    let end = position + ref_allele.len().saturating_sub(1);
    let region = build_region(chrom, position, end)?;
    let mut query = bam_reader.query(header, &region).map_err(|e| {
        AppError::validation(format!(
            "BAM query failed for {chrom}:{position}-{end}: {e}"
        ))
    })?;

    let mut unique_total: HashSet<Rc<MoleculeKey>> = HashSet::new();
    let mut unique_total_forward: HashSet<Rc<MoleculeKey>> = HashSet::new();
    let mut unique_total_reverse: HashSet<Rc<MoleculeKey>> = HashSet::new();
    let mut unique_alt: HashSet<Rc<MoleculeKey>> = HashSet::new();
    let mut unique_alt_forward: HashSet<Rc<MoleculeKey>> = HashSet::new();
    let mut unique_alt_reverse: HashSet<Rc<MoleculeKey>> = HashSet::new();
    // Molecules where one mate spans the locus and reads a different allele
    // from the one that supported the ALT. One of the two is wrong and there is
    // no telling which, so the molecule does not get to support the call; the
    // substitution counting takes the same view of a contradicted overlap.
    // Without this, a single mate carrying an alignment artefact was enough to
    // claim the whole molecule, and mate overlap is exactly where indel
    // realignment artefacts appear.
    let mut contradicted: HashSet<Rc<MoleculeKey>> = HashSet::new();

    let mut record = bam::Record::default();
    while query
        .read_record(&mut record)
        .map_err(|e| AppError::validation(format!("BAM read error: {e}")))?
        != 0
    {
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

        let observed = observed_allele_for_ref_span(&record, position, ref_allele.len());

        // Depth (denominator) eligibility. By default a read must fully span the
        // REF allele with every base either observed or deleted. With
        // `anchor_depth` we instead count any read that observes the anchor base
        // at sufficient quality, so reads that only partially overlap a
        // multi-base deletion still contribute to the locus depth (and thus to a
        // less biased EFREQ).
        let counts_for_depth = if anchor_depth {
            anchor_base_quality(&record, position).is_some_and(|q| q >= min_phred_quality)
        } else {
            observed
                .as_ref()
                .is_some_and(|o| o.min_quality >= min_phred_quality)
        };

        if !counts_for_depth {
            // In the default mode this also means the read cannot support the
            // ALT allele (no full-span observation), so skipping is safe. In
            // anchor-depth mode a read that fails the anchor-quality gate is not
            // counted toward depth, and ALT support below is gated on its own
            // full-span observation.
            if !anchor_depth {
                continue;
            }
        }

        let key = Rc::new(build_molecule_key(&record, pair_aware));
        let is_reverse = flags.is_reverse_complemented();
        if counts_for_depth {
            unique_total.insert(key.clone());
            if is_reverse {
                unique_total_reverse.insert(key.clone());
            } else {
                unique_total_forward.insert(key.clone());
            }
        }

        if let Some(observed) = observed.as_ref() {
            if observed.min_quality >= min_phred_quality {
                if observed.allele.eq_ignore_ascii_case(alt_allele)
                    && observed_supports_components(observed, required_components)
                {
                    unique_alt.insert(key.clone());
                    if is_reverse {
                        unique_alt_reverse.insert(key);
                    } else {
                        unique_alt_forward.insert(key);
                    }
                } else {
                    contradicted.insert(key);
                }
            }
        }
    }

    for key in &contradicted {
        unique_alt.remove(key);
        unique_alt_forward.remove(key);
        unique_alt_reverse.remove(key);
    }

    let alt_count = unique_alt.len();
    let alt_forward = unique_alt_forward.len();
    let alt_reverse = unique_alt_reverse.len();
    let total = unique_total.len();
    let total_forward = unique_total_forward.len();
    let total_reverse = unique_total_reverse.len();

    Ok(ReadCountSummary {
        snp_counts: vec![alt_count],
        mnv_count: alt_count,
        total_reads: vec![total],
        total_forward_reads: vec![total_forward],
        total_reverse_reads: vec![total_reverse],
        snp_forward_counts: vec![alt_forward],
        snp_reverse_counts: vec![alt_reverse],
        mnv_forward_count: alt_forward,
        mnv_reverse_count: alt_reverse,
        mnv_total_reads: total,
        mnv_total_forward_reads: total_forward,
        mnv_total_reverse_reads: total_reverse,
        // An indel row is a single allele, not a multi-position haplotype, so
        // there is no linkage question to answer here.
        snp_only_informative_counts: vec![0],
    })
}
