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

    let observation_span = crate::variants::observation_ref_len(ref_allele, required_components);
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

    let mut record_index = 0usize;
    let mut record = bam::Record::default();
    while query
        .read_record(&mut record)
        .map_err(|e| AppError::validation(format!("BAM read error: {e}")))?
        != 0
    {
        record_index += 1;
        let flags = record.flags();
        if flags.is_duplicate()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_qc_fail()
        {
            continue;
        }
        let mapq = record
            .mapping_quality()
            .map(|q: noodles::sam::alignment::record::MappingQuality| q.get())
            .unwrap_or(255);
        if mapq < min_mapq {
            continue;
        }

        // The allele is read over the REF span, because exact support means
        // reproducing that allele and nothing else.
        let observed = observed_allele_for_ref_span(&record, position, ref_allele.len());
        // Whether the read is entitled to an opinion at all. For an insertion
        // that needs the base after the anchor too: a read stopping on the
        // anchor covers the REF span and has still seen nothing of the
        // junction, so its "reads the reference" observation is not a
        // disagreement and must not strip a partner mate's support.
        let extended = if observation_span == ref_allele.len() {
            None
        } else {
            observed_allele_for_ref_span(&record, position, observation_span)
        };
        let saw_the_whole_allele = observation_span == ref_allele.len() || extended.is_some();
        // The base just past the REF span still deleted means this read carries
        // a longer deletion than the allele being counted. Inside the span the
        // two are indistinguishable, which is how a read deleting 4-5 came to
        // be exact support for a deletion of 4 alone.
        let deletion_runs_on = extended.as_ref().is_some_and(|seen| {
            seen.deleted_positions
                .contains(&(position + ref_allele.len()))
        });

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

        let key = Rc::new(build_molecule_key(&record, pair_aware, record_index));
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
                    && !deletion_runs_on
                {
                    unique_alt.insert(key.clone());
                    if is_reverse {
                        unique_alt_reverse.insert(key);
                    } else {
                        unique_alt_forward.insert(key);
                    }
                } else if saw_the_whole_allele {
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
        // A single allele has no pair, so no combination to distribute.
        haplotype_patterns: Vec::new(),
    })
}
