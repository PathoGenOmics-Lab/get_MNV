//! Per-contig parallel processing: variant annotation, BAM read counting,
//! LRU cache, and output emission.
//!
//! Input parsing and per-gene read support live in submodules; this file owns
//! the worker state types and the per-contig orchestration (`process_contig`)
//! plus output emission.

use super::summary::{summarize_contig_variants, ContigSummary};
use crate::cli::Args;
use crate::error::{AppError, AppResult, ErrorCode};
use crate::io::{self, ReferenceMap, VcfPosition};
use crate::output;
use crate::read_count::{self, RegionObservationCache};
use crate::variants::{self, Gene, VariantInfo, VariantType};
use log::info;
use lru::LruCache;
use noodles::bam;
use noodles::sam::Header;
use rayon::prelude::*;
use std::collections::HashMap;
use std::num::NonZeroUsize;

mod inputs;
mod read_support;

pub(crate) use inputs::{parse_inputs, resolve_variant_input_format};
use read_support::{
    annotate_declared_phase, append_supported_phased_indel_haplotypes, apply_read_summary,
    compute_frameshift_phasing, count_gene_variant_reads, count_indel_reads_into,
};

#[derive(Debug)]
pub(crate) struct ParsedInputs {
    pub(crate) base_name: String,
    pub(crate) references: ReferenceMap,
    pub(crate) snp_by_contig: HashMap<String, Vec<VcfPosition>>,
    pub(crate) contigs: Vec<String>,
    pub(crate) command_line: String,
    pub(crate) preloaded_gff: Option<HashMap<String, Vec<Gene>>>,
    pub(crate) original_info_headers: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct RegionCacheKey {
    contig: String,
    start: usize,
    end: usize,
    positions: Vec<usize>,
    min_mapq: u8,
}

struct WorkerState {
    bam: Option<bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>>,
    bam_header: Option<Header>,
    region_cache: LruCache<RegionCacheKey, RegionObservationCache>,
}

#[derive(Default)]
struct WorkerResult {
    variants: Vec<VariantInfo>,
    region_cache_hits: usize,
    region_cache_misses: usize,
}

struct PhasedIndelInputs<'a, 'r> {
    contig: &'a str,
    gene: &'a Gene,
    reference: &'a io::Reference<'r>,
    snp_list: &'a [VcfPosition],
    genetic_code: crate::genetic_code::GeneticCode,
}

pub(crate) fn sort_variants(variants: &mut [VariantInfo]) {
    variants.sort_by(|a, b| {
        let a_min = a.positions.iter().copied().min().unwrap_or(usize::MAX);
        let b_min = b.positions.iter().copied().min().unwrap_or(usize::MAX);
        a.chrom
            .cmp(&b.chrom)
            .then_with(|| a_min.cmp(&b_min))
            .then_with(|| a.positions.cmp(&b.positions))
            .then_with(|| a.variant_type.cmp(&b.variant_type))
            .then_with(|| a.gene.cmp(&b.gene))
    });
}

pub(crate) fn reclassify_generic_as_validation(error: AppError) -> AppError {
    if error.code == ErrorCode::Generic {
        AppError::validation(error.message)
    } else {
        error
    }
}

fn annotate_variants_for_gene(
    gene: &Gene,
    snp_list: &[VcfPosition],
    reference: &io::Reference,
    contig: &str,
    genetic_code: crate::genetic_code::GeneticCode,
    indel_config: &variants::IndelAnnotationConfig,
    phasing: &variants::FrameshiftPhasing,
) -> Vec<VariantInfo> {
    // A non-coding transcript has no reading frame, so translating it would
    // invent amino-acid changes for something that is never translated. Report
    // its variants against the gene as non_coding_transcript_exon_variant.
    if !gene.biotype.is_coding() {
        return snp_list
            .iter()
            .filter(|snp| io::gene_overlaps_variant(gene, snp))
            .map(|snp| {
                variants::build_non_coding_variant(
                    contig,
                    snp,
                    &gene.name,
                    variants::NonCodingReason::NonCodingTranscript,
                )
            })
            .collect();
    }
    variants::get_mnv_variants_for_gene_with_config(
        gene,
        snp_list,
        reference,
        contig,
        genetic_code,
        indel_config,
        phasing,
    )
}

/// Count BAM read support for intergenic variants. These fall outside every
/// gene, so the per-gene counting loop never reaches them. Intergenic SNPs are
/// counted so they can be filtered by their real support like any other SNP;
/// intergenic indels are left uncounted (as elsewhere).
///
/// Nearby intergenic SNPs are grouped into bounded windows and counted from a
/// single region observation cache per window, so clustered sites share one BAM
/// query instead of one query each. The window span and position count are
/// capped so memory stays bounded when sites are dense; distant sites simply
/// fall into separate windows.
/// First gene whose intron places this variant in a splice region, with the
/// splice consequence. The variant's anchor position is tested; only multi-exon
/// transcript models have internal junctions, so single-feature (prokaryotic)
/// annotations never match.
fn splice_site_for_variant(
    genes: &[Gene],
    positions: &[usize],
) -> Option<(String, crate::variants::SpliceConsequence)> {
    genes.iter().find_map(|gene| {
        positions
            .iter()
            .filter_map(|&position| {
                variants::splice::splice_consequence_for_position(gene, position)
            })
            .max_by_key(|consequence| consequence.severity())
            .map(|consequence| (gene.name.clone(), consequence))
    })
}

/// First gene with a real intron containing this variant. Only reached when the
/// variant is not in any CDS and not at a splice site, so it describes a variant
/// sitting deep inside an intron.
fn intron_for_variant(genes: &[Gene], positions: &[usize]) -> Option<String> {
    genes
        .iter()
        .find(|gene| {
            positions
                .iter()
                .any(|&position| variants::splice::is_intronic_position(gene, position))
        })
        .map(|gene| gene.name.clone())
}

/// How a single base that no gene annotated is accounted for.
#[derive(Debug, Clone, PartialEq, Eq)]
enum FallbackClass {
    Splice(String, crate::variants::SpliceConsequence),
    Intron(String),
    NonCoding(String),
    Intergenic,
}

/// Whether the gene's coding span covers this single base.
fn gene_contains_position(gene: &Gene, position: usize) -> bool {
    if gene.cds_segments.is_empty() {
        return position >= gene.start && position <= gene.end;
    }
    gene.cds_segments
        .iter()
        .any(|segment| position >= segment.start && position <= segment.end)
}

/// Drop the entries of a fallback row whose change some gene already annotated.
///
/// The fallback runs for a record no gene fully covered, and it used to rebuild
/// the whole record. A record straddling a gene's edge, `74 AAC>GGT` on a gene
/// ending at 75, then produced a second row repeating 74 and 75 with a
/// contradictory "Unknown" call and claiming base 76 for a gene that ends before
/// it. Only the parts nobody annotated belong in this row.
fn retain_unannotated_entries(
    variant: &mut VariantInfo,
    annotated: &std::collections::HashSet<&str>,
) -> bool {
    retain_entries(variant, |position, reference, alternate| {
        !annotated.contains(format!("SNV:{position}:{reference}>{alternate}").as_str())
    })
}

/// Keep only the entries of a row the predicate accepts, and say whether any
/// survived.
fn retain_entries(
    variant: &mut VariantInfo,
    keep_entry: impl Fn(usize, &str, &str) -> bool,
) -> bool {
    if variant.positions.len() != variant.ref_bases.len()
        || variant.positions.len() != variant.base_changes.len()
        || variant.positions.len() < 2
    {
        return true;
    }
    let keep: Vec<bool> = (0..variant.positions.len())
        .map(|index| {
            keep_entry(
                variant.positions[index],
                &variant.ref_bases[index],
                &variant.base_changes[index],
            )
        })
        .collect();
    if keep.iter().all(|kept| *kept) {
        return true;
    }
    if keep.iter().all(|kept| !*kept) {
        return false;
    }
    let mut index = 0;
    variant.positions.retain(|_| {
        index += 1;
        keep[index - 1]
    });
    let mut index = 0;
    variant.ref_bases.retain(|_| {
        index += 1;
        keep[index - 1]
    });
    let mut index = 0;
    variant.base_changes.retain(|_| {
        index += 1;
        keep[index - 1]
    });
    if variant.positions.len() == 1 {
        variant.variant_type = VariantType::Snp;
    }
    let kept: std::collections::BTreeSet<usize> = variant.positions.iter().copied().collect();
    variant.event_components.retain(|label| {
        label
            .split(':')
            .nth(1)
            .and_then(|value| value.split('-').next())
            .and_then(|value| value.parse::<usize>().ok())
            .is_none_or(|position| kept.contains(&position))
    });
    variant.snp_aa_changes = vec!["-".to_string(); variant.positions.len()];
    variant.snp_aa_changes_local = vec!["-".to_string(); variant.positions.len()];
    true
}

fn count_intergenic_variant_reads(
    args: &Args,
    contig: &str,
    variants: &mut [VariantInfo],
) -> AppResult<()> {
    const MAX_WINDOW_SPAN: usize = 50_000;
    const MAX_WINDOW_POSITIONS: usize = 256;

    let bam_path = match args.bam_file.as_ref() {
        Some(path) => path,
        None => return Ok(()),
    };

    // Where each substitution row's change actually happens, which is not always
    // where its record starts. A left-padded record such as `350 AA>AG` changes
    // base 351, and counting at 350 with the raw ALT took the padding base as
    // the alternate: that base is the reference, so every read scored as support
    // and the row reported full frequency for a change most reads do not carry.
    // The gene path never had this because its positions come from the
    // decomposed components; these rows now use the same source.
    let substitution_sites: Vec<(usize, Vec<usize>, Vec<String>)> = variants
        .iter()
        .enumerate()
        .filter(|(_, v)| matches!(v.variant_type, VariantType::Snp | VariantType::Mnv))
        .filter_map(|(i, v)| {
            let components: Vec<_> = v
                .event_components
                .iter()
                .filter_map(|label| variants::parse_component_label(label))
                .filter(|c| c.kind == variants::AlleleComponentKind::Snp)
                .collect();
            // A multi-base substitution outside a CDS used to match neither this
            // selector nor the indel one, so it reached no counter at all and a
            // threshold then read its absent counts as zero and deleted the row.
            // Guard on the shape the counts will be attached to rather than on
            // how many bases changed.
            if components.is_empty() || components.len() != v.positions.len() {
                return None;
            }
            Some((
                i,
                components.iter().map(|c| c.position).collect(),
                components.iter().map(|c| c.alt_allele.clone()).collect(),
            ))
        })
        .collect();
    let mut snp_indices: Vec<usize> = substitution_sites.iter().map(|(i, _, _)| *i).collect();
    let site_of: std::collections::HashMap<usize, (Vec<usize>, Vec<String>)> = substitution_sites
        .into_iter()
        .map(|(i, positions, alts)| (i, (positions, alts)))
        .collect();
    // Indels outside every feature reached no counter at all, so the row went
    // out claiming `Event Reads = 0` at `Event Depth = 0` for an allele every
    // read might carry. They are counted here, one at a time, because the
    // window cache only understands single-position substitutions.
    let indel_indices: Vec<usize> = variants
        .iter()
        .enumerate()
        .filter(|(_, v)| v.variant_type == VariantType::Indel && !v.positions.is_empty())
        .map(|(i, _)| i)
        .collect();
    if snp_indices.is_empty() && indel_indices.is_empty() {
        return Ok(());
    }
    snp_indices.sort_by_key(|i| site_of[i].0[0]);

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| {
            AppError::validation(format!(
                "Failed to open BAM '{bam_path}' for intergenic read counting: {e}"
            ))
        })?;
    let header = reader.read_header().map_err(|e| {
        AppError::validation(format!("Failed to read BAM header from '{bam_path}': {e}"))
    })?;

    let mut start = 0usize;
    while start < snp_indices.len() {
        let window_start_pos = site_of[&snp_indices[start]].0[0];
        let mut end = start + 1;
        while end < snp_indices.len()
            && end - start < MAX_WINDOW_POSITIONS
            && site_of[&snp_indices[end]].0[0].saturating_sub(window_start_pos) <= MAX_WINDOW_SPAN
        {
            end += 1;
        }

        let window = &snp_indices[start..end];
        // Every position the window's rows will be counted at, not just the one
        // each row starts on: a multi-base substitution needs all of its bases
        // in the cache or the count fails outright.
        let mut positions: Vec<usize> = window
            .iter()
            .flat_map(|i| site_of[i].0.iter().copied())
            .collect();
        positions.sort_unstable();
        positions.dedup();
        let region_end_pos = positions.last().copied().unwrap_or(window_start_pos);

        let cache = read_count::build_region_observation_cache(
            &mut reader,
            &header,
            contig,
            window_start_pos,
            region_end_pos,
            &positions,
            args.min_mapq,
        )
        .map_err(|e| {
            AppError::validation(format!(
                "Failed building read cache for intergenic window {contig}:{window_start_pos}-{region_end_pos}: {e}"
            ))
        })?;

        for &i in window {
            let (site_positions, site_alts) = &site_of[&i];
            let summary = read_count::count_reads_from_cache(
                &cache,
                site_positions,
                site_alts,
                args.min_quality,
                !args.count_mates_separately,
            )
            .map_err(|e| {
                AppError::validation(format!(
                    "Failed counting reads for intergenic position {contig}:{}: {e}",
                    variants[i].positions[0]
                ))
            })?;
            apply_read_summary(&mut variants[i], summary);
        }

        start = end;
    }

    for &i in &indel_indices {
        count_indel_reads_into(
            &mut reader,
            &header,
            args,
            contig,
            "outside any annotated feature",
            &mut variants[i],
        )?;
    }
    Ok(())
}

pub(crate) fn process_contig(
    args: &Args,
    contig: &str,
    references: &ReferenceMap,
    snp_by_contig: &HashMap<String, Vec<VcfPosition>>,
    preloaded_gff: Option<&HashMap<String, Vec<Gene>>>,
    genetic_code: crate::genetic_code::GeneticCode,
) -> AppResult<(Vec<VariantInfo>, ContigSummary)> {
    let reference =
        io::reference_for_chrom(references, contig).map_err(reclassify_generic_as_validation)?;
    let snp_list = snp_by_contig
        .get(contig)
        .ok_or_else(|| AppError::validation(format!("Missing VCF data for contig '{contig}'")))?;
    io::validate_vcf_reference_alleles(contig, snp_list, &reference)
        .map_err(reclassify_generic_as_validation)?;

    // An insertion written against the base after it anchors before its own POS,
    // which every observer here reads forward from. Re-anchor those records once
    // so the rest of the pipeline sees the spelling it can count.
    let (snp_list, re_anchored) = io::left_anchor_insertions(snp_list, &reference);
    let (snp_list, reference_calls) = io::drop_reference_calls(snp_list);
    let snp_list = &snp_list;
    if reference_calls > 0 {
        info!(
            "Contig '{contig}' -> {reference_calls} record(s) whose ALT repeats REF describe no \
             change and were skipped"
        );
    }
    if re_anchored > 0 {
        info!(
            "Contig '{contig}' -> {re_anchored} insertion(s) re-anchored onto the base they \
             follow, so their read support can be counted"
        );
    }

    let genes = if let Some(gff_genes) = preloaded_gff {
        let all_contig_genes = gff_genes.get(contig).cloned().unwrap_or_default();
        let filtered = io::filter_genes_with_snps(&all_contig_genes, snp_list);
        log::info!(
            "GFF/GFF3 contig '{}': {} gene entries, {} mapped to SNPs, {} without SNPs",
            contig,
            all_contig_genes.len(),
            filtered.len(),
            all_contig_genes.len() - filtered.len()
        );
        filtered
    } else {
        io::load_genes(
            args.genes_file(),
            snp_list,
            Some(contig),
            &args.gff_features(),
        )
        .map_err(reclassify_generic_as_validation)?
    };
    info!(
        "Contig '{}' -> {} SNP/variant records in VCF, {} mapped genes",
        contig,
        snp_list.len(),
        genes.len()
    );

    let bam_path = args.bam_file.clone();
    let should_count_reads = args.bam_file.is_some() && !args.dry_run;
    let worker_results: Result<Vec<WorkerResult>, AppError> = genes
        .par_iter()
        .map_init(
            || -> AppResult<WorkerState> {
                let (bam, bam_header) = if should_count_reads {
                    if let Some(path) = bam_path.as_ref() {
                        let mut reader = bam::io::indexed_reader::Builder::default()
                            .build_from_path(path)
                            .map_err(|e| {
                                AppError::validation(format!(
                                    "Failed to open BAM '{path}' in worker thread: {e}"
                                ))
                            })?;
                        let header = reader.read_header().map_err(|e| {
                            AppError::validation(format!(
                                "Failed to read BAM header from '{path}': {e}"
                            ))
                        })?;
                        (Some(reader), Some(header))
                    } else {
                        (None, None)
                    }
                } else {
                    (None, None)
                };
                Ok(WorkerState {
                    bam,
                    bam_header,
                    region_cache: LruCache::new(NonZeroUsize::new(64).expect("64 > 0")),
                })
            },
            |state_result, gene| -> AppResult<WorkerResult> {
                let state = state_result
                    .as_mut()
                    .map_err(|err| AppError::validation(err.to_string()))?;
                // BAM-derived cis/trans phasing so an upstream indel does not
                // frameshift a downstream codon it is not on the same molecule
                // as, and the frequency the reads give each upstream indel,
                // which is what the frameshift gate weighs when a BAM was given.
                let evidence = compute_frameshift_phasing(state, args, contig, gene, snp_list)?;
                let indel_config = variants::IndelAnnotationConfig {
                    frameshift_min_freq: args.frameshift_min_freq,
                    observed_indel_freq: evidence.observed_indel_freq,
                };
                let phasing = evidence.phasing;
                let mut variants = annotate_variants_for_gene(
                    gene,
                    snp_list,
                    &reference,
                    contig,
                    genetic_code,
                    &indel_config,
                    &phasing,
                );
                if variants.is_empty() {
                    return Ok(WorkerResult::default());
                }
                let (cache_hits, cache_misses) =
                    count_gene_variant_reads(state, args, contig, gene, &mut variants)?;
                append_supported_phased_indel_haplotypes(
                    state,
                    args,
                    PhasedIndelInputs {
                        contig,
                        gene,
                        reference: &reference,
                        snp_list,
                        genetic_code,
                    },
                    &mut variants,
                )?;
                annotate_declared_phase(&mut variants, snp_list);
                Ok(WorkerResult {
                    variants,
                    region_cache_hits: cache_hits,
                    region_cache_misses: cache_misses,
                })
            },
        )
        .collect();

    let mut all_variants = Vec::new();
    let mut cache_hits = 0usize;
    let mut cache_misses = 0usize;
    for result in worker_results.map_err(|e| {
        AppError::validation(format!("Error while processing contig '{contig}': {e}"))
    })? {
        cache_hits += result.region_cache_hits;
        cache_misses += result.region_cache_misses;
        all_variants.extend(result.variants);
    }

    // Anything the gene path did not turn into a row.
    //
    // Overlapping a gene is not the same as having been annotated by it: a
    // substitution in a trailing partial codon (a feature whose length is not a
    // multiple of three) is dropped there, and treating "a gene overlaps this"
    // as "this was annotated" removed such a variant from every output with no
    // warning at all. Ask the rows that were actually produced instead.
    {
        // Keyed by the event each row describes, not by the coordinate it sits
        // at. Two records can share a position, and a position-keyed test then
        // reads "some row exists here" as "this record was annotated": a
        // substitution at 1524 made an insertion anchored on 1524 vanish from
        // every output. The component label identifies the change itself.
        let annotated: std::collections::HashSet<&str> = all_variants
            .iter()
            .flat_map(|variant| variant.event_components.iter().map(String::as_str))
            .collect();
        // A record was annotated when every part of the change it describes
        // appears among the rows that were produced. Asking about its POS alone
        // got this wrong twice: a left-padded record such as `28 GC>GA` is
        // annotated at 29, and an insertion re-anchored onto a coordinate a
        // substitution already occupies is not the substitution.
        let produced_a_row = |snp: &crate::io::VcfPosition| {
            let labels = snp.event().component_labels();
            // Nothing to annotate means nothing is missing.
            labels.is_empty()
                || labels
                    .iter()
                    .all(|label| annotated.contains(label.as_str()))
        };
        let mut intergenic: Vec<VariantInfo> = Vec::new();
        for snp in snp_list.iter() {
            if !produced_a_row(snp) {
                // A variant not in any CDS may still sit inside a gene: at a
                // splice site, or deeper in the intron. Both are annotated
                // against their gene rather than called intergenic, which a
                // position inside a transcript is not. Splice sites are checked
                // first because they are the more specific consequence. These
                // stay in this vector so their read support is counted with the
                // intergenic SNPs.
                // The bases of this record nobody annotated. When that is the
                // whole record the tests below see all of it, as before; when a
                // gene covered part of it, they see only the rest, so a record
                // straddling a gene's edge no longer claims that gene for the
                // bases lying outside it.
                let unannotated: Vec<usize> = {
                    let mut positions = Vec::new();
                    for component in snp.event().components {
                        if annotated.contains(component.label().as_str()) {
                            continue;
                        }
                        let end = component.position + component.ref_allele.len().saturating_sub(1);
                        match component.kind {
                            crate::variants::AlleleComponentKind::Snp
                            | crate::variants::AlleleComponentKind::Insertion => {
                                positions.push(component.position)
                            }
                            _ => positions.extend(component.position..=end),
                        }
                    }
                    positions.sort_unstable();
                    positions.dedup();
                    if positions.is_empty() {
                        snp.changed_positions()
                    } else {
                        positions
                    }
                };
                // A record of substitutions can straddle a boundary, and each of
                // its bases then belongs somewhere different: one inside the gene,
                // the next past its end. Rebuilding the whole record under a
                // single answer put a base outside the gene into that gene's row.
                // Each base is accounted for on its own and the row is built per
                // group. An indel is one indivisible change, so it keeps the
                // single-row path below.
                let all_substitutions =
                    snp.event().components.iter().all(|component| {
                        component.kind == crate::variants::AlleleComponentKind::Snp
                    });
                if all_substitutions && !unannotated.is_empty() {
                    let mut groups: Vec<(FallbackClass, std::collections::BTreeSet<usize>)> =
                        Vec::new();
                    for &position in &unannotated {
                        let class = if let Some((gene_name, splice)) =
                            splice_site_for_variant(&genes, &[position])
                        {
                            FallbackClass::Splice(gene_name, splice)
                        } else if let Some(gene_name) = intron_for_variant(&genes, &[position]) {
                            FallbackClass::Intron(gene_name)
                        } else if let Some(gene) = genes
                            .iter()
                            .find(|gene| gene_contains_position(gene, position))
                        {
                            FallbackClass::NonCoding(gene.name.clone())
                        } else {
                            FallbackClass::Intergenic
                        };
                        match groups.iter_mut().find(|(existing, _)| *existing == class) {
                            Some((_, members)) => {
                                members.insert(position);
                            }
                            None => {
                                groups.push((class, std::collections::BTreeSet::from([position])))
                            }
                        }
                    }
                    for (class, members) in groups {
                        if matches!(class, FallbackClass::Intergenic) && args.exclude_intergenic {
                            continue;
                        }
                        let mut variant = match &class {
                            FallbackClass::Splice(gene_name, splice) => {
                                variants::build_splice_variant(contig, snp, gene_name, *splice)
                            }
                            FallbackClass::Intron(gene_name) => {
                                variants::build_intron_variant(contig, snp, gene_name)
                            }
                            FallbackClass::NonCoding(gene_name) => {
                                variants::build_non_coding_variant(
                                    contig,
                                    snp,
                                    gene_name,
                                    variants::NonCodingReason::IncompleteCodon,
                                )
                            }
                            FallbackClass::Intergenic => {
                                variants::build_intergenic_variant(contig, snp)
                            }
                        };
                        if !retain_entries(&mut variant, |position, _, _| {
                            members.contains(&position)
                        }) {
                            continue;
                        }
                        intergenic.push(variant);
                    }
                    continue;
                }

                let fully_unannotated = unannotated.len() == snp.changed_positions().len();
                let variant = match splice_site_for_variant(&genes, &unannotated) {
                    Some((gene_name, splice)) => {
                        variants::build_splice_variant(contig, snp, &gene_name, splice)
                    }
                    None => match intron_for_variant(&genes, &unannotated) {
                        Some(gene_name) => variants::build_intron_variant(contig, snp, &gene_name),
                        None => match genes.iter().find(|gene| {
                            if fully_unannotated {
                                io::gene_overlaps_variant(gene, snp)
                            } else {
                                unannotated
                                    .iter()
                                    .any(|&position| gene_contains_position(gene, position))
                            }
                        }) {
                            // Inside a gene, but the gene path could not build a
                            // codon around it. Keep the row, name the gene, and
                            // say the amino-acid effect is unknown.
                            Some(gene) => variants::build_non_coding_variant(
                                contig,
                                snp,
                                &gene.name,
                                variants::NonCodingReason::IncompleteCodon,
                            ),
                            None => {
                                // Genuinely outside every feature: this is the
                                // only kind `--exclude-intergenic` removes.
                                if args.exclude_intergenic {
                                    continue;
                                }
                                variants::build_intergenic_variant(contig, snp)
                            }
                        },
                    },
                };
                let mut variant = variant;
                if !retain_unannotated_entries(&mut variant, &annotated) {
                    continue;
                }
                intergenic.push(variant);
            }
        }
        // A base can be coding in one gene and a splice site of another that
        // overlaps it. The fallback above only runs for a record no gene
        // annotated at all, so once any gene emitted a codon row the splice
        // consequence of every other gene at that base was lost: a HIGH
        // splice_donor_variant disappeared because an unrelated overlapping gene
        // happened to cover the base. Each gene answers for itself.
        {
            let named: std::collections::HashSet<(&str, &str)> = all_variants
                .iter()
                .chain(intergenic.iter())
                .flat_map(|variant| {
                    variant
                        .event_components
                        .iter()
                        .map(move |label| (variant.gene.as_str(), label.as_str()))
                })
                .collect();
            let mut missed = Vec::new();
            for snp in snp_list.iter() {
                let positions = snp.changed_positions();
                let labels = snp.event().component_labels();
                for gene in &genes {
                    if labels
                        .iter()
                        .any(|label| named.contains(&(gene.name.as_str(), label.as_str())))
                    {
                        continue;
                    }
                    // Whatever this gene has to say about these bases: a splice
                    // site first, since it is the more specific consequence, and
                    // otherwise its intron. Asking only about splice sites left a
                    // gene whose intron holds the base with no row at all as soon
                    // as another gene annotated it, which is the same silent loss
                    // one consequence further in.
                    let splice = positions
                        .iter()
                        .filter_map(|&position| {
                            variants::splice::splice_consequence_for_position(gene, position)
                        })
                        .max_by_key(|consequence| consequence.severity());
                    let inside: std::collections::BTreeSet<usize> = positions
                        .iter()
                        .copied()
                        .filter(|&position| {
                            variants::splice::splice_consequence_for_position(gene, position)
                                .is_some()
                                || variants::splice::is_intronic_position(gene, position)
                        })
                        .collect();
                    if inside.is_empty() {
                        continue;
                    }
                    let mut variant = match splice {
                        Some(consequence) => {
                            variants::build_splice_variant(contig, snp, &gene.name, consequence)
                        }
                        None => variants::build_intron_variant(contig, snp, &gene.name),
                    };
                    if !retain_entries(&mut variant, |position, _, _| inside.contains(&position)) {
                        continue;
                    }
                    missed.push(variant);
                }
            }
            intergenic.extend(missed);
        }

        if !intergenic.is_empty() {
            let intergenic_count = intergenic.len();
            // These rows carry a declared phase like any other. The claim was
            // only read inside the per-gene fold, so an MNP with the same
            // `GT`/`PS` reported cis inside a gene and nothing at all one base
            // outside it: the same input claim, kept or lost by which code path
            // happened to build the row.
            annotate_declared_phase(&mut intergenic, snp_list);
            if should_count_reads {
                count_intergenic_variant_reads(args, contig, &mut intergenic)?;
            }
            all_variants.extend(intergenic);
            info!("Contig '{contig}' -> {intergenic_count} intergenic variant(s) added");
        }
    }

    sort_variants(&mut all_variants);
    let summary = summarize_contig_variants(
        contig,
        snp_list.len(),
        genes.len(),
        &all_variants,
        cache_hits,
        cache_misses,
    );
    Ok((all_variants, summary))
}

pub(crate) fn emit_contig_variants(
    tsv_writer: Option<&mut output::TsvWriter>,
    vcf_writer: Option<&mut output::VcfWriter>,
    variants: &[VariantInfo],
    references: &ReferenceMap,
) -> AppResult<()> {
    if let Some(writer) = tsv_writer {
        writer.write_variants(variants)?;
    }
    if let Some(writer) = vcf_writer {
        writer.write_variants(variants, references)?;
    }
    Ok(())
}
