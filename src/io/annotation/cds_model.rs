//! GFF gene records and spliced transcript-CDS model construction.

use crate::error::AppResult;
use crate::variants::{CdsSegment, Gene};
use std::fs::File;
use std::io::{BufRead, BufReader};

use super::gff_parse::{
    gene_name_from_gff, parse_gff_attributes, parse_gff_phase, parse_interval, parse_strand,
};

/// A raw parsed gene record from a GFF line, before any filtering.
pub(crate) struct GffGeneRecord {
    pub(crate) contig: String,
    pub(crate) gene: Gene,
    /// Feature type from GFF column 3 (e.g. "gene", "CDS"). Needed to decide
    /// whether a row participates in the per-transcript CDS aggregation.
    pub(crate) feature_type: String,
    /// Transcript identifier (`transcript_id` attribute preferred, falling
    /// back to `Parent`). Used to group CDS exons belonging to the same
    /// transcript so we can compute full-protein amino-acid positions.
    /// `None` means the row does not participate in aggregation.
    pub(crate) transcript_id: Option<String>,
}

/// Parse a GFF/GFF3 file, yielding one GffGeneRecord per matching feature type.
pub(crate) fn parse_gff_gene_records(
    genes_file: &str,
    feature_types: &[String],
) -> AppResult<Vec<GffGeneRecord>> {
    let file = File::open(genes_file)?;
    parse_gff_gene_records_from_reader(BufReader::new(file), feature_types)
}

/// The parsing itself, over any reader. Split out from the file entry point so
/// annotation can be parsed from memory, which is how test fixtures are built:
/// a fixture that goes through the real parser cannot encode a gene model the
/// parser would never produce.
pub(crate) fn parse_gff_gene_records_from_reader<R: BufRead>(
    reader: R,
    feature_types: &[String],
) -> AppResult<Vec<GffGeneRecord>> {
    let mut records = Vec::new();

    for (line_idx, line) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let entry = line.map_err(|e| {
            format!("Failed to read line {line_number} in GFF/GFF3 annotation file: {e}")
        })?;
        let trimmed = entry.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() != 9 {
            return Err(format!(
                "Invalid GFF/GFF3 entry at line {}: expected 9 tab-separated fields, got {}. Content: '{}'",
                line_number,
                fields.len(),
                entry
            )
            .into());
        }

        if !feature_types.iter().any(|ft| ft == fields[2]) {
            continue;
        }

        let (start, end) = parse_interval(fields[3], fields[4], line_number)?;
        let strand = parse_strand(fields[6], line_number)?;
        let phase = parse_gff_phase(fields[7], line_number)?;
        let attrs = parse_gff_attributes(fields[8]);
        let gene_name = gene_name_from_gff(&attrs);
        let feature_type = fields[2].to_string();
        let transcript_id = attrs
            .get("transcript_id")
            .or_else(|| attrs.get("Parent"))
            .cloned();

        records.push(GffGeneRecord {
            contig: fields[0].to_string(),
            gene: Gene {
                name: gene_name,
                start,
                end,
                strand,
                phase,
                protein_offset: 0,
                transcript_id: transcript_id.clone(),
                cds_segments: Vec::new(),
                // The GFF path selects CDS features, so it is coding by construction.
                biotype: crate::variants::Biotype::ProteinCoding,
            },
            feature_type,
            transcript_id,
        });
    }

    Ok(records)
}

/// Walk all parsed CDS rows of the same transcript and assign each one its
/// cumulative `protein_offset` such that
///
/// ```text
/// reported_aa = protein_offset + local_aa_pos
/// ```
///
/// where `local_aa_pos == 1` corresponds to the first **complete** codon
/// inside the current exon (i.e. the codon that starts at exon position
/// `phase + 1`). Per the GFF3 specification, if `S` is the sum of the
/// lengths of all prior exons of the same transcript and `phase_i` is the
/// phase of exon `i`, then `S + phase_i` is divisible by 3 and the first
/// complete codon of exon `i` is codon
///
/// ```text
/// ((S + phase_i) / 3) + 1
/// ```
///
/// so the offset is `(S + phase_i) / 3`. Counting `(len - phase)/3` per exon
/// would under-count by exactly one for every exon-crossing split codon
/// (every non-zero-phase exon).
///
/// Rows that are not CDS, or that have no transcript identifier, keep
/// `protein_offset = 0` and the historical per-feature numbering is
/// preserved.
pub(crate) fn assign_cds_protein_offsets(records: &mut [GffGeneRecord]) {
    use std::collections::BTreeMap;
    // Group row indices by (contig, transcript_id). Only CDS rows participate.
    let mut groups: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (idx, rec) in records.iter().enumerate() {
        if rec.feature_type != "CDS" {
            continue;
        }
        let Some(tid) = rec.transcript_id.clone() else {
            continue;
        };
        groups
            .entry((rec.contig.clone(), tid))
            .or_default()
            .push(idx);
    }

    for ((_contig, _tid), mut indices) in groups {
        // Sort exons in transcript order: ascending for plus strand, descending
        // for minus strand. We take strand from the first exon (all exons of a
        // transcript share the same strand).
        let strand = records[indices[0]].gene.strand;
        indices.sort_by(|&a, &b| {
            let sa = records[a].gene.start;
            let sb = records[b].gene.start;
            match strand {
                crate::variants::Strand::Plus => sa.cmp(&sb),
                crate::variants::Strand::Minus => sb.cmp(&sa),
            }
        });
        let mut sum_prior_lengths: usize = 0;
        for idx in indices {
            let gene = &mut records[idx].gene;
            let phase = gene.phase as usize;
            gene.protein_offset = (sum_prior_lengths + phase) / 3;
            let len = gene.end.saturating_sub(gene.start) + 1;
            sum_prior_lengths = sum_prior_lengths.saturating_add(len);
        }
    }
}

pub(crate) fn count_multi_exon_cds_transcripts(records: &[GffGeneRecord]) -> usize {
    use std::collections::BTreeMap;

    let mut counts: BTreeMap<(String, String), usize> = BTreeMap::new();
    for rec in records {
        if rec.feature_type != "CDS" {
            continue;
        }
        let Some(transcript_id) = rec.transcript_id.as_ref() else {
            continue;
        };
        *counts
            .entry((rec.contig.clone(), transcript_id.clone()))
            .or_default() += 1;
    }

    counts.values().filter(|count| **count > 1).count()
}

pub(super) fn warn_if_multi_exon_cds_detected(records: &[GffGeneRecord]) {
    let count = count_multi_exon_cds_transcripts(records);
    if count > 0 {
        log::info!(
            "Detected {count} multi-exon CDS transcript(s). get_MNV will build \
             transcript-aware CDS models for CDS rows with transcript identifiers; \
             CDS rows without transcript identifiers keep per-feature annotation."
        );
    }
}

pub(crate) fn build_transcript_cds_records(records: Vec<GffGeneRecord>) -> Vec<GffGeneRecord> {
    use std::collections::BTreeMap;

    let mut groups: BTreeMap<(String, String), Vec<GffGeneRecord>> = BTreeMap::new();
    let mut out = Vec::new();

    for rec in records {
        if rec.feature_type == "CDS" {
            if let Some(transcript_id) = rec.transcript_id.clone() {
                groups
                    .entry((rec.contig.clone(), transcript_id))
                    .or_default()
                    .push(rec);
                continue;
            }
        }
        out.push(rec);
    }

    let mut spliced_transcripts = 0usize;
    let mut joined_transcripts = 0usize;

    for ((contig, transcript_id), mut group) in groups {
        let strand = group[0].gene.strand;

        // Segments that do not share an orientation are not one transcript. The
        // sort below, and everything downstream that reads `cds_segments` in
        // transcript order, would otherwise be built on the first row's strand
        // and silently produce a nonsense reading frame.
        if group.iter().any(|rec| rec.gene.strand != strand) {
            log::warn!(
                "CDS transcript '{transcript_id}' on contig '{contig}' has rows on both \
                 strands; keeping per-feature annotation because a spliced CDS cannot be \
                 built from segments that disagree about orientation."
            );
            out.extend(group);
            continue;
        }

        group.sort_by(|a, b| match strand {
            crate::variants::Strand::Plus => a.gene.start.cmp(&b.gene.start),
            crate::variants::Strand::Minus => b.gene.start.cmp(&a.gene.start),
        });

        // A CDS row listed twice (concatenated annotation files, or a GFF that
        // repeats a segment) would count those bases twice in the spliced CDS
        // and shift every downstream codon. Sorting puts identical rows next to
        // each other, so dropping the repeats here is enough.
        let before_dedup = group.len();
        group.dedup_by(|a, b| a.gene.start == b.gene.start && a.gene.end == b.gene.end);
        if group.len() != before_dedup {
            log::warn!(
                "CDS transcript '{transcript_id}' on contig '{contig}' lists {} duplicate \
                 CDS row(s); ignoring the repeats, which would otherwise be counted twice \
                 in the spliced coding sequence.",
                before_dedup - group.len()
            );
        }

        if group.len() == 1 && group[0].gene.phase != 0 {
            log::warn!(
                "CDS transcript '{transcript_id}' on contig '{contig}' has a single \
                 CDS row with non-zero phase {}; keeping per-feature annotation because \
                 the full spliced CDS cannot be reconstructed from the selected rows.",
                group[0].gene.phase
            );
            out.extend(group);
            continue;
        }

        if group.first().is_some_and(|rec| rec.gene.phase != 0) {
            log::warn!(
                "CDS transcript '{transcript_id}' on contig '{contig}' starts with \
                 non-zero phase {}; keeping per-feature annotation because the selected \
                 CDS rows do not appear to contain the full coding sequence.",
                group[0].gene.phase
            );
            out.extend(group);
            continue;
        }

        let start = group
            .iter()
            .map(|rec| rec.gene.start)
            .min()
            .unwrap_or_default();
        let end = group
            .iter()
            .map(|rec| rec.gene.end)
            .max()
            .unwrap_or_default();
        let mut gene = group[0].gene.clone();
        gene.start = start;
        gene.end = end;
        gene.phase = 0;
        gene.protein_offset = 0;
        gene.transcript_id = Some(transcript_id.clone());
        gene.cds_segments = group
            .iter()
            .map(|rec| CdsSegment {
                start: rec.gene.start,
                end: rec.gene.end,
            })
            .collect();

        // Record what the model actually came out as, so the summary below can
        // tell the user which transcripts get_MNV read as spliced and which as
        // one continuous reading frame. The two are annotated very differently.
        if gene.cds_segments.len() > 1 {
            let intron_count = gene.introns().count();
            if intron_count == 0 {
                joined_transcripts += 1;
            } else {
                spliced_transcripts += 1;
            }
        }

        out.push(GffGeneRecord {
            contig,
            gene,
            feature_type: "CDS".to_string(),
            transcript_id: Some(transcript_id),
        });
    }

    if spliced_transcripts > 0 || joined_transcripts > 0 {
        log::info!(
            "Built {spliced_transcripts} spliced CDS model(s) and {joined_transcripts} \
             multi-segment model(s) whose segments abut or overlap. The latter are one \
             continuous reading frame (a ribosomal-slippage join such as SARS-CoV-2 \
             ORF1ab), so they carry no splice sites and no exon-exon junction."
        );
    }

    out
}
