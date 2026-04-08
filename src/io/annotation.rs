//! Gene annotation parsing: GFF/GFF3 and TSV formats, format detection, and
//! gene-to-SNP interval filtering.

use super::vcf::VcfPosition;
use crate::error::AppResult;
use crate::variants::Gene;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnnotationFormat {
    Tsv,
    Gff,
}

pub fn detect_annotation_format(genes_file: &str) -> AppResult<AnnotationFormat> {
    if let Some(ext) = Path::new(genes_file).extension().and_then(|e| e.to_str()) {
        if ext.eq_ignore_ascii_case("gff") || ext.eq_ignore_ascii_case("gff3") {
            return Ok(AnnotationFormat::Gff);
        }
        if ext.eq_ignore_ascii_case("txt") || ext.eq_ignore_ascii_case("tsv") {
            return Ok(AnnotationFormat::Tsv);
        }
    }

    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let entry = line?;
        let trimmed = entry.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() >= 9 {
            return Ok(AnnotationFormat::Gff);
        }
        return Ok(AnnotationFormat::Tsv);
    }

    Ok(AnnotationFormat::Tsv)
}

pub(crate) fn decode_percent_encoded(value: &str) -> String {
    let mut buf: Vec<u8> = Vec::with_capacity(value.len());
    let mut bytes = value.as_bytes().iter().copied();
    while let Some(current) = bytes.next() {
        if current == b'%' {
            let hi = bytes.next();
            let lo = bytes.next();
            if let (Some(hi), Some(lo)) = (hi, lo) {
                let hex_pair = [hi, lo];
                if let Ok(hex) = std::str::from_utf8(&hex_pair) {
                    if let Ok(decoded) = u8::from_str_radix(hex, 16) {
                        buf.push(decoded);
                        continue;
                    }
                }
                buf.extend_from_slice(&[b'%', hi, lo]);
                continue;
            }
            buf.push(b'%');
            if let Some(hi) = hi {
                buf.push(hi);
            }
            if let Some(lo) = lo {
                buf.push(lo);
            }
            continue;
        }
        buf.push(current);
    }
    String::from_utf8(buf).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned())
}

pub(crate) fn split_attribute_fields(attributes: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut start = 0usize;
    let mut in_quotes = false;
    let mut escaped = false;
    for (idx, ch) in attributes.char_indices() {
        match ch {
            '"' if !escaped => in_quotes = !in_quotes,
            ';' if !in_quotes => {
                let field = attributes[start..idx].trim();
                if !field.is_empty() {
                    fields.push(field);
                }
                start = idx + 1;
            }
            _ => {}
        }
        escaped = ch == '\\' && !escaped;
    }
    let trailing = attributes[start..].trim();
    if !trailing.is_empty() {
        fields.push(trailing);
    }
    fields
}

pub(crate) fn parse_gff_attributes(attributes: &str) -> HashMap<String, String> {
    let mut parsed = HashMap::new();
    for field in split_attribute_fields(attributes) {
        let (key_raw, value_raw) = if let Some((key, value)) = field.split_once('=') {
            (key.trim(), value.trim())
        } else if let Some((key, value)) = field.split_once(' ') {
            (key.trim(), value.trim().trim_matches('"'))
        } else {
            continue;
        };
        if key_raw.is_empty() || value_raw.is_empty() {
            continue;
        }
        parsed.insert(
            key_raw.to_string(),
            decode_percent_encoded(value_raw.trim_matches('"').trim()),
        );
    }
    parsed
}

/// Requires `snp_list` to be sorted by position (guaranteed by `load_vcf_positions_by_contig`).
pub(crate) fn has_snp_in_interval(snp_list: &[VcfPosition], start: usize, end: usize) -> bool {
    let idx = snp_list.partition_point(|snp| snp.position < start);
    idx < snp_list.len() && snp_list[idx].position <= end
}

fn parse_strand(raw: &str, line_number: usize) -> AppResult<crate::variants::Strand> {
    raw.parse::<crate::variants::Strand>().map_err(|_| {
        format!("Invalid strand at line {line_number} ('{raw}'). Expected '+' or '-'").into()
    })
}

/// Parse the GFF phase column (field 8, 0-indexed 7).
///
/// Per GFF3 spec, valid values are `0`, `1`, `2` (required for CDS features)
/// or `.` when not applicable (e.g. gene, exon, UTR). Any `.` or empty value
/// is normalised to 0 so features that do not carry phase information keep the
/// historical behaviour of the tool.
fn parse_gff_phase(raw: &str, line_number: usize) -> AppResult<u8> {
    match raw.trim() {
        "." | "" => Ok(0),
        "0" => Ok(0),
        "1" => Ok(1),
        "2" => Ok(2),
        other => Err(format!(
            "Invalid GFF phase at line {line_number} ('{other}'). Expected '0', '1', '2' or '.'"
        )
        .into()),
    }
}

fn parse_interval(start_raw: &str, end_raw: &str, line_number: usize) -> AppResult<(usize, usize)> {
    let start = start_raw.parse::<usize>().map_err(|e| {
        format!("Invalid start coordinate at line {line_number} ('{start_raw}'): {e}")
    })?;
    let end = end_raw
        .parse::<usize>()
        .map_err(|e| format!("Invalid end coordinate at line {line_number} ('{end_raw}'): {e}"))?;

    if start > end {
        return Err(format!(
            "Invalid gene interval at line {line_number}: start ({start}) is greater than end ({end})"
        )
        .into());
    }
    if start == 0 || end == 0 {
        return Err(format!(
            "Invalid gene interval at line {line_number}: coordinates must be 1-based positive integers (start={start}, end={end})"
        )
        .into());
    }

    Ok((start, end))
}

pub(crate) fn gene_name_from_gff(attrs: &HashMap<String, String>) -> String {
    // Prefer human-readable names common in eukaryotic GTF/GFF3
    // (gene_name, gene) before falling back to Name/locus_tag/ID.
    let primary = attrs
        .get("gene_name")
        .or_else(|| attrs.get("gene"))
        .or_else(|| attrs.get("Name"))
        .or_else(|| attrs.get("locus_tag"))
        .or_else(|| attrs.get("gene_id"))
        .or_else(|| attrs.get("ID"))
        .map(|value| value.trim_start_matches("gene-").to_string())
        .unwrap_or_else(|| "unknown_gene".to_string());

    // Only append a locus_tag suffix when it is actually present and
    // different from the primary name. This avoids the historical
    // `primary_primary` duplication when no locus_tag exists.
    match attrs.get("locus_tag") {
        Some(locus) if locus != &primary => format!("{primary}_{locus}"),
        _ => primary,
    }
}

fn load_genes_from_tsv(genes_file: &str, snp_list: &[VcfPosition]) -> AppResult<Vec<Gene>> {
    log::info!("Loading TSV gene file: {genes_file}");
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    let mut genes: Vec<Gene> = Vec::new();
    let mut parsed_entries = 0usize;
    let mut genes_with_snps = 0usize;
    let mut genes_without_snps = 0usize;

    for (line_idx, line) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let entry =
            line.map_err(|e| format!("Failed to read line {line_number} in genes file: {e}"))?;
        let trimmed = entry.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        parsed_entries += 1;

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 4 {
            return Err(format!(
                "Invalid genes entry at line {}: expected 4 tab-separated fields, got {}. Content: '{}'",
                line_number,
                fields.len(),
                entry
            ).into());
        }

        let (start, end) = parse_interval(fields[1], fields[2], line_number)?;
        let strand = parse_strand(fields[3], line_number)?;
        // Optional 5th column: phase (0|1|2|.). Defaults to 0 (prokaryote-style)
        // when omitted, preserving the historical 4-column TSV format.
        let phase = if fields.len() >= 5 {
            parse_gff_phase(fields[4], line_number)?
        } else {
            0
        };

        if has_snp_in_interval(snp_list, start, end) {
            genes_with_snps += 1;
            genes.push(crate::variants::Gene {
                name: fields[0].to_string(),
                start,
                end,
                strand,
                phase,
                protein_offset: 0,
            });
        } else {
            genes_without_snps += 1;
        }
    }

    log::info!(
        "TSV gene entries parsed: {parsed_entries} | mapped to SNPs: {genes_with_snps} | without SNPs: {genes_without_snps}"
    );

    Ok(genes)
}

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
    let reader = BufReader::new(file);
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

/// Scan a GFF file once and report whether it contains any CDS row with a
/// non-zero phase. Used to warn users that they should pass
/// `--gff-features CDS` when working with eukaryotic annotations.
pub(crate) fn gff_has_non_zero_phase_cds(genes_file: &str) -> AppResult<bool> {
    let file = File::open(genes_file)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let entry = line?;
        let trimmed = entry.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() != 9 {
            continue;
        }
        if fields[2] != "CDS" {
            continue;
        }
        match fields[7].trim() {
            "1" | "2" => return Ok(true),
            _ => {}
        }
    }
    Ok(false)
}

pub(crate) fn load_genes_from_gff(
    genes_file: &str,
    snp_list: &[VcfPosition],
    expected_chrom: Option<&str>,
    feature_types: &[String],
) -> AppResult<Vec<Gene>> {
    log::info!("Loading GFF/GFF3 annotation file: {genes_file}");
    warn_if_cds_phase_ignored(genes_file, feature_types)?;
    let mut records = parse_gff_gene_records(genes_file, feature_types)?;
    assign_cds_protein_offsets(&mut records);

    let mut genes: Vec<Gene> = Vec::new();
    let mut parsed_entries = 0usize;
    let mut genes_with_snps = 0usize;
    let mut genes_without_snps = 0usize;
    let mut genes_other_contigs = 0usize;
    let mut other_contigs: HashSet<String> = HashSet::new();

    for rec in records {
        if let Some(chrom) = expected_chrom {
            if rec.contig != chrom {
                genes_other_contigs += 1;
                other_contigs.insert(rec.contig);
                continue;
            }
        }

        parsed_entries += 1;

        if has_snp_in_interval(snp_list, rec.gene.start, rec.gene.end) {
            genes_with_snps += 1;
            genes.push(rec.gene);
        } else {
            genes_without_snps += 1;
        }
    }

    log::info!(
        "GFF/GFF3 gene entries parsed: {} | mapped to SNPs: {} | without SNPs: {} | skipped other contigs: {}{}",
        parsed_entries,
        genes_with_snps,
        genes_without_snps,
        genes_other_contigs,
        if other_contigs.is_empty() {
            String::new()
        } else {
            let mut values = other_contigs.into_iter().collect::<Vec<_>>();
            values.sort();
            format!(" ({})", values.join(", "))
        }
    );

    Ok(genes)
}

pub fn load_genes(
    genes_file: &str,
    snp_list: &[VcfPosition],
    expected_chrom: Option<&str>,
    feature_types: &[String],
) -> AppResult<Vec<Gene>> {
    match detect_annotation_format(genes_file)? {
        AnnotationFormat::Tsv => load_genes_from_tsv(genes_file, snp_list),
        AnnotationFormat::Gff => {
            load_genes_from_gff(genes_file, snp_list, expected_chrom, feature_types)
        }
    }
}

/// Pre-load all gene entries from a GFF/GFF3 file, grouped by contig.
/// No SNP filtering is applied; use `filter_genes_with_snps` afterwards.
pub fn preload_gff_genes(
    genes_file: &str,
    feature_types: &[String],
) -> AppResult<HashMap<String, Vec<Gene>>> {
    log::info!("Pre-loading GFF/GFF3 annotation file: {genes_file}");
    warn_if_cds_phase_ignored(genes_file, feature_types)?;
    let mut records = parse_gff_gene_records(genes_file, feature_types)?;
    assign_cds_protein_offsets(&mut records);
    let total_entries = records.len();

    let mut genes_by_contig: HashMap<String, Vec<Gene>> = HashMap::new();
    for rec in records {
        genes_by_contig
            .entry(rec.contig)
            .or_default()
            .push(rec.gene);
    }

    let contig_count = genes_by_contig.len();
    log::info!("GFF/GFF3 pre-loaded: {total_entries} gene entries across {contig_count} contigs");
    Ok(genes_by_contig)
}

/// Emit a single warning line if the GFF contains CDS rows with non-zero
/// phase but the user did not select `CDS` in `--gff-features`. In that
/// configuration the phase-aware codon grouping introduced in 1.1.2 is
/// effectively disabled and codon boundaries may be wrong for eukaryotic
/// annotations.
fn warn_if_cds_phase_ignored(genes_file: &str, feature_types: &[String]) -> AppResult<()> {
    let selects_cds = feature_types.iter().any(|f| f == "CDS");
    if selects_cds {
        return Ok(());
    }
    if gff_has_non_zero_phase_cds(genes_file)? {
        log::warn!(
            "GFF/GFF3 file '{genes_file}' contains CDS rows with non-zero phase, \
             but --gff-features does not include 'CDS' (currently: {}). \
             Codon grouping will ignore the phase information and may be incorrect \
             for eukaryotic annotations. Re-run with '--gff-features CDS' to use \
             phase-aware codon boundaries.",
            feature_types.join(",")
        );
    }
    Ok(())
}

/// Filter genes that overlap with at least one SNP position.
/// Requires `snp_list` to be sorted by position.
pub fn filter_genes_with_snps(genes: &[Gene], snp_list: &[VcfPosition]) -> Vec<Gene> {
    genes
        .iter()
        .filter(|gene| has_snp_in_interval(snp_list, gene.start, gene.end))
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- decode_percent_encoded ----

    #[test]
    fn test_decode_percent_plain() {
        assert_eq!(decode_percent_encoded("hello"), "hello");
    }

    #[test]
    fn test_decode_percent_space() {
        assert_eq!(decode_percent_encoded("hello%20world"), "hello world");
    }

    #[test]
    fn test_decode_percent_multiple() {
        assert_eq!(decode_percent_encoded("%2C%3B"), ",;");
    }

    #[test]
    fn test_decode_percent_trailing_incomplete() {
        // Incomplete percent sequence at end preserved literally
        assert_eq!(decode_percent_encoded("a%2"), "a%2");
        assert_eq!(decode_percent_encoded("a%"), "a%");
    }

    #[test]
    fn test_decode_percent_invalid_hex() {
        // Invalid hex digits preserved literally
        assert_eq!(decode_percent_encoded("%ZZ"), "%ZZ");
    }

    // ---- split_attribute_fields ----

    #[test]
    fn test_split_simple() {
        let fields = split_attribute_fields("ID=gene1;Name=rpoB");
        assert_eq!(fields, vec!["ID=gene1", "Name=rpoB"]);
    }

    #[test]
    fn test_split_quoted_semicolon() {
        // Semicolons inside quotes should not split
        let fields = split_attribute_fields("Note=\"a;b\";ID=gene1");
        assert_eq!(fields, vec!["Note=\"a;b\"", "ID=gene1"]);
    }

    #[test]
    fn test_split_trailing_semicolon() {
        let fields = split_attribute_fields("ID=gene1;");
        assert_eq!(fields, vec!["ID=gene1"]);
    }

    #[test]
    fn test_split_empty() {
        let fields = split_attribute_fields("");
        assert!(fields.is_empty());
    }

    // ---- parse_gff_attributes ----

    #[test]
    fn test_parse_gff_key_value() {
        let attrs = parse_gff_attributes("ID=gene1;Name=rpoB;locus_tag=Rv0667");
        assert_eq!(attrs.get("ID").unwrap(), "gene1");
        assert_eq!(attrs.get("Name").unwrap(), "rpoB");
        assert_eq!(attrs.get("locus_tag").unwrap(), "Rv0667");
    }

    #[test]
    fn test_parse_gff_space_separated() {
        // GTF-style: key "value"
        let attrs = parse_gff_attributes("gene_id \"ENSG001\"; gene_name \"TP53\"");
        assert_eq!(attrs.get("gene_id").unwrap(), "ENSG001");
        assert_eq!(attrs.get("gene_name").unwrap(), "TP53");
    }

    #[test]
    fn test_parse_gff_percent_encoded() {
        let attrs = parse_gff_attributes("Name=rpoB%20gene");
        assert_eq!(attrs.get("Name").unwrap(), "rpoB gene");
    }

    // ---- has_snp_in_interval ----

    #[test]
    fn test_has_snp_in_interval_found() {
        let snps = vec![
            VcfPosition { position: 100, ref_allele: "A".into(), alt_allele: "T".into(), original_dp: None, original_freq: None, original_info: None },
            VcfPosition { position: 200, ref_allele: "G".into(), alt_allele: "C".into(), original_dp: None, original_freq: None, original_info: None },
        ];
        assert!(has_snp_in_interval(&snps, 50, 150));
        assert!(has_snp_in_interval(&snps, 100, 100)); // exact match
        assert!(has_snp_in_interval(&snps, 150, 250));
    }

    #[test]
    fn test_has_snp_in_interval_not_found() {
        let snps = vec![
            VcfPosition { position: 100, ref_allele: "A".into(), alt_allele: "T".into(), original_dp: None, original_freq: None, original_info: None },
        ];
        assert!(!has_snp_in_interval(&snps, 200, 300));
        assert!(!has_snp_in_interval(&snps, 50, 99));
    }

    #[test]
    fn test_has_snp_in_interval_empty() {
        assert!(!has_snp_in_interval(&[], 1, 1000));
    }

    // ---- gene_name_from_gff ----

    #[test]
    fn test_gene_name_from_gff_full() {
        let mut attrs = HashMap::new();
        attrs.insert("gene".to_string(), "rpoB".to_string());
        attrs.insert("locus_tag".to_string(), "Rv0667".to_string());
        assert_eq!(gene_name_from_gff(&attrs), "rpoB_Rv0667");
    }

    #[test]
    fn test_gene_name_from_gff_no_gene_uses_name() {
        let mut attrs = HashMap::new();
        attrs.insert("Name".to_string(), "gene-katG".to_string());
        // "gene-" prefix should be stripped; no locus_tag → no duplicated suffix
        assert_eq!(gene_name_from_gff(&attrs), "katG");
    }

    #[test]
    fn test_gene_name_from_gff_empty() {
        let attrs = HashMap::new();
        assert_eq!(gene_name_from_gff(&attrs), "unknown_gene");
    }

    #[test]
    fn test_gene_name_from_gff_prefers_gene_name_over_id() {
        // Regression: eukaryotic GTF/GFF3 rows often carry gene_name (human
        // readable) plus an ID coming from annotation tools like agat. We
        // should prefer gene_name over ID so the TSV shows GNAQ instead of
        // agat-cds-37838.
        let mut attrs = HashMap::new();
        attrs.insert("ID".to_string(), "agat-cds-37838".to_string());
        attrs.insert("gene_id".to_string(), "ENSG00000156052.11".to_string());
        attrs.insert("gene_name".to_string(), "GNAQ".to_string());
        assert_eq!(gene_name_from_gff(&attrs), "GNAQ");
    }

    #[test]
    fn test_gene_name_from_gff_uses_gene_id_when_no_name() {
        let mut attrs = HashMap::new();
        attrs.insert("ID".to_string(), "agat-cds-1".to_string());
        attrs.insert("gene_id".to_string(), "ENSG00000156052.11".to_string());
        assert_eq!(gene_name_from_gff(&attrs), "ENSG00000156052.11");
    }

    #[test]
    fn test_gene_name_from_gff_locus_tag_equal_to_primary_not_duplicated() {
        // When locus_tag is the ONLY name available, primary and locus_tag
        // are equal so we must not emit `primary_primary`.
        let mut attrs = HashMap::new();
        attrs.insert("locus_tag".to_string(), "Rv0007".to_string());
        assert_eq!(gene_name_from_gff(&attrs), "Rv0007");
    }

    // ---- filter_genes_with_snps ----

    // ---- assign_cds_protein_offsets / multi-exon aggregation ----

    fn record(contig: &str, start: usize, end: usize, strand: crate::variants::Strand, phase: u8, ft: &str, tid: Option<&str>) -> GffGeneRecord {
        GffGeneRecord {
            contig: contig.to_string(),
            gene: Gene {
                name: "x".into(),
                start,
                end,
                strand,
                phase,
                protein_offset: 0,
            },
            feature_type: ft.to_string(),
            transcript_id: tid.map(|s| s.to_string()),
        }
    }

    #[test]
    fn test_assign_protein_offsets_gnaq_minus_strand() {
        use crate::variants::Strand;
        // Seven CDS rows of ENST00000286548.9 (GNAQ, minus strand) copied from
        // the minimal example. On minus strand, the first exon of the
        // transcript is the one with the HIGHEST genomic coordinate.
        let mut recs = vec![
            // (genomic order, not transcript order)
            record("chr9", 77_721_323, 77_721_513, Strand::Minus, 2, "CDS", Some("ENST00000286548.9")), // exon 7
            record("chr9", 77_728_514, 77_728_667, Strand::Minus, 0, "CDS", Some("ENST00000286548.9")), // exon 6
            record("chr9", 77_794_463, 77_794_592, Strand::Minus, 1, "CDS", Some("ENST00000286548.9")), // exon 5
            record("chr9", 77_797_520, 77_797_648, Strand::Minus, 1, "CDS", Some("ENST00000286548.9")), // exon 4
            record("chr9", 77_815_616, 77_815_770, Strand::Minus, 0, "CDS", Some("ENST00000286548.9")), // exon 3
            record("chr9", 77_922_161, 77_922_345, Strand::Minus, 2, "CDS", Some("ENST00000286548.9")), // exon 2
            record("chr9", 78_031_100, 78_031_235, Strand::Minus, 0, "CDS", Some("ENST00000286548.9")), // exon 1
        ];
        assign_cds_protein_offsets(&mut recs);
        // Build a map (start -> offset) so we don't depend on ordering.
        let by_start: std::collections::HashMap<usize, usize> = recs
            .iter()
            .map(|r| (r.gene.start, r.gene.protein_offset))
            .collect();
        // Expected cumulative codon offset before each exon, derived from the
        // GFF spec formula `(sum_prior_lengths + phase_i) / 3`. Hand-checked
        // against the canonical Ensembl GNAQ-201 protein numbering: the
        // famous oncogenic mutation Q209L (chr9:77794572) must end up at
        // codon 209, which only happens with these offsets.
        //   exon 1: S=0,    phase=0  → offset = 0      (codons 1..45)
        //   exon 2: S=136,  phase=2  → offset = 46     (codons 47..107; 46 is split 1↔2)
        //   exon 3: S=321,  phase=0  → offset = 107    (codons 108..158)
        //   exon 4: S=476,  phase=1  → offset = 159    (codons 160..201; 159 is split 3↔4)
        //   exon 5: S=605,  phase=1  → offset = 202    (codons 203..245; 202 is split 4↔5)
        //   exon 6: S=735,  phase=0  → offset = 245    (codons 246..296)
        //   exon 7: S=889,  phase=2  → offset = 297    (codons 298..360; 297 is split 6↔7)
        assert_eq!(by_start[&78_031_100], 0, "exon 1");
        assert_eq!(by_start[&77_922_161], 46, "exon 2");
        assert_eq!(by_start[&77_815_616], 107, "exon 3");
        assert_eq!(by_start[&77_797_520], 159, "exon 4");
        assert_eq!(by_start[&77_794_463], 202, "exon 5");
        assert_eq!(by_start[&77_728_514], 245, "exon 6");
        assert_eq!(by_start[&77_721_323], 297, "exon 7");
    }

    #[test]
    fn test_assign_protein_offsets_gnaq_q209l_hotspot() {
        use crate::variants::Strand;
        // Hard-coded sanity check: GNAQ Q209L canonical mutation lives at
        // chr9:77794572 (middle base of codon 209). Verify the chain
        //   exon-5 effective end → local codon index → +protein_offset
        // produces codon 209 exactly.
        let mut recs = vec![
            record("chr9", 78_031_100, 78_031_235, Strand::Minus, 0, "CDS", Some("ENST00000286548.9")),
            record("chr9", 77_922_161, 77_922_345, Strand::Minus, 2, "CDS", Some("ENST00000286548.9")),
            record("chr9", 77_815_616, 77_815_770, Strand::Minus, 0, "CDS", Some("ENST00000286548.9")),
            record("chr9", 77_797_520, 77_797_648, Strand::Minus, 1, "CDS", Some("ENST00000286548.9")),
            record("chr9", 77_794_463, 77_794_592, Strand::Minus, 1, "CDS", Some("ENST00000286548.9")),
            record("chr9", 77_728_514, 77_728_667, Strand::Minus, 0, "CDS", Some("ENST00000286548.9")),
            record("chr9", 77_721_323, 77_721_513, Strand::Minus, 2, "CDS", Some("ENST00000286548.9")),
        ];
        assign_cds_protein_offsets(&mut recs);
        let exon5 = recs.iter().find(|r| r.gene.start == 77_794_463).unwrap();
        // The exon-5 effective end (after subtracting phase=1) is 77794591.
        // For position 77794572: offset = 77794591 - 77794572 = 19,
        // local codon index = 19/3 = 6, local_aa_pos = 7,
        // total = protein_offset(202) + 7 = 209.
        let pos = 77_794_572usize;
        let eff_end = exon5.gene.end - exon5.gene.phase as usize;
        let local_aa = (eff_end - pos) / 3 + 1;
        let aa = exon5.gene.protein_offset + local_aa;
        assert_eq!(aa, 209, "GNAQ Q209L canonical position must map to codon 209");
    }

    #[test]
    fn test_assign_protein_offsets_plus_strand_two_exons() {
        use crate::variants::Strand;
        // Plus strand, two exons. exon1 length 12, phase 0 → 4 complete codons
        // contributed; exon2 starts at protein_offset=4.
        let mut recs = vec![
            record("chrX", 100, 111, Strand::Plus, 0, "CDS", Some("T1")),
            record("chrX", 200, 217, Strand::Plus, 0, "CDS", Some("T1")),
        ];
        assign_cds_protein_offsets(&mut recs);
        let by_start: std::collections::HashMap<usize, usize> = recs
            .iter()
            .map(|r| (r.gene.start, r.gene.protein_offset))
            .collect();
        assert_eq!(by_start[&100], 0);
        assert_eq!(by_start[&200], 4);
    }

    #[test]
    fn test_assign_protein_offsets_non_cds_untouched() {
        use crate::variants::Strand;
        let mut recs = vec![
            record("chr1", 10, 30, Strand::Plus, 0, "gene", None),
            record("chr1", 10, 30, Strand::Plus, 0, "exon", Some("T1")),
        ];
        assign_cds_protein_offsets(&mut recs);
        for rec in &recs {
            assert_eq!(rec.gene.protein_offset, 0);
        }
    }

    // ---- filter_genes_with_snps ----

    #[test]
    fn test_filter_genes_with_snps() {
        let genes = vec![
            Gene { name: "gene1".into(), start: 100, end: 200, strand: crate::variants::Strand::Plus, phase: 0, protein_offset: 0 },
            Gene { name: "gene2".into(), start: 500, end: 600, strand: crate::variants::Strand::Minus, phase: 0, protein_offset: 0 },
        ];
        let snps = vec![
            VcfPosition { position: 150, ref_allele: "A".into(), alt_allele: "T".into(), original_dp: None, original_freq: None, original_info: None },
        ];
        let filtered = filter_genes_with_snps(&genes, &snps);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].name, "gene1");
    }
}
