//! Self-contained interactive HTML report.
//!
//! The report is always built from emitted get_MNV TSV files, so the three entry
//! points share one code path: a normal run (`--report`), a `--sample all` run
//! (one report covering every sample), and standalone aggregation of TSVs from
//! separate runs (`--report-from`, typical for per-sample pipelines).
//!
//! Output is a single HTML file with the data embedded as compact JSON and no
//! external requests, so it opens offline by double-clicking. Repeated fields
//! (sample, contig, gene, variant type, change type, SO term, impact) are
//! dictionary-encoded to keep the file small for cohort-sized inputs.

use crate::error::{AppError, AppResult};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

const TEMPLATE: &str = include_str!("report_template.html");
const DATA_PLACEHOLDER: &str = "__GET_MNV_REPORT_DATA__";

/// String interner backing the report's dictionary encoding.
#[derive(Default)]
struct Interner {
    values: Vec<String>,
    index: HashMap<String, usize>,
}

impl Interner {
    fn intern(&mut self, value: &str) -> usize {
        if let Some(existing) = self.index.get(value) {
            return *existing;
        }
        let id = self.values.len();
        self.values.push(value.to_string());
        self.index.insert(value.to_string(), id);
        id
    }

    fn ensure(&mut self, values: &[&str]) {
        for value in values {
            self.intern(value);
        }
    }
}

#[derive(Default)]
struct Row {
    sample: usize,
    chrom: usize,
    gene: usize,
    positions: String,
    ref_bases: String,
    alt_bases: String,
    aa_change: String,
    variant_type: usize,
    change_type: usize,
    so_term: usize,
    impact: usize,
    grantham: Option<u32>,
    shift: Option<String>,
    dbs: Option<String>,
    nmd: Option<String>,
    hgvs_g: Option<String>,
    hgvs_c: Option<String>,
    freq: Option<f64>,
    depth: Option<u64>,
    ref_codon: Option<String>,
    alt_codon: Option<String>,
    event_class: Option<String>,
    event_components: Option<String>,
    reads: Option<u64>,
    phasing: Option<f64>,
    phasing_reads: Option<u64>,
}

#[derive(Default)]
pub struct ReportBuilder {
    samples: Interner,
    chroms: Interner,
    genes: Interner,
    types: Interner,
    changes: Interner,
    so_terms: Interner,
    impacts: Interner,
    rows: Vec<Row>,
}

/// A TSV placeholder meaning "not applicable"; treated as an absent value.
fn optional(value: &str) -> Option<String> {
    let trimmed = value.trim();
    if trimmed.is_empty() || trimmed == "-" {
        None
    } else {
        Some(trimmed.to_string())
    }
}

/// Leading integer of a cell such as `177 (radical)`.
fn leading_u32(value: &str) -> Option<u32> {
    let digits: String = value
        .trim()
        .chars()
        .take_while(char::is_ascii_digit)
        .collect();
    digits.parse().ok()
}

/// First entry of a comma-separated numeric cell (per-position columns hold one
/// value per constituent SNV).
fn first_number<T: std::str::FromStr>(value: &str) -> Option<T> {
    value
        .split(',')
        .next()
        .map(str::trim)
        .filter(|part| !part.is_empty() && *part != "-")
        .and_then(|part| part.parse().ok())
}

/// Sample label for a TSV path: the file stem with the `.MNV` suffix removed,
/// which is how the pipeline names per-sample outputs.
fn sample_label(path: &str) -> String {
    let stem = Path::new(path)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string());
    let stem = stem.strip_suffix(".gz").unwrap_or(&stem);
    let stem = stem.strip_suffix(".tsv").unwrap_or(stem);
    let stem = stem.strip_suffix(".MNV").unwrap_or(stem);
    stem.to_string()
}

/// Whether a file begins with the gzip magic number (which also covers BGZF).
/// A file that cannot be opened or read is left to the reader, which reports the
/// I/O error with its own path context.
fn is_gzip_compressed(path: &str) -> bool {
    use std::io::Read;
    let mut magic = [0u8; 2];
    std::fs::File::open(path)
        .and_then(|mut file| file.read_exact(&mut magic))
        .is_ok()
        && magic == [0x1f, 0x8b]
}

/// Label candidates for a TSV path, from the bare file stem to progressively
/// more qualified forms that prepend one more parent directory at a time, ending
/// at the path itself.
fn label_candidates(path: &str) -> Vec<String> {
    let stem = sample_label(path);
    let mut candidates = vec![stem.clone()];
    let parents: Vec<String> = Path::new(path)
        .parent()
        .map(|parent| {
            parent
                .components()
                .filter_map(|component| match component {
                    std::path::Component::Normal(name) => Some(name.to_string_lossy().into_owned()),
                    _ => None,
                })
                .collect()
        })
        .unwrap_or_default();
    for take in 1..=parents.len() {
        let prefix = parents[parents.len() - take..].join("/");
        candidates.push(format!("{prefix}/{stem}"));
    }
    candidates
}

/// Distinct sample labels for a set of TSV paths, in input order.
///
/// The label is normally just the file stem, but per-sample pipelines routinely
/// write identically named files into per-sample directories
/// (`results/<sample>/variants.MNV.tsv`), so a bare stem is not unique in
/// general. Colliding paths are qualified with as many trailing directory
/// components as it takes to tell them apart; paths that are still
/// indistinguishable (the same file listed twice) get a numeric suffix. Without
/// this the report interned one shared label and silently merged distinct
/// samples into a single row of the matrix.
pub fn sample_labels(paths: &[String]) -> Vec<String> {
    let ladders: Vec<Vec<String>> = paths.iter().map(|path| label_candidates(path)).collect();
    let mut levels = vec![0usize; paths.len()];

    // Qualify every member of a colliding group until the collision clears or no
    // member has any more parent directories to add.
    loop {
        let mut groups: HashMap<&str, Vec<usize>> = HashMap::new();
        for (idx, level) in levels.iter().enumerate() {
            groups
                .entry(ladders[idx][*level].as_str())
                .or_default()
                .push(idx);
        }
        let colliding: Vec<usize> = groups
            .into_values()
            .filter(|members| members.len() > 1)
            .flatten()
            .collect();
        let mut progressed = false;
        for idx in colliding {
            if levels[idx] + 1 < ladders[idx].len() {
                levels[idx] += 1;
                progressed = true;
            }
        }
        if !progressed {
            break;
        }
    }

    // Anything still sharing a label is the same path more than once.
    let mut used: HashMap<String, usize> = HashMap::new();
    levels
        .iter()
        .enumerate()
        .map(|(idx, level)| {
            let base = ladders[idx][*level].clone();
            let seen = used.entry(base.clone()).or_insert(0);
            *seen += 1;
            if *seen == 1 {
                base
            } else {
                format!("{base} #{seen}")
            }
        })
        .collect()
}

impl ReportBuilder {
    pub fn new() -> Self {
        let mut builder = Self::default();
        // Seed the fixed vocabularies so the report's chips and legends keep a
        // stable order regardless of which values a given run happens to emit.
        builder.types.ensure(&["SNP", "MNV", "SNP/MNV", "INDEL"]);
        builder
            .impacts
            .ensure(&["HIGH", "MODERATE", "LOW", "MODIFIER"]);
        builder
    }

    /// Number of variant rows collected so far.
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Read one get_MNV TSV, labelling its rows with `sample` (defaults to the
    /// file stem when `None`).
    pub fn add_tsv(&mut self, path: &str, sample: Option<&str>) -> AppResult<()> {
        // get_MNV writes plain TSVs, but a compressed one is an easy mistake to
        // make and the CSV reader reports it as "invalid UTF-8 near byte index 1",
        // which says nothing about compression.
        if is_gzip_compressed(path) {
            return Err(AppError::validation(format!(
                "'{path}' is compressed; the report reads plain get_MNV TSV files. \
                 Decompress it first, for example with `gunzip -c '{path}' > out.MNV.tsv`"
            )));
        }
        let label = match sample {
            Some(name) => name.to_string(),
            None => sample_label(path),
        };
        let sample_id = self.samples.intern(&label);

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .from_path(path)
            .map_err(|e| {
                AppError::validation(format!("Cannot read TSV '{path}' for the report: {e}"))
            })?;

        let headers = reader
            .headers()
            .map_err(|e| AppError::validation(format!("Cannot read TSV header of '{path}': {e}")))?
            .clone();
        let column: HashMap<&str, usize> = headers
            .iter()
            .enumerate()
            .map(|(idx, name)| (name, idx))
            .collect();

        for required in ["Chromosome", "Gene", "Positions", "Variant Type"] {
            if !column.contains_key(required) {
                return Err(AppError::validation(format!(
                    "'{path}' does not look like a get_MNV TSV: missing the '{required}' column"
                )));
            }
        }

        let get = |record: &csv::StringRecord, name: &str| -> String {
            column
                .get(name)
                .and_then(|idx| record.get(*idx))
                .unwrap_or("")
                .to_string()
        };

        for record in reader.records() {
            let record = record.map_err(|e| {
                AppError::validation(format!("Failed reading a row of TSV '{path}': {e}"))
            })?;

            let variant_type = get(&record, "Variant Type");
            let is_indel = variant_type == "INDEL";
            let freq = if is_indel {
                first_number::<f64>(&get(&record, "Event Frequency"))
            } else {
                first_number::<f64>(&get(&record, "MNV Frequencies"))
                    .or_else(|| first_number::<f64>(&get(&record, "SNP Frequencies")))
            };
            let depth = if is_indel {
                first_number::<u64>(&get(&record, "Event Depth"))
            } else {
                first_number::<u64>(&get(&record, "Total Reads"))
            };
            let reads = if is_indel {
                first_number::<u64>(&get(&record, "Event Reads"))
            } else {
                first_number::<u64>(&get(&record, "MNV Reads"))
                    .or_else(|| first_number::<u64>(&get(&record, "SNP Reads")))
            };

            let row = Row {
                sample: sample_id,
                chrom: self.chroms.intern(&get(&record, "Chromosome")),
                gene: self.genes.intern(&get(&record, "Gene")),
                positions: get(&record, "Positions"),
                ref_bases: get(&record, "Reference Bases"),
                alt_bases: get(&record, "Base Changes"),
                aa_change: get(&record, "AA Changes"),
                variant_type: self.types.intern(&variant_type),
                change_type: self.changes.intern(&get(&record, "Change Type")),
                so_term: self.so_terms.intern(&get(&record, "SO Term")),
                impact: self.impacts.intern(&get(&record, "Impact")),
                grantham: leading_u32(&get(&record, "Grantham")),
                shift: optional(&get(&record, "MNV Consequence Shift")),
                dbs: optional(&get(&record, "DBS Class")),
                nmd: optional(&get(&record, "NMD Prediction")),
                hgvs_g: optional(&get(&record, "HGVS g.")),
                hgvs_c: optional(&get(&record, "HGVS c.")),
                freq,
                depth,
                ref_codon: optional(&get(&record, "Reference Codon")),
                alt_codon: optional(&get(&record, "MNV Codon")),
                event_class: optional(&get(&record, "Event Class")),
                event_components: optional(&get(&record, "Event Components")),
                reads,
                phasing: first_number::<f64>(&get(&record, "MNV Phasing Support")),
                phasing_reads: first_number::<u64>(&get(&record, "MNV Phasing Reads")),
            };
            self.rows.push(row);
        }

        Ok(())
    }

    /// Render the report to `path`. `command` is the sanitized command line shown
    /// in the header; `generated` is a human-readable timestamp.
    pub fn write(&self, path: &str, command: &str, generated: &str) -> AppResult<()> {
        let json = self.to_json(command, generated);
        let html = TEMPLATE.replace(DATA_PLACEHOLDER, &json);
        let mut file = std::fs::File::create(path).map_err(|e| {
            AppError::validation(format!("Cannot create report file '{path}': {e}"))
        })?;
        file.write_all(html.as_bytes())
            .map_err(|e| AppError::validation(format!("Cannot write report file '{path}': {e}")))?;
        Ok(())
    }

    fn to_json(&self, command: &str, generated: &str) -> String {
        let mut out = String::with_capacity(self.rows.len() * 160 + 1024);
        out.push_str("{\"meta\":{\"version\":");
        push_json_str(&mut out, env!("CARGO_PKG_VERSION"));
        out.push_str(",\"generated\":");
        push_json_str(&mut out, generated);
        out.push_str(",\"command\":");
        push_json_str(&mut out, command);
        out.push_str("},\"dict\":{");
        push_dict(&mut out, "samples", &self.samples.values, true);
        push_dict(&mut out, "chroms", &self.chroms.values, false);
        push_dict(&mut out, "genes", &self.genes.values, false);
        push_dict(&mut out, "types", &self.types.values, false);
        push_dict(&mut out, "changes", &self.changes.values, false);
        push_dict(&mut out, "so", &self.so_terms.values, false);
        push_dict(&mut out, "impacts", &self.impacts.values, false);
        out.push_str("},\"rows\":[");
        for (idx, row) in self.rows.iter().enumerate() {
            if idx > 0 {
                out.push(',');
            }
            push_row(&mut out, row);
        }
        out.push_str("]}");
        out
    }
}

fn push_dict(out: &mut String, key: &str, values: &[String], first: bool) {
    if !first {
        out.push(',');
    }
    out.push('"');
    out.push_str(key);
    out.push_str("\":[");
    for (idx, value) in values.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        push_json_str(out, value);
    }
    out.push(']');
}

fn push_row(out: &mut String, row: &Row) {
    out.push('[');
    push_usize(out, row.sample);
    out.push(',');
    push_usize(out, row.chrom);
    out.push(',');
    push_usize(out, row.gene);
    out.push(',');
    push_json_str(out, &row.positions);
    out.push(',');
    push_json_str(out, &row.ref_bases);
    out.push(',');
    push_json_str(out, &row.alt_bases);
    out.push(',');
    push_json_str(out, &row.aa_change);
    out.push(',');
    push_usize(out, row.variant_type);
    out.push(',');
    push_usize(out, row.change_type);
    out.push(',');
    push_usize(out, row.so_term);
    out.push(',');
    push_usize(out, row.impact);
    out.push(',');
    push_opt_num(out, row.grantham);
    out.push(',');
    push_opt_str(out, row.shift.as_deref());
    out.push(',');
    push_opt_str(out, row.dbs.as_deref());
    out.push(',');
    push_opt_str(out, row.nmd.as_deref());
    out.push(',');
    push_opt_str(out, row.hgvs_g.as_deref());
    out.push(',');
    push_opt_str(out, row.hgvs_c.as_deref());
    out.push(',');
    push_opt_freq(out, row.freq);
    out.push(',');
    push_opt_num(out, row.depth);
    out.push(',');
    push_opt_str(out, row.ref_codon.as_deref());
    out.push(',');
    push_opt_str(out, row.alt_codon.as_deref());
    out.push(',');
    push_opt_str(out, row.event_class.as_deref());
    out.push(',');
    push_opt_str(out, row.event_components.as_deref());
    out.push(',');
    push_opt_num(out, row.reads);
    out.push(',');
    push_opt_freq(out, row.phasing);
    out.push(',');
    push_opt_num(out, row.phasing_reads);
    out.push(']');
}

fn push_usize(out: &mut String, value: usize) {
    out.push_str(&value.to_string());
}

fn push_opt_num<T: std::fmt::Display>(out: &mut String, value: Option<T>) {
    match value {
        Some(v) => out.push_str(&v.to_string()),
        None => out.push_str("null"),
    }
}

fn push_opt_freq(out: &mut String, value: Option<f64>) {
    match value {
        Some(v) if v.is_finite() => out.push_str(&format!("{v:.4}")),
        _ => out.push_str("null"),
    }
}

fn push_opt_str(out: &mut String, value: Option<&str>) {
    match value {
        Some(v) => push_json_str(out, v),
        None => out.push_str("null"),
    }
}

/// Append a JSON string literal. Also escapes `<` and `/` inside `</` so the
/// payload can never terminate the surrounding `<script>` element.
fn push_json_str(out: &mut String, value: &str) {
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '<' => out.push_str("\\u003c"),
            '>' => out.push_str("\\u003e"),
            '&' => out.push_str("\\u0026"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_label_strips_mnv_suffixes() {
        assert_eq!(sample_label("/data/TB-001.MNV.tsv"), "TB-001");
        assert_eq!(sample_label("TB-002.MNV.tsv"), "TB-002");
        assert_eq!(sample_label("plain.tsv"), "plain");
    }

    #[test]
    fn sample_labels_disambiguate_identical_file_names() {
        // The canonical per-sample pipeline layout: same file name, different
        // directory. These must stay distinct samples.
        let paths = vec![
            "results/TB-001/variants.MNV.tsv".to_string(),
            "results/TB-002/variants.MNV.tsv".to_string(),
        ];
        assert_eq!(
            sample_labels(&paths),
            vec!["TB-001/variants", "TB-002/variants"]
        );
    }

    #[test]
    fn sample_labels_keep_the_bare_stem_when_it_is_unique() {
        let paths = vec![
            "a/TB-001.MNV.tsv".to_string(),
            "b/TB-002.MNV.tsv".to_string(),
        ];
        assert_eq!(sample_labels(&paths), vec!["TB-001", "TB-002"]);
    }

    #[test]
    fn sample_labels_qualify_only_the_colliding_group() {
        let paths = vec![
            "run/x/out.MNV.tsv".to_string(),
            "run/y/out.MNV.tsv".to_string(),
            "run/TB-003.MNV.tsv".to_string(),
        ];
        assert_eq!(sample_labels(&paths), vec!["x/out", "y/out", "TB-003"]);
    }

    #[test]
    fn sample_labels_walk_up_until_the_paths_differ() {
        let paths = vec![
            "runA/s/out.MNV.tsv".to_string(),
            "runB/s/out.MNV.tsv".to_string(),
        ];
        assert_eq!(sample_labels(&paths), vec!["runA/s/out", "runB/s/out"]);
    }

    #[test]
    fn sample_labels_number_repeats_of_the_same_path() {
        // The same file twice is indistinguishable by path, so both exhaust the
        // ladder and the numeric suffix keeps them apart.
        let paths = vec!["x/out.MNV.tsv".to_string(), "x/out.MNV.tsv".to_string()];
        assert_eq!(sample_labels(&paths), vec!["x/out", "x/out #2"]);
    }

    #[test]
    fn every_input_file_keeps_its_own_sample_row() {
        let dir = std::env::temp_dir().join("get_mnv_report_label_test");
        let header = "Chromosome\tGene\tPositions\tReference Bases\tBase Changes\t\
             AA Changes\tVariant Type\tChange Type\tSO Term\tImpact\n";
        let mut paths = Vec::new();
        for (sample, pos) in [("S1", 500), ("S2", 900)] {
            let sub = dir.join(sample);
            std::fs::create_dir_all(&sub).expect("temp dir");
            let path = sub.join("out.MNV.tsv");
            std::fs::write(
                &path,
                format!(
                    "{header}chr1\tgeneA\t{pos}\tA\tG\t-\tSNP\tNon-synonymous\t\
                     missense_variant\tMODERATE\n"
                ),
            )
            .expect("write tsv");
            paths.push(path.to_string_lossy().into_owned());
        }

        let labels = sample_labels(&paths);
        let mut builder = ReportBuilder::new();
        for (path, label) in paths.iter().zip(labels.iter()) {
            builder
                .add_tsv(path, Some(label))
                .expect("tsv should parse");
        }
        assert_eq!(builder.len(), 2);
        assert_eq!(
            builder.samples.values.len(),
            2,
            "identically named files must not merge into one sample: {:?}",
            builder.samples.values
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn compressed_input_is_reported_as_such() {
        let dir = std::env::temp_dir().join("get_mnv_report_gz_test");
        std::fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("S1.MNV.tsv.gz");
        // A gzip header is enough; the check never gets as far as the payload.
        std::fs::write(&path, [0x1f, 0x8b, 0x08, 0x00]).expect("write file");

        let mut builder = ReportBuilder::new();
        let err = builder
            .add_tsv(path.to_str().expect("utf8 path"), None)
            .expect_err("compressed input must be rejected");
        let message = err.to_string();
        assert!(
            message.contains("compressed") && message.contains("gunzip"),
            "error should name the cause and the fix: {message}"
        );

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn leading_u32_parses_grantham_cell() {
        assert_eq!(leading_u32("177 (radical)"), Some(177));
        assert_eq!(leading_u32("5 (conservative)"), Some(5));
        assert_eq!(leading_u32("-"), None);
        assert_eq!(leading_u32(""), None);
    }

    #[test]
    fn first_number_takes_the_leading_entry() {
        assert_eq!(first_number::<f64>("0.2500, 0.3000"), Some(0.25));
        assert_eq!(first_number::<u64>("12, 30"), Some(12));
        assert_eq!(first_number::<u64>("-"), None);
        assert_eq!(first_number::<f64>(""), None);
    }

    #[test]
    fn optional_treats_dash_as_absent() {
        assert_eq!(optional("-"), None);
        assert_eq!(optional("  "), None);
        assert_eq!(optional("CC>TT"), Some("CC>TT".to_string()));
    }

    #[test]
    fn json_string_escaping_cannot_break_out_of_the_script_tag() {
        let mut out = String::new();
        push_json_str(&mut out, "</script><img src=x>");
        assert!(!out.contains("</script>"), "got: {out}");
        assert!(out.contains("\\u003c"), "got: {out}");
    }

    #[test]
    fn interner_reuses_ids() {
        let mut interner = Interner::default();
        assert_eq!(interner.intern("a"), 0);
        assert_eq!(interner.intern("b"), 1);
        assert_eq!(interner.intern("a"), 0);
        assert_eq!(interner.values.len(), 2);
    }

    #[test]
    fn builder_reads_a_tsv_and_emits_embeddable_json() {
        let dir = std::env::temp_dir().join("get_mnv_report_test");
        std::fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("S1.MNV.tsv");
        std::fs::write(
            &path,
            "Chromosome\tGene\tPositions\tReference Bases\tBase Changes\tAA Changes\t\
             Variant Type\tChange Type\tSO Term\tImpact\tGrantham\tMNV Consequence Shift\t\
             DBS Class\tNMD Prediction\tHGVS g.\tHGVS c.\tReference Codon\tMNV Codon\t\
             Event Class\tEvent Components\n\
             chr1\tgeneA\t28, 30\tG, T\tT, A\tAla10Ser\tSNP/MNV\tNon-synonymous\t\
             missense_variant\tMODERATE\t99 (moderately conservative)\tMNV-gained\t\
             CC>TT\t-\tg.[28G>T;30T>A]\tc.[28G>T;30T>A]\tGCT\tTCA\tmnv\tSNV:28:G>T\n",
        )
        .expect("write tsv");

        let mut builder = ReportBuilder::new();
        builder
            .add_tsv(path.to_str().expect("utf8 path"), None)
            .expect("tsv should parse");
        assert_eq!(builder.len(), 1);

        let json = builder.to_json("get_mnv --report r.html", "2026-07-25");
        assert!(json.contains("\"S1\""), "sample label missing: {json}");
        assert!(json.contains("geneA"));
        assert!(
            json.contains("CC\\u003eTT"),
            "DBS should be escaped: {json}"
        );
        // The payload must parse as JSON and must not close the script element.
        assert!(!json.contains("</script>"));
        let _: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");

        let _ = std::fs::remove_file(&path);
    }
}
