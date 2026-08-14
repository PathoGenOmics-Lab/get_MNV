//! GFF/GFF3 attribute and field parsing primitives.

use crate::error::AppResult;
use std::collections::HashMap;

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
    // A quote delimits a value only where a value begins: right after the `=`
    // of a GFF3 pair or the space of a GTF one. GFF3 reserves `%`, `;`, `=`,
    // `&`, `,` and the control characters and nothing else, so a quote inside a
    // value is an ordinary character, and treating it as a delimiter made every
    // following `;` stop separating fields: the rest of the column became one
    // field and every attribute after the quote was lost, `Parent` among them.
    // The exon then left its transcript, was annotated on its own and restarted
    // the amino-acid numbering at 1, so one character in a Note moved a residue
    // by seventeen positions. A value that genuinely opens with a quote still
    // has its semicolons protected.
    let mut previous: Option<char> = None;
    for (idx, ch) in attributes.char_indices() {
        match ch {
            '"' if !escaped && (in_quotes || matches!(previous, Some('=') | Some(' '))) => {
                in_quotes = !in_quotes;
            }
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
        previous = Some(ch);
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
pub(super) fn parse_strand(raw: &str, line_number: usize) -> AppResult<crate::variants::Strand> {
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
pub(super) fn parse_gff_phase(raw: &str, line_number: usize) -> AppResult<u8> {
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

pub(super) fn parse_interval(
    start_raw: &str,
    end_raw: &str,
    line_number: usize,
) -> AppResult<(usize, usize)> {
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

    // Only append a locus_tag suffix when it is actually present, non-empty
    // and different from the primary name. This avoids both the historical
    // `primary_primary` duplication when no locus_tag exists and a trailing
    // underscore (`primary_`) when locus_tag is the empty string.
    match attrs.get("locus_tag") {
        Some(locus) if !locus.is_empty() && locus != &primary => {
            format!("{primary}_{locus}")
        }
        _ => primary,
    }
}
