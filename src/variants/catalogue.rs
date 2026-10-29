//! Matching calls against a catalogue of known variants.
//!
//! A catalogue is the user's own table of variants that mean something: drug
//! resistance mutations, lineage markers, known pathogenic changes. Each entry
//! names a gene and a variant, at the protein (`p.Ser450Leu`), coding
//! (`c.1349C>T`, `c.-15C>T`) or genomic (`g.761155C>T`) level, or `LoF` for any
//! loss of function in the gene, and carries a free-text label.
//!
//! Protein entries are where codon-level annotation matters. A codon with two
//! substitutions has two readings, the residue each base gives on its own and
//! the one the whole codon gives, and a catalogue built from base-by-base
//! annotation can be right about one and wrong about the other. A match is
//! reported with which reading made it: the codon only (a catalogued change
//! that no single base gives, which a base-by-base annotator misses) or the
//! single bases only (a change such an annotator would report, which the codon
//! does not make).

use crate::error::{AppError, AppResult};
use crate::variants::{parse_component_label, AlleleComponentKind, Gene, Strand, VariantInfo};
use std::collections::HashMap;

/// What an entry describes, normalized for comparison.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CatalogueVariant {
    /// A protein change in get_MNV's own three-letter form (`Ser450Leu`,
    /// `Met1?`), or a frameshift reduced to its first residue (`Ala11fs`).
    Protein(String),
    /// A single-base substitution at a coding position, in coding-strand bases.
    CodingSubstitution {
        position: CodingAnchor,
        reference: char,
        alternate: char,
    },
    /// Any other coding descriptor, compared as written (`c.517_519dup`).
    CodingOther(String),
    /// A single-base substitution at a genomic position.
    GenomicSubstitution {
        position: usize,
        reference: char,
        alternate: char,
    },
    /// Any other genomic descriptor, compared as written.
    GenomicOther(String),
    /// A frameshift, premature stop or start loss anywhere in the gene.
    LossOfFunction,
}

/// Where a coding position lies relative to the CDS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CodingAnchor {
    /// `c.N`, inside the CDS.
    Cds(usize),
    /// `c.-N`, N bases upstream of the first coding base.
    Upstream(usize),
    /// `c.*N`, N bases downstream of the last coding base.
    Downstream(usize),
}

#[derive(Debug, Clone)]
pub struct CatalogueEntry {
    pub gene: String,
    /// The variant as the catalogue wrote it, for the output.
    pub written: String,
    pub variant: CatalogueVariant,
    pub label: String,
}

#[derive(Debug, Clone, Default)]
pub struct Catalogue {
    pub entries: Vec<CatalogueEntry>,
}

fn three_letter(residue: char) -> Option<&'static str> {
    let name = crate::utils::aa_three_letter(residue.to_ascii_uppercase());
    (name != "X").then_some(name)
}

/// A protein change as get_MNV writes it: three-letter residues, `*` for a
/// stop, the reference residue repeated for a synonymous change, and a
/// frameshift reduced to its first residue and position.
pub fn normalize_protein(change: &str) -> Option<String> {
    let mut text = change.trim().trim_start_matches("p.").to_string();
    if text.starts_with('(') && text.ends_with(')') {
        text = text[1..text.len() - 1].to_string();
    }
    let text = text.replace("Ter", "*");
    // One-letter form: S450L, M1?, *100Q.
    let one_letter = regex_one_letter(&text);
    let text = match one_letter {
        Some((reference, position, rest)) => {
            let first = if reference == '*' {
                "*".to_string()
            } else {
                three_letter(reference)?.to_string()
            };
            let rest = match rest.as_str() {
                "=" => first.clone(),
                "?" => "?".to_string(),
                "fs" => "fs".to_string(),
                r if r.chars().count() == 1 => {
                    let c = r.chars().next()?;
                    if c == '*' {
                        "*".to_string()
                    } else {
                        three_letter(c)?.to_string()
                    }
                }
                r => r.to_string(),
            };
            format!("{first}{position}{rest}")
        }
        None => text,
    };
    // Synonymous in the three-letter form: Ser450= is Ser450Ser.
    let text = match text.strip_suffix('=') {
        Some(head) => {
            let residue: String = head.chars().take_while(|c| !c.is_ascii_digit()).collect();
            format!("{head}{residue}")
        }
        None => text,
    };
    Some(reduce_frameshift(&text))
}

/// `Ala11LeufsTer5`, `Ala11Leufs*5` and `Ala11fs` all reduce to `Ala11fs`.
fn reduce_frameshift(change: &str) -> String {
    let Some(at) = change.find("fs") else {
        return change.to_string();
    };
    let head = &change[..at];
    let residue: String = head.chars().take_while(|c| !c.is_ascii_digit()).collect();
    let digits: String = head[residue.len()..]
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    if residue.is_empty() || digits.is_empty() {
        return change.to_string();
    }
    format!("{residue}{digits}fs")
}

/// `S450L` as (S, 450, "L"); `M1?` as (M, 1, "?").
fn regex_one_letter(text: &str) -> Option<(char, String, String)> {
    let mut chars = text.chars();
    let reference = chars.next()?;
    if !(reference.is_ascii_uppercase() || reference == '*') {
        return None;
    }
    let rest: String = chars.collect();
    let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
    if digits.is_empty() {
        return None;
    }
    let tail = rest[digits.len()..].to_string();
    // A three-letter change has lowercase letters right after the first one.
    if text.chars().nth(1).is_some_and(|c| c.is_ascii_lowercase()) {
        return None;
    }
    let single = tail.chars().count() == 1 || tail == "fs" || tail == "=" || tail == "?";
    single.then_some((reference, digits, tail))
}

fn parse_substitution(text: &str) -> Option<(String, char, char)> {
    let (position, change) = text.split_at(
        text.find(|c: char| c.is_ascii_alphabetic() && c != '-' && c.is_ascii_uppercase())?,
    );
    let mut parts = change.split('>');
    let (r, a) = (parts.next()?, parts.next()?);
    if parts.next().is_some() || r.len() != 1 || a.len() != 1 {
        return None;
    }
    let (r, a) = (
        r.chars().next()?.to_ascii_uppercase(),
        a.chars().next()?.to_ascii_uppercase(),
    );
    ("ACGT".contains(r) && "ACGT".contains(a)).then(|| (position.to_string(), r, a))
}

/// Read one entry's variant column.
pub fn parse_variant(text: &str) -> Option<CatalogueVariant> {
    let text = text.trim();
    if text.eq_ignore_ascii_case("lof") {
        return Some(CatalogueVariant::LossOfFunction);
    }
    if let Some(rest) = text.strip_prefix("c.") {
        if let Some((position, reference, alternate)) = parse_substitution(rest) {
            let anchor = if let Some(n) = position.strip_prefix('-') {
                CodingAnchor::Upstream(n.parse().ok()?)
            } else if let Some(n) = position.strip_prefix('*') {
                CodingAnchor::Downstream(n.parse().ok()?)
            } else {
                CodingAnchor::Cds(position.parse().ok()?)
            };
            return Some(CatalogueVariant::CodingSubstitution {
                position: anchor,
                reference,
                alternate,
            });
        }
        return Some(CatalogueVariant::CodingOther(text.to_string()));
    }
    if let Some(rest) = text.strip_prefix("g.") {
        if let Some((position, reference, alternate)) = parse_substitution(rest) {
            return Some(CatalogueVariant::GenomicSubstitution {
                position: position.parse().ok()?,
                reference,
                alternate,
            });
        }
        return Some(CatalogueVariant::GenomicOther(text.to_string()));
    }
    normalize_protein(text).map(CatalogueVariant::Protein)
}

impl Catalogue {
    /// Read a tab-separated catalogue with a header naming at least `gene` and
    /// `variant` columns, and optionally `label`. Lines starting with `#` are
    /// comments.
    pub fn load(path: &str) -> AppResult<Self> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| AppError::validation(format!("Cannot read catalogue '{path}': {e}")))?;
        let mut lines = text
            .lines()
            .enumerate()
            .filter(|(_, line)| !line.trim().is_empty() && !line.starts_with('#'));
        let (_, header) = lines
            .next()
            .ok_or_else(|| AppError::validation(format!("Catalogue '{path}' is empty")))?;
        let columns: Vec<String> = header
            .split('\t')
            .map(|c| c.trim().to_ascii_lowercase())
            .collect();
        let index = |name: &str| columns.iter().position(|c| c == name);
        let (gene_at, variant_at) = match (index("gene"), index("variant")) {
            (Some(g), Some(v)) => (g, v),
            _ => {
                return Err(AppError::validation(format!(
                    "Catalogue '{path}' needs a header with 'gene' and 'variant' columns, got: {header}"
                )))
            }
        };
        let label_at = index("label");
        let mut entries = Vec::new();
        for (number, line) in lines {
            let fields: Vec<&str> = line.split('\t').collect();
            let cell = |i: usize| fields.get(i).map(|f| f.trim()).unwrap_or("");
            let (gene, written) = (cell(gene_at), cell(variant_at));
            if gene.is_empty() || written.is_empty() {
                return Err(AppError::validation(format!(
                    "Catalogue '{path}' line {}: a gene and a variant are both needed",
                    number + 1
                )));
            }
            let variant = parse_variant(written).ok_or_else(|| {
                AppError::validation(format!(
                    "Catalogue '{path}' line {}: cannot read the variant '{written}' \
                     (expected p., c., g. notation or LoF)",
                    number + 1
                ))
            })?;
            entries.push(CatalogueEntry {
                gene: gene.to_string(),
                written: written.to_string(),
                variant,
                label: label_at.map(cell).unwrap_or("").to_string(),
            });
        }
        log::info!("Catalogue '{path}': {} entries", entries.len());
        Ok(Self { entries })
    }
}

/// Whether a row's gene is the gene a catalogue names: the same name, or one
/// of the `_`-separated parts of a combined name such as `rpoB_Rv0667`.
pub fn gene_matches(row_gene: &str, entry_gene: &str) -> bool {
    row_gene == entry_gene || row_gene.split('_').any(|part| part == entry_gene)
}

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        other => other,
    }
}

/// The genomic position and genomic bases of a coding substitution, from the
/// gene objects that carry the gene's name on this contig.
fn place_coding(
    genes: &[&Gene],
    anchor: CodingAnchor,
    reference: char,
    alternate: char,
) -> Option<(usize, char, char)> {
    let strand = genes.first()?.strand;
    let position = match anchor {
        CodingAnchor::Cds(n) => genes
            .iter()
            .find_map(|g| crate::variants::genomic_for_cds(g, n))?,
        CodingAnchor::Upstream(n) => {
            let first = genes
                .iter()
                .find_map(|g| crate::variants::genomic_for_cds(g, 1))?;
            match strand {
                Strand::Plus => first.checked_sub(n)?,
                Strand::Minus => first + n,
            }
        }
        CodingAnchor::Downstream(n) => {
            let (gene, last) = genes
                .iter()
                .filter_map(|g| Some((*g, crate::variants::last_cds_position(g)?)))
                .max_by_key(|(_, last)| *last)?;
            let end = crate::variants::genomic_for_cds(gene, last)?;
            match strand {
                Strand::Plus => end + n,
                Strand::Minus => end.checked_sub(n)?,
            }
        }
    };
    Some(match strand {
        Strand::Plus => (position, reference, alternate),
        Strand::Minus => (position, complement(reference), complement(alternate)),
    })
}

fn substitutions_of(variant: &VariantInfo) -> Vec<(usize, char)> {
    variant
        .event_components
        .iter()
        .filter_map(|label| parse_component_label(label))
        .filter(|c| c.kind == AlleleComponentKind::Snp)
        .filter_map(|c| {
            Some((
                c.position,
                c.alt_allele.chars().next()?.to_ascii_uppercase(),
            ))
        })
        .collect()
}

fn describe(entry: &CatalogueEntry, note: Option<&str>) -> String {
    let what = format!("{} {}", entry.gene, entry.written);
    let detail = match note {
        Some(note) => format!("{what}, {note}"),
        None => what,
    };
    if entry.label.is_empty() {
        detail
    } else {
        format!("{} ({detail})", entry.label)
    }
}

const LOSS_OF_FUNCTION: [&str; 3] = ["frameshift_variant", "stop_gained", "start_lost"];

/// Annotate every row of a contig with the catalogue entries it matches.
///
/// `genes` are all the annotation's genes on this contig, with or without a
/// variant: an entry upstream of a gene is placed from the gene's start. An
/// entry whose reference base is not the one this FASTA has at that position
/// was written against another reference and matches nothing; how many there
/// were is returned so the run can say so.
pub fn annotate_catalogue(
    catalogue: &Catalogue,
    variants: &mut [VariantInfo],
    genes: &[Gene],
    sequence: &str,
) -> usize {
    let mut other_reference = 0usize;
    // Resolve the positional entries once per contig.
    let mut by_position: HashMap<(usize, char), Vec<usize>> = HashMap::new();
    for (index, entry) in catalogue.entries.iter().enumerate() {
        let placed = match &entry.variant {
            CatalogueVariant::CodingSubstitution {
                position,
                reference,
                alternate,
            } => {
                let named: Vec<&Gene> = genes
                    .iter()
                    .filter(|g| gene_matches(&g.name, &entry.gene))
                    .collect();
                if named.is_empty() {
                    continue;
                }
                place_coding(&named, *position, *reference, *alternate)
            }
            CatalogueVariant::GenomicSubstitution {
                position,
                reference,
                alternate,
            } => Some((*position, *reference, *alternate)),
            _ => continue,
        };
        let Some((position, reference, alternate)) = placed else {
            continue;
        };
        let on_this_fasta = sequence
            .as_bytes()
            .get(position.wrapping_sub(1))
            .is_some_and(|b| b.eq_ignore_ascii_case(&(reference as u8)));
        if !on_this_fasta {
            other_reference += 1;
            continue;
        }
        by_position
            .entry((position, alternate))
            .or_default()
            .push(index);
    }

    for variant in variants.iter_mut() {
        let mut hits: Vec<String> = Vec::new();
        for (position, alternate) in substitutions_of(variant) {
            for &index in by_position
                .get(&(position, alternate))
                .into_iter()
                .flatten()
            {
                hits.push(describe(&catalogue.entries[index], None));
            }
        }
        let (so_term, _) = crate::output::so_consequence(variant);
        let genomic = variant
            .annotations
            .hgvs_g
            .clone()
            .or_else(|| crate::variants::hgvs::genomic(variant));
        for entry in &catalogue.entries {
            if !gene_matches(&variant.gene, &entry.gene) {
                continue;
            }
            match &entry.variant {
                CatalogueVariant::Protein(change) => {
                    let codon = variant
                        .aa_changes
                        .iter()
                        .filter_map(|a| normalize_protein(a))
                        .any(|a| &a == change);
                    let single = variant.snp_aa_changes.len() > 1
                        && variant
                            .snp_aa_changes
                            .iter()
                            .filter_map(|a| normalize_protein(a))
                            .any(|a| &a == change);
                    let several = variant.snp_aa_changes.len() > 1;
                    let note = match (codon, single) {
                        (true, false) if several => Some("codon reading only"),
                        (false, true) => Some("per-SNV reading only"),
                        _ => None,
                    };
                    if codon || single {
                        hits.push(describe(entry, note));
                    }
                }
                CatalogueVariant::CodingOther(written)
                    if variant.annotations.hgvs_c.as_deref() == Some(written.as_str()) =>
                {
                    hits.push(describe(entry, None));
                }
                CatalogueVariant::GenomicOther(written)
                    if genomic.as_deref() == Some(written.as_str()) =>
                {
                    hits.push(describe(entry, None));
                }
                CatalogueVariant::LossOfFunction
                    if so_term
                        .split('&')
                        .any(|term| LOSS_OF_FUNCTION.contains(&term)) =>
                {
                    hits.push(describe(entry, None));
                }
                _ => {}
            }
        }
        hits.dedup();
        variant.annotations.catalogue = hits;
    }
    other_reference
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_changes_are_read_in_every_usual_spelling() {
        for written in [
            "p.Ser450Leu",
            "Ser450Leu",
            "p.S450L",
            "S450L",
            "p.(Ser450Leu)",
        ] {
            assert_eq!(
                normalize_protein(written).as_deref(),
                Some("Ser450Leu"),
                "{written}"
            );
        }
        assert_eq!(normalize_protein("p.Gln10Ter").as_deref(), Some("Gln10*"));
        assert_eq!(normalize_protein("p.Q10*").as_deref(), Some("Gln10*"));
        assert_eq!(normalize_protein("p.Ser450=").as_deref(), Some("Ser450Ser"));
        assert_eq!(normalize_protein("p.M1?").as_deref(), Some("Met1?"));
        for written in [
            "p.Ala11fs",
            "p.Ala11LeufsTer5",
            "p.Ala11Leufs*5",
            "Ala11Leufs",
        ] {
            assert_eq!(
                normalize_protein(written).as_deref(),
                Some("Ala11fs"),
                "{written}"
            );
        }
    }

    #[test]
    fn coding_and_genomic_entries_are_read_by_position() {
        assert_eq!(
            parse_variant("c.-15C>T"),
            Some(CatalogueVariant::CodingSubstitution {
                position: CodingAnchor::Upstream(15),
                reference: 'C',
                alternate: 'T'
            })
        );
        assert_eq!(
            parse_variant("c.1349C>T"),
            Some(CatalogueVariant::CodingSubstitution {
                position: CodingAnchor::Cds(1349),
                reference: 'C',
                alternate: 'T'
            })
        );
        assert_eq!(
            parse_variant("c.*3G>A"),
            Some(CatalogueVariant::CodingSubstitution {
                position: CodingAnchor::Downstream(3),
                reference: 'G',
                alternate: 'A'
            })
        );
        assert_eq!(
            parse_variant("g.761155C>T"),
            Some(CatalogueVariant::GenomicSubstitution {
                position: 761155,
                reference: 'C',
                alternate: 'T'
            })
        );
        assert_eq!(
            parse_variant("c.517_519dup"),
            Some(CatalogueVariant::CodingOther("c.517_519dup".into()))
        );
        assert_eq!(parse_variant("LoF"), Some(CatalogueVariant::LossOfFunction));
    }

    #[test]
    fn a_combined_gene_name_answers_to_each_of_its_parts() {
        assert!(gene_matches("rpoB_Rv0667", "rpoB"));
        assert!(gene_matches("rpoB_Rv0667", "Rv0667"));
        assert!(gene_matches("rpoB", "rpoB"));
        assert!(!gene_matches("rpoBC", "rpoB"));
    }
}
