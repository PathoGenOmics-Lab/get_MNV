//! Matching calls against a catalogue of known variants.
//!
//! A catalogue is the user's own table of variants that mean something: drug
//! resistance mutations, lineage markers, known pathogenic changes. Each entry
//! names a gene and a variant, at the protein (`p.Ser450Leu`), coding
//! (`c.1349C>T`, `c.-15C>T`), non-coding (`n.1401A>G`, for a gene that is not
//! translated) or genomic (`g.761155C>T`) level, or `LoF` for any loss of
//! function in the gene, and carries a free-text label.
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
use crate::variants::CodingNumbering;
use crate::variants::{parse_component_label, AlleleComponentKind, Gene, Strand, VariantInfo};
use std::collections::HashMap;

/// What an entry describes, normalized for comparison.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CatalogueVariant {
    /// A protein change in get_MNV's own three-letter form (`Ser450Leu`,
    /// `Met1?`), or a frameshift reduced to its first residue (`Ala11fs`).
    Protein(String),
    /// A single-base substitution at a coding position, in coding-strand bases:
    /// `c.N`, or `c.N+d` / `c.N-d` for a base `d` into the intron or flank
    /// beside coding base `N` (`offset` is then `d` or `-d`). `n.N` numbers a
    /// gene that is not translated from its first base, as `c.N` does there.
    CodingSubstitution {
        position: CodingAnchor,
        offset: isize,
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
    // An rRNA or other untranslated gene is numbered `n.` from its first base,
    // the base a `c.` number starts from in such a gene: the WHO catalogue
    // writes rrs and rrl changes so, and they were read as protein changes that
    // matched nothing.
    if let Some(rest) = text.strip_prefix("c.").or_else(|| text.strip_prefix("n.")) {
        if let Some((position, reference, alternate)) = parse_substitution(rest) {
            // An intronic or flanking base is written from the coding base
            // beside it: c.30+20, c.31-11, c.-15-3. The sign that opens the
            // number belongs to the anchor; the one after a digit, to the offset.
            let split = position
                .char_indices()
                .skip(1)
                .find(|&(i, c)| {
                    (c == '+' || c == '-') && position[..i].ends_with(|d: char| d.is_ascii_digit())
                })
                .map(|(i, _)| i);
            let (anchor, offset) = match split {
                Some(i) => (&position[..i], position[i..].parse::<isize>().ok()?),
                None => (position.as_str(), 0),
            };
            let anchor = if let Some(n) = anchor.strip_prefix('-') {
                CodingAnchor::Upstream(n.parse().ok()?)
            } else if let Some(n) = anchor.strip_prefix('*') {
                CodingAnchor::Downstream(n.parse().ok()?)
            } else {
                CodingAnchor::Cds(anchor.parse().ok()?)
            };
            return Some(CatalogueVariant::CodingSubstitution {
                position: anchor,
                offset,
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
/// gene objects that carry the gene's name on this contig, each with its coding
/// numbering (`c.` number to genomic position).
fn place_coding(
    genes: &[(&Gene, &CodingNumbering)],
    anchor: CodingAnchor,
    offset: isize,
    reference: char,
    alternate: char,
) -> Option<(usize, char, char)> {
    let strand = genes.first()?.0.strand;
    // Along the transcript: downstream is up the genome on the plus strand and
    // down it on the minus.
    let along = |position: usize, by: isize| -> Option<usize> {
        let by = match strand {
            Strand::Plus => by,
            Strand::Minus => -by,
        };
        position.checked_add_signed(by).filter(|&p| p >= 1)
    };
    let position = match anchor {
        CodingAnchor::Cds(n) => genes
            .iter()
            .find_map(|(_, numbering)| numbering.genomic(n))?,
        CodingAnchor::Upstream(n) => {
            let first = genes
                .iter()
                .find_map(|(_, numbering)| numbering.genomic(1))?;
            along(first, -(n as isize))?
        }
        CodingAnchor::Downstream(n) => {
            let (_, end) = genes
                .iter()
                .filter_map(|(_, numbering)| numbering.last())
                .max_by_key(|(last, _)| *last)?;
            along(end, n as isize)?
        }
    };
    let position = along(position, offset)?;
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
    // Each gene's coding numbering, built once per contig however many entries
    // name the gene: building it per entry made a catalogue of tens of
    // thousands of entries cost a pass over the gene each.
    let mut numberings: HashMap<usize, CodingNumbering> = HashMap::new();
    // The genes by every name they answer to (see `gene_matches`), so an entry
    // finds its gene without comparing names with every gene of the contig.
    let mut genes_by_name: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, gene) in genes.iter().enumerate() {
        for name in std::iter::once(gene.name.as_str()).chain(gene.name.split('_')) {
            let list = genes_by_name.entry(name).or_default();
            if list.last() != Some(&i) {
                list.push(i);
            }
        }
    }
    // Resolve the positional entries once per contig.
    let mut by_position: HashMap<(usize, char), Vec<usize>> = HashMap::new();
    // A genomic descriptor names its own place too, so it is looked up for
    // every row of the contig rather than among the rows of its gene: the row
    // of a deletion upstream of a gene is intergenic, and an entry met only
    // through the gene's name could never match it.
    let mut by_descriptor: HashMap<&str, Vec<usize>> = HashMap::new();
    for (index, entry) in catalogue.entries.iter().enumerate() {
        let named: Vec<usize> = genes_by_name
            .get(entry.gene.as_str())
            .cloned()
            .unwrap_or_default();
        let placed = match &entry.variant {
            CatalogueVariant::CodingSubstitution {
                position,
                offset,
                reference,
                alternate,
            } => {
                if named.is_empty() {
                    continue;
                }
                for &i in &named {
                    numberings
                        .entry(i)
                        .or_insert_with(|| CodingNumbering::of(&genes[i]));
                }
                let named: Vec<(&Gene, &CodingNumbering)> = named
                    .iter()
                    .map(|&i| (&genes[i], &numberings[&i]))
                    .collect();
                place_coding(&named, *position, *offset, *reference, *alternate)
            }
            // A genomic position means nothing without its contig, and the
            // entry's gene is what says which one: on any other contig the same
            // number is another base.
            CatalogueVariant::GenomicSubstitution {
                position,
                reference,
                alternate,
            } if !named.is_empty() => Some((*position, *reference, *alternate)),
            CatalogueVariant::GenomicOther(written) if !named.is_empty() => {
                by_descriptor
                    .entry(written.as_str())
                    .or_default()
                    .push(index);
                continue;
            }
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

    // The entries each row is compared with, by the gene name they give: a row
    // looks up its own name and each `_`-separated part of it, and meets them
    // in the catalogue's own order.
    let mut by_gene: HashMap<&str, Vec<usize>> = HashMap::new();
    for (index, entry) in catalogue.entries.iter().enumerate() {
        by_gene.entry(entry.gene.as_str()).or_default().push(index);
    }
    // get_MNV's own labels are not protein changes: a codon read in a shifted
    // frame (`Ala20Thr (fs)`) has no frameshift of its own at that codon, and
    // reading the label as one matched every frameshift entry there.
    let protein_changes = |changes: &[String]| -> Vec<String> {
        changes
            .iter()
            .filter(|change| !change.contains(char::is_whitespace))
            .filter_map(|change| normalize_protein(change))
            .collect()
    };

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
        let genomic = variant
            .annotations
            .hgvs_g
            .clone()
            .or_else(|| crate::variants::hgvs::genomic(variant));
        for &index in genomic
            .as_deref()
            .and_then(|written| by_descriptor.get(written))
            .into_iter()
            .flatten()
        {
            hits.push(describe(&catalogue.entries[index], None));
        }
        let (so_term, _) = crate::output::so_consequence(variant);
        let codon_changes = protein_changes(&variant.aa_changes);
        let single_changes = protein_changes(&variant.snp_aa_changes);
        let mut indices: Vec<usize> = std::iter::once(variant.gene.as_str())
            .chain(variant.gene.split('_'))
            .filter_map(|name| by_gene.get(name))
            .flatten()
            .copied()
            .collect();
        indices.sort_unstable();
        indices.dedup();
        for entry in indices.iter().map(|&index| &catalogue.entries[index]) {
            match &entry.variant {
                CatalogueVariant::Protein(change) => {
                    let codon = codon_changes.iter().any(|a| a == change);
                    let single = variant.snp_aa_changes.len() > 1
                        && single_changes.iter().any(|a| a == change);
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

        // Each base of a row of several, on its own: its position, its own
        // residue and consequence, its own descriptors.
        if variant.positions.len() > 1 {
            let mut bases: Vec<Vec<String>> = vec![Vec::new(); variant.positions.len()];
            for (i, base_hits) in bases.iter_mut().enumerate() {
                let position = variant.positions[i];
                let reference = variant.ref_bases.get(i).map(String::as_str).unwrap_or("");
                let alternate = variant
                    .base_changes
                    .get(i)
                    .map(String::as_str)
                    .unwrap_or("");
                if let Some(alt) = alternate.chars().next().filter(|_| alternate.len() == 1) {
                    for &index in by_position
                        .get(&(position, alt.to_ascii_uppercase()))
                        .into_iter()
                        .flatten()
                    {
                        base_hits.push(describe(&catalogue.entries[index], None));
                    }
                }
                let single = variant.annotations.single_bases.get(i);
                let own_change =
                    protein_changes(variant.snp_aa_changes.get(i..=i).unwrap_or_default());
                let own_term = single.map(|effect| {
                    crate::output::so_consequence_for(
                        variant,
                        effect.in_context(variant.change_type),
                        effect.splice,
                    )
                    .0
                });
                let own_genomic = format!(
                    "g.{position}{}>{}",
                    reference.to_ascii_uppercase(),
                    alternate.to_ascii_uppercase()
                );
                for &index in by_descriptor
                    .get(own_genomic.as_str())
                    .into_iter()
                    .flatten()
                {
                    base_hits.push(describe(&catalogue.entries[index], None));
                }
                for entry in indices.iter().map(|&index| &catalogue.entries[index]) {
                    let hit = match &entry.variant {
                        CatalogueVariant::Protein(change) => own_change.iter().any(|a| a == change),
                        CatalogueVariant::CodingOther(written) => {
                            single.and_then(|e| e.hgvs_c.as_deref()) == Some(written.as_str())
                        }
                        CatalogueVariant::LossOfFunction => own_term.as_deref().is_some_and(|t| {
                            t.split('&').any(|term| LOSS_OF_FUNCTION.contains(&term))
                        }),
                        _ => false,
                    };
                    if hit {
                        base_hits.push(describe(entry, None));
                    }
                }
                base_hits.dedup();
            }
            variant.annotations.catalogue_bases = bases;
        }
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
                offset: 0,
                reference: 'C',
                alternate: 'T'
            })
        );
        assert_eq!(
            parse_variant("c.1349C>T"),
            Some(CatalogueVariant::CodingSubstitution {
                position: CodingAnchor::Cds(1349),
                offset: 0,
                reference: 'C',
                alternate: 'T'
            })
        );
        assert_eq!(
            parse_variant("c.*3G>A"),
            Some(CatalogueVariant::CodingSubstitution {
                position: CodingAnchor::Downstream(3),
                offset: 0,
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

    /// An untranslated gene is numbered `n.`, as the WHO catalogue numbers rrs
    /// and rrl; such an entry was taken for a protein change and matched
    /// nothing.
    #[test]
    fn non_coding_entries_are_read_by_position_like_coding_ones() {
        assert_eq!(
            parse_variant("n.1401A>G"),
            Some(CatalogueVariant::CodingSubstitution {
                position: CodingAnchor::Cds(1401),
                offset: 0,
                reference: 'A',
                alternate: 'G'
            })
        );
        assert_eq!(
            parse_variant("n.-119C>T"),
            Some(CatalogueVariant::CodingSubstitution {
                position: CodingAnchor::Upstream(119),
                offset: 0,
                reference: 'C',
                alternate: 'T'
            })
        );
        assert_eq!(
            parse_variant("n.2269_2270insT"),
            Some(CatalogueVariant::CodingOther("n.2269_2270insT".into()))
        );
    }

    /// A base in an intron or a flank is written from the coding base beside it,
    /// as get_MNV writes it itself; such an entry stopped the whole run.
    #[test]
    fn intronic_and_flanking_coding_entries_carry_their_offset() {
        for (written, position, offset) in [
            ("c.30+20A>G", CodingAnchor::Cds(30), 20),
            ("c.31-11A>G", CodingAnchor::Cds(31), -11),
            ("c.-15-3A>G", CodingAnchor::Upstream(15), -3),
            ("c.*2+1A>G", CodingAnchor::Downstream(2), 1),
        ] {
            assert_eq!(
                parse_variant(written),
                Some(CatalogueVariant::CodingSubstitution {
                    position,
                    offset,
                    reference: 'A',
                    alternate: 'G'
                }),
                "{written}"
            );
        }
    }

    #[test]
    fn a_combined_gene_name_answers_to_each_of_its_parts() {
        assert!(gene_matches("rpoB_Rv0667", "rpoB"));
        assert!(gene_matches("rpoB_Rv0667", "Rv0667"));
        assert!(gene_matches("rpoB", "rpoB"));
        assert!(!gene_matches("rpoBC", "rpoB"));
    }
}
