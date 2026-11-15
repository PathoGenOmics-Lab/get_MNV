//! Catalogue codons read from the molecules (`--read-catalogue-codons`).
//!
//! A catalogue match on a VCF row needs the caller to have reported the
//! variant. A resistant population below the caller's threshold, or two
//! minority changes in one codon merged into an amino acid nobody carries,
//! never reach the catalogue that way. This module asks the reads directly:
//! every codon a catalogue protein entry names is read molecule by molecule
//! (see [`crate::read_count::read_codon_triplets`]), and each codon allele the
//! molecules carry is reported with its count, its frequency and the entries it
//! matches. A triplet is one molecule's codon, so the amino acid reported is
//! one that a molecule carries.

use crate::genetic_code::GeneticCode;
use crate::read_count::CodonTriplets;
use crate::variants::catalogue::{gene_matches, normalize_protein, Catalogue, CatalogueVariant};
use crate::variants::{CodingNumbering, Gene, Strand};
use std::collections::BTreeMap;
use std::io::Write;

/// One codon a catalogue entry names, placed on this contig.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CatalogueCodon {
    /// Index of the gene in the list the codons were found in.
    pub gene: usize,
    /// Protein position, 1-based.
    pub codon: usize,
    /// Genomic positions of the codon's three bases, in the order the ribosome
    /// reads them.
    pub positions: [usize; 3],
    /// Catalogue entries naming this codon, by index.
    pub entries: Vec<usize>,
}

/// The protein position a change names, when it names exactly one residue
/// replaced by another (or itself, or a stop): `Asp94Gly`, `Gln445*`,
/// `Ser450Ser`. `None` for frameshifts, start loss, indels and ranges, which a
/// single codon cannot answer.
pub fn single_residue_position(change: &str) -> Option<usize> {
    let residue: String = change
        .chars()
        .take_while(char::is_ascii_alphabetic)
        .collect();
    if residue.len() != 3 {
        return None;
    }
    let rest = &change[residue.len()..];
    let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
    let tail = &rest[digits.len()..];
    let replaced = tail == "*"
        || (tail.len() == 3
            && tail.chars().next().is_some_and(|c| c.is_ascii_uppercase())
            && tail.chars().skip(1).all(|c| c.is_ascii_lowercase()));
    if !replaced {
        return None;
    }
    digits.parse().ok().filter(|&n| n > 0)
}

/// Every codon the catalogue's protein entries name in these genes, once each,
/// in gene and codon order.
pub fn catalogue_codons(catalogue: &Catalogue, genes: &[Gene]) -> Vec<CatalogueCodon> {
    let mut found: BTreeMap<(usize, usize), CatalogueCodon> = BTreeMap::new();
    let mut numberings: BTreeMap<usize, CodingNumbering> = BTreeMap::new();
    for (index, entry) in catalogue.entries.iter().enumerate() {
        let CatalogueVariant::Protein(change) = &entry.variant else {
            continue;
        };
        let Some(codon) = single_residue_position(change) else {
            continue;
        };
        for (gene_index, gene) in genes.iter().enumerate() {
            if !gene.biotype.is_coding() || !gene_matches(&gene.name, &entry.gene) {
                continue;
            }
            let numbering = numberings
                .entry(gene_index)
                .or_insert_with(|| CodingNumbering::of(gene));
            let first = codon * 3 - 2;
            let (Some(a), Some(b), Some(c)) = (
                numbering.genomic(first),
                numbering.genomic(first + 1),
                numbering.genomic(first + 2),
            ) else {
                continue;
            };
            let slot = found
                .entry((gene_index, codon))
                .or_insert_with(|| CatalogueCodon {
                    gene: gene_index,
                    codon,
                    positions: [a, b, c],
                    entries: Vec::new(),
                });
            if !slot.entries.contains(&index) {
                slot.entries.push(index);
            }
        }
    }
    found.into_values().collect()
}

/// One codon allele the molecules carry.
#[derive(Debug, Clone, PartialEq)]
pub struct CodonAllele {
    pub chrom: String,
    pub gene: String,
    pub codon: usize,
    pub positions: [usize; 3],
    /// Reference codon, coding strand.
    pub reference_codon: String,
    /// The molecules' codon, coding strand.
    pub codon_allele: String,
    /// `Asp94Gly`, `Asp94Asp` for a synonymous codon, `Gln445Ter` for a stop.
    pub aa_change: String,
    pub molecules: usize,
    pub spanning: usize,
    pub partial: usize,
    pub frequency: f64,
    /// Catalogue entries this allele matches, described as in the `Catalogue`
    /// column.
    pub catalogue: Vec<String>,
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn residue_name(residue: char) -> &'static str {
    if residue == '*' {
        "Ter"
    } else {
        crate::utils::aa_three_letter(residue)
    }
}

/// Thresholds a codon allele must meet to be reported.
#[derive(Debug, Clone, Copy)]
pub struct AlleleFloor {
    pub min_frequency: f64,
    pub min_molecules: usize,
}

/// The non-reference codon alleles of one catalogue codon that meet `floor`,
/// most frequent first.
#[allow(clippy::too_many_arguments)]
pub fn codon_alleles(
    chrom: &str,
    gene: &Gene,
    codon: &CatalogueCodon,
    triplets: &CodonTriplets,
    reference: &str,
    genetic_code: GeneticCode,
    catalogue: &Catalogue,
    floor: AlleleFloor,
) -> Vec<CodonAllele> {
    let coding = |genomic: [u8; 3]| -> [u8; 3] {
        match gene.strand {
            Strand::Plus => genomic,
            Strand::Minus => genomic.map(complement),
        }
    };
    let reference_bytes = reference.as_bytes();
    let Some(reference_genomic) = codon
        .positions
        .iter()
        .map(|&p| {
            reference_bytes
                .get(p.wrapping_sub(1))
                .map(u8::to_ascii_uppercase)
        })
        .collect::<Option<Vec<u8>>>()
    else {
        return Vec::new();
    };
    let reference_codon = coding([
        reference_genomic[0],
        reference_genomic[1],
        reference_genomic[2],
    ]);
    let reference_residue = genetic_code.translate(&reference_codon);
    let mut alleles: Vec<CodonAllele> = Vec::new();
    for (genomic, &molecules) in &triplets.counts {
        let allele = coding(*genomic);
        if allele == reference_codon || triplets.spanning == 0 {
            continue;
        }
        let frequency = molecules as f64 / triplets.spanning as f64;
        if molecules < floor.min_molecules || frequency < floor.min_frequency {
            continue;
        }
        let residue = genetic_code.translate(&allele);
        let aa_change = format!(
            "{}{}{}",
            residue_name(reference_residue),
            codon.codon,
            residue_name(residue)
        );
        let normalized = normalize_protein(&aa_change);
        let premature_stop = residue == '*' && reference_residue != '*';
        let matched = codon
            .entries
            .iter()
            .map(|&index| &catalogue.entries[index])
            .filter(|entry| match &entry.variant {
                CatalogueVariant::Protein(change) => Some(change) == normalized.as_ref(),
                _ => false,
            })
            .chain(catalogue.entries.iter().filter(|entry| {
                premature_stop
                    && entry.variant == CatalogueVariant::LossOfFunction
                    && gene_matches(&gene.name, &entry.gene)
            }))
            .map(|entry| {
                let what = format!("{} {}", entry.gene, entry.written);
                if entry.label.is_empty() {
                    what
                } else {
                    format!("{} ({what})", entry.label)
                }
            })
            .collect();
        alleles.push(CodonAllele {
            chrom: chrom.to_string(),
            gene: gene.name.clone(),
            codon: codon.codon,
            positions: codon.positions,
            reference_codon: String::from_utf8_lossy(&reference_codon).into_owned(),
            codon_allele: String::from_utf8_lossy(&allele).into_owned(),
            aa_change,
            molecules,
            spanning: triplets.spanning,
            partial: triplets.partial,
            frequency,
            catalogue: matched,
        });
    }
    alleles.sort_by(|a, b| {
        b.molecules
            .cmp(&a.molecules)
            .then(a.codon_allele.cmp(&b.codon_allele))
    });
    alleles
}

pub const CATALOGUE_CODON_COLUMNS: [&str; 12] = [
    "Chromosome",
    "Gene",
    "Codon",
    "Positions",
    "Reference Codon",
    "Codon Allele",
    "AA Change",
    "Molecules",
    "Spanning Molecules",
    "Partial Molecules",
    "Frequency",
    "Catalogue",
];

pub fn write_header(out: &mut impl Write) -> std::io::Result<()> {
    writeln!(out, "{}", CATALOGUE_CODON_COLUMNS.join("\t"))
}

pub fn write_allele(out: &mut impl Write, allele: &CodonAllele) -> std::io::Result<()> {
    let catalogue = if allele.catalogue.is_empty() {
        "-".to_string()
    } else {
        allele.catalogue.join("; ")
    };
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
        allele.chrom,
        allele.gene,
        allele.codon,
        allele
            .positions
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", "),
        allele.reference_codon,
        allele.codon_allele,
        allele.aa_change,
        allele.molecules,
        allele.spanning,
        allele.partial,
        allele.frequency,
        catalogue
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::catalogue::CatalogueEntry;

    fn gene(name: &str, start: usize, end: usize, strand: Strand) -> Gene {
        Gene {
            name: name.to_string(),
            start,
            end,
            strand,
            phase: 0,
            protein_offset: 0,
            transcript_id: None,
            cds_segments: Vec::new(),
            biotype: crate::variants::Biotype::ProteinCoding,
        }
    }

    fn entry(gene: &str, written: &str, label: &str) -> CatalogueEntry {
        CatalogueEntry {
            gene: gene.to_string(),
            written: written.to_string(),
            variant: crate::variants::catalogue::parse_variant(written).expect("parses"),
            label: label.to_string(),
        }
    }

    #[test]
    fn only_single_residue_changes_name_a_codon() {
        assert_eq!(single_residue_position("Asp94Gly"), Some(94));
        assert_eq!(single_residue_position("Gln445*"), Some(445));
        assert_eq!(single_residue_position("Ser450Ser"), Some(450));
        assert_eq!(single_residue_position("Ala11fs"), None);
        assert_eq!(single_residue_position("Met1?"), None);
        assert_eq!(single_residue_position("Asp94"), None);
    }

    #[test]
    fn codons_are_placed_in_ribosome_order_on_both_strands() {
        let catalogue = Catalogue {
            entries: vec![
                entry("plus", "p.Asp2Gly", "R"),
                entry("plus", "p.Asp2Asn", "R"),
                entry("minus", "p.Ser1Leu", "R"),
            ],
        };
        let genes = [
            gene("plus", 11, 19, Strand::Plus),
            gene("minus", 31, 39, Strand::Minus),
        ];
        let codons = catalogue_codons(&catalogue, &genes);
        assert_eq!(codons.len(), 2, "two entries at one codon make one codon");
        assert_eq!(codons[0].positions, [14, 15, 16]);
        assert_eq!(codons[0].entries, vec![0, 1]);
        assert_eq!(codons[1].positions, [39, 38, 37]);
    }

    /// The case the reader exists for: two minority changes in one codon on
    /// different molecules. Each triplet is reported as the amino acid its
    /// molecules carry, and both match the catalogue; the merged codon, which
    /// no molecule carries, is not reported at all.
    #[test]
    fn a_trans_codon_gives_each_population_its_own_amino_acid() {
        // plus gene at 1..9: codon 2 is positions 4..6, reference GAC (Asp).
        let reference = "ATGGACTAA";
        let catalogue = Catalogue {
            entries: vec![
                entry("g", "p.Asp2Gly", "LEV R"),
                entry("g", "p.Asp2Asn", "LEV R"),
                entry("g", "p.Asp2Ser", "made up"),
            ],
        };
        let genes = [gene("g", 1, 9, Strand::Plus)];
        let codon = &catalogue_codons(&catalogue, &genes)[0];
        let triplets = CodonTriplets {
            counts: BTreeMap::from([(*b"GAC", 80), (*b"GGC", 12), (*b"AAC", 8)]),
            spanning: 100,
            partial: 3,
        };
        let floor = AlleleFloor {
            min_frequency: 0.05,
            min_molecules: 3,
        };
        let code = GeneticCode::new(11).expect("table 11");
        let alleles = codon_alleles(
            "chr", &genes[0], codon, &triplets, reference, code, &catalogue, floor,
        );
        let got: Vec<(&str, &str, usize)> = alleles
            .iter()
            .map(|a| (a.codon_allele.as_str(), a.aa_change.as_str(), a.molecules))
            .collect();
        assert_eq!(got, vec![("GGC", "Asp2Gly", 12), ("AAC", "Asp2Asn", 8)]);
        assert_eq!(
            alleles[0].catalogue,
            vec!["LEV R (g p.Asp2Gly)".to_string()]
        );
        assert_eq!(
            alleles[1].catalogue,
            vec!["LEV R (g p.Asp2Asn)".to_string()]
        );
        assert!(alleles.iter().all(|a| a.aa_change != "Asp2Ser"));
    }

    #[test]
    fn minus_strand_triplets_are_complemented_and_floors_apply() {
        // minus gene at 1..6: codon 1 reads 6,5,4 complemented.
        // Genomic 4..6 = "GTC" -> reading 6,5,4 = C,T,G -> complement GAC (Asp).
        let reference = "AAAGTC";
        let catalogue = Catalogue {
            entries: vec![entry("m", "p.Asp1Gly", "R")],
        };
        let genes = [gene("m", 1, 6, Strand::Minus)];
        let codon = &catalogue_codons(&catalogue, &genes)[0];
        assert_eq!(codon.positions, [6, 5, 4]);
        // Molecules read genomic bases in positions order: C,T,G is reference;
        // C,C,G complements to GGC (Gly).
        let triplets = CodonTriplets {
            counts: BTreeMap::from([(*b"CTG", 95), (*b"CCG", 4), (*b"CAG", 1)]),
            spanning: 100,
            partial: 0,
        };
        let code = GeneticCode::new(11).expect("table 11");
        let floor = AlleleFloor {
            min_frequency: 0.02,
            min_molecules: 2,
        };
        let alleles = codon_alleles(
            "chr", &genes[0], codon, &triplets, reference, code, &catalogue, floor,
        );
        assert_eq!(
            alleles.len(),
            1,
            "the single-molecule triplet is below the floor"
        );
        assert_eq!(alleles[0].reference_codon, "GAC");
        assert_eq!(alleles[0].codon_allele, "GGC");
        assert_eq!(alleles[0].aa_change, "Asp1Gly");
        assert_eq!(alleles[0].catalogue.len(), 1);
    }
}
