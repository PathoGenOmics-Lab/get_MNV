use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum VariantType {
    Snp,
    Mnv,
    SnpMnv,
    Indel,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum ChangeType {
    Unknown,
    Synonymous,
    NonSynonymous,
    StopGained,
    StopLost,
    IndelOverlap,
    FrameshiftSynonymous,
    FrameshiftNonSynonymous,
    FrameshiftStopGained,
    FrameshiftStopLost,
    FrameshiftUnknown,
    FrameshiftIndel,
    InFrameIndel,
}

impl ChangeType {
    pub fn from_label(label: &str) -> Self {
        match label {
            "Synonymous" => ChangeType::Synonymous,
            "Non-synonymous" => ChangeType::NonSynonymous,
            "Stop gained" => ChangeType::StopGained,
            "Stop lost" => ChangeType::StopLost,
            "Unknown" => ChangeType::Unknown,
            "Indel overlap" => ChangeType::IndelOverlap,
            "Synonymous (frameshift)" => ChangeType::FrameshiftSynonymous,
            "Non-synonymous (frameshift)" => ChangeType::FrameshiftNonSynonymous,
            "Stop gained (frameshift)" => ChangeType::FrameshiftStopGained,
            "Stop lost (frameshift)" => ChangeType::FrameshiftStopLost,
            "Unknown (frameshift)" => ChangeType::FrameshiftUnknown,
            "Frameshift Indel" => ChangeType::FrameshiftIndel,
            "In-frame Indel" => ChangeType::InFrameIndel,
            _ => ChangeType::Unknown,
        }
    }

    pub fn with_frameshift(self) -> Self {
        match self {
            ChangeType::Synonymous => ChangeType::FrameshiftSynonymous,
            ChangeType::NonSynonymous => ChangeType::FrameshiftNonSynonymous,
            ChangeType::StopGained => ChangeType::FrameshiftStopGained,
            ChangeType::StopLost => ChangeType::FrameshiftStopLost,
            ChangeType::Unknown => ChangeType::FrameshiftUnknown,
            ChangeType::IndelOverlap => ChangeType::IndelOverlap,
            ChangeType::FrameshiftSynonymous => ChangeType::FrameshiftSynonymous,
            ChangeType::FrameshiftNonSynonymous => ChangeType::FrameshiftNonSynonymous,
            ChangeType::FrameshiftStopGained => ChangeType::FrameshiftStopGained,
            ChangeType::FrameshiftStopLost => ChangeType::FrameshiftStopLost,
            ChangeType::FrameshiftUnknown => ChangeType::FrameshiftUnknown,
            ChangeType::FrameshiftIndel => ChangeType::FrameshiftIndel,
            ChangeType::InFrameIndel => ChangeType::InFrameIndel,
        }
    }
}

impl VariantType {
    pub fn as_str(self) -> &'static str {
        match self {
            VariantType::Snp => "SNP",
            VariantType::Mnv => "MNV",
            VariantType::SnpMnv => "SNP/MNV",
            VariantType::Indel => "INDEL",
        }
    }
}

impl Display for VariantType {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.write_str(self.as_str())
    }
}

impl Display for ChangeType {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        let value = match self {
            ChangeType::Unknown => "Unknown",
            ChangeType::Synonymous => "Synonymous",
            ChangeType::NonSynonymous => "Non-synonymous",
            ChangeType::StopGained => "Stop gained",
            ChangeType::StopLost => "Stop lost",
            ChangeType::IndelOverlap => "Indel overlap",
            ChangeType::FrameshiftSynonymous => "Synonymous (frameshift)",
            ChangeType::FrameshiftNonSynonymous => "Non-synonymous (frameshift)",
            ChangeType::FrameshiftStopGained => "Stop gained (frameshift)",
            ChangeType::FrameshiftStopLost => "Stop lost (frameshift)",
            ChangeType::FrameshiftUnknown => "Unknown (frameshift)",
            ChangeType::FrameshiftIndel => "Frameshift Indel",
            ChangeType::InFrameIndel => "In-frame Indel",
        };
        f.write_str(value)
    }
}

impl FromStr for Strand {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Plus),
            "-" => Ok(Strand::Minus),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Gene {
    pub name: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Snp {
    pub index: usize,
    pub position: usize,
    pub ref_base: String,
    pub base: String,
    pub original_dp: Option<usize>,
    pub original_freq: Option<f64>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CodonInfo {
    pub codon_list: Vec<Snp>,
    pub original_codon: String,
    pub gene_name: String,
    pub gene_start: usize,
    pub gene_end: usize,
    pub codon_start: usize,
    pub codon_end: usize,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct VariantInfo {
    pub chrom: String,
    pub gene: String,
    pub positions: Vec<usize>,
    pub ref_bases: Vec<String>,
    pub base_changes: Vec<String>,
    pub aa_changes: Vec<String>,
    pub snp_aa_changes: Vec<String>,
    pub variant_type: VariantType,
    pub change_type: ChangeType,
    pub snp_reads: Option<Vec<usize>>,
    pub snp_forward_reads: Option<Vec<usize>>,
    pub snp_reverse_reads: Option<Vec<usize>>,
    pub mnv_reads: Option<usize>,
    pub mnv_forward_reads: Option<usize>,
    pub mnv_reverse_reads: Option<usize>,
    pub mnv_total_reads: Option<usize>,
    pub total_reads: Option<Vec<usize>>,
    pub total_forward_reads: Option<Vec<usize>>,
    pub total_reverse_reads: Option<Vec<usize>>,
    pub mnv_total_forward_reads: Option<usize>,
    pub mnv_total_reverse_reads: Option<usize>,
    pub ref_codon: Option<String>,
    pub snp_codon: Option<String>,
    pub mnv_codon: Option<String>,
    pub original_dp: Option<Vec<usize>>,
    pub original_freq: Option<Vec<f64>>,
}

use crate::utils::{determine_change_type, iupac_aa, process_translate, reverse_complement};

fn collect_all_usize(values: impl Iterator<Item = Option<usize>>) -> Option<Vec<usize>> {
    let mut out = Vec::new();
    for value in values {
        out.push(value?);
    }
    Some(out)
}

fn collect_all_f64(values: impl Iterator<Item = Option<f64>>) -> Option<Vec<f64>> {
    let mut out = Vec::new();
    for value in values {
        out.push(value?);
    }
    Some(out)
}

fn construct_codon(codon_info: &CodonInfo, target_snps: &[&Snp]) -> String {
    let mut codon = String::new();
    for i in 0..3 {
        let current_pos = codon_info.codon_start + i;
        if let Some(snp) = target_snps.iter().find(|&&s| s.position == current_pos) {
            codon.push_str(&snp.base);
        } else {
            codon.push(codon_info.original_codon.chars().nth(i).unwrap_or('N'));
        }
    }
    codon
}

pub fn process_codon(codon_info: CodonInfo, strand: Strand, chrom: &str) -> VariantInfo {
    let ref_codon = codon_info.original_codon.clone();

    let mnv_snps: Vec<&Snp> = codon_info.codon_list.iter().collect();
    let mnv_codon = construct_codon(&codon_info, &mnv_snps);

    let snp_codons: Vec<String> = codon_info
        .codon_list
        .iter()
        .map(|snp| construct_codon(&codon_info, &[snp]))
        .collect();
    let snp_codon = snp_codons.join(" ; ");

    let (orig_aa, mut_aa) = match strand {
        Strand::Minus => (
            process_translate(reverse_complement(&ref_codon).as_bytes()),
            process_translate(reverse_complement(&mnv_codon).as_bytes()),
        ),
        Strand::Plus => (
            process_translate(ref_codon.as_bytes()),
            process_translate(mnv_codon.as_bytes()),
        ),
    };

    let aa_pos = match strand {
        Strand::Plus => (codon_info.codon_list[0].position - codon_info.gene_start) / 3 + 1,
        Strand::Minus => (codon_info.gene_end - codon_info.codon_list[0].position) / 3 + 1,
    };

    let combined_change = format!("{}{}{}", orig_aa, aa_pos, mut_aa);
    let combined_aa = iupac_aa(&combined_change);
    let change_type = ChangeType::from_label(&determine_change_type(&combined_change));

    let snp_changes: Vec<String> = codon_info
        .codon_list
        .iter()
        .map(|snp| {
            let single_codon = construct_codon(&codon_info, &[snp]);
            let single = match strand {
                Strand::Minus => reverse_complement(&single_codon),
                Strand::Plus => single_codon,
            };
            let single_aa = process_translate(single.as_bytes());
            iupac_aa(&format!("{}{}{}", orig_aa, aa_pos, single_aa))
        })
        .collect();

    VariantInfo {
        chrom: chrom.to_string(),
        gene: codon_info.gene_name,
        positions: codon_info.codon_list.iter().map(|s| s.position).collect(),
        ref_bases: codon_info
            .codon_list
            .iter()
            .map(|s| s.ref_base.clone())
            .collect(),
        base_changes: codon_info
            .codon_list
            .iter()
            .map(|s| s.base.clone())
            .collect(),
        aa_changes: vec![combined_aa],
        snp_aa_changes: snp_changes,
        variant_type: if codon_info.codon_list.len() == 1 {
            VariantType::Snp
        } else {
            VariantType::SnpMnv
        },
        change_type,
        snp_reads: None,
        snp_forward_reads: None,
        snp_reverse_reads: None,
        mnv_reads: None,
        mnv_forward_reads: None,
        mnv_reverse_reads: None,
        mnv_total_reads: None,
        total_reads: None,
        total_forward_reads: None,
        total_reverse_reads: None,
        mnv_total_forward_reads: None,
        mnv_total_reverse_reads: None,
        ref_codon: Some(ref_codon),
        snp_codon: Some(snp_codon),
        mnv_codon: Some(mnv_codon),
        original_dp: collect_all_usize(codon_info.codon_list.iter().map(|s| s.original_dp)),
        original_freq: collect_all_f64(codon_info.codon_list.iter().map(|s| s.original_freq)),
    }
}

// Process one gene at a time to keep memory usage low.
fn codon_bounds_for_position(gene: &Gene, position: usize) -> Option<(usize, usize)> {
    if position < gene.start || position > gene.end {
        return None;
    }
    match gene.strand {
        Strand::Plus => {
            let offset = position - gene.start;
            let codon_start = gene.start + (offset / 3) * 3;
            let codon_end = codon_start + 2;
            if codon_end <= gene.end {
                Some((codon_start, codon_end))
            } else {
                None
            }
        }
        Strand::Minus => {
            let offset = gene.end - position;
            let codon_end = gene.end.saturating_sub((offset / 3) * 3);
            if codon_end < gene.start + 2 {
                return None;
            }
            let codon_start = codon_end - 2;
            if codon_start >= gene.start {
                Some((codon_start, codon_end))
            } else {
                None
            }
        }
    }
}

pub fn get_mnv_variants_for_gene(
    gene: &Gene,
    snp_list: &[crate::io::VcfPosition],
    reference: &crate::io::Reference,
    chrom: &str,
) -> Vec<VariantInfo> {
    let mut variants = Vec::new();

    let mut codon_to_snps: BTreeMap<usize, Vec<crate::variants::Snp>> = BTreeMap::new();
    let mut indels: Vec<crate::io::VcfPosition> = Vec::new();

    for variant in snp_list
        .iter()
        .filter(|snp| snp.position >= gene.start && snp.position <= gene.end)
    {
        let is_snp = variant.ref_allele.len() == 1
            && variant.alt_allele.len() == 1
            && !variant.alt_allele.starts_with('<');
        if is_snp {
            if let Some((codon_start, _)) = codon_bounds_for_position(gene, variant.position) {
                codon_to_snps
                    .entry(codon_start)
                    .or_default()
                    .push(crate::variants::Snp {
                        index: variant.position,
                        position: variant.position,
                        ref_base: variant.ref_allele.clone(),
                        base: variant.alt_allele.clone(),
                        original_dp: variant.original_dp,
                        original_freq: variant.original_freq,
                    });
            }
        } else {
            indels.push(variant.clone());
        }
    }

    if codon_to_snps.is_empty() && indels.is_empty() {
        return variants;
    }

    let mut codon_starts: Vec<usize> = codon_to_snps.keys().copied().collect();
    match gene.strand {
        Strand::Plus => codon_starts.sort_unstable(),
        Strand::Minus => codon_starts.sort_unstable_by(|a, b| b.cmp(a)),
    }

    for codon_start in codon_starts {
        let codon_snps = codon_to_snps.get(&codon_start).cloned().unwrap_or_default();
        if codon_snps.is_empty() {
            continue;
        }
        let codon_end = codon_start + 2;
        if codon_end > reference.sequence.len() {
            continue;
        }
        let mut codon_snps = codon_snps;
        codon_snps.sort_by_key(|s| s.position);

        let codon_seq = &reference.sequence[(codon_start - 1)..codon_end];

        let overlaps_indel = indels.iter().any(|indel| {
            let indel_end = indel.position + indel.ref_allele.len() - 1;
            indel.position <= codon_end && indel_end >= codon_start
        });

        let mut upstream_shift: isize = 0;
        let mut has_symbolic_sv = false;

        for indel in &indels {
            let is_upstream = match gene.strand {
                Strand::Plus => indel.position < codon_start,
                Strand::Minus => indel.position > codon_end,
            };

            if is_upstream {
                if indel.alt_allele.starts_with('<') {
                    has_symbolic_sv = true;
                } else {
                    upstream_shift +=
                        (indel.alt_allele.len() as isize) - (indel.ref_allele.len() as isize);
                }
            }
        }

        let is_frameshifted = has_symbolic_sv || upstream_shift % 3 != 0;

        let codon_info = CodonInfo {
            codon_list: codon_snps,
            original_codon: codon_seq.to_string(),
            gene_name: gene.name.clone(),
            gene_start: gene.start,
            gene_end: gene.end,
            codon_start,
            codon_end,
        };

        let mut var_info = process_codon(codon_info, gene.strand, chrom);

        if overlaps_indel {
            var_info.change_type = ChangeType::IndelOverlap;
            var_info.aa_changes = vec!["Unknown".to_string()];
            var_info.snp_aa_changes = vec!["Unknown".to_string(); var_info.snp_aa_changes.len()];
        } else if is_frameshifted {
            var_info.change_type = var_info.change_type.with_frameshift();
            var_info.aa_changes = var_info
                .aa_changes
                .into_iter()
                .map(|s| format!("{} (fs)", s))
                .collect();
            var_info.snp_aa_changes = var_info
                .snp_aa_changes
                .into_iter()
                .map(|s| format!("{} (fs)", s))
                .collect();
        }

        variants.push(var_info);
    }

    for indel in indels {
        let is_frameshift = if indel.alt_allele.starts_with('<') {
            true
        } else {
            (indel.ref_allele.len() as isize - indel.alt_allele.len() as isize) % 3 != 0
        };

        let change_type = if is_frameshift {
            ChangeType::FrameshiftIndel
        } else {
            ChangeType::InFrameIndel
        };

        variants.push(VariantInfo {
            chrom: chrom.to_string(),
            gene: gene.name.clone(),
            positions: vec![indel.position],
            ref_bases: vec![indel.ref_allele.clone()],
            base_changes: vec![indel.alt_allele.clone()],
            aa_changes: vec!["-".to_string()],
            snp_aa_changes: vec!["-".to_string()],
            variant_type: VariantType::Indel,
            change_type,
            snp_reads: None,
            snp_forward_reads: None,
            snp_reverse_reads: None,
            mnv_reads: None,
            mnv_forward_reads: None,
            mnv_reverse_reads: None,
            mnv_total_reads: None,
            total_reads: None,
            total_forward_reads: None,
            total_reverse_reads: None,
            mnv_total_forward_reads: None,
            mnv_total_reverse_reads: None,
            ref_codon: None,
            snp_codon: None,
            mnv_codon: None,
            original_dp: indel.original_dp.map(|value| vec![value]),
            original_freq: indel.original_freq.map(|value| vec![value]),
        });
    }

    variants
}

#[cfg(test)]
mod tests {
    use super::{
        codon_bounds_for_position, process_codon, ChangeType, CodonInfo, Gene, Snp, Strand,
        VariantType,
    };
    use crate::utils::reverse_complement;

    fn next_u64(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        *seed
    }

    fn random_base(seed: &mut u64) -> char {
        match next_u64(seed) % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        }
    }

    #[test]
    fn test_property_reverse_complement_involution_for_random_sequences() {
        let mut seed = 123456789u64;
        for _ in 0..200 {
            let mut seq = String::new();
            for _ in 0..30 {
                seq.push(random_base(&mut seed));
            }
            let twice = reverse_complement(&reverse_complement(&seq));
            assert_eq!(seq, twice);
        }
    }

    #[test]
    fn test_property_codon_bounds_plus_strand_cover_position() {
        let gene = Gene {
            name: "gene_plus".to_string(),
            start: 100,
            end: 399,
            strand: Strand::Plus,
        };
        for pos in gene.start..=gene.end {
            let bounds = codon_bounds_for_position(&gene, pos).expect("expected codon bounds");
            assert_eq!(bounds.1 - bounds.0 + 1, 3);
            assert!(bounds.0 <= pos && pos <= bounds.1);
            assert!(bounds.0 >= gene.start && bounds.1 <= gene.end);
        }
    }

    #[test]
    fn test_property_codon_bounds_minus_strand_cover_position() {
        let gene = Gene {
            name: "gene_minus".to_string(),
            start: 100,
            end: 399,
            strand: Strand::Minus,
        };
        for pos in gene.start..=gene.end {
            let bounds = codon_bounds_for_position(&gene, pos).expect("expected codon bounds");
            assert_eq!(bounds.1 - bounds.0 + 1, 3);
            assert!(bounds.0 <= pos && pos <= bounds.1);
            assert!(bounds.0 >= gene.start && bounds.1 <= gene.end);
        }
    }

    #[test]
    fn test_process_codon_emits_expected_variant_type_for_single_snp() {
        let codon_info = CodonInfo {
            codon_list: vec![Snp {
                index: 101,
                position: 101,
                ref_base: "T".to_string(),
                base: "C".to_string(),
                original_dp: Some(20),
                original_freq: Some(0.5),
            }],
            original_codon: "ATG".to_string(),
            gene_name: "gene1".to_string(),
            gene_start: 100,
            gene_end: 399,
            codon_start: 100,
            codon_end: 102,
        };

        let result = process_codon(codon_info, Strand::Plus, "chr1");
        assert_eq!(result.variant_type, VariantType::Snp);
        assert!(matches!(
            result.change_type,
            ChangeType::Synonymous
                | ChangeType::NonSynonymous
                | ChangeType::StopGained
                | ChangeType::StopLost
                | ChangeType::Unknown
        ));
    }
}
