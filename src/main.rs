use bio::io::fasta;
use clap::{Arg, Command, value_parser};
use env_logger::Env;
use log::{info, warn};
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::bam::{IndexedReader, Read as BamReadTrait};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::error::Error;
use std::collections::HashSet;
use protein_translate::translate;
use csv::WriterBuilder;

#[derive(Debug, Serialize, Deserialize, Clone)]
struct CodonInfo {
    codon_list: Vec<Snp>,
    original_codon: String,
    gene_name: String,
    gene_start: usize,
    gene_end: usize,
    codon_start: usize,
    codon_end: usize,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Snp {
    index: usize,
    position: usize,
    base: String,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct VariantInfo {
    gene: String,
    positions: Vec<usize>,
    base_changes: Vec<String>,
    aa_changes: Vec<String>,
    snp_aa_changes: Vec<String>,  
    variant_type: String, 
    change_type: String,  
    snp_reads: Option<Vec<usize>>,
    mnv_reads: Option<usize>,
}

fn reference_fasta(fasta_file: &str) -> Result<(String, String), Box<dyn std::error::Error>> {
    info!("Reading FASTA file: {}", fasta_file);
    let reader = fasta::Reader::new(File::open(fasta_file)?);
    let record = reader.records().next().ok_or("No records found")??;
    let chrom_name = record.id().to_string(); 
    let sequence = String::from_utf8(record.seq().to_vec())?;
    Ok((chrom_name, sequence))
}

fn getseq_posbase(vcf_file: &str) -> Result<Vec<(usize, String)>, Box<dyn std::error::Error>> {
    info!("Reading VCF file: {}", vcf_file);
    let mut vcf = bcf::Reader::from_path(vcf_file)?;
    let mut positions = Vec::new();

    for result in vcf.records() {
        let record = result?;
        let pos = record.pos() as usize + 1;
        let alt_bases: Vec<String> = record
            .alleles()
            .iter()
            .skip(1)
            .map(|allele| String::from_utf8(allele.to_vec()).unwrap())
            .collect();
        
        if let Some(alt_base) = alt_bases.get(0) {
            positions.push((pos, alt_base.clone()));
        }
    }

    Ok(positions)
}

fn check_genes(list_snp: &[(usize, String)], gene_file: &str) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    info!("Checking genes in the file: {}", gene_file);
    let file = File::open(gene_file)?;
    let reader = BufReader::new(file);
    let mut genes_of_interest = Vec::new();

    for gene_entry in reader.lines() {
        let gene_entry = gene_entry?;
        let fields: Vec<&str> = gene_entry.split('\t').collect();
        if fields.len() < 4 {
            continue; 
        }
        let start_gene: usize = fields[1].parse()?;
        let end_gene: usize = fields[2].parse()?;

        let relevant_snps: Vec<&(usize, String)> = list_snp
            .iter()
            .filter(|&&(pos, _)| pos >= start_gene && pos <= end_gene)
            .collect();

        if !relevant_snps.is_empty() {
            genes_of_interest.push(gene_entry);
        }
    }

    Ok(genes_of_interest)
}

fn determine_change_type(aa_change: &str) -> String {

    let original_aa = &aa_change[..3];
    let mutated_aa = &aa_change[aa_change.len()-3..];
    
    if original_aa == mutated_aa {
        "Synonymous".to_string()
    } else if mutated_aa == "*" || original_aa == "*" {
        "Stop gained".to_string()
    } else {
        "Non-synonymous".to_string()
    }
}

fn iupac_aa(codon: &str) -> String {
    let aa_code: HashMap<char, &str> = [
        ('A', "Ala"), ('C', "Cys"), ('D', "Asp"), ('E', "Glu"), ('F', "Phe"),
        ('G', "Gly"), ('H', "His"), ('I', "Ile"), ('K', "Lys"), ('L', "Leu"),
        ('M', "Met"), ('N', "Asn"), ('P', "Pro"), ('Q', "Gln"), ('R', "Arg"),
        ('S', "Ser"), ('T', "Thr"), ('V', "Val"), ('W', "Trp"), ('Y', "Tyr"),
        ('*', "*")
    ].iter().cloned().collect();

    if codon.len() < 3 {
        return "Invalid_codon".to_string();
    }

    let ref_aa = aa_code.get(&codon.chars().next().unwrap()).unwrap_or(&"X");
    let chg_aa = aa_code.get(&codon.chars().last().unwrap()).unwrap_or(&"X");
    
    format!("{}{}{}", ref_aa, &codon[1..codon.len()-1], chg_aa)
}

fn reverse_complement(seq: &str) -> String {
    let revcomp = bio::alphabets::dna::revcomp(seq.as_bytes());
    String::from_utf8(revcomp).unwrap()
}

fn process_translate(seq: &[u8]) -> String {
    translate(seq)
}

fn process_codons_based_on_strand(codon_info: CodonInfo, strand: &str) -> Vec<VariantInfo> {
    let mut results = Vec::new();
    let mut mutable_codon: Vec<char> = codon_info.original_codon.chars().collect();

    for snp in &codon_info.codon_list {
        if snp.position < codon_info.codon_start || snp.position > codon_info.codon_end {
            continue; 
        }
        mutable_codon[snp.position - codon_info.codon_start] = snp.base.chars().next().unwrap();
    }

    let original_aa = if strand == "-" {
        process_translate(reverse_complement(&codon_info.original_codon).as_bytes())
    } else {
        process_translate(codon_info.original_codon.as_bytes())
    };

    let mutated_codon: String = if strand == "-" {
        reverse_complement(&mutable_codon.iter().collect::<String>())
    } else {
        mutable_codon.iter().collect()
    };

    let mutated_aa = process_translate(mutated_codon.as_bytes());

    let mut positions = Vec::new();
    let mut base_changes = Vec::new();

    for snp in &codon_info.codon_list {
        positions.push(snp.position);
        base_changes.push(snp.base.clone());
    }

    let aa_position = if strand == "+" {
        (positions[0] - codon_info.gene_start) / 3 + 1
    } else {
        (codon_info.gene_end - positions[0]) / 3 + 1
    };
    let chg_aa_combined = format!("{}{}{}", original_aa, aa_position, mutated_aa);
    let chg_aa_combined = iupac_aa(&chg_aa_combined);

    let change_type = determine_change_type(&chg_aa_combined);

    let mut snp_aa_changes = Vec::new();
    for snp in &codon_info.codon_list {
        let mut single_mutable_codon = codon_info.original_codon.chars().collect::<Vec<char>>();
        single_mutable_codon[snp.position - codon_info.codon_start] = snp.base.chars().next().unwrap();
        let single_mutated_codon: String = if strand == "-" {
            reverse_complement(&single_mutable_codon.iter().collect::<String>())
        } else {
            single_mutable_codon.iter().collect()
        };
        let single_mutated_aa = process_translate(single_mutated_codon.as_bytes());
        let single_chg_aa = format!("{}{}{}", original_aa, aa_position, single_mutated_aa);
        let single_chg_aa = iupac_aa(&single_chg_aa);
        snp_aa_changes.push(single_chg_aa);
    }

    let variant_type = if codon_info.codon_list.len() == 1 {
        "SNP".to_string()
    } else {
        "MNV".to_string()
    };

    results.push(VariantInfo {
        gene: codon_info.gene_name.clone(),
        positions,
        base_changes,
        aa_changes: vec![chg_aa_combined], 
        snp_aa_changes,
        variant_type, 
        change_type,  
        snp_reads: None, 
        mnv_reads: None, 
    });

    results
}

fn get_mnv_variants(gene_list: &[String], snp_list: &[(usize, String)], sequence: &str) -> Vec<VariantInfo> {
    let mut results = Vec::new();

    for gene_info in gene_list {
        let details: Vec<&str> = gene_info.split('\t').collect();
        if details.len() < 4 {
            continue; 
        }
        let gene_name = details[0].to_string();
        let strand = details[3];
        let gene_start: usize = match details[1].parse() {
            Ok(val) => val,
            Err(_) => continue,
        };
        let gene_end: usize = match details[2].parse() {
            Ok(val) => val,
            Err(_) => continue,
        };

        let relevant_snps: Vec<(usize, String)> = snp_list
            .iter()
            .filter(|&&(pos, _)| pos >= gene_start && pos <= gene_end)
            .cloned()
            .collect();

        
        for position in (gene_start..=gene_end).step_by(3) {
            let codon_start = position;
            let codon_end = position + 2;
            let codon_snps: Vec<Snp> = relevant_snps
                .iter()
                .filter(|&&(pos, _)| pos >= codon_start && pos <= codon_end)
                .map(|&(pos, ref base)| Snp {
                    index: pos, 
                    position: pos,
                    base: base.clone(),
                })
                .collect();

            if !codon_snps.is_empty() {
                if codon_end <= sequence.len() {
                    let codon_sequence = &sequence[codon_start - 1..codon_end];
                    let info = CodonInfo {
                        codon_list: codon_snps,
                        original_codon: codon_sequence.to_string(),
                        gene_name: gene_name.clone(),
                        gene_start,
                        gene_end,
                        codon_start,
                        codon_end,
                    };
                    results.extend(process_codons_based_on_strand(info, strand));
                } else {
                    warn!("The gene's codon {} exceeds the length of the reference sequence.", gene_name);
                }
            }
        }
    }

    results
}

fn count_reads_per_position(
    bam_file: &str,
    chrom: &str,
    variant_positions: &[usize],
    alt_bases: &[String],
    min_phred_quality: u8,
) -> Result<(Vec<usize>, usize), Box<dyn Error>> {
    let mut snp_counts = vec![0; variant_positions.len()]; 
    let mut mnv_count = 0;
    let mut bam = IndexedReader::from_path(bam_file)?;
    let mut unique_mnv_reads = HashSet::new();
    let mut unique_snp_reads_per_position: Vec<HashSet<Vec<u8>>> = vec![HashSet::new(); variant_positions.len()];

    let region_start = (variant_positions[0] - 1) as i64; 
    let region_end = variant_positions[variant_positions.len() - 1] as i64;

    if let Err(e) = bam.fetch((chrom, region_start, region_end)) {
        warn!("Error when searching for the region in the BAM: {}", e);
        return Ok((snp_counts, mnv_count));
    }

    for result in bam.records() {
        let record = result?;
        let read_id = record.qname().to_vec(); 
        let read_pos = record.pos() as usize + 1; 
        let read_length = record.seq().len();
        let read_end = read_pos + read_length - 1;


        if record.is_duplicate() {
            continue;
        }

        if record.mapq() < min_phred_quality {
            continue;
        }

        if variant_positions.iter().all(|&pos| pos >= read_pos && pos <= read_end) {
            let mut matches_mnv = true;

            for (&pos, alt_base) in variant_positions.iter().zip(alt_bases.iter()) {
                let read_base_index = pos - read_pos;

                if read_base_index >= read_length {
                    matches_mnv = false;
                    break;
                }

                let read_base = record.seq().as_bytes()[read_base_index] as char;

                let base_quality = record.qual()[read_base_index];
                if base_quality < min_phred_quality {
                    matches_mnv = false;
                    break;
                }

                if read_base.to_ascii_uppercase() != alt_base.chars().next().unwrap().to_ascii_uppercase() {
                    matches_mnv = false;
                    break;
                }
            }

            if matches_mnv {
                if !unique_mnv_reads.contains(&read_id) {
                    mnv_count += 1;
                    unique_mnv_reads.insert(read_id.clone());
                }
                continue; 
            }
        }

        let mut valid_snp_read = true;

        for (&pos, alt_base) in variant_positions.iter().zip(alt_bases.iter()) {
            if pos >= read_pos && pos <= read_end {
                let read_base_index = pos - read_pos;

                if read_base_index < read_length {
                    let read_base = record.seq().as_bytes()[read_base_index] as char;
                    let base_quality = record.qual()[read_base_index];

                    if base_quality < min_phred_quality || read_base.to_ascii_uppercase() != alt_base.chars().next().unwrap().to_ascii_uppercase() {
                        valid_snp_read = false;
                        break;
                    }
                }
            }
        }

        if valid_snp_read {
            for (i, &pos) in variant_positions.iter().enumerate() {
                if pos >= read_pos && pos <= read_end {
                    if !unique_snp_reads_per_position[i].contains(&read_id) {
                        snp_counts[i] += 1;
                        unique_snp_reads_per_position[i].insert(read_id.clone());
                    }
                }
            }
        }
    }


    Ok((snp_counts, mnv_count))
}

fn write_to_tsv_with_reads(
    data_list: Vec<VariantInfo>,
    filename: &str,
    bam_file: Option<&String>,
    chrom_name: &str,
    min_phred_quality: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(format!("{}.MNV.tsv", filename))?;

    // Escribir encabezado condicionalmente
    if bam_file.is_some() {
        writer.write_record(&["Gene", "Positions", "Base Changes", "AA Changes", "SNP AA Changes", "Variant Type", "Change Type", "SNP Reads", "MNV Reads"])?;
    } else {
        writer.write_record(&["Gene", "Positions", "Base Changes", "AA Changes", "SNP AA Changes", "Variant Type", "Change Type"])?;
    }

    for variant in data_list {
        if let Some(bam_path) = bam_file {
            let (snp_counts, mnv_reads) = count_reads_per_position(
                bam_path,
                chrom_name,
                &variant.positions,
                &variant.base_changes,
                min_phred_quality
            )?;

            let mut variant = variant.clone();
            variant.snp_reads = Some(snp_counts.clone());
            variant.mnv_reads = Some(mnv_reads);

            if variant.snp_reads.as_ref().unwrap().iter().any(|&count| count > 0) && variant.mnv_reads.unwrap() > 0 {
                variant.variant_type = "SNP/MNV".to_string();
            }

            let positions_str = variant.positions.iter().map(|pos| pos.to_string()).collect::<Vec<String>>().join(",");
            let base_changes_str = variant.base_changes.join(","); 
            let aa_changes_str = variant.aa_changes.join("; ");
            let snp_aa_changes_str = variant.snp_aa_changes.join("; ");
            let change_type_str = variant.change_type;
            let snp_reads_output = variant.snp_reads.as_ref().unwrap().iter().map(|&count| count.to_string()).collect::<Vec<String>>().join(",");
            let mnv_reads_output = variant.mnv_reads.unwrap().to_string();

            writer.write_record(&[
                &variant.gene,
                &positions_str,
                &base_changes_str,
                &aa_changes_str,
                &snp_aa_changes_str,
                &variant.variant_type,
                &change_type_str,
                &snp_reads_output,
                &mnv_reads_output,
            ])?;
        } else {
            let positions_str = variant.positions.iter().map(|pos| pos.to_string()).collect::<Vec<String>>().join(",");
            let base_changes_str = variant.base_changes.join(","); 
            let aa_changes_str = variant.aa_changes.join("; ");
            let snp_aa_changes_str = variant.snp_aa_changes.join("; ");
            let change_type_str = variant.change_type;

            writer.write_record(&[
                &variant.gene,
                &positions_str,
                &base_changes_str,
                &aa_changes_str,
                &snp_aa_changes_str,
                &variant.variant_type,
                &change_type_str,
            ])?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn get_base_name(file_path: &str) -> String {
    Path::new(file_path)
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    let matches = Command::new("get_mnv")
        .version("1.0.0")
        .author("Paula Ruiz Rodriguez")
        .about("Identifies multiple SNVs within the same codon, reclassifies them as MNVs, and accurately computes resulting amino acid changes from genomic reads")
        .arg(Arg::new("vcf")
            .short('v')
            .long("vcf")
            .value_name("VCF_FILE")
            .help("VCF file with SNPs")
            .required(true))
        .arg(Arg::new("bam")
            .short('b')
            .long("bam")
            .value_name("BAM_FILE")
            .help("BAM file with aligned reads (optional)")
            .required(false))
        .arg(Arg::new("fasta")
            .short('f')
            .long("fasta")
            .value_name("FASTA_FILE")
            .help("FASTA file with reference sequence")
            .required(true))
        .arg(Arg::new("genes")
            .short('g')
            .long("genes")
            .value_name("GENES_FILE")
            .help("File with gene information")
            .required(true))
        .arg(Arg::new("quality")
            .short('q')
            .long("quality")
            .value_name("QUALITY")
            .help("Minimum Phred quality score (default: 20)")
            .num_args(1)
            .default_value("20") 
            .value_parser(value_parser!(u8))) 
        .get_matches();

    let vcf_file = matches.get_one::<String>("vcf").unwrap();
    let bam_file = matches.get_one::<String>("bam").cloned(); 
    let fasta_file = matches.get_one::<String>("fasta").unwrap();
    let genes_file = matches.get_one::<String>("genes").unwrap();
    
    let min_phred_quality: u8 = *matches.get_one::<u8>("quality").unwrap();

    info!("VCF file: {}", vcf_file);
    if let Some(ref bam) = bam_file {
        info!("BAM file: {}", bam);
    } else {
        info!("No BAM file provided. Read counts will be omitted.");
    }
    info!("FASTA file: {}", fasta_file);
    info!("Gene file: {}", genes_file);
    info!("Minimum Phred quality: {}", min_phred_quality);

    let base_name = get_base_name(vcf_file);
    
    let (chrom_name, sequence) = reference_fasta(fasta_file)?;

    info!("Extracted chromosome name: {}", chrom_name);

    let snp_list = getseq_posbase(vcf_file)?;
    let gene_list = check_genes(&snp_list, genes_file)?;

    let mut mnv_variants = get_mnv_variants(&gene_list, &snp_list, &sequence);

    if let Some(ref bam) = bam_file {
        for variant in &mut mnv_variants {
            let (snp_counts, mnv_reads) = count_reads_per_position(
                bam,
                &chrom_name,
                &variant.positions,
                &variant.base_changes,
                min_phred_quality
            )?;
            variant.snp_reads = Some(snp_counts);
            variant.mnv_reads = Some(mnv_reads);
        }
    }

    write_to_tsv_with_reads(mnv_variants, &base_name, bam_file.as_ref(), &chrom_name, min_phred_quality)?;
    
    Ok(())
}
