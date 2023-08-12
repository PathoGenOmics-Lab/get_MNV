import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from typing import List, Tuple, Union
from collections import namedtuple

def reference_fasta(fasta_file: str = 'MTB_ancestor.fas') -> str:
    '''
    Parse a FASTA file and return the sequence of the first record.

    Parameters:
    - fasta_file (str): Path to the FASTA file. Defaults to 'MTB_ancestor.fas'.

    Returns:
    - str: Sequence of the first record in the FASTA file.
    '''
    with open(fasta_file, 'r') as file:
        fasta_sequence = next(SeqIO.parse(file, 'fasta'))
        return str(fasta_sequence.seq)
    
def read_genes_names(genes_file: str) -> list[str]:
    '''
    Extract gene names from a tab-separated file.
    
    Parameters:
    - genes_file (str): Path to the text file with gene names in the first column.
    
    Returns:
    - list[str]: List containing gene names.
    '''
    with open(genes_file, 'r') as in_file:
        return [line.split('\t')[0].strip() for line in in_file]

def getseq_posbase(vcf_file: str = 'G35894.var.snp.vcf') -> List[List[str]]:
    '''
    Extract positions and alternative base from an annotated VCF file.

    Parameters:
    - vcf_file (str): Path to the annotated VCF file with snpEff. Default is 'G35894.var.snp.vcf'.

    Returns:
    - List[List[str]]: List of lists containing positions and the corresponding alternative base.

    Notes:
    - Ignores intergenic positions.
    - Assumes the VCF format where the second column is the position and the fifth column is the alternative base.
    '''
    
    with open(vcf_file, 'r') as in_file:
        return [
            [fields[1], fields[4]]
            for fields in (line.strip().split('\t') for line in in_file if not line.startswith('#') and 'intergenic' not in line)
        ]

def check_genes(list_snp: list, gene_file: str) -> list:
    '''
    Function to get genes containing SNPs in the analyzed vcf.
        Input  -> list_snps: list with SNPs of vcf file, gene_file: 
                  genes with coordinates.
        Output -> list with genes of interest containing SNPs of vcf file.
    '''
    genes_of_interest = set()

    with open(gene_file, 'r') as in_file:
        for gene_entry in in_file:
            chromosome, start_gene, end_gene, *rest = gene_entry.strip('\n').split('\t')
            start_gene, end_gene = int(start_gene), int(end_gene)

            relevant_snps = [snp for snp in list_snp if start_gene <= int(snp[0]) <= end_gene]

            if relevant_snps:
                genes_of_interest.add(gene_entry)

    return list(genes_of_interest)

# Define a structure to group related parameters
CodonInfo = namedtuple("CodonInfo", ["codon_list", "new_codon", "original_codon", "amino_acid", "gene"])
SNP = namedtuple("SNP", ["index", "position", "base"])

def iupac_aa(codon: str):
    '''
    Function to translate iupac code to three aa code letter  
        Input  -> aa string with 1 letter     -->   refPOSalt  
        Output -> aa string with 3 letters    -->   refPOSalt  
    '''
    aa_code = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*':'*'
    }

    ref = aa_code.get(codon[0])
    chg = aa_code.get(codon[-1])
    if ref == chg:
        status = 'synonymous_variant'
    elif chg == '*':
        status = 'stop_gained'
    else:
        status = 'missense_variant'
    cod_codon = ''.join([aa_code.get(codon[0]), codon[1:-1], aa_code.get(codon[-1])]) + '\t' + status
    return cod_codon

def process_codons_based_on_strand(codon_info: CodonInfo, strand: str = '+') -> List[List[str]]:
    output_list = []

    for snp in codon_info.codon_list:
        print(snp)
        codon_info.new_codon[snp.index] = snp.base

    if strand == '-':
        translated_codon = Seq(''.join(codon_info.new_codon)).reverse_complement()
    else:
        translated_codon = Seq(''.join(codon_info.new_codon))

    my_newaa = translated_codon.translate()

    chg_aa = ''.join([str(codon_info.amino_acid), str(codon_info.original_codon), str(my_newaa)])
    chg_aa = iupac_aa(chg_aa)

    for snp in codon_info.codon_list:
        sentence = '\t'.join([codon_info.gene, snp.position, snp.base, chg_aa]) + '\n'
        
        for entry in output_list:
            pos = entry[0].split('\t')[1]
            if snp.position == pos:
                entry.append(sentence)
                break
        else:
            output_list.append([sentence])

    return output_list

def extract_codon_snp(codon_start: int, codon_end: int, snp_list: list):
    """Extract SNPs that fall within a specific codon."""
    return [snp for snp in snp_list if codon_start <= int(snp[0]) <= codon_end]

def get_mnv_variants(gene_list: list, snp_list: list, sequence: str):
    """Identify and process multiple nucleotide variants (MNVs) from a given SNP list."""
    
    results = []
    for gene_info in gene_list:
        details = gene_info.strip('\n').split('\t')
        gene_name, strand, gene_start, gene_end = details[0], details[3], int(details[1]), int(details[2])
        print(gene_name, strand, gene_start, gene_end)
        relevant_snps = [SNP(i, *snp) for i, snp in enumerate(snp_list) if gene_start <= int(snp[0]) <= gene_end]
        print(relevant_snps)
        for position in range(gene_start, gene_end, 3):
            codon_start, codon_end = position, position + 2
            codon_snps = [snp for snp in relevant_snps if codon_start <= int(snp.position) <= codon_end]
            if len(codon_snps) > 1:  # If multiple SNPs are in the same codon
                codon_sequence = Seq(sequence[codon_start-1:codon_end])
                
                info = CodonInfo(codon_snps, codon_sequence, None, None, gene_name)
                #print(info,strand)
                #results.extend(process_codons_based_on_strand(info, strand))

    return results

parser = argparse.ArgumentParser(description = 'script to annotate MNV') 
parser.add_argument('-v', dest = 'vcf', required =True, help = 'Vcf file with snps')
parser.add_argument('-f', dest = 'fasta', required =True, help = 'Name of the reference fasta') 
parser.add_argument('-g', dest = 'genes', required = True, help = 'File with gene info')
args = parser.parse_args()

sequence = reference_fasta(args.fasta)
lista_snp = getseq_posbase(args.vcf)
gene_list = check_genes(lista_snp,args.genes)
print(len(gene_list))
#get_mnv_variants(gene_list, lista_snp, sequence)

#python3 new.py -v G35894.var.snp.vcf -f MTB_ancestor.fas -g anot_genes.txt