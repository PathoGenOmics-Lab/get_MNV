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

def check_genes(list_snp: List[Tuple[int, str]], gene_file: str) -> List[str]:
    '''
    Extract genes containing SNPs from an analyzed VCF.

    Parameters:
    - list_snp (List[Tuple[int, str]]): List of SNPs from a VCF file, where each SNP is represented by (position, base).
    - gene_file (str): File with genes and their coordinates.

    Returns:
    - List[str]: List of genes of interest containing SNPs from the VCF file.
    '''
    
    analyze_genelist = set()
    list_snp.sort(key=lambda x: x[0])  # Sort SNP list by position
    snp_index = 0
    snp_count = len(list_snp)
    
    with open(gene_file, 'r') as in_file:
        for line in in_file:
            fields = line.strip().split('\t')
            start_gene, end_gene = int(fields[1]), int(fields[2])
            
            while snp_index < snp_count:
                snp_position = int(list_snp[snp_index][0])
                if snp_position > end_gene:
                    break
                elif start_gene <= snp_position <= end_gene:
                    analyze_genelist.add(line)
                snp_index += 1

    return list(analyze_genelist)

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

def process_codons_based_on_strand(codon_info: CodonInfo, strand: str = 'positive') -> List[List[str]]:
    output_list = []

    for snp in codon_info.codon_list:
        codon_info.new_codon[snp.index] = snp.base

    if strand == 'negative':
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