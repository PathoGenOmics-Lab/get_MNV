import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from typing import List

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