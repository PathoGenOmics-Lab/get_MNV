import argparse
from Bio import SeqIO
from Bio.Seq import Seq

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
