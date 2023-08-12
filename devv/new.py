import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from typing import List
from collections import namedtuple
import vcf
import pandas as pd

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
CodonInfo = namedtuple("CodonInfo", ["codon_list", "new_codon", "original_codon", "gene_name", "gene_start", "gene_end", "codon_start", "codon_end"])
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

def process_codons_based_on_strand(codon_info: CodonInfo, strand: str = '+') -> List[str]:
    """
    Process codons based on the provided information and strand orientation.

    Args:
        codon_info (CodonInfo): Information about the codon and associated SNPs.
        strand (str, optional): Strand orientation ('+' for positive strand, '-' for negative strand).
            Defaults to '+'

    Returns:
        List[str]: A list of strings containing processed information for each SNP in the codon.
    """
    output_list = []
    # Create a mutable version of the codon for mutation
    mutable_codon = list(codon_info.original_codon)
    for snp in codon_info.codon_list:
        position_in_codon = snp['codon']
        mutable_codon[position_in_codon] = snp['base']

    # Translate original and mutated codons
    original_aa = Seq(''.join(codon_info.original_codon)).translate()
    
    if strand == '-':
        mutated_aa = Seq(''.join(mutable_codon)).reverse_complement().translate()
        original_aa = Seq(''.join(codon_info.original_codon)).reverse_complement().translate()
    else:
        mutated_aa = Seq(''.join(mutable_codon)).translate()

    # Construct the output list for each SNP in the codon
    for snp in codon_info.codon_list:
        aa_position = (codon_info.codon_start - codon_info.gene_start + int(snp['codon'])) // 3 + 1
        if strand == '-':
            aa_position = (codon_info.gene_end - codon_info.codon_end + int(snp['codon'])) // 3 + 1
        chg_aa = ''.join([str(original_aa), str(aa_position), str(mutated_aa)])
        chg_aa = iupac_aa(chg_aa)
        
        sentence = '\t'.join([codon_info.gene_name, snp['position'], snp['base'], chg_aa])
        output_list.append(sentence)

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
        relevant_snps = [(i, *snp) for i, snp in enumerate(snp_list) if gene_start <= int(snp[0]) <= gene_end]
        
        for position in range(gene_start, gene_end, 3):
            codon_start, codon_end = position, position + 2
            # We calculate the SNP position within the codon while constructing this list
            codon_snps = [{
                'index': index,
                'position': snp_position,
                'base': snp_base,
                'codon': int(snp_position) - codon_start
            } for index, snp_position, snp_base in relevant_snps if codon_start <= int(snp_position) <= codon_end]
            
            if len(codon_snps) > 1:  # If multiple SNPs are in the same codon
                codon_sequence = Seq(sequence[codon_start-1:codon_end])
                
                info = CodonInfo(
                    codon_list=codon_snps,
                    original_codon=codon_sequence, 
                    new_codon=None,  
                    gene_name=gene_name, 
                    gene_start=gene_start,
                    gene_end=gene_end,
                    codon_start=codon_start,
                    codon_end=codon_end
                )
                  
                results.extend(process_codons_based_on_strand(info, strand))
    return results


def vcf_to_dataframe(vcf_filename):
    """
    Convert a VCF file to a pandas DataFrame.

    Args:
    - vcf_filename (str): Path to the VCF file.

    Returns:
    - pd.DataFrame: DataFrame representation of the VCF file, indexed by position number.
    """
    
    vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
    
    records = []
    for record in vcf_reader:
        record_data = {
            'CHROM': record.CHROM,
            'POS': record.POS,
            'ID': record.ID,
            'REF': record.REF,
            'ALT': [str(alt) for alt in record.ALT],
            'QUAL': record.QUAL,
            'FILTER': record.FILTER,
            'INFO': record.INFO,
            'FORMAT': record.FORMAT
        }
        if record.samples:
            sample_data = record.samples[0].data
            for key in sample_data._fields:
                record_data[key] = getattr(sample_data, key)  # Use getattr to access named tuple fields
        
        records.append(record_data)

    df = pd.DataFrame.from_records(records)
    df.set_index('POS', inplace=True, drop=False)  # Set 'POS' as index but still retain it in the DataFrame
    return df

def arguments():
    parser = argparse.ArgumentParser(description = 'script to annotate MNV') 
    parser.add_argument('-v', dest = 'vcf', required =True, help = 'Vcf file with snps')
    parser.add_argument('-f', dest = 'fasta', required =True, help = 'Name of the reference fasta') 
    parser.add_argument('-g', dest = 'genes', required = True, help = 'File with gene info')
    args = parser.parse_args()
    return args

def main():
    args = arguments()
    sequence = reference_fasta(args.fasta)
    lista_snp = getseq_posbase(args.vcf)
    gene_list = check_genes(lista_snp,args.genes)
    get_mnv_variants(gene_list, lista_snp, sequence)
    df = vcf_to_dataframe(args.vcf)
    print(df.head())
    print(df['FORMAT'])
    print(df.iloc[0])
    
    df['ANN'] = df['INFO'].apply(lambda x: x.get('ANN') if isinstance(x, dict) else None)
    df['ANN'].to_csv('column_ANN.txt', index=False, header=True)
#python3 new.py -v G35894.var.snp.vcf -f MTB_ancestor.fas -g anot_genes.txt

main()