#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                           .            .       :            ...   ::            ::.  
                                                                                           =:          :=       *+:          **-   %@*         .%@#   
                                                                                         .===.    .   -==.    .:***=.    .  :**:   -@@-.   .. :@@#    
                                                                                         -===-      .====.     =*++*+:      =**     #@#   .  -@@+.    
                                                                                    .   .==:==:    .==-==:  . .**=.:**+.   .**=    .-@@-  ..=@@=..    
                                                                                      ..==-.-== ..:==-.-=-   ::**....+**-. .**:.    .*@@. .*@@:  ..   
                                                                   .****:              -==   ==- :==:  :==    -**     :**=.-**    .  .@@* #@%:        
                                                                 .**.  =%             :==.   .==-==.   .==    +*=   .  .=****=        =@@%@#. .       
                                                                +#:  .**.     #+    .:==: . . -===. ...:==. :.**:..   ...:***-.    ... %@@#.   ..   :@
                                                         .-----%=   +#:       #*     -==       ==       ==:  :**.          =*          :@+          :@
                            :.                           .:::-:    =-         #*               .              ..            .           :           :@
               :=******: -#*-*%           -+*****+.        .%=  .******:      ##====================================================================*@
            .+#+:  :=: :#+   *+        -**-. -+: =%       -%:  =%:            ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          :**:  =*+-..+#.   :%       =#=  =*+-   =#     .#*  :#=                                     .                                                
        .#*.  =%-   .%=     %-     =%- .*#-     =%.    =%: .*#.                            .::::   -.                                                 
       +#.  .%=    -%: -%  +#    .#+  *#:     -#+    :%+  -%-                            :-:.:::-.                                                    
      +-   =%:    **   %= :@.   =%: :%-:..:=**=     +*.  *#         -%:             .:..:-                                                            
          -%    :%-   =#  +- .=#-  :@. :==-:     .+#:   +#        :#*            ::-:::-:                                                             
      .   %-   +#. : :@.  :+*+:    +#         .=#*:     %=      :#*:               :             -::                                                  
     -@   +*-+#- :#=:@: .%=.    =+ .#*=--=++**=:     *: :%+--=**=.        :::::::::::::::::::::. : .                                                  
     .@:   .:. :** :%:  %=      .%=   .::.  .=**:    *#   .::..-*#:      -::::::::::::::::::: .=                                                      
      :**=--=**=. -%.  #+         -***+++***+:        =#*+=++*+-                                                                                      
         .::.:   +#.  #=                                  .                                                                                           
          -*#= .#+  .%=                                                                                                                               
       -*#=   =%:  :%:                                                                                                                                
    .*#=    :#=   +#.                                                                                                                                 
  .**.    :#+   :%=                                                                                                                                   
  @:   .=#+   :#*                                                                                                                                     
  *+++*=:   -*+.                                                                                                                                      
@-      .-*#=                                                                                                                                         
:**++***+-                                                                                                                                            

Script for annotate vcf file with multiple-nucleotide variants

**IMPORTANT these script works with SNV against a reference, 
insertions and deletions modifiying reading frame are not currently supported**

Command:
    python3 get_mnv.py -v G35894.var.snp.vcf -f MTB_ancestor.fas -g anot_genes.txt
Outputs:
    vcf_name + .MNV.vcf  -> annotated vcfs
    vcf_name + .MNV.tsv  -> tsv file with MNV variants
Contact -> @paururo
"""

__author__ = 'Paula Ruiz Rodriguez'
__credits__ = ['Paula Ruiz Rodriguez', 'PathoGenOmics']
__license__ = 'GPL'
__version__ = '0.1.2'
__maintainer__ = 'Paula Ruiz Rodriguez: @paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'developing'
__lab__ = 'PathoGenOmics, I2SysBio'

import vcf
import pandas as pd
import ast
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from typing import List
from collections import namedtuple
from io import StringIO

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
CodonInfo = namedtuple('CodonInfo', ['codon_list', 'new_codon', 'original_codon', 'gene_name', 'gene_start', 'gene_end', 'codon_start', 'codon_end'])
SNP = namedtuple('SNP', ['index', 'position', 'base'])

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
    '''
    Process codons based on the provided information and strand orientation.

    Args:
        codon_info (CodonInfo): Information about the codon and associated SNPs.
        strand (str, optional): Strand orientation ('+' for positive strand, '-' for negative strand).
            Defaults to '+'

    Returns:
        List[str]: A list of strings containing processed information for each SNP in the codon.
    '''
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
    '''Extract SNPs that fall within a specific codon.'''
    return [snp for snp in snp_list if codon_start <= int(snp[0]) <= codon_end]

def get_mnv_variants(gene_list: list, snp_list: list, sequence: str) -> list:
    '''
    Identify codons with multiple SNPs within genes from the provided gene and SNP lists.
    
    Parameters:
    - gene_list: Contains gene name, start, end, and strand data in tab-separated format.
    - snp_list: Contains SNP position and base mutation information.
    - sequence: DNA sequence for deriving codons.
    
    Returns:
    - List of processed codons containing multiple SNPs based on gene strand.
    '''
    
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


def vcf_to_dataframe(vcf_filename: str):
    '''
    Convert a VCF file to a pandas DataFrame.

    Args:
    - vcf_filename (str): Path to the VCF file.

    Returns:
    - pd.DataFrame: DataFrame representation of the VCF file, indexed by position number.
    '''
    
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

def convert_to_dict(info_data):
    '''
    Converts a string representation of a dictionary to an actual dictionary.
    
    Args:
    - info_data (str or dict): Data to be converted into a dictionary.

    Returns:
    - dict: Converted dictionary or the original data if already a dictionary.
    - None: If the conversion fails or data type is unexpected.
    '''
    if isinstance(info_data, str):
        try:
            return ast.literal_eval(info_data)
        except ValueError:
            print(f'Failed to convert {info_data} to dictionary.')
            return None
    elif isinstance(info_data, dict):
        return info_data
    else:
        print(f'Unexpected data type for {info_data}.')
        return None

def modify_info_dict(info_dict, chg, aa):
    '''
    Modifies a dictionary based on the provided gene, change, and amino acid info.

    Args:
    - info_dict (dict): Dictionary to be modified.
    - gene (str): Gene information.
    - chg (str): Change information.
    - aa (str): Amino acid information.

    Returns:
    - dict: Modified dictionary.
    '''
    segments = info_dict['ANN'][0].split('|')
    segments.insert(1, 'MNV')
    segments[2] = chg
    segments[11] = f'p.{aa}'
    
    info_dict['ANN'][0] = '|'.join(segments)
    return info_dict

def change_vcf(df, mnv):
    '''
    Modifies a DataFrame based on a list of MNVs.

    Args:
    - df (pd.DataFrame): DataFrame with VCF data.
    - mnv (list): List of MNVs.

    Returns:
    - pd.DataFrame: Modified DataFrame.
    - mnv_position (set): Set of MNV positions.
    '''
    mnv_position = set()
    for snp in mnv:
        list_snp = snp.split('\t')
        gene = list_snp[0]
        position = int(list_snp[1])
        mnv_position.add(position)
        aa = list_snp[3]
        chg = list_snp[4]

        info_data = df.loc[position, 'INFO']
        info_dict = convert_to_dict(info_data)
        
        if info_dict is None:
            continue

        info_dict = modify_info_dict(info_dict, chg, aa)
        df.loc[position, 'INFO'] = str(info_dict) 

    return df, mnv_position

def modify_info_for_snv(info_dict):
    '''
    Modifies a dictionary to add 'SNV' for a specific segment.

    Args:
    - info_dict (dict): Dictionary to be modified.

    Returns:
    - dict: Modified dictionary.
    '''
    segments = info_dict['ANN'][0].split('|')
    segments.insert(1, 'SNV')
    info_dict['ANN'][0] = '|'.join(segments)
    return info_dict

def update_vcf(df, exclude_set):
    '''
    Updates a DataFrame's 'ANN' information by adding 'SNV' to specific segments.
    Rows with index in the 'exclude_set' are not updated.

    Args:
    - df (pd.DataFrame): DataFrame with VCF data.
    - exclude_set (set): Set of row indices to be excluded from the update.

    Returns:
    - pd.DataFrame: Updated DataFrame.
    '''
    for index, row in df.iterrows():
        if index in exclude_set:
            continue
        
        info_data = row['INFO']
        info_dict = convert_to_dict(info_data)
        
        if info_dict is None:
            continue

        info_dict = modify_info_for_snv(info_dict)
        df.at[index, 'INFO'] = str(info_dict)

    return df

def convert_info(row):
    if isinstance(row, str):
        try:
            info_dict = eval(row)
        except:
            return row  # if conversion fails, return the original string

        # Convert dictionary to ';' separated string format
        return ';'.join([f'{k}={v}' for k, v in info_dict.items()])
    return row

def process_info_string(s):
    replacements = [('{', ''), ('}', ''),
                    ('[', ''), (']', ''),
                    ("'", ''), (': ', '='),
                    (',', ';'), ('; ', ',')
    ]
    for old, new in replacements:
        s = s.replace(old, new)
    return s

def convert_to_vcf_format(df):
    # Convert list values in ALT column to string
    df['ALT'] = df['ALT'].str[0]
    
    # Convert list values in ANN column to string
    df['INFO'] = df['INFO'].apply(lambda x: x.replace('['', '').replace('']', ''))
    
    # Set FILTER column to 'PASS'
    df['FILTER'] = 'PASS'
    
    # Replace missing values with VCF's '.' representation
    df = df.fillna('.')
    
    # Convert the 'INFO' column's dictionary format to VCF's ';' separated format
    df['INFO'] = df['INFO'].apply(convert_info)
    df['INFO'] = df['INFO'].apply(process_info_string)
    # Drop unnecessary columns and format the DataFrame to resemble VCF
    columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    df = df[columns]

    return df

def write_to_tsv(data_list:List[str] ,filename:str):
    '''
    Convert a list of tab-separated values to a DataFrame, sort it by the second column, 
    and save it to a TSV file.
    
    Parameters:
    - data_list (List[str]): A list of strings containing tab-separated values.
    - filename (str): The base name of the file to save the sorted DataFrame.
    
    Returns:
    None. Writes sorted data to a TSV file named '{filename}.mnv.tsv'.
    '''
    data_str = "\n".join(data_list)
    df = pd.read_csv(StringIO(data_str), sep='\t', header=None)
    df[1] = pd.to_numeric(df[1])
    df = df.sort_values(by=1)
    df.to_csv(filename+'.MNV.tsv', sep='\t', index=False, header=False)

def write_to_vcf(df, filename):
    """Writes a DataFrame to a VCF format file."""
    with open(filename, 'w') as f:
        # Add VCF header
        f.write('##fileformat=VCFv4.2\n')
        f.write('#' + '\t'.join(df.columns) + '\n')
        df.to_csv(f, sep='\t', index=False, header=False)

def get_snpeff(vcf, fasta, genes):
    name = '.'.join(vcf.split('.')[:-1])
    # Obtain the sequence from a reference FASTA
    sequence = reference_fasta(fasta)
    # Get SNP list from VCF
    lista_snp = getseq_posbase(vcf)
    # Check which genes have been mentioned in the provided list
    gene_list = check_genes(lista_snp, genes)
    mnv = get_mnv_variants(gene_list, lista_snp, sequence)# Obtain MNV variants
    write_to_tsv(mnv, name)# Write MNV variants to TSV
    df = vcf_to_dataframe(vcf) # Convert VCF to DataFrame
    updated_df, mnv_position = change_vcf(df, mnv)# Change VCF entries based on MNV information
    ultimate_df = update_vcf(updated_df, mnv_position)# Update VCF with positions
    converted_df = convert_to_vcf_format(ultimate_df)# Convert the updated DataFrame back to VCF format
    
    write_to_vcf(converted_df, name + '.MNV.vcf')