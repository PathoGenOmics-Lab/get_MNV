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
 
Requires biopython

Contact -> @paururo
"""

__author__ = 'Paula Ruiz Rodriguez'
__credits__ = ['Paula Ruiz Rodriguez', 'PathoGenOmics']
__license__ = 'GPL'
__version__ = '0.1.1'
__maintainer__ = 'Paula Ruiz Rodriguez: @paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'developing'
__lab__ = 'PathoGenOmics, I2SysBio'

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

#constant variable 
SNPEFF = ['chromosome_number_variation', 'exon_loss_variant', 
'frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost', 
'splice_acceptor_variant', 'splice_donor_variant', 'rare_amino_acid_variant', 
'missense_variant', 'disruptive_inframe_insertion', 'conservative_inframe_insertion', 
'disruptive_inframe_deletion', 'conservative_inframe_deletion', 
'5_prime_UTR_truncation+exon_loss_variant', '3_prime_UTR_truncation+exon_loss', 
'splice_branch_variant', 'splice_region_variant', 'stop_retained_variant', 
'initiator_codon_variant', 'synonymous_variant', 
'initiator_codon_variant+non_canonical_start_codon', 'stop_retained_variant', 
'coding_sequence_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 
'5_prime_UTR_premature_start_codon_gain_variant', 'upstream_gene_variant', 
'downstream_gene_variant', 'TF_binding_site_variant', 
'regulatory_region_variant', 'miRNA', 'custom', 'sequence_feature', 
'conserved_intron_variant', 'intron_variant', 'intragenic_variant', 
'conserved_intergenic_variant', 'intergenic_region', 'coding_sequence_variant', 
'non_coding_exon_variant', 'nc_transcript_variant', 'gene_variant', 'chromosome']

def reference_fasta(fasta_file: str = 'MTB_ancestor.fas'):
    '''
    Function to parse and get reference sequence, requires SeqIO
    Input  -> fasta file with id and sequence
    Output -> only sequence in str()
    '''
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence

def read_genes_names(genes_file: str):
    '''
    Function to get gene names  
    Input  -> text file with gene names in first column (tab sep)  
    Output -> list with gene names  
    '''
    list_name_genes = list()    
    with open(genes_file,'r') as in_file:
        for line in in_file:
            gene_name = line.strip('\n').split('\t')[0]
            list_name_genes.append(gene_name) 
    return list_name_genes  

def getseq_posbase(vcf_file: str = 'G35894.var.snp.vcf'):
    '''
    Function to get positions and alternative base of annotated vcf file 
    ignoring intergenic regions
    Input  -> annotated vcf file with snpEff
    Output -> list of list (positions with alternative base)
    '''
    list_snp = list()    
    with open(vcf_file,'r') as in_file:
        for line in in_file:
            if '#' not in line: #ignore header
                if 'intergenic' not in line: #ignore intergenic pos
                    l = line.strip('\n').split('\t')
                    list_snp.append([l[1],l[4]]) #[1:pos,4:alt_base]
    return list_snp

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

def check_genes(list_snp: list, gene_file: str):
    '''
    Function to get genes containing snps in the analyzed vcf  
        Input  -> list_snps: list with snps of vcf file, gene_file: 
                  genes with coord  
        Output -> list with genes of interest containing snps of vcf file
    '''
    analyze_genelist = list()
    with open(gene_file,'r') as in_file:
        for line in in_file:
            l = line.strip('\n').split('\t')
            start_gene, end_gene = int(l[1]), int(l[2])
            for elemento in list_snp:
                if int(elemento[0]) > end_gene:
                    break
                elif start_gene <= int(elemento[0]) <= end_gene:
                    if line not in analyze_genelist:
                        analyze_genelist.append(line)
    return analyze_genelist

def process_listcodon(lista_codon: list, new_codon, codon, my_aa, gene, lista_salida: list):
    '''
    Function for process snps in positive strand, gets old codon, and translates with new snps
    '''
    for i in lista_codon:
        new_codon[int(i[0])] = str(i[2])
    new_codon = Seq(''.join(new_codon))    
    my_newaa = new_codon.translate()
    
    chg_aa = ''.join([str(my_aa),str(codon),str(my_newaa)])
    chg_aa = iupac_aa(chg_aa)
    esta = False
    for i in lista_codon:
        sentence = '\t'.join ([gene, i[1], i[2], chg_aa]) + "\n"
        for elemento in lista_salida:
            pos = elemento[0].split("\t")[1]
            if int(i[1]) == int(pos):
                elemento.append(sentence)
                esta = True
        if not esta:
            lista_salida.append([sentence])