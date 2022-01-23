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
__version__ = '0.1.0'
__maintainer__ = 'Paula Ruiz Rodriguez: @paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'developing'
__lab__ = 'PathoGenOmics, I2SysBio'

import argparse
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def install(package = 'biopython'):
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])

def reference_fasta(fasta_file = 'MTB_ancestor.fas'):
    '''
    Function to parse and get reference sequence, requires SeqIO
    Input  -> fasta file with id and sequence
    Output -> only sequence in str()
    '''
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence

def getseq_posbase(vcf_file = 'G35894.var.snp.vcf'):
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

def process_listcodon(lista_codon, new_codon, codon, my_aa, gene, out_write):
    for i in lista_codon:
        new_codon[int(i[0])] = str(i[2])

    new_codon = Seq(''.join(new_codon))
    try:
        my_newaa = new_codon.translate()
    except:
        if "-" in new_codon:
            my_newaa = "-"
    print(my_newaa)
    #chg_aa = str(my_aa) + str(codon) + str(my_newaa)
    chg_aa = ''.join([str(my_aa),str(codon),str(my_newaa)])
    chg_aa = iupac_aa(chg_aa)
    for i in lista_codon:
        sentence = '\t'.join ([gene, i[1], i[2], chg_aa]) + "\n"
        out_write.write(sentence)

def process_listcodonN(lista_codon, new_codon, codon, my_aa, gene, out_write):

    for i in lista_codon:

        new_codon[int(i[0])] = str(i[2])

    new_codon2 = Seq(''.join(new_codon))
    new2 = new_codon2.reverse_complement()
    try:
        my_newaa = new2.translate()
    except:
        if '-' in new2:
            my_newaa = '-'

    #chg_aa = str(my_aa) + str(codon) + str(my_newaa)
    chg_aa = ''.join([str(my_aa),str(codon),str(my_newaa)])
    chg_aa = iupac_aa(chg_aa)
    for i in lista_codon:
        sentence = '\t'.join ([gene, i[1], i[2], chg_aa]) + '\n'
        out_write.write(sentence)

def check_genes(lista_snp, gene_file):
    '''
    Function to clean genes containing snps in the analyzed vcf
    '''
    analyze_genelist = list()
    with open(gene_file,'r') as in_file:
        for line in in_file:
            l = line.strip('\n').split('\t')
            start_gene, end_gene = int(l[1]), int(l[2])
            for elemento in lista_snp:
                if int(elemento[0]) > end_gene:
                    break
                elif start_gene <= int(elemento[0]) <= end_gene:
                    if line not in analyze_genelist:
                        analyze_genelist.append(line)
    return analyze_genelist

def getMNV(analyze_genelist, lista_snp, sequence, out_write):                
    '''
    Function to distinct SNPs being MNV and calculate new aa change 
    '''
    for line in analyze_genelist:
        l = line.strip('\n').split('\t')
        gene, orientation = l[0], l[3]   
        codon = 1
        start_gene, end_gene = int(l[1]), int(l[2])
        for ele in range(start_gene, end_gene, 3):
            new_cod = (end_gene - start_gene + 1) / 3
            start_codon = ele
            end_codon = ele + 2 

            lista_codon = list()
            for elemento in lista_snp:
                if int(elemento[0]) > end_gene:
                    break

                elif start_gene <= int(elemento[0]) <= end_gene:
                    if start_codon <= int(elemento[0]) <= end_codon:
                        place = (int(elemento[0]) - start_gene) / 3
                        if '.33' in str(place): #second place in codon
                            p_name = 1
                        elif '.66' in str(place): #third place in codon
                            p_name = 2
                        else: #first place in codon
                            p_name = 0
                        lista_codon.append([p_name,elemento[0], elemento[1]])
            
            if len(lista_codon) > 1: #if exists consecutive SNPs in same codon            
                my_codon = Seq(sequence[start_codon-1:end_codon])
                if orientation == '-':
                    my_aa = my_codon.reverse_complement().translate()
                    codon = int(new_cod - codon + 1)
                    new_codon = list(my_codon)

                    process_listcodonN(lista_codon, new_codon, codon, my_aa,gene, out_write)
                    
                else:
                    my_aa = my_codon.translate()
                
                    new_codon = list(my_codon)
                    process_listcodon(lista_codon, new_codon, codon, my_aa,gene, out_write)

            codon += 1

def iupac_aa(codon):
    '''
    Function to translate iupac code to three code letter
    '''
    aa_code = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*':'*', '-':'Deletion'
    }

    ref = aa_code.get(codon[0])
    chg = aa_code.get(codon[-1])
    if ref == chg:
        status = 'synonymous_variant'
    elif chg == '*':
        status = 'stop_gained'
    elif chg == 'Deletion':
        status = 'deleted_codon'
    else:
        status = 'missense_variant'
    cod_codon = ''.join([aa_code.get(codon[0]),codon[1:-1],aa_code.get(codon[-1])]) + '\t' + status
    return cod_codon

def get_results(in_results):
    
    
    with open(in_results,'r') as in_input:
        list_results_pos = [line.strip('\n').split('\t')[1] for line in in_input]
    with open(in_results,'r') as in_input:
        list_results = [line.strip('\n').split('\t') for line in in_input]
    return list_results_pos, list_results

def list_snpeff():
    '''
    Function containing a list with functional status of aa chg from snpEff
    '''
    
    list_snpeff = ['chromosome_number_variation', 'exon_loss_variant', 
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
    
    return list_snpeff

def write_vcf(in_vcf, outfile, list_MNV,list_results):
    '''
    
    '''
    with open(outfile,'w') as o_sfile:
        with open(in_vcf, 'r') as i_efile:
            for line in i_efile:
                if '#' in line:
                    o_sfile.write(line)
                else:
                    l = line.strip('\n').split('\t')
                    if l[1] in list_MNV:
                        for element in list_results:
                            if element[1] == l[1]:
                               aa = element[3]
                               chg = element[4]
                               break
                        info = l[7].split('|')
                        info2 ='|'.join([info[0],'MNV', chg, aa]) + '|' + '|'.join(info[1:])
                        sentence = '\t'.join(l[0:6]) + '\t' + info2 + '\n' 
                        o_sfile.write(sentence)
                    else:
                        aa = 'aa_NA' #intergenic region or without aa info
                        chg = 'chg_Custom' 
                        info = l[7].split('|')
                        for e in info:
                            if 'p.' in e:
                                aa = e.strip('p.')
                                break
                        if info[1] == 'custom':
                            for e in info:
                                lis_snp = list_snpeff()
                                if e in lis_snp:
                                    chg = e

                        else:
                            chg = info[1]
                        info2 = info[0] + '|' + 'SNP' + '|' + chg + '|' + aa + '|' + '|'.join(info[1:])
                        sentence = '\t'.join(l[0:6]) + '\t' + info2 + '\n'

                        o_sfile.write(sentence)

def main():
    parser = argparse.ArgumentParser(description = 'script to annotate MNV') 
    parser.add_argument('-v', dest = 'vcf', required =True, help = 'Vcf file with snps')
    parser.add_argument('-f', dest = 'fasta', required =True, help = 'Name of the reference fasta') 
    parser.add_argument('-g', dest = 'genes', required = True, help = 'File with gene info')
    parser.add_argument('-install', dest = 'install', required = False, help = '(optional) Install module with pip')
    args = parser.parse_args()
    
    print(__doc__)
    if args.install:
        install(args.install)

    sequence = reference_fasta(args.fasta)
    lista_snp = getseq_posbase(args.vcf)
    outtxtfile = args.vcf.strip('.vcf') + '.MNV.tsv'
    with open(outtxtfile,'w') as out_write:
        analyze_genelist = check_genes(lista_snp, args.genes)
        getMNV(analyze_genelist, lista_snp, sequence, out_write)
    list_results_pos, list_results = get_results(outtxtfile)
    outvcffile = args.vcf.strip('.vcf') + '.MNV.vcf'
    write_vcf(args.vcf, outvcffile, list_results_pos, list_results)

if __name__ == '__main__':
    main()