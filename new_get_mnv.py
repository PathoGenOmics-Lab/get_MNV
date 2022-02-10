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
    Works inside get_MNV function
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

def process_listcodonN(lista_codon: list, new_codon, codon, my_aa, gene, lista_salida: list):
    '''
    Function for process snps in negative strand, gets old codon, and translates with new snps
    Works inside get_MNV function
    '''
    for i in lista_codon:
        new_codon[int(i[0])] = str(i[2])

    new_codon2 = Seq(''.join(new_codon))
    new2 = new_codon2.reverse_complement()

    my_newaa = new2.translate()

    chg_aa = ''.join([str(my_aa),str(codon),str(my_newaa)])
    chg_aa = iupac_aa(chg_aa)
    esta = False
    for i in lista_codon:
        sentence = '\t'.join([gene, i[1], i[2], chg_aa]) + '\n'
        for elemento in lista_salida:
            pos = elemento[0].split('\t')[1]
            if int(i[1]) == int(pos):
                elemento.append(sentence)
                esta = True
        if not esta:
            lista_salida.append([sentence])

def get_place(place: float):
    '''
    Function to get position of snp in a codon
    '''
    
    if '.33' in str(place): #second place in codon
        p_name = 1
    elif '.66' in str(place): #third place in codon
        p_name = 2
    else: #first place in codon
        p_name = 0
    return p_name

def getMNV(analyze_genelist: list, lista_snp: list, sequence):                
    '''
    Function to distinct SNPs being MNV and calculate new aa change 
    '''
    lista_salida = []
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
                        p_name = get_place(place)
                        lista_codon.append([p_name,elemento[0], elemento[1]])
            
            if len(lista_codon) > 1: #if exists consecutive SNPs in same codon            
                my_codon = Seq(sequence[start_codon-1:end_codon])
                if orientation == '-':
                    my_aa = my_codon.reverse_complement().translate()
                    p_codon = int(new_cod - codon + 1)
                    new_codon = list(my_codon)
                    process_listcodonN(lista_codon, new_codon, p_codon, my_aa, gene, lista_salida)
                    
                else:
                    my_aa = my_codon.translate()
                    new_codon = list(my_codon)
                    process_listcodon(lista_codon, new_codon, codon, my_aa, gene, lista_salida)

            codon += 1
    return lista_salida

def write_vcf(in_vcf, outfile, list_MNV, list_name_genes):
    '''
    
    '''
    a_a=False
    with open(outfile,'w') as o_sfile:
        with open(in_vcf, 'r') as i_efile:
            for line in i_efile:
                lista_cambios = []                
                mnv = False
                if '#' in line:
                    o_sfile.write(line)
                else:
                    lista_genes_cambios = []
                    l = line.strip('\n').split('\t')
                    for element1 in list_MNV:
                        for element in element1:
                            pos = element.strip("\n").split("\t")[1]
                            if int(pos) == int(l[1]):
                                aa = element.strip("\n").split("\t")[3]
                                chg = element.strip("\n").split("\t")[4]
                                gene = element.strip("\n").split("\t")[0]
                                lista_genes_cambios.append(gene)
                                sent = '-'.join(['MNV', gene, aa, chg])
                                if sent not in lista_cambios:
                                    lista_cambios.append(sent)
                                mnv = True
                                 
                    info = l[7].split('|')        
                    lista_snp = []
                    
                    for elemento in info:
                        if elemento in list_name_genes:
                            if mnv:
                                if elemento not in lista_genes_cambios:
                                    gene = elemento
                                    if 'p.' in elemento:
                                        aa = elemento.strip('p.')
                                    elif elemento in SNPEFF:
                                        chg = elemento
                                    sent = '-'.join(['SNP',gene,aa,chg])

                                    lista_cambios.append(sent)
                            elif "intragenic" not in line:
                                gene = elemento
                                if 'p.' in elemento:
                                    aa = elemento.strip('p.')
                                    a_a = True
                                elif "custom" not in elemento: 
                                    if elemento in SNPEFF:
                                        chg = elemento
                    if a_a:
                        sent = '-'.join(['SNP',gene,aa,chg])
                        lista_snp.append(sent)
                        a_a=False
                                                                     
                    if len(lista_snp) != 0:
                        for ele in lista_snp:
                            lista_cambios.append(ele)

                    if "intergenic" not in line:
                        res = ";".join(lista_cambios)
                        if "MNV" in res:
                            inf = "MNV"
                            if "SNP" in res:
                                inf = "MNV+SNP"
                        else:
                            inf = "SNP"

                            g = False
                            p = False
                            t = False            
                            for elemento in info:
                                if 'p.' in elemento:
                                    aa = elemento.strip('p.')
                                    g = True
                                elif "gene" in elemento:
                                    gene = elemento.strip('Transcript_gene-').strip('gene-')
                                    p = True
                                elif "custom" not in elemento:
                                    if elemento in SNPEFF:
                                        chg = elemento
                                    t = True
                                if g and p and t:

                                    res = 'SNP-' + gene + '-' + aa + '-' + chg
                                    g = False
                                    p = False
                                    t = False 

                        sen1 = inf + '|' + res
                        sentence = '\t'.join(l[0:6]) + '\t' + info[0]+ '|' + sen1 + '|' + '|'.join(info[1:])+"\n"
                        o_sfile.write(sentence)
                    else:
                        sen1 = "intergenic" +"|"+"intergenic"
                        sentence = '\t'.join(l[0:6]) + '\t' + info[0]+ '|' + sen1 + '|' + '|'.join(info[1:])+"\n"
                        o_sfile.write(sentence)


def main():
    parser = argparse.ArgumentParser(description = 'script to annotate MNV') 
    parser.add_argument('-v', dest = 'vcf', required =True, help = 'Vcf file with snps')
    parser.add_argument('-f', dest = 'fasta', required =True, help = 'Name of the reference fasta') 
    parser.add_argument('-g', dest = 'genes', required = True, help = 'File with gene info')
    args = parser.parse_args()
    
    print(__doc__)

    sequence = reference_fasta(args.fasta)
    lista_snp = getseq_posbase(args.vcf)
    outtxtfile = args.vcf.strip('.vcf') + '.MNV.tsv'
    with open(outtxtfile,'w') as out_write:
        analyze_genelist = check_genes(lista_snp, args.genes)
        list_results = getMNV(analyze_genelist, lista_snp, sequence)
        for element in list_results:
            for i in element:
                out_write.write(i)
        list_name_genes = read_genes_names(args.genes)
    
    outvcffile = args.vcf.strip('.vcf') + '.MNV.vcf'
    write_vcf(args.vcf, outvcffile, list_results, list_name_genes)

if __name__ == '__main__':
    main()