[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/PathoGenOmics-Lab/get_MNV/blob/main/LICENSE)
<p align="center">
  <img src="https://github.com/Pathogenomics-Lab/get_MNV/blob/main/images/get_mnv2.png" height="230" alt="get_MNV">
</p>

__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>  



# get_Multi-NucleotideVariants 
<p align="justify">Single Nucleotide Variants (SNVs) represent one of the most common types of genetic mutations. However, a situation may arise where multiple SNVs occur within the same codon, leading to the translation of a different amino acid. This is referred to as a Multi-Nucleotide Variant (MNV).

Current annotation programs, such as ANNOVAR or SnpEff, are designed to work predominantly with SNVs, which implies a potential gap in their functionality. Consequently, we may overlook the actual amino acid changes that result from multiple SNVs within the same codon. 

**get_MNV** seeks to address this issue, enhancing the comprehensiveness of genetic variant interpretation.</p>

<p align="center"><img src="https://github.com/Pathogenomics-Lab/get_MNV/blob/main/images/get_mnv_aa.png" height="350" alt="get_MNV"></p>


**IMPORTANT this script works with SNV against a reference, insertions and deletions modifiying reading frame are not currently supported**

## Installation

This script is written in _python 3_, so you will need this version of _python_ and the following modules, if they are not present it can be installed with `pip` command.


```bash
module argparse
module subprocess
module sys
module biopython
```

## Usage
**Execution command:**
```bash
python3 get_mnv.py -v input.var.snp.vcf -f reference_fasta.fas -g annotated_genes.txt
```
Arguments of the `get_mnv.py` script
```bash
script to annotate MNV

arguments:
  -h, --help        show this help message and exit
  -v VCF            Vcf file with snps
  -f FASTA          Name of the reference fasta
  -g GENES          File with gene info
  -install INSTALL  (optional) Install module with pip
```

The following input files are required:


Vcf file with SNV from snpeff annotation `-v input.var.snp.vcf`, for each line ir must appear one position:



```bash
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	
MTB_anc	1977	.	G	A	.	PASS	ADP=11;WT=0;HET=0;HOM=1;NC=0;ANN=A|intergenic_region|MODIFIER|dnaA_Rv0001-dnaN_Rv0002|gene0-gene1|intergenic_region|gene0-gene1|||n.1977G>A||||||	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:58:12:11:0:11:100%:1.4176E-6:0:44:0:0:5:6
MTB_anc	2532	.	C	T	.	PASS	ADP=8;WT=0;HET=0;HOM=1;NC=0;ANN=T|synonymous_variant|LOW|dnaN_Rv0002|gene1|transcript|Transcript_gene1|Coding|1/1|c.481C>T|p.Leu161Leu|481/1209|481/1209|161/402||,T|custom|MODIFIER|||CUSTOM&additionnal_annotations|Essential_in_vitro|||n.2532C>T||||||,T|custom|MODIFIER|||CUSTOM&additionnal_annotations|information_pathways|||n.2532C>T||||||	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	1/1:41:9:8:0:8:100%:7.77E-5:0:38:0:0:4:4
```

Reference sequence in fasta format `-f reference_fasta.fas`:

```bash
>Id
Sequence
```
Annotation file `-g annotated_genes.txt` with codificant genes of the organism of study separated with tabs = "\t", first row is for the name of the gene, second need to have start coordinate in the referece, third is the ending coordinate, and fourth column is the orientation of the gene:

```bash
Rv0007_Rv0007	9914	10828	+
ileT_Rvnt01	10887	10960	+
alaT_Rvnt02	11112	11184	+
Rv0008c_Rv0008c	11874	12311	-
ppiA_Rv0009	12468	13016	+
Rv0010c_Rv0010c	13133	13558	-
```

Optional command:
`-install biopython` installs biopython with pip command and runs the script

## Output

Generates two files:
First endend in MNV.tsv, is a tab separated file with only SNV annotated as MNV.
```bash
Rv0095c_Rv0095c	104941	T	Gly92Gln	missense_variant
Rv0095c_Rv0095c	104942	G	Gly92Gln	missense_variant
Rv0278c_Rv0278c	336081	G	Val77Thr	missense_variant
Rv0278c_Rv0278c	336082	T	Val77Thr	missense_variant
```

Second file is MNV.vcf, the original vcf annotated with SNV or MNV effect and in the cases of MNVs it will appear the AA annotated correctly
```bash
MTB_anc	104824	.	A	C	.	ADP=17;WT=0;HET=0;HOM=1;NC=0;ANN=C|SNP|missense_variant|Ile131Ser|missense_variant|MODERATE|Rv0095c_Rv0095c|gene100|transcript|Transcript_gene100|Coding|1/1|c.392T>G|p.Ile131Ser|392/411|392/411|131/136||,C|custom|MODIFIER|||CUSTOM&additionnal_annotations|Non_essential_in_vitro|||n.104824A>C||||||,C|custom|MODIFIER|||CUSTOM&additionnal_annotations|insertion_seqs_and_phages|||n.104824A>C||||||
MTB_anc	104838	.	G	T	.	ADP=16;WT=0;HET=1;HOM=0;NC=0;ANN=T|SNP|missense_variant|Asp126Glu|missense_variant|MODERATE|Rv0095c_Rv0095c|gene100|transcript|Transcript_gene100|Coding|1/1|c.378C>A|p.Asp126Glu|378/411|378/411|126/136||,T|custom|MODIFIER|||CUSTOM&additionnal_annotations|Non_essential_in_vitro|||n.104838G>T||||||,T|custom|MODIFIER|||CUSTOM&additionnal_annotations|insertion_seqs_and_phages|||n.104838G>T||||||
MTB_anc	104941	.	C	T	.	ADP=22;WT=0;HET=0;HOM=1;NC=0;ANN=T|MNV|missense_variant|Gly92Gln|missense_variant|MODERATE|Rv0095c_Rv0095c|gene100|transcript|Transcript_gene100|Coding|1/1|c.275G>A|p.Gly92Glu|275/411|275/411|92/136||,T|custom|MODIFIER|||CUSTOM&additionnal_annotations|Non_essential_in_vitro|||n.104941C>T||||||,T|custom|MODIFIER|||CUSTOM&additionnal_annotations|insertion_seqs_and_phages|||n.104941C>T||||||
MTB_anc	104942	.	C	G	.	ADP=25;WT=0;HET=0;HOM=1;NC=0;ANN=G|MNV|missense_variant|Gly92Gln|missense_variant|MODERATE|Rv0095c_Rv0095c|gene100|transcript|Transcript_gene100|Coding|1/1|c.274G>C|p.Gly92Arg|274/411|274/411|92/136||,G|custom|MODIFIER|||CUSTOM&additionnal_annotations|Non_essential_in_vitro|||n.104942C>G||||||,G|custom|MODIFIER|||CUSTOM&additionnal_annotations|insertion_seqs_and_phages|||n.104942C>G||||||
```
## References

SnpEff program: http://pcingola.github.io/SnpEff/

```bash
Cingolani, P., Platts, A., Wang, l., Coon, M., Nguyen, T., Wang, L., Land, S. J., Lu, X., & Ruden, D. M. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly, 6(2), 80–92. https://doi.org/10.4161/fly.19695
```

MTB ancestor fasta from Comas et al 2010: https://zenodo.org/record/3497110

```bash
 Human T cell epitopes of Mycobacterium tuberculosis are evolutionarily hyperconserved. Comas et al 2010; Nature Genetics (doi:10.1038/ng.590).
```
