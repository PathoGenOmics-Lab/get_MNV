<div align="center">
  
[![get_mnv](https://img.shields.io/badge/get_mnv-rust-%23ff8000?style=flat-square)](https://github.com/PathoGenOmics-Lab/get_MNV)
[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-%23af64d1?style=flat-square)](https://github.com/PathoGenOmics-Lab/get_MNV/blob/main/LICENSE) 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/get_mnv/badges/downloads.svg)](https://anaconda.org/bioconda/get_mnv)
[![Anaconda-Version Badge](https://anaconda.org/bioconda/get_mnv/badges/version.svg)](https://anaconda.org/bioconda/get_mnv)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.13907423-%23ff0077?style=flat-square)](https://doi.org/10.5281/zenodo.13907423)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-red?style=flat-square)](https://github.com/PathoGenOmics-Lab)

</div>

<p align="center">
  <a href="https://github.com/PathoGenOmics-Lab/get_MNV">
    <img src="https://github.com/Pathogenomics-Lab/get_MNV/blob/main/images/get_mnv.png" height="350" alt="get_MNV">
  </a>
</p>

__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>  



# get Multi-Nucleotide Variants 
<p align="justify">get_MNV is a tool designed to identify Multi-Nucleotide Variants (MNVs) within the same codon in genomic sequences. MNVs occur when multiple Single Nucleotide Variants (SNVs) are present within the same codon, leading to the translation of a different amino acid. This tool addresses limitations in current annotation programs like ANNOVAR or SnpEff, which are primarily designed to work with individual SNVs and might overlook the actual amino acid changes resulting from MNVs.

**get_MNV** seeks to address this issue, enhancing the comprehensiveness of genetic variant interpretation.</p>

<p align="center"><img src="https://github.com/Pathogenomics-Lab/get_MNV/blob/main/images/get_mnv_aa.png" height="350" alt="get_MNV"></p>


**IMPORTANT this script works with SNV against a reference, insertions and deletions modifiying reading frame are not currently supported**

## Features
- MNV Identification: Detects SNVs occurring within the same codon and reclassifies them as MNVs.
- Accurate Amino Acid Change Calculation: Computes the resulting amino acid changes based on genomic reads.
- Integration with BAM and VCF Files: Supports input from VCF files for variants and optional BAM files for aligned reads.
- Quality Analysis: Allows setting a minimum Phred quality threshold to filter out low-quality reads.

## Installation
You can install get_MNV via conda, mamba or downloading [the binary file](https://github.com/PathoGenOmics-Lab/get_MNV/releases/download/1.0.0/get_mnv):
### Using conda
```
conda install -c bioconda get_mnv
```
### Using mamba
```
mamba install -c bioconda get_mnv
```
### Using binary
```
wget https://github.com/PathoGenOmics-Lab/get_MNV/releases/download/1.0.0/get_mnv
```
# Usage
```
get_mnv [OPTIONS] --vcf <VCF_FILE> --fasta <FASTA_FILE> --genes <GENES_FILE>
```
## Options:
- -v, --vcf <VCF_FILE>: VCF file containing the SNVs. (Required)
- -b, --bam <BAM_FILE>: BAM file with aligned reads. (Optional)
- -f, --fasta <FASTA_FILE>: FASTA file with the reference sequence. (Required)
- -g, --genes <GENES_FILE>: File containing gene information. (Required)
- -q, --quality <QUALITY>: Minimum Phred quality score (default: 20).
## Example:
```
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --genes genes.txt \
  --quality 30
```
## Input File Formats

- VCF File: Should contain the identified SNVs.
- BAM File: (Optional) Genomic reads aligned to the reference sequence.
- FASTA File: Reference genomic sequence.
- Gene File: A tab-delimited text file with the following structure per line (GeneName,GeneStart,GeneEnd,Strand):
```bash
Rv0007_Rv0007	9914	10828	+
ileT_Rvnt01	10887	10960	+
alaT_Rvnt02	11112	11184	+
Rv0008c_Rv0008c	11874	12311	-
ppiA_Rv0009	12468	13016	+
Rv0010c_Rv0010c	13133	13558	-
```
## Output
The program generates a TSV file named <vcf_filename>.MNV.tsv containing the following information:
- Gene: Name of the gene.
- Positions: Positions of the variants.
- Base Changes: Nucleotide base changes.
- AA Changes: Resulting amino acid changes.
- SNP AA Changes: Amino acid changes if considering individual SNVs.
- Variant Type: Type of variant (SNP, MNV, or SNP/MNV).
- Change Type: Type of change at the protein level (Synonymous, Non-synonymous, Stop gained).
- SNP Reads: (If BAM provided) Count of reads supporting each SNP.
- MNV Reads: (If BAM provided) Count of reads supporting the MNV.
    
# Limitations
- The script currently works only with SNVs compared against a reference sequence.
- Insertions and deletions that modify the reading frame are not supported in this version.

# Fun
## 3D model logo
Click for the stl file
<p align="center">
  <a href="https://www.printables.com/model/1030383-get_mnv-logo" target="_blank">
    <img src="https://media.printables.com/media/prints/1030383/images/7820375_62fee28c-1ef3-446a-9187-3a74e3912c09_7526c3fd-6f35-4ec1-ab2c-a8022ac592e9/thumbs/inside/1920x1440/jpg/img_3773.webp" height="350" alt="get_MNV">
  </a>
</p>

