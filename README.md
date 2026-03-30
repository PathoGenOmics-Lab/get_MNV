
<p align="center">
  <a href="https://github.com/PathoGenOmics-Lab/get_MNV">
    <img src="https://github.com/Pathogenomics-Lab/get_MNV/blob/main/images/get_mnv.png" height="350" alt="get_MNV">
  </a>
</p>
<div align="center">

[![get_mnv](https://img.shields.io/badge/get_mnv-rust-%23ff8000?style=flat-square)](https://github.com/PathoGenOmics-Lab/get_MNV)
[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-%23af64d1?style=flat-square)](https://github.com/PathoGenOmics-Lab/get_MNV/blob/main/LICENSE)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/get_mnv.svg?style=flat-square)](https://anaconda.org/bioconda/get_mnv)
[![Anaconda-Version Badge](https://anaconda.org/bioconda/get_mnv/badges/version.svg)](https://anaconda.org/bioconda/get_mnv)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.13907423-%23ff0077?style=flat-square)](https://doi.org/10.5281/zenodo.13907423)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-red?style=flat-square)](https://github.com/PathoGenOmics-Lab)

</div>

# get_MNV — Multi-Nucleotide Variant detection

**get_MNV** identifies Multi-Nucleotide Variants (MNVs) within the same codon in genomic sequences. When multiple SNVs co-occur in one codon, the resulting amino acid change can differ from what individual SNV annotation predicts. Standard tools like ANNOVAR or SnpEff annotate SNVs independently and may miss these compound effects.

get_MNV addresses this by:

- **Detecting codon-level MNVs** from VCF + reference + gene annotation
- **Recalculating amino acid changes** considering all SNVs in each codon together
- **Quantifying read support** from BAM files (optional) with strand-bias statistics
- **Outputting TSV, VCF, and BCF** with full reproducibility metadata

<p align="center"><img src="https://github.com/Pathogenomics-Lab/get_MNV/blob/main/images/get_mnv_aa.png" height="350" alt="get_MNV"></p>

## Installation

### Conda / Mamba (recommended)

```bash
conda install -c bioconda get_mnv
# or
mamba install -c bioconda get_mnv
```

### Pre-built binary

```bash
wget https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest/download/get_mnv
chmod +x get_mnv
```

### From source

```bash
git clone https://github.com/PathoGenOmics-Lab/get_MNV.git
cd get_MNV
cargo install --path .
```

## Quick start

```bash
# Basic: TSV output
get_mnv --vcf variants.vcf --fasta reference.fasta --gff genes.gff3

# With BAM reads and quality filters
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --quality 30 \
  --mapq 20

# Both TSV + VCF output with strand-bias filtering
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both --vcf-gz --emit-filtered \
  --min-strand-bias-p 0.05

# Mitochondrial genome with vertebrate genetic code
get_mnv \
  --vcf mito.vcf \
  --fasta mito.fasta \
  --gff mito.gff3 \
  --translation-table 2
```

## Features

| Feature | Description |
|---------|-------------|
| 🧬 MNV detection | Identifies SNVs in the same codon and reclassifies as MNVs |
| 🔬 Accurate AA changes | Computes amino acid changes from the full codon haplotype |
| 📊 Read support | BAM-based SNP/MNV read counts with strand-specific metrics |
| 📁 Multiple outputs | TSV, VCF (plain/BGZF), BCF, JSON summary, reproducibility manifest |
| ⚡ Parallel | Multi-threaded processing with Rayon |
| 🔍 Strand bias | Fisher exact test with configurable filtering |
| 🧩 Flexible input | GFF3 or TSV annotations, multi-contig, multi-sample |
| 🧪 Genetic codes | 9 NCBI translation tables (bacterial, mitochondrial, etc.) |
| ✅ Validation | Dry-run mode, strict metrics, input checksums |

## Documentation

| Document | Description |
|----------|-------------|
| [Usage](docs/usage.md) | Full CLI reference and examples |
| [Input formats](docs/input-formats.md) | VCF, FASTA, GFF, TSV, BAM specifications |
| [Output formats](docs/output-formats.md) | TSV, VCF, BCF, JSON output details |
| [Troubleshooting](docs/troubleshooting.md) | Common errors and exit codes |
| [Benchmarking](docs/benchmarking.md) | Performance testing |
| [Changelog](CHANGELOG.md) | Version history |

## Example output

```
Chromosome  Gene      Positions       Base Changes  AA Changes  Variant Type  Change Type
MTB_anc     Rv0095c   104838          T             Asp126Glu   SNP           Non-synonymous
MTB_anc     Rv0095c   104941,104942   T,G           Gly92Gln    SNP/MNV       Non-synonymous
MTB_anc     esxL      1341102,1341103 T,C           Arg33Ser    MNV           Non-synonymous
```

## Limitations

- Works with SNVs against a reference sequence
- Insertions and deletions that modify the reading frame are detected but not fully annotated
- Multiallelic VCF records require `--split-multiallelic` or pre-splitting (`bcftools norm -m -`)
- VCF contig names must match FASTA and GFF exactly

## Citation

If you use get_MNV in your research, please cite:

> Ruiz-Rodriguez P, Coscolla M. get_MNV: Multi-Nucleotide Variant detection tool. Zenodo. https://doi.org/10.5281/zenodo.13907423

## Contributors

<div align="center">
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      💻 🔬 🤔 🔣 🎨 🔧
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      🔍 🤔 🧑‍🏫 🔬 📓
    </td>
  </tr>
</table>
</div>

## Fun

Click for the 3D printable logo:

<p align="center">
  <a href="https://www.printables.com/model/1030383-get_mnv-logo" target="_blank">
    <img src="https://media.printables.com/media/prints/1030383/images/7820375_62fee28c-1ef3-446a-9187-3a74e3912c09_7526c3fd-6f35-4ec1-ab2c-a8022ac592e9/thumbs/inside/1920x1440/jpg/img_3773.webp" height="350" alt="get_MNV 3D logo">
  </a>
</p>
