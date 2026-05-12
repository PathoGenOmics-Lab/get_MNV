<p align="center">
  <img src="images/get_mnv.png" alt="get_MNV logo" width="800" />
</p>

<div align="center">

[![License: AGPL v3](https://img.shields.io/badge/license-AGPL%20v3-%23af64d1?style=flat-square)](LICENSE)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/get_mnv.svg?style=flat-square&label=bioconda)](https://anaconda.org/bioconda/get_mnv)
[![Version](https://img.shields.io/badge/version-1.1.3-%23149389?style=flat-square)](https://github.com/PathoGenOmics-Lab/get_MNV/releases)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.13907423-%23ff0077?style=flat-square)](https://doi.org/10.5281/zenodo.13907423)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square)](https://github.com/PathoGenOmics-Lab)

**Multi-Nucleotide Variant detection - codon-level annotation from VCF or iVar TSV.**
**Pure Rust · no C dependencies · cross-platform (macOS, Linux, Windows)**

[Quick Start](#quick-start) · [GUI](#desktop-gui) · [Features](#features) · [Docs](docs/) · [Citation](#citation)

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>

---

## What is get_MNV?

get_MNV finds cases where two or more SNVs fall in the same codon and should be interpreted together. These combined changes can produce a different amino acid effect than the individual SNVs alone.

The tool takes:

- Variant calls: VCF or iVar `variants.tsv`
- Reference sequence: FASTA
- Gene annotation: GFF/GFF3/GTF or a simple TSV file
- Optional aligned reads: BAM, used to count SNP and MNV read support

It writes annotated variants as TSV, VCF, or both.

<p align="center">
  <img src="images/get_mnv_aa.png" alt="MNV amino acid reclassification" width="650" />
</p>

**Main features:**

- Groups SNVs by codon and reports SNP, MNV, or SNP/MNV calls
- Recalculates amino acid changes from the full codon haplotype
- Reads VCF and iVar TSV variant calls
- Uses BAM reads, when provided, to count SNP/MNV support and strand bias
- Supports 9 NCBI genetic code tables
- Includes a desktop GUI for drag-and-drop analysis

## Installation

### Desktop GUI

Download the latest release for your platform:

| Platform | Download |
|---|---|
| 🍎 macOS (Apple Silicon) | [**get_MNV_1.1.3_aarch64.dmg**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🍎 macOS (Intel) | [**get_MNV_1.1.3_x64.dmg**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🐧 Linux | [**Releases page**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🪟 Windows | [**Releases page**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |

> [!NOTE]
> **macOS users**: The app is not signed with an Apple Developer certificate. On first launch, right-click the app → **Open** → click **Open** in the dialog. See [Apple support](https://support.apple.com/en-us/HT202491) for details.

All releases are available on the [Releases page](https://github.com/PathoGenOmics-Lab/get_MNV/releases).

### Command line

```bash
conda install -c bioconda get_mnv
```

or download a pre-built binary:

```bash
wget https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest/download/get_mnv
chmod +x get_mnv
./get_mnv --help
```

or build from source:

```bash
git clone https://github.com/PathoGenOmics-Lab/get_MNV.git
cd get_MNV
cargo install --path .
```

## Quick Start

### VCF input

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3
```

### iVar TSV input

```bash
get_mnv \
  --vcf sample_variants.tsv \
  --input-format ivar \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3
```

`--input-format ivar` is optional when the iVar header can be auto-detected, but adding it makes the command explicit.

### With BAM read support

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3
```

### TSV and VCF output

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both \
  --summary-json run.summary.json \
  --run-manifest run.manifest.json
```

Run `get_mnv --help` for the full list of options.

## Common Arguments

| Argument | What it does |
|---|---|
| `--vcf <FILE>` | Variant input file. It can be VCF or iVar TSV. |
| `--input-format auto|vcf|ivar` | Selects the variant parser. Default: `auto`. |
| `--bam <FILE>` | Optional sorted and indexed BAM for read support. |
| `--fasta <FILE>` | Reference FASTA. Contig names must match the variant file. |
| `--gff <FILE>` | Gene annotation in GFF/GFF3/GTF format. |
| `--genes <FILE>` | Simple gene annotation TSV. Use instead of `--gff`. |
| `--gff-features <LIST>` | Feature types to analyze, for example `CDS` or `gene,pseudogene`. |
| `--quality <N>` | Minimum variant quality. Default: `20`. |
| `--min-mapq <N>` | Minimum read mapping quality when using BAM. Default: `0`. |
| `--snp <N>` | Minimum SNP-supporting reads. Default: `0`. |
| `--min-snp-frequency <F>` | Minimum BAM-derived SNP frequency, from `0` to `1`. Default: `0`. |
| `--min-snp-strand <N>` | Minimum SNP-supporting reads required on each strand. Default: `0`. |
| `--mnv <N>` | Minimum MNV-supporting reads. Default: `0`. |
| `--min-mnv-frequency <F>` | Minimum BAM-derived MNV haplotype frequency, from `0` to `1`. Default: `0`. |
| `--min-mnv-strand <N>` | Minimum MNV-supporting reads required on each strand. Default: `0`. |
| `--both` | Write both TSV and VCF outputs. |
| `--summary-json <FILE>` | Write a machine-readable run summary. |
| `--run-manifest <FILE>` | Write command, version, inputs, outputs, and checksums. |

Frequency filters use read support recalculated from `--bam`, not the original
`OFREQ` value from VCF/iVar input. Use values such as `0.05` for 5% or `0.20`
for 20%. When VCF output is requested, low-frequency records are skipped by
default or marked with `FILTER=LowFrequency` when `--emit-filtered` is enabled.
SNP and MNV frequency filters are independent: `--min-snp-frequency` applies to
individual SNP observations, while `--min-mnv-frequency` applies to the phased
MNV haplotype. In mixed `SNP/MNV` calls, a strong MNV haplotype is not removed
just because the individual SNP observations are below the SNP threshold.
The read-count and strand-support filters are independent in the same way:
`--snp` and `--min-snp-strand` apply to SNP observations, while `--mnv` and
`--min-mnv-strand` apply to the MNV haplotype.

## Outputs

By default, get_MNV writes:

```text
<input_name>.MNV.tsv
```

With `--convert` or `--both`, it also writes:

```text
<input_name>.MNV.vcf
```

The most important output fields are:

| Column | Meaning |
|---|---|
| `Chromosome` | Contig name |
| `Gene` | Gene or feature name |
| `Positions` | One position for SNPs, multiple positions for MNVs |
| `Base Changes` | Alternative bases |
| `AA Changes` | Amino acid change after combining SNVs in the codon |
| `Variant Type` | `SNP`, `MNV`, `SNP/MNV`, or `INDEL` |
| `Change Type` | Synonymous, non-synonymous, stop gained/lost, unknown, etc. |

When a BAM is provided, extra columns report read depth, SNP support, MNV support, frequency, and strand counts.

## Features

| Feature | Description |
|---|---|
| 🧬 MNV detection | Groups SNVs in the same codon and reclassifies as MNVs |
| 🔬 Accurate AA changes | Computes amino acid changes from the full codon haplotype |
| 📊 Read support | BAM-based SNP/MNV read counts with strand-specific metrics |
| 🔍 Strand bias | Fisher exact test (SB, FS, SOR) with configurable filtering |
| 📁 Multiple outputs | TSV, VCF (plain/BGZF+Tabix), BCF, JSON summary, run manifest |
| ⚡ Parallel | Multi-threaded contig processing with Rayon |
| 🧪 Genetic codes | 9 NCBI translation tables (1, 2, 3, 4, 5, 6, 11, 12, 25) |
| 🧩 Flexible input | VCF or iVar TSV variant calls; GFF3/GTF or TSV annotations; multi-contig and multi-sample VCFs |
| ✅ Validation | Dry-run mode, strict metrics, input checksums, error JSON |
| 🖥️ Desktop GUI | Native Tauri app with drag-and-drop, genomic track viewer, dark mode |

## Desktop GUI

The desktop app gives the same analysis workflow in a visual interface:

- Drop VCF or iVar TSV variant files
- Drop FASTA, GFF/GTF/GFF3, and optional BAM files
- Choose common parameters from the form
- Run one sample or multiple matched samples
- Inspect, filter, and export results

```bash
bash scripts/dev.sh   # development
bash scripts/build_gui_bundle.sh  # production .app + .dmg bundle
```

## Example Output

```
Chromosome  Gene      Positions       Base Changes  AA Changes  Variant Type  Change Type
MTB_anc     Rv0095c   104838          T             Asp126Glu   SNP           Non-synonymous
MTB_anc     Rv0095c   104941,104942   T,G           Gly92Gln    SNP/MNV       Non-synonymous
MTB_anc     esxL      1341102,1341103 T,C           Arg33Ser    MNV           Non-synonymous
```

**Variant types:**
- **SNP**: single nucleotide change, one SNV per codon
- **MNV**: all reads carry multiple SNVs together (Multi-Nucleotide Variant)
- **SNP/MNV**: some reads carry individual SNVs, others carry the MNV combination

## Documentation

| Document | Description |
|---|---|
| [Usage](docs/usage.md) | Full CLI reference and examples |
| [Input formats](docs/input-formats.md) | VCF, FASTA, GFF, TSV, BAM specifications |
| [Output formats](docs/output-formats.md) | TSV, VCF, BCF, JSON output details |
| [Troubleshooting](docs/troubleshooting.md) | Common errors and solutions |
| [Benchmarking](docs/benchmarking.md) | Performance testing |
| [Changelog](CHANGELOG.md) | Version history |

## For Developers

The core CLI and library live in `src/`. The desktop app uses Tauri in
`src-tauri/` and React/TypeScript in `frontend/`.

Useful commands:

```bash
cargo test --workspace
npm run build --prefix frontend
bash scripts/build_get_mnv.sh
bash scripts/build_gui_bundle.sh
```

## Limitations

- Designed for SNVs against a reference sequence
- Insertions and deletions are detected but not fully codon-annotated
- Multiallelic VCF records require `--split-multiallelic` or pre-splitting (`bcftools norm -m -`)
- Variant contig names must match FASTA and GFF exactly
- **Multiple transcripts per gene**: when using `--gff-features CDS` with a GFF file that contains multiple transcripts for the same gene, each transcript is annotated independently, producing one output line per transcript per variant. If you want a single line per variant, filter your GFF to keep only the canonical transcript before running get_MNV (e.g., using [AGAT](https://github.com/NBISweden/AGAT) `agat_sp_keep_longest_isoform.pl` or a similar tool)

## Citation

If you use get_MNV in your research, please cite:

> Ruiz-Rodriguez P, Coscolla M. **get_MNV: Multi-Nucleotide Variant detection tool.** Zenodo. doi: [10.5281/zenodo.13907423](https://doi.org/10.5281/zenodo.13907423)

```bibtex
@software{ruiz-rodriguez_get_mnv_2026,
  title     = {get\_MNV: Multi-Nucleotide Variant detection tool},
  author    = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
  year      = {2026},
  doi       = {10.5281/zenodo.13907423},
  url       = {https://github.com/PathoGenOmics-Lab/get_MNV},
  version   = {1.1.3},
  license   = {AGPL-3.0}
}
```

## License

[GNU Affero General Public License v3.0](LICENSE)

## Fun

Click for the 3D printable logo:

<p align="center">
  <a href="https://www.printables.com/model/1030383-get_mnv-logo" target="_blank">
    <img src="https://media.printables.com/media/prints/1030383/images/7820375_62fee28c-1ef3-446a-9187-3a74e3912c09_7526c3fd-6f35-4ec1-ab2c-a8022ac592e9/thumbs/inside/1920x1440/jpg/img_3773.webp" height="350" alt="get_MNV 3D logo">
  </a>
</p>

---

<h2 id="contributors" align="center">

✨ Contributors
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START -->
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
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Design">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td>
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).
</div>
<!-- ALL-CONTRIBUTORS-LIST:END -->
