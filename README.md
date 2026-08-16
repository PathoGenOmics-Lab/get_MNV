<p align="center">
  <img src="docs/assets/logo.svg#gh-light-mode-only" alt="get_MNV logo" width="640">
  <img src="docs/assets/logo-dark.svg#gh-dark-mode-only" alt="get_MNV logo" width="640">
</p>

<div align="center">

<a href="LICENSE"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/license-AGPL%20v3-%23af64d1?style=flat-square&amp;labelColor=21262d"><img alt="License: AGPL v3" src="https://img.shields.io/badge/license-AGPL%20v3-%23af64d1?style=flat-square"></picture></a>
<a href="https://anaconda.org/bioconda/get_mnv"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/conda/dn/bioconda/get_mnv.svg?style=flat-square&amp;label=bioconda&amp;labelColor=21262d"><img alt="Bioconda" src="https://img.shields.io/conda/dn/bioconda/get_mnv.svg?style=flat-square&amp;label=bioconda"></picture></a>
<a href="https://github.com/PathoGenOmics-Lab/get_MNV/releases"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/version-1.1.5-%23149389?style=flat-square&amp;labelColor=21262d"><img alt="Version" src="https://img.shields.io/badge/version-1.1.5-%23149389?style=flat-square"></picture></a>
<a href="https://doi.org/10.5281/zenodo.13907422"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.13907422-%23ff0077?style=flat-square&amp;labelColor=21262d"><img alt="DOI" src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.13907422-%23ff0077?style=flat-square"></picture></a>
<a href="https://github.com/PathoGenOmics-Lab"><picture><source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square&amp;labelColor=21262d"><img alt="PGO" src="https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square"></picture></a>

**Multi-Nucleotide Variant detection - codon-level annotation from VCF or iVar TSV.**
**Pure Rust · no C dependencies · cross-platform (macOS, Linux, Windows)**

[Documentation](https://pathogenomics-lab.github.io/get_MNV/) · [Install](#installation) · [Quick start](#quick-start) · [Citation](#citation)

**English** · [Español](README.es.md)

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>

---

## What is get_MNV?

get_MNV finds cases where two or more SNVs fall in the same codon and should be
interpreted together. These combined changes can produce a different amino acid
effect than the individual SNVs alone.

It takes variant calls (VCF or iVar `variants.tsv`), a reference FASTA and a
gene annotation (GFF/GFF3/GTF or a simple TSV), optionally the aligned reads,
and writes annotated variants as TSV, VCF, or both, plus a self-contained
interactive HTML report.

<p align="center">
  <img src="docs/assets/get_mnv_aa.png#gh-light-mode-only" alt="MNV amino acid reclassification" width="650">
  <img src="docs/assets/get_mnv_aa-dark.png#gh-dark-mode-only" alt="MNV amino acid reclassification" width="650">
</p>

- Groups SNVs by codon and reports SNP, MNV, or SNP/MNV calls
- Recalculates amino acid changes from the full codon haplotype
- Decomposes `REF/ALT` alleles into SNV, MNV, insertion, deletion, delins and
  complex indel components
- Counts SNP, MNV and exact indel support from a BAM, with strand bias
- Supports 9 NCBI genetic code tables
- Ships a desktop GUI with drag and drop and a genomic track viewer

## Installation

### Desktop GUI

Download the latest release for your platform:

| Platform | Download |
|---|---|
| 🍎 macOS (Apple Silicon) | [**get_MNV_1.1.5_aarch64.dmg**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🍎 macOS (Intel) | [**get_MNV_1.1.5_x64.dmg**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🐧 Linux | [**Releases page**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |
| 🪟 Windows | [**Releases page**](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest) |

> [!NOTE]
> **macOS users**: The app is not signed with an Apple Developer certificate. On first launch, right-click the app → **Open** → click **Open** in the dialog. See [Apple support](https://support.apple.com/en-us/HT202491) for details.

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

## Quick start

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both \
  --report run.html
```

`--bam` is optional and is what turns "these two changes are in one codon" into
"these two changes are on one molecule". Use `--tsv` instead of `--vcf` for an
iVar `variants.tsv`, and `--genes` instead of `--gff` for a four-column
annotation. `get_mnv --help` lists every option, and the
[CLI reference](https://pathogenomics-lab.github.io/get_MNV/cli-reference/)
explains what each one changes about the answer.

The output looks like this:

```text
Chromosome  Gene      Positions       Base Changes  AA Changes  Variant Type  Change Type
MTB_anc     Rv0095c   104838          T             Asp126Glu   SNP           Non-synonymous
MTB_anc     Rv0095c   104941,104942   T,G           Gly92Gln    SNP/MNV       Non-synonymous
MTB_anc     esxL      1341102,1341103 T,C           Arg33Ser    SNP/MNV       Non-synonymous
```

A ready-to-run *M. tuberculosis* dataset (reference, genes, VCF, and a tiny demo
BAM for the read viewer) lives in [`example/`](example/README.md). The
[tutorial](https://pathogenomics-lab.github.io/get_MNV/getting-started/) walks
through it.

## Documentation

**<https://pathogenomics-lab.github.io/get_MNV/>** is the manual, in English and
Spanish: every option with its default, what each output column means, and what
the tool does and does not take on. The same pages are in [`docs/`](docs/) in
this repository.

| Start here | |
|---|---|
| [Command line tutorial](https://pathogenomics-lab.github.io/get_MNV/getting-started/) | A first run end to end on the bundled data, with the output explained line by line |
| [Common recipes](https://pathogenomics-lab.github.io/get_MNV/usage/) | Ready-to-run commands for the usual jobs |
| [Desktop GUI tutorial](https://pathogenomics-lab.github.io/get_MNV/gui-tutorial/) | The same run in the app, screen by screen |

| Reference | |
|---|---|
| [CLI reference](https://pathogenomics-lab.github.io/get_MNV/cli-reference/) | Every option, with its default and what it changes |
| [Input formats](https://pathogenomics-lab.github.io/get_MNV/input-formats/) | What the VCF, FASTA, annotation and BAM have to look like |
| [Output formats](https://pathogenomics-lab.github.io/get_MNV/output-formats/) | Every TSV column, VCF INFO key and JSON field |
| [Example report](https://pathogenomics-lab.github.io/get_MNV/assets/example-report.html) | A real HTML report, open it and click around |

| How it works | |
|---|---|
| [Scope and compatibility](https://pathogenomics-lab.github.io/get_MNV/indel-mnv-semantics/) | What get_MNV takes on, what it leaves to your caller, and where its limits are |
| [Indels and local haplotypes](https://pathogenomics-lab.github.io/get_MNV/indel-haplotypes/) | How an indel is read off the alignments and what each number counts |
| [Linkage](https://pathogenomics-lab.github.io/get_MNV/linkage/) | Telling a real haplotype from two variants that merely share a codon |
| [Troubleshooting](https://pathogenomics-lab.github.io/get_MNV/troubleshooting/) | The errors that stop a run, and what each warning is telling you |

Version history is in [CHANGELOG.md](CHANGELOG.md).

## For developers

The core CLI and library live in `src/`. The desktop app uses Tauri in
`src-tauri/` and React/TypeScript in `frontend/`.

```bash
cargo test --workspace
npm run build --prefix frontend
bash scripts/build_get_mnv.sh
bash scripts/build_gui_bundle.sh
```

`tests/scenarios/` is a Python harness that builds synthetic FASTA, GFF, VCF and
BAM inputs from declarative scenarios, runs the compiled binary and checks each
output row by row; it needs `samtools` on `PATH`. See
[tests/scenarios/README.md](tests/scenarios/README.md) for the cases it covers
and how to add one.

```bash
cargo build                                # produces target/debug/get_mnv
python3 tests/scenarios/run.py             # run every scenario
```

## Citation

If you use get_MNV in your research, please cite:

> Ruiz-Rodriguez P, Coscolla M. **get_MNV: Multi-Nucleotide Variant detection tool.** Zenodo. doi: [10.5281/zenodo.13907422](https://doi.org/10.5281/zenodo.13907422)

```bibtex
@software{ruiz-rodriguez_get_mnv_2026,
  title     = {get\_MNV: Multi-Nucleotide Variant detection tool},
  author    = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
  year      = {2026},
  doi       = {10.5281/zenodo.13907422},
  url       = {https://github.com/PathoGenOmics-Lab/get_MNV},
  version   = {1.1.5},
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
