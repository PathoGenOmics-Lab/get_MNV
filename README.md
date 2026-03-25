

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

__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>  



# get Multi-Nucleotide Variants 
<p align="justify">get_MNV is a tool designed to identify Multi-Nucleotide Variants (MNVs) within the same codon in genomic sequences. MNVs occur when multiple Single Nucleotide Variants (SNVs) are present within the same codon, leading to the translation of a different amino acid. This tool addresses limitations in current annotation programs like ANNOVAR or SnpEff, which are primarily designed to work with individual SNVs and might overlook the actual amino acid changes resulting from MNVs.

**get_MNV** seeks to address this issue, enhancing the comprehensiveness of genetic variant interpretation.</p>

<p align="center"><img src="https://github.com/Pathogenomics-Lab/get_MNV/blob/main/images/get_mnv_aa.png" height="350" alt="get_MNV"></p>


**IMPORTANT this script works with SNV against a reference, insertions and deletions modifiying reading frame are not currently supported**

## 💾 Features
- MNV Identification: Detects SNVs occurring within the same codon and reclassifies them as MNVs.
- Accurate Amino Acid Change Calculation: Computes the resulting amino acid changes based on genomic reads.
- Integration with BAM and VCF Files: Supports input from VCF files for variants and optional BAM files for aligned reads.
- Quality Analysis: Allows setting a minimum Phred quality threshold to filter out low-quality reads.

## 🛠️ Installation
You can install get_MNV via conda, mamba (for unix/mac) or downloading [the binary file](https://github.com/PathoGenOmics-Lab/get_MNV/releases/download/1.0.0/get_mnv) (unix):
### 🐍 Using conda
```
conda install -c bioconda get_mnv
```
### 🐍 Using mamba
```
mamba install -c bioconda get_mnv
```
### 📨 Using binary
```
wget https://github.com/PathoGenOmics-Lab/get_MNV/releases/download/1.0.0/get_mnv
```
# 📎 Usage
```
get_mnv [OPTIONS] --vcf <VCF_FILE> --fasta <FASTA_FILE> (--genes <GENES_FILE> | --gff <GFF_FILE>)
```
## 🗃️ Options:
- -v, --vcf <VCF_FILE>: VCF file containing the SNVs. (Required)
- -b, --bam <BAM_FILE>: BAM file with aligned reads. (Optional)
- -f, --fasta <FASTA_FILE>: FASTA file with the reference sequence. (Required)
- -g, --genes <GENES_FILE>: Gene annotation file in TSV format. (Required if `--gff` is not used)
- --gff <GFF_FILE>: Gene annotation file in GFF/GFF3 format. (Required if `--genes` is not used)
- --gff-features <FEATURES>: Comma-separated GFF feature types to analyze (default: `gene,pseudogene`). Example: `--gff-features CDS,tRNA`.
- --sample <SAMPLE>: Sample name used to read original FORMAT metrics (`DP`/`AF`/`FREQ`) from multi-sample VCF. Default: first sample.
  - Use `--sample all` to process every sample in the VCF and write one output per sample.
- --chrom <CHROM>: Optional contig selection. Accepts one or multiple comma-separated contigs (e.g. `--chrom chr1,chr2`). Default: all contigs in VCF.
- -q, --quality <QUALITY>: Minimum Phred quality score (default: 20).
- --mapq <MAPQ>: Minimum mapping quality MAPQ (default: 0).
- --threads <N>: Number of worker threads for parallel processing (default: Rayon auto).
- --normalize-alleles: Normalize REF/ALT alleles (trim shared prefix/suffix) before processing.
- --split-multiallelic: Split multiallelic VCF records (`ALT=A,C,...`) into independent ALT alleles.
- --dry-run: Validate inputs and print per-contig summary without writing output files.
- --strict: Fail if original input metrics (`ODP`/`OFREQ`) cannot be recovered from VCF FORMAT/INFO fields.
- --vcf-gz: Write VCF output as BGZF-compressed `.MNV.vcf.gz` (requires `--convert` or `--both`).
- --index-vcf-gz: Build Tabix index (`.tbi`) for `.MNV.vcf.gz` output.
- --bcf: Also write BCF output (`.MNV.bcf`) converted from generated VCF (requires `--convert` or `--both`).
- --min-snp-strand <N>: Require at least `N` supporting reads in forward and reverse strand for SNP records.
- --min-mnv-strand <N>: Require at least `N` supporting reads in forward and reverse strand for MNV records.
- --min-strand-bias-p <P>: Filter variants with Fisher exact strand-bias p-value `< P` (`0 <= P <= 1`).
- --emit-filtered: Keep threshold-failing variants in VCF and mark them with FILTER tags instead of skipping them.
- --strand-bias-info: Add Fisher exact strand-bias p-values (`SBP` for SNP, `MSBP` for MNV) to VCF INFO.
- --keep-original-info: Preserve all original INFO fields from the input VCF in the output VCF (e.g. SnpEff `ANN`, VEP `CSQ`). Requires `--convert` or `--both`.
- --exclude-intergenic: Exclude intergenic SNPs (variants outside annotated genes) from the output. By default, intergenic variants are included with gene = `intergenic` and change type = `Unknown`.
- --summary-json <JSON_FILE>: Write structured run summary in JSON format.
- --error-json <JSON_FILE>: Write structured JSON error details when the command fails.
- --run-manifest <JSON_FILE>: Write reproducibility manifest with command line, summary and checksums.
- --convert: Output only VCF (`.MNV.vcf` or `.MNV.vcf.gz` with `--vcf-gz`) instead of TSV.
- --both: Output both TSV and VCF in one run.
## Example:
```
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --quality 30 \
  --mapq 20
```
## Input File Formats

- VCF File: Should contain the identified SNVs.
- BAM File: (Optional) Genomic reads aligned to the reference sequence.
- FASTA File: Reference genomic sequence.
- Gene annotation file:
  - TSV: tab-delimited text file with one entry per line `(GeneName, GeneStart, GeneEnd, Strand)`.
  - GFF/GFF3: features of type `gene` and `pseudogene` are used by default (columns 4/5/7 for start/end/strand). Use `--gff-features` to customize which feature types are analyzed.

Contig naming contract:
- VCF contigs, FASTA record IDs, and GFF sequence IDs must match exactly (case-sensitive).
- Input VCF should declare all used contigs in the header (`##contig=<ID=...>`).
- If contig names do not match, normalize them before running `get_mnv`.

TSV example:
```bash
Rv0007_Rv0007	9914	10828	+
ileT_Rvnt01	10887	10960	+
alaT_Rvnt02	11112	11184	+
Rv0008c_Rv0008c	11874	12311	-
ppiA_Rv0009	12468	13016	+
Rv0010c_Rv0010c	13133	13558	-
```
## 🎴Output
The program generates a TSV file named <vcf_filename>.MNV.tsv containing the following information:
- Chromosome: Contig/chromosome name.
- Gene: Name of the gene.
- Positions: Positions of the variants.
- Reference Bases: Reference nucleotide bases at variant positions.
- Base Changes: Nucleotide base changes.
- AA Changes: Resulting amino acid changes.
- SNP AA Changes: Amino acid changes if considering individual SNVs.
- Variant Type: Type of variant (SNP, MNV, or SNP/MNV).
- Change Type: Type of change at the protein level (Synonymous, Non-synonymous, Stop gained).
- SNP Reads: (If BAM provided) Count of reads supporting each SNP.
- SNP Forward Reads / SNP Reverse Reads: (If BAM provided) Strand-specific SNP support counts.
- MNV Reads: (If BAM provided) Count of reads supporting the MNV.
- MNV Forward Reads / MNV Reverse Reads: (If BAM provided) Strand-specific MNV support counts.

When `--convert` or `--both` is used, the VCF INFO field includes:
- `ODP` / `OFREQ`: Original depth and frequency from the input VCF (parsed from FORMAT/INFO `DP` + `AF/FREQ`, with caller-compatible fallback from `AD` or `AO/RO` when needed).
- `DP` / `FREQ`: Recalculated depth and frequency from BAM reads.
- `SRF` / `SRR`: SNP supporting reads on forward/reverse strand.
- `MRF` / `MRR`: MNV supporting reads on forward/reverse strand.

When `--emit-filtered` is used, VCF FILTER can contain:
- `LowSupport`: fails minimum read support thresholds (`--snp` / `--mnv`).
- `StrandSupport`: fails minimum per-strand support (`--min-snp-strand` / `--min-mnv-strand`).
- `StrandBias`: fails Fisher exact p-value threshold (`--min-strand-bias-p`).

When `--bcf` is enabled, an additional `.MNV.bcf` output is written.

VCF header metadata includes:
- `##get_mnv_version`
- `##get_mnv_command`
- `##get_mnv_min_quality`
- `##get_mnv_min_mapq`
- `##get_mnv_min_snp_reads`
- `##get_mnv_min_mnv_reads`

Notes:
- For MNV records, recalculated `DP/FREQ` use the depth of reads that span all SNP positions in the MNV haplotype.
- Recalculated frequencies (`FREQ`, SNP/MNV frequencies in TSV) are printed with 4 decimal places.
- `--dry-run` prints parsing/annotation summaries and skips TSV/VCF creation.
- `--dry-run` also skips `.bcf` creation.
- `--summary-json` stores per-contig and global counters (variants, genes, cache hits/misses).
- `summary-json` includes `schema_version` for stable downstream parsing.
- `summary-json` includes per-phase timings in milliseconds (`parse_inputs_ms`, `process_ms`, `emit_ms`, `total_ms`).
- `summary-json` includes SHA-256 checksums for VCF/FASTA/annotation/BAM inputs.
- `--run-manifest` includes summary plus output checksums (VCF/TSV/BCF when present).
- `--sample all` writes one output set per sample using suffix `.sample_<sample_name>`.
- If `--chrom` is provided, GFF/GFF3 gene features are filtered to those contigs.
- Multi-contig processing without `--chrom` is supported with GFF/GFF3 annotation files.
- TSV annotation files do not contain contig information, so for multi-contig VCF inputs use `--gff` or restrict with `--chrom`.
- `--strict` is useful for QC pipelines that require original depth/frequency traceability in every output record.
- Multiallelic records can be processed with `--split-multiallelic`; default behavior remains strict fail.
- Intergenic variants (positions outside all annotated genes) are included by default with gene = `intergenic`, no codon/AA data, and change type = `Unknown`. Use `--exclude-intergenic` to omit them.
- `--keep-original-info` carries through all non-get_mnv INFO fields from the input VCF header and records.

Example:
```
Chromosome	Gene	Positions	Base Changes	AA Changes	SNP AA Changes	Variant Type	Change Type	SNP Reads	SNP Forward Reads	SNP Reverse Reads	MNV Reads	MNV Forward Reads	MNV Reverse Reads
MTB_anc	Rv0095c_Rv0095c	104838	T	Asp126Glu	Asp126Glu	SNP	Non-synonymous	16	8	8	0	0	0
MTB_anc	Rv0095c_Rv0095c	104941,104942	T,G	Gly92Gln	Gly92Glu; Gly92Arg	SNP/MNV	Non-synonymous	0,3	0,1	0,2	22	11	11
MTB_anc	esxL_Rv1198	1341102,1341103	T,C	Arg33Ser	Arg33Cys; Arg33Pro	MNV	Non-synonymous	0,0	0,0	0,0	9	5	4
```

## Benchmark (development)
Run a baseline benchmark for codon variant generation:
```bash
cargo run --release --bin bench_variants -- --warmup 5 --iters 30
```
Benchmark against the bundled example dataset and export CSV:
```bash
cargo run --release --bin bench_variants -- --dataset example --warmup 5 --iters 20 --threads 4 --csv benchmark.csv --max-avg-ms 200
```
This appends one line per run to `benchmark.csv`.

Larger synthetic stress benchmark:
```bash
cargo run --release --bin bench_variants -- --warmup 3 --iters 20 --threads 4 --synthetic-scale 4
```

For reproducible benchmark batches and example-run artifact generation, see scripts in `/analysis`:
- `/analysis/run_reproducible_benchmark.sh`
- `/analysis/reproduce_example_run.sh`

## Common Errors
- `[E002] Contig validation failed ... Missing in VCF/FASTA`:
  contig names do not match between files. Normalize names in VCF/FASTA/GFF and re-run.
- `[E002] VCF REF/FASTA mismatch at <contig>:<pos>`:
  REF in VCF does not match the FASTA base. Normalize/fix the VCF against the same reference.
- `[E002] Invalid base 'X' in REF/ALT allele ...`:
  non-IUPAC base detected. Clean invalid alleles before running.
- `Sample '<name>' not found in VCF header`:
  check sample name spelling or omit `--sample` to use the first sample.
- `TSV annotation does not include contig names`:
  for multi-contig VCF use `--gff` or restrict with `--chrom`.
- `--strict enabled, but original VCF metrics are missing ...`:
  ensure each selected record has recoverable `DP` and `AF/FREQ` in FORMAT or INFO.
- `--index-vcf-gz requires --vcf-gz`:
  enable BGZF VCF output before requesting Tabix index generation.
- `--bcf requires --convert or --both`:
  enable VCF output mode before requesting BCF conversion.
- `--min-strand-bias-p must be between 0 and 1`:
  provide a probability threshold inside `[0, 1]`.
- `Requested --sample all but input VCF has no sample columns`:
  use a VCF with FORMAT/sample columns or select a single sample mode.
- `Multiallelic VCF record ... is not supported`:
  split multiallelic sites first (`bcftools norm -m -`) or run with `--split-multiallelic`.

Exit status codes:
- `1`: generic error (`E000`)
- `2`: config/CLI error (`E001`)
- `3`: input/validation error (`E002`)
- `10-14`: I/O/CSV/HTSlib/UTF-8/parse specific errors
  
# 📉 Limitations
- The script currently works only with SNVs compared against a reference sequence.
- Insertions and deletions that modify the reading frame are not supported in this version.
- Multiallelic VCF records require either pre-splitting (e.g. `bcftools norm -m -`) or enabling `--split-multiallelic`.
- VCF contig names must be properly declared in the VCF header and present in the FASTA file.
---
<h2 id="contributors" align="center">

✨ [Contributors]((https://github.com/PathoGenOmics-Lab/AMAP/graphs/contributors))
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
get_MNV is developed with ❤️ by:
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
      <a href="" title="Desing">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td> 
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td> 
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
---  
# Fun
## 3D model logo
Click for the stl file
<p align="center">
  <a href="https://www.printables.com/model/1030383-get_mnv-logo" target="_blank">
    <img src="https://media.printables.com/media/prints/1030383/images/7820375_62fee28c-1ef3-446a-9187-3a74e3912c09_7526c3fd-6f35-4ec1-ab2c-a8022ac592e9/thumbs/inside/1920x1440/jpg/img_3773.webp" height="350" alt="get_MNV">
  </a>
</p>

