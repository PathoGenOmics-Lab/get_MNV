# CLI Reference

Full reference for the `get_mnv` command-line options (version 1.1.5). Run
`get_mnv --help` for the same list in your terminal.

## Synopsis

```text
get_mnv [OPTIONS] --fasta <FASTA_FILE> <--vcf <VCF_FILE>|--tsv <TSV_FILE>>
```

You must provide a reference (`--fasta`) and exactly one variant source
(`--vcf` or `--tsv`), plus a gene annotation (`--gff` or `--genes`).

## Input

| Option | Description |
|---|---|
| `-v, --vcf <FILE>` | Variant input in plain or BGZF-compressed VCF (SNVs/MNVs and indels). |
| `--tsv <FILE>` | iVar `variants.tsv` input. |
| `-b, --bam <FILE>` | Optional aligned reads for read support. Must be coordinate-sorted and indexed. |
| `-f, --fasta <FILE>` | Reference FASTA (required). |
| `--sample <NAME>` | Sample for original FORMAT metrics in a multi-sample VCF (default: first sample; `all` for all). |
| `--chrom <NAME>` | Restrict processing to one contig (default: all contigs in the input). |

## Annotation

| Option | Description |
|---|---|
| `--gff <FILE>` | Gene annotation in GFF/GFF3 format. |
| `-g, --genes <FILE>` | Simple gene table TSV: `gene,start,end,strand`. Use instead of `--gff`. |
| `--gff-features <LIST>` | Comma-separated GFF feature types to analyze (default: `gene,pseudogene`). Use `CDS` for spliced transcripts. |
| `--translation-table <N>` | NCBI genetic code (default: `11`, bacterial). Supported: 1, 2, 3, 4, 5, 6, 11, 12, 25. |
| `--exclude-intergenic` | Drop variants outside annotated genes. |

## Read support and quality

These apply only when `--bam` is provided.

| Option | Default | Description |
|---|---|---|
| `-q, --quality <N>` | `20` | Minimum base Phred quality. |
| `--min-mapq <N>` | `0` | Minimum mapping quality (MAPQ). |
| `-s, --snp <N>` | `0` | Minimum SNP-supporting reads. |
| `--min-snp-frequency <F>` | `0.0` | Minimum BAM-derived SNP allele frequency (0.0–1.0). |
| `-m, --mnv <N>` | `0` | Minimum MNV-supporting reads. |
| `--min-mnv-frequency <F>` | `0.0` | Minimum BAM-derived MNV haplotype frequency (0.0–1.0). |
| `--min-snp-strand <N>` | `0` | Minimum SNP-supporting reads on each strand. |
| `--min-mnv-strand <N>` | `0` | Minimum MNV-supporting reads on each strand. |
| `--min-strand-bias-p <F>` | `0.0` | Minimum Fisher exact p-value accepted for strand-bias metrics. |

!!! note
    Frequency and read-count filters use support recalculated from `--bam`, not
    the original `OFREQ`/`ODP` from the input. SNP and MNV filters are
    independent: a strong MNV haplotype is kept even when its individual SNVs
    fall below the SNP threshold.

## Indel tuning

| Option | Default | Description |
|---|---|---|
| `--frameshift-min-freq <F>` | `0.0` | Minimum frequency an upstream indel must reach to mark downstream SNV/MNV codons as frameshifted. Raise it on intra-host data to avoid relabelling high-frequency substitutions because of a low-frequency upstream indel on a different molecule. |
| `--indel-anchor-depth` | off | Count indel-locus depth (the `EFREQ` denominator) from reads observing the anchor base, not only reads spanning the full REF allele. Reduces depth under-counting for multi-base deletions. |
| `--phased-indel-min-reads <N>` | `1` | Minimum BAM-supporting reads to emit a phased indel/complex haplotype row. |
| `--phased-indel-min-freq <F>` | `0.0` | Minimum BAM-derived frequency to emit a phased indel/complex haplotype row. |
| `--normalize-alleles` | off | Trim shared REF/ALT prefix/suffix before processing. |
| `--split-multiallelic` | off | Split multiallelic VCF records into independent ALT alleles instead of failing. |

## Output

| Option | Description |
|---|---|
| `--convert` | Write VCF output (`.MNV.vcf`) instead of TSV. |
| `--both` | Write both TSV and VCF in one run. |
| `--vcf-gz` | Write BGZF-compressed `.vcf.gz` (VCF output mode). |
| `--index-vcf-gz` | Build a Tabix `.tbi` index (requires `--vcf-gz`). |
| `--bcf` | Also write a BCF converted from the generated VCF (requires `--convert`/`--both`). |
| `--strand-bias-info` | Add Fisher exact strand-bias p-values to VCF INFO (`SBP`/`MSBP`). |
| `--keep-original-info` | Preserve original VCF INFO fields in the output (requires `--convert`/`--both`). |
| `--emit-filtered` | Emit records that fail thresholds with `FILTER` tags instead of skipping them. |

## Validation and metadata

| Option | Description |
|---|---|
| `--dry-run` | Validate inputs and print a per-contig summary without writing outputs. |
| `--strict` | Fail if original depth/frequency metrics (`ODP`/`OFREQ`) are missing in the input. |
| `--summary-json <FILE>` | Write a machine-readable run summary. |
| `--run-manifest <FILE>` | Write a reproducibility manifest (inputs, outputs, checksums, runtime metadata). |
| `--error-json <FILE>` | Write structured error details as JSON when the command fails. |
| `--threads <N>` | Number of worker threads (default: Rayon auto). |
| `-h, --help` | Print help. |
| `-V, --version` | Print version. |
