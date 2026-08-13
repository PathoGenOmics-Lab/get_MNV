# CLI Reference

Full reference for the `get_mnv` command-line options (version 1.1.5). Run
`get_mnv --help` for the same list in your terminal.

## Synopsis

```text
get_mnv [OPTIONS] --fasta <FASTA_FILE> <--vcf <VCF_FILE>|--tsv <TSV_FILE>>
```

You must provide a reference (`--fasta`) and exactly one variant source
(`--vcf` or `--tsv`), plus a gene annotation (`--gff` or `--genes`). The one
exception is `--report-from`, which builds a report out of existing TSV outputs
and needs none of them.

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

!!! note "How the SNP and MNV thresholds combine"
    Frequency and read-count filters use support recalculated from `--bam`, not
    the original `OFREQ`/`ODP` from the input.

    A codon-level (`SNP/MNV`) row is kept when **either** side clears its bar:
    its individual SNVs pass the SNP thresholds, **or** its haplotype passes the
    MNV thresholds. That keeps a well-supported haplotype whose individual SNVs
    are weak, and it also works the other way: with the SNP thresholds left at
    their default of `0` the SNP side always passes, so raising `--mnv` alone
    removes nothing. Raise both, or neither.

    Which thresholds govern which row:

    | Row | Judged by |
    |---|---|
    | `SNP` | the SNP thresholds |
    | `MNV`, `SNP/MNV` | either side, as above |
    | `INDEL` | the **MNV** thresholds, measured against the indel's own event support (`Event Reads`, `Event Forward/Reverse Reads`, `Event Depth`), not against any SNP column |

    So an indel is filtered with `--mnv`, `--min-mnv-frequency` and
    `--min-mnv-strand`; `--snp` never touches it. An intergenic indel carries no
    recomputed support at all, so it cannot satisfy a read-based filter and is
    dropped whenever any MNV threshold is active.

## Indel tuning

| Option | Default | Description |
|---|---|---|
| `--frameshift-min-freq <F>` | `0.5` | Minimum frequency an upstream indel must reach to mark downstream SNV/MNV codons as frameshifted. The default propagates only from a majority upstream indel; set `0.0` to propagate from every one. Indels with no known frequency always propagate. |
| `--legacy-indel-depth` | off | Restrict indel-locus depth (the `EFREQ` denominator) to reads spanning the whole REF allele. By default it is counted from reads observing the anchor base, which avoids under-counting depth on multi-base deletions; this flag restores the older, narrower denominator. |
| `--phased-indel-min-reads <N>` | `2` | Minimum BAM-supporting reads to emit a phased indel/complex haplotype row. One read is not evidence of a haplotype. |
| `--count-mates-separately` | off | Count the two mates of a paired-end fragment as two observations instead of one molecule. |
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

## Interactive HTML report

| Option | Description |
|---|---|
| `--report <HTML_FILE>` | Write a self-contained interactive HTML report of the called variants. Needs the TSV output (the default; with `--convert` add `--both`). With `--sample all` the report covers every sample. |
| `--report-from <TSV>...` | Build the report from existing get_MNV TSV files instead of running the pipeline, for cohorts processed one sample per run. Each file becomes one sample, labelled by its file name. Requires `--report` for the output path. |

See [Output formats](output-formats.md#interactive-html-report) for what the report contains.

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
