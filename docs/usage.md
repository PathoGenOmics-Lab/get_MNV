# Usage

This page shows the most common commands and what the main arguments mean.

## Basic Command

```bash
get_mnv \
  (--vcf <VCF_FILE> | --tsv <IVAR_TSV_FILE>) \
  --fasta <REFERENCE_FASTA> \
  (--gff <ANNOTATION_GFF> | --genes <ANNOTATION_TSV>)
```

Use `--vcf` for VCF/BCF input and `--tsv` for the `variants.tsv` file produced
by `ivar variants`.

## Common Recipes

### VCF Input

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3
```

### iVar TSV Input

```bash
get_mnv \
  --tsv sample_variants.tsv \
  --fasta reference.fasta \
  --gff genes.gff3
```

### Add BAM Read Support

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3
```

### Write Both TSV and VCF

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both
```

### Analyze CDS Features in a GFF

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --gff-features CDS
```

Use `--gff-features CDS` when you want codon-aware protein annotation from CDS
features, especially for eukaryotic GFF/GTF files.

## Required Arguments

| Argument | Meaning |
|---|---|
| `--vcf <FILE>` | Variant calls in VCF/BCF format. |
| `--tsv <FILE>` | iVar `variants.tsv` calls. |
| `--fasta <FILE>` | Reference FASTA used to call the variants. |
| `--gff <FILE>` | Gene annotation in GFF/GFF3/GTF format. |
| `--genes <FILE>` | Simple gene annotation TSV. Use this instead of `--gff`. |

You must provide either `--gff` or `--genes`.

## Input Arguments

| Argument | Default | Meaning |
|---|---:|---|
| `--bam <FILE>` | none | Sorted and indexed BAM used to count read support. |
| `--sample <NAME>` | first sample | Sample to read from a multi-sample VCF. Use `all` for every sample. |
| `--chrom <NAME>` | all contigs | Restrict the run to one contig. |
| `--gff-features <LIST>` | `gene,pseudogene` | Feature types to analyze from GFF/GTF. |
| `--translation-table <N>` | `11` | NCBI genetic code table. Supported: `1,2,3,4,5,6,11,12,25`. |

## Filter Arguments

| Argument | Default | Meaning |
|---|---:|---|
| `--quality <N>` | `20` | Minimum variant quality. |
| `--min-mapq <N>` | `0` | Minimum mapping quality for BAM reads. |
| `--snp <N>` | `0` | Minimum SNP-supporting reads. |
| `--mnv <N>` | `0` | Minimum MNV-supporting reads. |
| `--min-snp-frequency <F>` | `0` | Minimum BAM-derived SNP allele frequency (`0` to `1`). |
| `--min-mnv-frequency <F>` | `0` | Minimum BAM-derived MNV haplotype frequency (`0` to `1`). |
| `--min-snp-strand <N>` | `0` | Minimum SNP reads on each strand. |
| `--min-mnv-strand <N>` | `0` | Minimum MNV reads on each strand. |
| `--min-strand-bias-p <P>` | `0` | Minimum Fisher exact p-value for strand-bias filtering. |
| `--strict` | off | Fail when original depth/frequency metrics are missing. |

Frequency filters require `--bam` because they use read support recalculated by
get_MNV. They do not filter on the original input `OFREQ` value. For example,
`--min-snp-frequency 0.05` keeps SNP records at 5% or higher, and
`--min-mnv-frequency 0.20` keeps MNV haplotypes at 20% or higher.
These two thresholds are independent: in mixed `SNP/MNV` calls, the SNP
threshold does not remove a strong MNV haplotype, and the MNV threshold does not
remove SNP observations that pass the SNP threshold.
Read-count and strand-support filters follow the same rule: `--snp` and
`--min-snp-strand` apply to SNP observations, while `--mnv` and
`--min-mnv-strand` apply to MNV haplotypes.

## Output Arguments

| Argument | Meaning |
|---|---|
| `--convert` | Write VCF instead of TSV. |
| `--both` | Write both TSV and VCF. |
| `--vcf-gz` | Write compressed `.MNV.vcf.gz` output. |
| `--index-vcf-gz` | Create a Tabix index for `.MNV.vcf.gz`. |
| `--bcf` | Also write BCF output. Requires VCF output. |
| `--emit-filtered` | In VCF output, keep records that fail filters and mark them in `FILTER`. TSV output still omits failed rows. |
| `--strand-bias-info` | Add strand-bias p-values to VCF INFO fields. |
| `--keep-original-info` | Preserve non-get_MNV INFO fields from the input VCF. Requires VCF output. |
| `--exclude-intergenic` | Skip variants outside annotated features. |
| `--summary-json <FILE>` | Write a JSON run summary. |
| `--error-json <FILE>` | Write JSON error details if the run fails. |
| `--run-manifest <FILE>` | Write command, version, inputs, outputs, and checksums. |

## Utility Arguments

| Argument | Meaning |
|---|---|
| `--dry-run` | Validate inputs without writing output files. |
| `--threads <N>` | Number of worker threads. Default: automatic. |
| `--normalize-alleles` | Trim shared REF/ALT context before processing. |
| `--split-multiallelic` | Split multiallelic VCF records inside get_MNV. |

## Notes

- Contig names must match exactly across the variant file, FASTA, GFF, and BAM.
- iVar TSV parsing keeps passing SNV rows and skips iVar indel notation such as
  `+A` or `-N`.
- If you use `--genes`, the annotation TSV has no contig column. For
  multi-contig data, prefer `--gff`.
