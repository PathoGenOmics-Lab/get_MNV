# Common Recipes

Ready-to-run commands for the jobs people ask for most. Every option, with its
default and its meaning, lives in one place: the
[CLI reference](cli-reference.md).

## Basic Command

```bash
get_mnv \
  (--vcf <VCF_FILE> | --tsv <IVAR_TSV_FILE>) \
  --fasta <REFERENCE_FASTA> \
  (--gff <ANNOTATION_GFF> | --genes <ANNOTATION_TSV>)
```

Use `--vcf` for plain `.vcf` or BGZF-compressed `.vcf.gz` input and `--tsv`
for the `variants.tsv` file produced by `ivar variants`. BCF input is not
accepted directly; convert it to VCF first with `bcftools view`.

## VCF Input

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3
```

## iVar TSV Input

```bash
get_mnv \
  --tsv sample_variants.tsv \
  --fasta reference.fasta \
  --gff genes.gff3
```

## Add BAM Read Support

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3
```

Without a BAM, get_MNV annotates what the caller reported. With one, it counts
the reads itself, and only then can it tell a real codon-level haplotype from
two substitutions that never shared a molecule. See
[Linkage](linkage.md) for what that buys you.

## Filter on Recounted Read Support

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --min-snp-frequency 0.05 \
  --min-mnv-frequency 0.20
```

Frequency and read-count filters use the support get_MNV recalculates from the
BAM, not the `OFREQ`/`ODP` values in the input, so they need `--bam`. The SNP
and MNV thresholds are independent: a strong MNV haplotype survives even when
its individual substitutions fall below the SNP threshold.

## Write Both TSV and VCF

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both
```

## Build the Interactive Report

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --report sample.html
```

The report is a single self-contained HTML file, so it opens with no server and
travels as one attachment. It needs the TSV output, which is the default. `--convert` writes the VCF
*instead* of the TSV, so reach for `--both` when you want a report and a VCF.

For a cohort already processed one sample per run, build the report from the
existing outputs instead of running the pipeline again:

```bash
get_mnv --report-from run1.MNV.tsv run2.MNV.tsv --report cohort.html
```

## Analyze CDS Features in a GFF

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --gff-features CDS
```

Use `--gff-features CDS` when you want codon-aware protein annotation from CDS
features, especially for eukaryotic GFF/GTF files. CDS rows with
`transcript_id` or `Parent` are reconstructed as spliced transcript CDS models.

## Notes

- Contig names must match exactly across the variant file, FASTA, GFF, and BAM.
- iVar TSV parsing keeps passing SNV and indel rows. Indel notation such as
  `+SEQ` or `-SEQ` is converted to VCF-like anchored `REF/ALT` alleles using
  the FASTA reference.
- If you use `--genes`, the annotation TSV has no contig column. For
  multi-contig data, prefer `--gff`.
- Indel and frameshift behaviour has its own knobs, and two of their defaults
  deliberately differ from the tool's original behaviour. The reasoning is in
  [Scope and Compatibility](indel-mnv-semantics.md#tuning-knobs).
