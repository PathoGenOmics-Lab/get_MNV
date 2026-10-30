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

Use `--vcf` for plain `.vcf`, BGZF-compressed `.vcf.gz` or `.bcf` input and
`--tsv` for the `variants.tsv` file produced by `ivar variants`.

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

## Flag Known Resistance Mutations in *M. tuberculosis*

`example/who_catalogue_2023.tsv` is the WHO catalogue of *M. tuberculosis*
resistance mutations (2023), its groups 1 and 2, in the form `--catalogue`
reads. On the example's reference and gene table:

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta example/MTB_ancestor.fas \
  --genes example/anot_genes.txt \
  --catalogue example/who_catalogue_2023.tsv
```

A call that matches gets the entry in its `Catalogue` column. A codon with two
substitutions is matched as the codon it makes, which base-by-base annotation
does not see: rpoB `p.Ser450Phe` needs two changes in codon 450, and its row
reads

```text
RIF: Assoc w R (rpoB p.Ser450Leu, per-SNV reading only); RIF: Assoc w R (rpoB p.Ser450Phe, codon reading only)
```

that is, the codon is Ser450Phe, and an annotator reading one base at a time
would have reported Ser450Leu. For `p.Ser441Gln` it would have reported
Ser441Pro and a premature stop. To keep other groups, or to place the table on
another reference, build it from WHO's files:

```bash
base="https://raw.githubusercontent.com/GTB-tbsequencing/mutation-catalogue-2023/main/Final%20Result%20Files"
curl -LO "$base/WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
curl -LO "$base/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
python3 tools/who_catalogue_to_tsv.py WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt \
  --grades 1,2,5 \
  --coordinates WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt \
  --fasta example/MTB_ancestor.fas \
  --unconverted who_unconverted.tsv \
  -o who_catalogue.tsv
```

What the converter rewrites, what it cannot convert, and what a gene table or
reference of your own needs are in
[Input formats](input-formats.md#the-who-catalogue-of-m-tuberculosis-resistance-mutations).

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
- If you use `--genes` with multi-contig data, name each gene's contig in the
  seventh column of the annotation TSV, or use `--gff`.
- Indel and frameshift behaviour has its own knobs, and two of their defaults
  deliberately differ from the tool's original behaviour. The reasoning is in
  [Scope and Compatibility](indel-mnv-semantics.md#tuning-knobs).
