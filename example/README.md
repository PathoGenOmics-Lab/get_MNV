# Example dataset

**English** · [Español](README.es.md)

A small *M. tuberculosis* dataset for trying get_MNV end to end, including the
GUI read viewer.

| File | What it is |
|---|---|
| `MTB_ancestor.fas` (`.fai`) | Reference, single contig `MTB_anc` (4,411,532 bp) |
| `anot_genes.txt` | Gene table: `name`, `start`, `end`, `strand` (tab-separated); its 70 RNA genes add `phase` and `biotype` |
| `G35894.var.snp.vcf` | VarScan-style VCF calls for sample G35894 |
| `G35894.var.snp.MNV.tsv` | Pre-computed get_MNV output for that VCF (no BAM) |
| `G35894.demo.bam` (`.bai`) | Tiny synthetic alignment for the read viewer (see below) |
| `make_demo_bam.py` | Regenerates `G35894.demo.bam` from the reference |
| `who_catalogue_2023.tsv` | WHO catalogue of resistance mutations (2023), groups 1 and 2, for `--catalogue` (see below) |

## Run the CLI

```bash
get_mnv \
  --vcf example/G35894.var.snp.vcf \
  --fasta example/MTB_ancestor.fas \
  --genes example/anot_genes.txt
```

This reproduces `G35894.var.snp.MNV.tsv`: 941 variants (797 SNP, 10 SNP/MNV,
134 intergenic) across 635 annotated genes.

## See the reads

`G35894.demo.bam` is a **tiny synthetic** alignment (24 reads, ~200 bp) over a
single real MNV in this sample so read support and the GUI pileup work out of
the box:

> gene `Rv2036` (+ strand), genomic `2282376` & `2282377`, codon `GTT → GCC`
> (Val93Ala). Every read carries both ALT bases, split evenly across strands.

Add `--bam` to count read support from it:

```bash
get_mnv \
  --vcf example/G35894.var.snp.vcf \
  --fasta example/MTB_ancestor.fas \
  --genes example/anot_genes.txt \
  --bam example/G35894.demo.bam --mnv 1
```

The `Rv2036` `2282376, 2282377` row then reports `MNV Reads=24` (12 forward /
12 reverse), `MNV Frequencies=1.0000`.

In the desktop GUI, load the VCF, FASTA, gene table and `G35894.demo.bam`, run,
then open that `SNP/MNV` row to see the 24 reads in the pileup.

> The demo BAM only covers that one locus, so reads render only there. For a
> real analysis, supply your own coordinate-sorted, indexed BAM.

## Regenerating the demo BAM

```bash
cd example && python3 make_demo_bam.py   # needs samtools on PATH (or $SAMTOOLS)
```

## Flag WHO resistance mutations

`who_catalogue_2023.tsv` holds the 1,379 entries of groups 1 and 2 (associated
with resistance) of the WHO catalogue of mutations in the *M. tuberculosis*
complex, second edition (2023), converted by `tools/who_catalogue_to_tsv.py`
and placed on `MTB_ancestor.fas`:

```bash
get_mnv \
  --vcf example/G35894.var.snp.vcf \
  --fasta example/MTB_ancestor.fas \
  --genes example/anot_genes.txt \
  --catalogue example/who_catalogue_2023.tsv
```

G35894 carries none of them, so every row's `Catalogue` stays `-`. How the
table is built, and what a sample that does carry them looks like, is in the
documentation's recipe and input-format pages.

The table contains information from the
[WHO mutation catalogue](https://github.com/GTB-tbsequencing/mutation-catalogue-2023),
which is made available under the
[ODC Attribution License](https://opendatacommons.org/licenses/by/1-0/).
