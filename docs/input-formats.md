# Input Formats

get_MNV needs three files: variant calls, a reference FASTA, and gene
annotation. A BAM file is optional.

## 1. Variant Calls

Pass variant calls with:

```bash
--vcf <VCF_FILE>
# or
--tsv <IVAR_TSV_FILE>
```

The supported variant-call inputs are:

- VCF (`.vcf` or `.vcf.gz`)
- iVar variants TSV (`.tsv`)

Use `--vcf` for plain `.vcf` or BGZF-compressed `.vcf.gz` files and `--tsv`
for iVar `variants.tsv` files. BCF input is not accepted directly; convert it
first, for example with `bcftools view input.bcf > input.vcf`.
Older commands that pass an iVar TSV through `--vcf` are still auto-detected
when the header has the standard iVar columns.

### VCF

Use a standard VCF file containing SNV/MNV calls, indels, or complex alleles.

Requirements:

- VCF contig names must match the FASTA and GFF/GTF contig names.
- REF alleles must match the FASTA sequence.
- Multiallelic records should be pre-split, or run with
  `--split-multiallelic`.

get_MNV can read original depth/frequency metrics from common INFO or FORMAT
fields, including `DP`, `AF`, `FREQ`, `AD`, `AO`, and `RO`.

These input frequency values are kept for reporting as `OFREQ`. Command-line
frequency filters (`--min-snp-frequency`, `--min-mnv-frequency`) use
BAM-derived read support instead, so they require `--bam`.

### iVar TSV

Use the TSV produced by `ivar variants`.

Required columns:

| Column | Meaning |
|---|---|
| `REGION` | Contig name |
| `POS` | 1-based position |
| `REF` | Reference base |
| `ALT` | Alternative base |

Optional columns used when present:

| Column | Used as |
|---|---|
| `TOTAL_DP` | Original depth (`ODP`) |
| `ALT_FREQ` | Original frequency (`OFREQ`) |
| `REF_DP`, `ALT_DP` | Used to infer depth/frequency if needed |
| `PASS` | Used to keep passing rows |

Filtering:

- If `PASS` exists, get_MNV keeps truthy values such as `TRUE`, `PASS`, `1`,
  or `YES`.
- Rows where `REF == ALT` are skipped.
- iVar indel notation such as `+SEQ` or `-SEQ` is converted to VCF-like
  anchored alleles using the FASTA reference, then analysed with the same
  allele-event model as VCF input.
- `ALT_FREQ` is reported as original frequency (`OFREQ`). It is separate from
  the BAM-derived frequency filters.

## 2. Reference FASTA

Pass the reference with:

```bash
--fasta reference.fasta
```

Requirements:

- FASTA record IDs must match the variant contig names.
- Bases must be valid IUPAC DNA bases.
- Duplicate contig names are not allowed.

## 3. Gene Annotation

Provide either `--gff` or `--genes`.

### GFF/GFF3/GTF

Recommended for most datasets:

```bash
--gff genes.gff3
```

By default, get_MNV analyzes `gene,pseudogene` features. For protein-coding
codon annotation, use:

```bash
--gff-features CDS
```

Important details:

- Coordinates are read from columns 4 and 5.
- Strand is read from column 7.
- For `CDS` features, phase from column 8 is used when present.
- For `CDS` rows with `transcript_id` or `Parent`, get_MNV builds the spliced
  CDS sequence for each transcript. Codon grouping, MNV amino-acid effects, and
  indel frameshift context are then evaluated on the full transcript CDS.
- If a GFF/GTF contains multiple transcripts for the same gene, one variant can
  produce one output line per overlapping transcript. Each is annotated on its
  own terms, which is the honest answer when the annotation offers several. For
  one line per variant, keep a single transcript per gene in the GFF before the
  run, for example with [AGAT](https://github.com/NBISweden/AGAT)'s
  `agat_sp_keep_longest_isoform.pl`.

Gene names are read from common attributes such as `gene_name`, `gene`, `Name`,
`locus_tag`, `gene_id`, and `ID`.

### Simple TSV Annotation

Use `--genes` for a small, simple annotation file:

```bash
--genes genes.tsv
```

Four-column format:

```text
GeneName	GeneStart	GeneEnd	Strand
```

Example:

```text
Rv0007_Rv0007	9914	10828	+
Rv0008c_Rv0008c	11874	12311	-
```

Optional five-column format with phase:

```text
GeneName	GeneStart	GeneEnd	Strand	Phase
```

Phase can be `0`, `1`, `2`, or `.`. If the phase column is omitted, it defaults
to `0`.

Optional six-column format with biotype:

```text
GeneName	GeneStart	GeneEnd	Strand	Phase	Biotype
```

```text
Rv0007_Rv0007	9914	10828	+	0	protein_coding
mcr11_RVnc0013	1413094	1413224	-	0	ncRNA
```

Accepted values are `protein_coding`, `coding`, `CDS` and `mRNA` for translated
features, and `ncRNA`, `rRNA`, `tRNA`, `tmRNA`, `miRNA`, `snRNA`, `snoRNA`,
`misc_RNA`, `antisense_RNA`, `SRP_RNA`, `RNase_P_RNA`, `non_coding` and
`pseudogene` for features that are not. Anything else is rejected with an error
rather than guessed, because guessing wrong either invents a protein or hides a
real one.

A non-coding feature is reported against its gene as
`non_coding_transcript_exon_variant` (`MODIFIER`) with no amino-acid change.
**When the biotype column is omitted every feature is assumed to be
protein-coding**, which is what four- and five-column files have always done: an
RNA gene is then translated as though it were a protein. Declare the biotype
when your annotation contains non-coding features.

GFF and GTF input says what it is on its own. A record's `gene_biotype`,
`biotype` or `transcript_biotype` attribute decides, read through the same
vocabulary as the column above; failing that, the feature type in column 3 does,
so `rRNA`, `tRNA`, `ncRNA` and `pseudogene` rows are read as non-coding without
being told. A feature that says nothing about itself is assumed to be
protein-coding, which is what a bacterial `gene` row is. Note that the default
`--gff-features` list is `gene,pseudogene`, so pseudogenes are selected by
default and are not translated.

Limitations of TSV annotation:

- It has no contig column.
- For multi-contig data, use GFF/GTF or restrict the run with `--chrom`.

## 4. BAM Reads (Optional)

Pass BAM reads with:

```bash
--bam reads.bam
```

When a BAM is provided, get_MNV calculates:

- SNP read support
- MNV haplotype read support
- Total depth and frequency
- Forward/reverse strand counts
- Optional strand-bias statistics

Requirements:

- BAM must be sorted.
- BAM must be indexed (`.bai`).
- BAM contig names must match the variant file and FASTA.
- Duplicate, secondary, supplementary and QC-fail reads are ignored, matching what `samtools mpileup` excludes by default.

## Contig Names

All input files must agree on contig names:

```text
variant contig == FASTA record ID == GFF sequence ID == BAM reference name
```

For example, `chr1` and `1` are different names. Rename or normalize inputs
before running get_MNV.
