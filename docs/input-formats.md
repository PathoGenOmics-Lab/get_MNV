# Input Formats

get_MNV needs three files: variant calls, a reference FASTA, and gene
annotation. A BAM file is optional.

## 1. Variant Calls

Pass variant calls with:

```bash
--vcf <FILE>
```

The option name is kept for backwards compatibility. The file can be:

- VCF (`.vcf` or `.vcf.gz`)
- iVar variants TSV (`.tsv`)

Use one of:

```bash
--input-format auto   # default
--input-format vcf
--input-format ivar
```

### VCF

Use a standard VCF file containing SNV calls.

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
- iVar indel notation such as `+A` or `-N` is skipped for now. get_MNV is
  SNV-based.
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
- For multi-exon CDS annotations, amino acid numbering is reported against the
  full transcript when transcript information is available.
- If a GFF/GTF contains multiple transcripts for the same gene, one variant can
  produce one output line per overlapping transcript.

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
- Duplicate, secondary, and supplementary reads are ignored.

## Contig Names

All input files must agree on contig names:

```text
variant contig == FASTA record ID == GFF sequence ID == BAM reference name
```

For example, `chr1` and `1` are different names. Rename or normalize inputs
before running get_MNV.
