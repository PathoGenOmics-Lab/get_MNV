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

- VCF (`.vcf` or `.vcf.gz`) and BCF (`.bcf`)
- iVar variants TSV (`.tsv`)

Use `--vcf` for plain `.vcf`, BGZF-compressed `.vcf.gz` or `.bcf` files and
`--tsv` for iVar `variants.tsv` files. A BCF is decoded to VCF text as it is
read, so it annotates exactly as the VCF it encodes; it is recognised by its
`.bcf` extension.
Older commands that pass an iVar TSV through `--vcf` are still auto-detected
when the header has the standard iVar columns.

### VCF

Use a standard VCF file containing SNV/MNV calls, indels, or complex alleles.

Requirements:

- VCF contig names must match the FASTA and GFF/GTF contig names.
- REF alleles must match the FASTA sequence.
- Multiallelic records are split into one row per ALT allele.

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

Optional seventh column with the contig:

```text
GeneName	GeneStart	GeneEnd	Strand	Phase	Biotype	Contig
```

```text
Rv0007_Rv0007	9914	10828	+	0	protein_coding	NC_000962.3
geneD	1	300	+	0	protein_coding	plasmid1
```

A row that names a contig is only applied to variants on that contig. A row
without one, as in the shorter formats, is applied to whichever contig is being
annotated, which is all a single-contig file needs. A VCF with more than one
contig therefore needs every row to name its contig; otherwise the run stops
rather than lay one contig's genes over another's variants. Using it means
writing the phase and biotype columns too: `.` is accepted for an unknown
phase.

## 4. BAM or CRAM Reads (Optional)

Pass aligned reads with:

```bash
--bam reads.bam
--bam reads.cram
```

When reads are provided, get_MNV calculates:

- SNP read support
- MNV haplotype read support
- Total depth and frequency
- Forward/reverse strand counts
- Optional strand-bias statistics

Requirements:

- The file must be coordinate-sorted.
- It must be indexed: `.bai` (or `.csi`) for a BAM, `.crai` for a CRAM.
  `samtools index` writes either.
- Contig names must match the variant file and FASTA.
- A CRAM stores reads as differences from a reference, so it is decoded
  against the `--fasta` of the run. That has to be the reference the CRAM was
  written against; a different one fails with a checksum error rather than
  giving wrong bases. A `.fai` next to the FASTA lets the sequence be read on
  demand; without one it is loaded whole.

A CRAM gives the same counts as the BAM it was made from: the test suite runs
every scenario both ways and compares the outputs byte for byte.
- Duplicate, secondary, supplementary and QC-fail reads are ignored, matching what `samtools mpileup` excludes by default.

## 5. Catalogue of Known Variants (Optional)

`--catalogue known.tsv` flags calls that match a table of variants you already
know something about: drug-resistance mutations, lineage markers, known
pathogenic changes. It is tab-separated, with a header naming a `gene` and a
`variant` column and, optionally, a `label`; other columns are ignored and
lines starting with `#` are comments.

```text
gene	variant	label
rpoB	p.Ser450Leu	RIF: Assoc w R
fabG1	c.-15C>T	INH: Assoc w R
katG	LoF	INH: Assoc w R
```

A variant can be written as:

| Form | Example | Matches |
|---|---|---|
| Protein | `p.Ser450Leu`, `p.S450L`, `p.Gln10Ter`, `p.Ala11fs` | rows of the gene with that protein change. A frameshift matches by its first residue and position, whatever follows `fs`. |
| Coding substitution | `c.1349C>T`, `c.-15C>T`, `c.*3G>A` | the substituted base, placed in the genome from the gene model (`c.-N` upstream of the first coding base, `c.*N` past the last), on any row that carries it, intergenic ones included |
| Non-coding substitution | `n.1401A>G`, `n.-119C>T` | the same, in a gene that is not translated, such as an rRNA, numbered from its first base |
| Genomic substitution | `g.761155C>T` | the substituted base at that position |
| Other HGVS | `c.516_518dup`, `g.1673425del` | rows of the gene whose `HGVS c.` is written exactly so, or rows anywhere on the gene's contig, intergenic ones included, whose `HGVS g.` is |
| `LoF` | `LoF` | any frameshift, premature stop or start loss in the gene |

The gene is matched by name, or by any `_`-separated part of a combined name, so
`rpoB` matches a gene named `rpoB_Rv0667`. A substitution whose reference base
is not the one the FASTA has at that position was written against another
reference; it matches nothing, and the run says how many there were.

A protein entry is compared with both readings of a codon row: the residue the
whole codon gives and the one each substitution would give alone. When only one
reading matches, the row says so, `codon reading only` for a catalogued change
that no single base makes, which base-by-base annotation misses, and
`per-SNV reading only` for one such annotation would report that the codon does
not make. Matches go to the `Catalogue` column (see
[Output formats](output-formats.md)).

### The WHO Catalogue of *M. tuberculosis* Resistance Mutations

The WHO catalogue of mutations in the *Mycobacterium tuberculosis* complex
(second edition, 2023) grades about 48,000 variant and drug pairs into five
groups, from `1) Assoc w R` to `5) Not assoc w R`. `example/who_catalogue_2023.tsv`
holds groups 1 and 2, the variants associated with resistance, ready for
`--catalogue`; [Common recipes](usage.md#flag-known-resistance-mutations-in-m-tuberculosis)
shows the run. Each entry is labelled with the drug's WHO abbreviation and the
grade, so a match reads `RIF: Assoc w R (rpoB p.Ser450Leu)`, and a variant
graded for several drugs is one entry per drug. Three more columns keep WHO's
drug name, grade and variant name (`rpoB_p.Ser450Leu`), for filtering the table
or joining it back to the catalogue; get_MNV ignores them.

`tools/who_catalogue_to_tsv.py` builds the table from WHO's own files, so it can
be rebuilt with other groups (`--grades 1,2,5` adds the variants shown not to be
associated, `--grades all` keeps every group) or placed on another reference.
It needs only Python's standard library. Most entries are written the way
get_MNV writes the same change, and are kept; the rest are rewritten:

| WHO writes | The table holds |
|---|---|
| `p.Ser450Leu`, `p.Trp122*`, `p.Asp379fs`, `p.Met1?`, an in-frame `del`, `ins`, `dup` or `delins` | the same |
| `c.-15C>T`, `c.102G>A`, `n.1401A>G`, `LoF` | the same |
| `p.Ter187Argext*?`, a lost stop | `p.Ter187Arg`, the residue that replaces it |
| `c.421_423delTTGinsCTC`, the same length in and out | `c.[421T>C;423G>C]`, the substitutions it makes |
| `c.-5delG` or `n.2269_2270insT`, an indel outside a coding sequence | the `HGVS g.` get_MNV gives the change, `g.2289247del`, from WHO's genomic coordinates and the run's FASTA (`--coordinates`, `--fasta`) |
| `deletion`, the whole gene | nothing: WHO gives it no genomic variant, and a VCF of changes inside genes has no row for a gene the sample lacks |

A rewritten entry keeps WHO's name in its label: `PZA: Assoc w R, WHO c.-5delG
(pncA g.2289247del)`. `--unconverted` writes what was not converted, with the
reason. For groups 1 and 2 that is the four whole-gene deletions, of katG, ethA,
fgd1 and pncA, which each have a `LoF` entry for the same drug; without
`--coordinates`, the three indels outside coding sequences (pncA `c.-5delG`, eis
`c.-8delC`, rrl `n.2269_2270insT`) are listed as well.

Matching protein changes, rather than genomic records, is what lets a codon with
two substitutions be found. 22 of the rpoB variants in groups 1 and 2, among
them `p.Ser450Phe`, `p.Ser450Tyr`, `p.Ser441Gln` and `p.His445Cys`, need two
changes in one codon; read base by base, `p.Ser450Phe` is `p.Ser450Leu` and a
synonymous change. get_MNV flags each of them with `codon reading only`.

WHO also lists, for every graded variant, the genomic changes that make it:
26,231 of them for groups 1 and 2. Each was annotated on the example's MTBC
ancestor with the shipped table. 26,223 come out flagged with their own entry,
198 of them through the codon reading only, and 2 more with another entry for
the same drug, the gene's `LoF`. The 6 left are 2 insertions just before a
start codon, which WHO counts as losses of function and get_MNV as intergenic;
2 changes whose consequence on the ancestor is not the one they have on H37Rv;
and 2 deletions WHO links to fbiC `p.Ala855fs` that the bases after them
repeat, which HGVS places past fbiC's stop codon, leaving its coding sequence
whole. Given the folder that holds WHO's two files, `tests/who_catalogue/run.py`
repeats the check.

With your own data:

- **Gene names.** WHO names genes by their H37Rv symbols and locus tags (`rpoB`,
  `katG`, `Rv0678`, `rrs`), and an entry matches a gene of that name or with it
  as one `_`-separated part. All 65 genes of the catalogue are found that way in
  `example/anot_genes.txt` and in NCBI's annotation of H37Rv (GCF_000195955.2),
  whose genes get_MNV names `rpoB_Rv0667`. With an annotation that names them
  otherwise, change the table's `gene` column to the annotation's names.
- **Reference.** The catalogue is written against H37Rv (NC_000962.3). Protein,
  `c.` and `n.` entries are numbered along their gene, so they hold on any
  reference whose genes match H37Rv's. The `g.` entries are H37Rv positions:
  build the table with your own FASTA in `--fasta` when it is on H37Rv
  coordinates, and leave `--coordinates` out when it is not. The MTBC ancestor
  of the example shares H37Rv's coordinates but not the alleles H37Rv's own
  lineage carries: 38 entries name a residue or base the ancestor does not have
  there, none in groups 1 and 2 (`embB p.Glu378Ala`, a lineage marker, is one).
  Against the ancestor they match nothing, and a sample carrying H37Rv's allele
  shows the reverse change, `Ala378Glu`.

The table and the files it comes from contain information from the
[WHO mutation catalogue](https://github.com/GTB-tbsequencing/mutation-catalogue-2023),
which is made available under the
[ODC Attribution License](https://opendatacommons.org/licenses/by/1-0/); cite
it as *Catalogue of mutations in Mycobacterium tuberculosis complex and their
association with drug resistance, second edition. Geneva: World Health
Organization; 2023*.

## Contig Names

All input files must agree on contig names:

```text
variant contig == FASTA record ID == GFF sequence ID == BAM reference name
```

For example, `chr1` and `1` are different names. Rename or normalize inputs
before running get_MNV.
