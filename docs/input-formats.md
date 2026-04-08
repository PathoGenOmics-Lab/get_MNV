# Input formats

## VCF file (required)

Standard VCF (v4.x) containing identified SNVs.

**Requirements:**
- Contigs must be declared in the header (`##contig=<ID=...>`)
- Contig names must match exactly (case-sensitive) across VCF, FASTA, and GFF
- Multiallelic records require `--split-multiallelic` or pre-splitting with `bcftools norm -m -`

**Supported metrics extraction** (for `ODP`/`OFREQ` fields):
- FORMAT: `DP`, `AF`, `FREQ`, `AD`, `AO`/`RO`
- INFO: `DP`, `AF`, `FREQ`, `AD`, `AO`/`RO`
- Automatic fallback chain for maximum caller compatibility

## FASTA file (required)

Reference genomic sequence. Record IDs must match VCF contig names exactly.

**Requirements:**
- Valid IUPAC DNA bases (A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N)
- No duplicate contig names
- UTF-8 encoding

## Gene annotation (required — one of the following)

### GFF/GFF3 format (recommended)

Standard GFF3 with features of type `gene` and `pseudogene` by default.

- Columns 4/5/7 are used for start, end, and strand
- **Column 8 (phase) is honoured for `CDS` features** (since v1.1.2). Valid
  values are `0`, `1`, `2` or `.`. The phase is the number of bases that must
  be skipped from the feature start (or from the feature end on the minus
  strand) to reach the first base of the first complete codon. This is
  required for eukaryotic annotations whose CDS exons begin mid-codon — for
  example GRCh38 Ensembl `GNAQ` exon 5 has `phase=1`, and without this
  correction SNVs were grouped into the wrong codons. Features without phase
  information (`gene`, `exon`, `UTR`, …) implicitly use `phase = 0`, which is
  the historical behaviour.
- Gene name is extracted from attributes: `gene` > `Name` > `locus_tag` > `ID`
- Use `--gff-features CDS,tRNA` to customize which feature types are analyzed
- Supports percent-encoded attribute values and quoted semicolons
- Multi-contig aware — genes are automatically filtered per contig

### TSV format

Tab-delimited file with one gene per line.

#### Prokaryote (4 columns)

The legacy 4-column format. Every CDS starts in frame 0, so no phase is
needed:

```
GeneName	GeneStart	GeneEnd	Strand
```

Example:
```
Rv0007_Rv0007	9914	10828	+
ileT_Rvnt01	10887	10960	+
Rv0008c_Rv0008c	11874	12311	-
```

#### Eukaryote (5 columns, since v1.1.2)

For eukaryotic CDS exons that do not start in frame 0, add a 5th column with
the GFF-style phase (`0`, `1`, `2` or `.`):

```
GeneName	GeneStart	GeneEnd	Strand	Phase
```

Example (GRCh38 GNAQ exon 5, minus strand, phase=1):
```
GNAQ_exon5	77794463	77794592	-	1
```

- The phase column is **optional**. When omitted (4-column rows), phase
  defaults to `0` — i.e. the prokaryote behaviour above.
- Valid values: `0`, `1`, `2` or `.` (dot is treated as `0`).
- 4-column and 5-column rows can be mixed in the same file.

**Limitations:**
- No contig information — for multi-contig VCF use `--gff` or restrict with `--chrom`
- Coordinates must be 1-based, positive integers
- Strand must be `+` or `-`

## BAM file (optional)

Indexed BAM file with aligned reads. When provided, get_MNV calculates:

- Per-position read support counts (SNP and MNV)
- Strand-specific counts (forward/reverse)
- Recalculated depth and frequency from reads

**Requirements:**
- Must be sorted and indexed (`.bai`)
- Contig names must match VCF
- Duplicate, secondary, and supplementary reads are automatically excluded

## Contig naming contract

All input files must use consistent contig names:

```
VCF contig == FASTA record ID == GFF sequence ID
```

If names don't match, normalize them before running get_MNV (e.g., `chr1` vs `1`).
