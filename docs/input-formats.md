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
- Gene name is extracted from attributes: `gene` > `Name` > `locus_tag` > `ID`
- Use `--gff-features CDS,tRNA` to customize which feature types are analyzed
- Supports percent-encoded attribute values and quoted semicolons
- Multi-contig aware — genes are automatically filtered per contig

### TSV format

Tab-delimited file with one gene per line:

```
GeneName	GeneStart	GeneEnd	Strand
```

Example:
```
Rv0007_Rv0007	9914	10828	+
ileT_Rvnt01	10887	10960	+
Rv0008c_Rv0008c	11874	12311	-
```

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
