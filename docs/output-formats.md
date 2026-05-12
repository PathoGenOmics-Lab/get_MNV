# Output Formats

get_MNV can write TSV, VCF, BCF, and JSON metadata files.

## Default TSV Output

Default file name:

```text
<input_name>.MNV.tsv
```

Use this format for spreadsheets, downstream parsing, and quick inspection.

Main columns:

| Column | Meaning |
|---|---|
| `Chromosome` | Contig name |
| `Gene` | Gene or feature name. Intergenic variants are marked as `intergenic`. |
| `Positions` | One position for SNPs, multiple comma-separated positions for MNVs. |
| `Reference Bases` | Reference bases at those positions. |
| `Base Changes` | Alternative bases. |
| `AA Changes` | Amino acid change after combining all SNVs in the codon. |
| `SNP AA Changes` | Amino acid change for each SNV considered separately. |
| `Local AA Changes` | Exon/local numbering, useful for older downstream workflows. |
| `Local SNP AA Changes` | Per-SNP amino acid changes in local numbering. |
| `Variant Type` | `SNP`, `MNV`, `SNP/MNV`, or `INDEL`. |
| `Change Type` | Synonymous, non-synonymous, stop gained/lost, unknown, etc. |
| `Reference Codon` | Original codon. |
| `SNP Codon` | Codon with individual SNP substitutions. |
| `MNV Codon` | Codon with all grouped substitutions. |

Extra columns when `--bam` is used:

| Column | Meaning |
|---|---|
| `SNP Reads` | Reads supporting each individual SNV. |
| `SNP Forward/Reverse Reads` | Strand-specific SNP support. |
| `MNV Reads` | Reads supporting the full MNV haplotype. |
| `MNV Forward/Reverse Reads` | Strand-specific MNV support. |
| `Total Reads` | Depth at the variant positions. |
| `SNP Frequencies` | Per-position SNP frequencies. |
| `MNV Frequencies` | MNV haplotype frequency. |

Frequency columns are calculated from BAM support. `--min-snp-frequency` and
`--min-mnv-frequency` use these same BAM-derived values. The filters are
independent: `--min-snp-frequency` applies to individual SNP observations, and
`--min-mnv-frequency` applies to phased MNV haplotypes. In mixed `SNP/MNV`
calls, a row or VCF record is kept when either component passes its own active
threshold.
Read-count and strand-support filters (`--snp`, `--mnv`, `--min-snp-strand`,
and `--min-mnv-strand`) follow the same independent SNP/MNV behavior.

Example:

```text
Chromosome	Gene	Positions	Base Changes	AA Changes	Variant Type	Change Type
MTB_anc	Rv0095c	104838	T	Asp126Glu	SNP	Non-synonymous
MTB_anc	Rv0095c	104941,104942	T,G	Gly92Gln	SNP/MNV	Non-synonymous
```

## VCF Output

Write VCF with:

```bash
--convert
```

or write both TSV and VCF with:

```bash
--both
```

Default file name:

```text
<input_name>.MNV.vcf
```

Use `--vcf-gz` for compressed output:

```text
<input_name>.MNV.vcf.gz
```

Common INFO fields:

| Field | Meaning |
|---|---|
| `GENE` | Gene or feature name |
| `AA` | Amino acid change |
| `CT` | Change type |
| `TYPE` | Variant type |
| `ODP` | Original depth from the input variant file |
| `OFREQ` | Original allele frequency from the input variant file |
| `SR`, `SRF`, `SRR` | SNP reads: total, forward, reverse |
| `MR`, `MRF`, `MRR` | MNV reads: total, forward, reverse |
| `DP` | Depth recalculated from BAM |
| `FREQ` | Frequency recalculated from BAM |
| `SBP` | SNP strand-bias p-value |
| `MSBP` | MNV strand-bias p-value |

The VCF header records the get_MNV version, command line, and thresholds used.
When `--emit-filtered` is enabled, VCF records below read-support, frequency,
strand-support, or strand-bias thresholds are written with FILTER tags such as
`LowSupport`, `LowFrequency`, `StrandSupport`, or `StrandBias`; otherwise they
are skipped.

## BCF Output

Write BCF with:

```bash
--bcf
```

BCF requires VCF output mode, so use it with `--convert` or `--both`.

Default file name:

```text
<input_name>.MNV.bcf
```

## JSON Files

### Summary JSON

Write with:

```bash
--summary-json run.summary.json
```

Includes:

- Input file checksums
- Per-contig variant counts
- Global variant counts
- Runtime timings
- Output paths

### Run Manifest

Write with:

```bash
--run-manifest run.manifest.json
```

Includes the summary plus:

- Command line
- Tool version
- Output file checksums
- Timestamp

### Error JSON

Write errors as JSON with:

```bash
--error-json run.error.json
```

This is useful in automated pipelines.

## Notes

- For MNV records, depth and frequency are calculated from reads spanning all
  positions in the grouped haplotype.
- Frequencies are printed with 4 decimal places.
- `--min-snp-frequency` and `--min-mnv-frequency` are values from `0` to `1`
  and require `--bam`.
- SNP and MNV frequency filters are independent, so a strong MNV haplotype is
  not removed by a stricter SNP-frequency threshold.
- SNP and MNV read-support and strand-support filters are also independent.
- `--sample all` writes one output set per VCF sample.
- `--keep-original-info` preserves non-get_MNV INFO fields from the input VCF.
