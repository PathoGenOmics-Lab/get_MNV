# Output formats

## TSV output (default)

File: `<vcf_filename>.MNV.tsv`

| Column | Description |
|--------|-------------|
| Chromosome | Contig/chromosome name |
| Gene | Gene name (or `intergenic`) |
| Positions | Variant positions (comma-separated for MNVs) |
| Reference Bases | Reference nucleotides |
| Base Changes | Alternative nucleotides |
| AA Changes | Amino acid changes (MNV haplotype), **full-protein numbering** (compatible with Ensembl VEP, ANNOVAR, SnpEff, UniProt) |
| SNP AA Changes | Amino acid changes per individual SNP, full-protein numbering |
| Local AA Changes | Same MNV haplotype change but using **per-feature (exon-local) numbering** — what versions ≤ 1.1.1 used to report. Identical to `AA Changes` for prokaryotes and single-exon features |
| Local SNP AA Changes | Per-SNP amino acid change in per-feature (exon-local) numbering |
| Variant Type | `SNP`, `MNV`, `SNP/MNV`, or `INDEL` |
| Change Type | `Synonymous`, `Non-synonymous`, `Stop gained`, `Stop lost`, `Unknown`, `Indel overlap`, frameshift variants, `In-frame Indel`, `Frameshift Indel` |
| Reference Codon | Original codon sequence |
| SNP Codon | Codon with individual SNP substitutions |
| MNV Codon | Codon with all MNV substitutions |

**With BAM (additional columns):**

| Column | Description |
|--------|-------------|
| SNP Reads | Supporting reads per SNP position |
| SNP Forward/Reverse Reads | Strand-specific SNP counts |
| MNV Reads | Reads supporting the full MNV haplotype |
| MNV Forward/Reverse Reads | Strand-specific MNV counts |
| Total Reads | Total depth at each position |
| SNP Frequencies | Per-position SNP frequencies |
| MNV Frequencies | MNV haplotype frequency |

### Example

```
Chromosome	Gene	Positions	Base Changes	AA Changes	Variant Type	Change Type
MTB_anc	Rv0095c	104838	T	Asp126Glu	SNP	Non-synonymous
MTB_anc	Rv0095c	104941,104942	T,G	Gly92Gln	SNP/MNV	Non-synonymous
MTB_anc	esxL	1341102,1341103	T,C	Arg33Ser	MNV	Non-synonymous
```

## VCF output (`--convert` or `--both`)

File: `<vcf_filename>.MNV.vcf` (or `.MNV.vcf.gz` with `--vcf-gz`)

### INFO fields

| Field | Description |
|-------|-------------|
| `GENE` | Gene name |
| `AA` | Amino acid change |
| `CT` | Change type |
| `TYPE` | Variant type (SNP/MNV/INDEL) |
| `ODP` | Original depth from input VCF |
| `OFREQ` | Original allele frequency from input VCF |
| `SR/SRF/SRR` | SNP reads (total/forward/reverse) — BAM only |
| `MR/MRF/MRR` | MNV reads (total/forward/reverse) — BAM only |
| `DP` | Recalculated depth from BAM |
| `FREQ` | Recalculated frequency from BAM |
| `SBP` | SNP strand-bias p-value (with `--strand-bias-info`) |
| `MSBP` | MNV strand-bias p-value (with `--strand-bias-info`) |

### FILTER tags (with `--emit-filtered`)

| Tag | Description |
|-----|-------------|
| `LowSupport` | Below minimum read support thresholds |
| `StrandSupport` | Below minimum per-strand support |
| `StrandBias` | Below Fisher exact p-value threshold |

### Header metadata

The VCF header includes `##get_mnv_version`, `##get_mnv_command`, and all parameter thresholds used.

## BCF output (`--bcf`)

Binary VCF converted from the generated VCF. File: `<vcf_filename>.MNV.bcf`

## JSON outputs

### Summary JSON (`--summary-json`)

Structured run summary with:
- Per-contig variant counts (SNP, MNV, SNP/MNV, INDEL, intergenic)
- Global aggregated counts
- Per-phase timings (ms): `parse_inputs_ms`, `process_ms`, `emit_ms`, `total_ms`
- SHA-256 checksums for all input files
- Schema version for stable parsing

### Run manifest (`--run-manifest`)

Reproducibility manifest including:
- Everything in summary JSON
- Output file checksums (VCF/TSV/BCF)
- Tool version and command line
- Unix timestamp

### Error JSON (`--error-json`)

Structured error details when the command fails, including error code and message.

## Notes

- For MNV records, `DP/FREQ` use the depth of reads spanning all SNP positions in the haplotype
- Frequencies are printed with 4 decimal places
- `--sample all` writes one output set per sample with suffix `.sample_<name>`
- `--keep-original-info` preserves all non-get_mnv INFO fields from the input VCF
