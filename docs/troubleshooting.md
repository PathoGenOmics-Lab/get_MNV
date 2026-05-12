# Troubleshooting

This page lists common errors and the quickest fixes.

## Contig Names Do Not Match

Example:

```text
[E002] Contig validation failed
```

Fix:

- Make sure the same contig names are used in the variant file, FASTA, GFF/GTF,
  and BAM.
- Names are case-sensitive.
- `chr1` and `1` are different names.

## VCF REF Does Not Match FASTA

Example:

```text
[E002] VCF REF/FASTA mismatch at <contig>:<pos>
```

Fix:

- Check that the VCF was called against the same reference FASTA.
- Check that coordinates and contig names were not changed after variant
  calling.

## iVar TSV Is Not Detected

Fix:

- Pass the file with `--tsv sample_variants.tsv`.
- Check that the TSV header contains at least `REGION`, `POS`, `REF`, and
  `ALT`.
- If the file is actually a gene annotation TSV, pass it with `--genes`, not
  `--tsv`.

## Invalid Bases

Example:

```text
[E002] Invalid base 'X' in REF/ALT allele
```

Fix:

- Remove invalid alleles from the variant file.
- get_MNV accepts IUPAC DNA bases, but not arbitrary symbols in REF/ALT.

## Multiallelic VCF Records

Example:

```text
Multiallelic VCF record is not supported
```

Fix:

```bash
get_mnv ... --split-multiallelic
```

or pre-split the VCF:

```bash
bcftools norm -m - input.vcf > split.vcf
```

## TSV Annotation with Multiple Contigs

Example:

```text
TSV annotation does not include contig names
```

Fix:

- Use `--gff` for multi-contig data.
- Or restrict the run to one contig with `--chrom`.

## Sample Name Not Found

Example:

```text
Sample '<name>' not found in VCF header
```

Fix:

- Check the sample name in the VCF header.
- Omit `--sample` to use the first sample.
- Use `--sample all` only when the VCF has sample columns.

## Strict Mode Fails

Example:

```text
--strict enabled, but original VCF metrics are missing
```

Fix:

- Disable `--strict`, or
- Make sure each variant has depth and frequency metrics that get_MNV can read.

## Filtering by Allele Frequency

Use `--min-snp-frequency <F>` for SNP records and `--min-mnv-frequency <F>`
for MNV haplotypes. Values are fractions from `0` to `1`, so `0.05` means 5%.

These filters require `--bam` because get_MNV calculates them from read
support. They do not use the original VCF/iVar `OFREQ` value.
The SNP and MNV thresholds are independent. In mixed `SNP/MNV` calls,
`--min-snp-frequency` filters SNP observations and `--min-mnv-frequency`
filters the phased MNV haplotype; a strong MNV should not disappear only
because individual SNP observations are below the SNP threshold.

Common fixes:

- If you see a `requires --bam` error, add a sorted/indexed BAM or remove the
  frequency filters.
- If you want to filter by the caller's original allele frequency (`OFREQ`),
  pre-filter the VCF or iVar TSV before running get_MNV.
- Combine frequency filters with read-support filters such as `--snp`, `--mnv`,
  `--min-snp-strand`, and `--min-mnv-strand` for stricter calls.

## Output Directory Is Not Writable

Example:

```text
Read-only file system
```

Fix:

- Run the command from a writable folder, or
- In the GUI, choose an output directory where you have write permission.

## Flag Conflicts

| Error | Fix |
|---|---|
| `--index-vcf-gz requires --vcf-gz` | Add `--vcf-gz`. |
| `--bcf requires --convert or --both` | Add `--convert` or `--both`. |
| `--keep-original-info requires --convert or --both` | Add `--convert` or `--both`. |
| `--min-strand-bias-p must be between 0 and 1` | Use a value from `0` to `1`. |

## Exit Codes

| Code | Meaning |
|---:|---|
| `0` | Success |
| `1` | Generic error |
| `2` | Configuration error |
| `3` | Input validation error |
| `10` | File read/write error |
| `11` | CSV/TSV parsing error |
| `12` | BAM/VCF parsing error |
| `13` | UTF-8 encoding error |
| `14` | Data parsing error |
