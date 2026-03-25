# Troubleshooting

## Common errors

### Contig mismatches

```
[E002] Contig validation failed ... Missing in VCF/FASTA
```

Contig names don't match between files. All inputs (VCF, FASTA, GFF) must use identical contig names (case-sensitive). Normalize names before running.

### Reference mismatch

```
[E002] VCF REF/FASTA mismatch at <contig>:<pos>
```

The REF allele in VCF doesn't match the FASTA base at that position. Ensure the VCF was called against the same reference.

### Invalid bases

```
[E002] Invalid base 'X' in REF/ALT allele ...
```

Non-IUPAC base detected. Clean invalid alleles from the VCF before running.

### Sample not found

```
Sample '<name>' not found in VCF header
```

Check sample name spelling, or omit `--sample` to use the first sample.

### Multi-contig with TSV annotation

```
TSV annotation does not include contig names
```

TSV format has no contig column. Use `--gff` for multi-contig VCFs, or restrict with `--chrom`.

### Strict mode failures

```
--strict enabled, but original VCF metrics are missing ...
```

Ensure each VCF record has recoverable `DP` and `AF/FREQ` in FORMAT or INFO fields.

### Multiallelic records

```
Multiallelic VCF record ... is not supported
```

Either pre-split with `bcftools norm -m -` or run with `--split-multiallelic`.

### Flag conflicts

| Error | Fix |
|-------|-----|
| `--index-vcf-gz requires --vcf-gz` | Enable BGZF output first |
| `--bcf requires --convert or --both` | Enable VCF output mode |
| `--keep-original-info requires --convert or --both` | Enable VCF output mode |
| `--min-strand-bias-p must be between 0 and 1` | Use a value in `[0, 1]` |

### No samples

```
Requested --sample all but input VCF has no sample columns
```

Use a VCF with FORMAT/sample columns, or omit `--sample all`.

## Exit codes

| Code | Name | Description |
|------|------|-------------|
| 0 | Success | |
| 1 | Generic | Unclassified error (`E000`) |
| 2 | Config | CLI/configuration error (`E001`) |
| 3 | Validation | Input/validation error (`E002`) |
| 10 | I/O | File read/write error |
| 11 | CSV | CSV parsing error |
| 12 | HTSlib | BAM/VCF library error |
| 13 | UTF-8 | Encoding error |
| 14 | Parse | Data parsing error |
