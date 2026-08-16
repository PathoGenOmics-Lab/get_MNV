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

`--split-multiallelic` emits one annotated row per ALT, even when multiple
alts share the same codon position. Each row reports its own amino-acid
effect, codon, and BAM-derived read support. Pre-splitting the VCF with
`bcftools norm -m -` produces the same output.

## BCF Input

Example:

```text
BCF input is not supported. Convert to VCF first
```

Fix:

```bash
bcftools view input.bcf > input.vcf
get_mnv --vcf input.vcf ...
```

`--bcf` is an output option only; it does not make BCF valid as input.

## TSV Annotation with Multiple Contigs

Example:

```text
TSV annotation does not include contig names
```

Fix:

- Use `--gff` for multi-contig data.
- Or restrict the run to one contig with `--chrom`.

## Fewer Rows Than Records in the VCF

Example:

```text
Skipped 148 ALT alleles the selected sample's genotype does not carry
```

Cause: a VCF record lists every ALT seen at that site across every sample, so the
selected sample's `GT` decides which of them belong to it. A genotype of `0/0`
carries none, and on a multiallelic record `1/1` carries only the first ALT.

Fix:

- Nothing, if the genotypes are right: the run annotated what the sample has.
- Pick the intended sample with `--sample`, or write one file per sample with
  `--sample all`.
- If the calls have no meaningful genotype, remove the `GT` field or set it to
  `./.`, which keeps every allele because unknown is not absence.

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

## Warnings

A warning never stops the run, so it is easy to scroll past. These are the ones
that change what ends up in the output, and none of them is visible in the TSV
itself: if the annotation of a variant is not what you expected, or its metrics
came from a sample you did not pick, the reason was on stderr.

None of them removes a variant. A base get_MNV cannot place in a codon is still
reported, as a row naming its gene with `-` in every amino acid column, so the
symptom to look for is an unannotated row rather than a missing one.

### The amino acid call is missing for a variant that is inside a gene

| Warning | What it means | What to do |
|---|---|---|
| `falls in the phase-skipped region of CDS ...` | The codon spans into a neighbouring exon, and a single GFF row cannot say what the neighbouring bases are, so no codon is built for that base. The variant is still reported: one row naming the gene, `-` in every amino acid column, and SO term `coding_sequence_variant` (`MODIFIER`). | Give the CDS rows a `transcript_id` or `Parent` so get_MNV can splice the transcript, and run with `--gff-features CDS`. |

### The annotation is coarser than the file could support

The first three fall back to per-feature annotation, so codons are grouped inside
one CDS row instead of across the spliced transcript, and a variant near an exon
boundary can then get a different amino acid call. The fourth changes nothing on
its own.

| Warning | What it means | What to do |
|---|---|---|
| `has rows on both strands` | Two rows of one transcript disagree about orientation, and a spliced CDS cannot be built from them. | Fix the strand column, or split the rows into separate transcripts. |
| `has a single CDS row with non-zero phase` / `starts with non-zero phase` | The selected rows do not look like a whole coding sequence, so the phase cannot be trusted to place the codons. | Include every CDS row of the transcript in the run. |
| `contains CDS rows with non-zero phase, but --gff-features does not include 'CDS'` | Codon grouping is ignoring phase the file provides. | Add `--gff-features CDS`. |
| `lists N duplicate CDS row(s)` | get_MNV drops the repeats itself and still builds the spliced transcript, so the annotation is the same one a clean file would give. | Nothing. Deduplicating the GFF only silences the warning. |

### The metrics come from a sample you did not choose

| Warning | What it means | What to do |
|---|---|---|
| `Multi-sample VCF detected (N samples). Using first sample` | `ODP` and `OFREQ` are read from whichever sample comes first in the header, which is rarely a deliberate choice. That sample's genotype also decides which ALT alleles are kept, so the choice changes the row count too (see [Fewer Rows Than Records in the VCF](#fewer-rows-than-records-in-the-vcf)). | Name one with `--sample`, or annotate every one with `--sample all`. |

### The read support is not what the input claimed

| Warning | What it means | What to do |
|---|---|---|
| `has N spanning read(s) but 0 with exact CIGAR support` | Reads cover the site but none of them carries that exact indel, which is what happens when the allele is written at a different position inside a homopolymer or repeat. | Left-align the input first, for example `bcftools norm -f ref.fa`. |
| `was declared X by the input VCF, but N of M reads spanning the site carry the whole haplotype` | The caller and the reads disagree about whether these changes travel together. | Look at the alignment before trusting either. |
| `shows N distinct combinations on the reads; reporting the K best supported` | A window with more local haplotypes than get_MNV reports. The dropped ones are the weakly supported tail. See [Indels and Local Haplotypes](indel-haplotypes.md). | Nothing, unless a rare haplotype matters; `--chrom` takes whole contig names, not regions, so inspect the window in the alignment instead. |

### A flag did nothing

| Warning | What it means |
|---|---|
| `--keep-original-info has no effect with iVar TSV input` | An iVar TSV has no INFO column to preserve. |
| `--gff-features is ignored when using TSV annotation format (--genes)` | Feature types are a GFF/GTF idea; the four-column TSV has none. |

### A threshold reached less than you meant

| Warning | What it means |
|---|---|
| `The MNV thresholds will not remove any SNP/MNV row` | The MNV thresholds are still filtering MNV, indel and intergenic rows. What they cannot reach is the SNP/MNV rows: one survives when *either* its SNVs or its haplotype clears its bar, so with the SNP thresholds at `0` the SNP side always passes. Raise `--snp`, `--min-snp-frequency` or `--min-snp-strand` as well to filter codon-level rows. |


### An optional external program was missing

| Warning | What it means |
|---|---|
| `bcftools not found in PATH. Skipping BCF output.` | `--bcf` was asked for and no BCF was written. The run continues and exits `0`; the TSV and VCF are unaffected, and the summary JSON reports no BCF. Install bcftools, or drop `--bcf`. |
| `tabix not found in PATH. Skipping .tbi index` | `--index-vcf-gz` was asked for and no `.tbi` was written. The `.vcf.gz` itself is complete. Install samtools/htslib, or index it later with `tabix -p vcf <file>`. |

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
| `14` | Integer parsing error |
| `15` | Floating-point parsing error |
| `16` | JSON error |


## Still stuck

Everything above is a message the tool can produce and what it is asking for. If
what you are seeing is not here, or the answer here did not fit your data, ask in
[Q&A](https://github.com/PathoGenOmics-Lab/get_MNV/discussions/categories/q-a). Answers there get marked as answers, so the next person with your
problem finds it before they have to ask.
