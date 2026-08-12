# Getting Started

This short tutorial runs get_MNV end to end on the bundled *M. tuberculosis*
example and shows how to inspect a real codon-level MNV, including its read
support.

## 1. Get the example data

The [`example/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/example)
folder of the repository ships everything you need:

| File | Role |
|---|---|
| `MTB_ancestor.fas` | Reference (single contig `MTB_anc`) |
| `anot_genes.txt` | Gene table (`name`, `start`, `end`, `strand`) |
| `G35894.var.snp.vcf` | VarScan-style variant calls |
| `G35894.demo.bam` | Tiny aligned-read subset for the read viewer |

```bash
git clone https://github.com/PathoGenOmics-Lab/get_MNV.git
cd get_MNV/example
```

## 2. Annotate the variants

```bash
get_mnv \
  --vcf G35894.var.snp.vcf \
  --fasta MTB_ancestor.fas \
  --genes anot_genes.txt
```

This writes `G35894.var.snp.MNV.tsv` next to the input. For this sample it
produces 941 variants: 797 SNP, 10 SNP/MNV and 134 intergenic.

!!! note "Genetic code"
    get_MNV defaults to NCBI translation table 11 (bacterial), which is correct
    for *M. tuberculosis*. Use `--translation-table` for other organisms.

## 3. Read the output

Open the TSV. One row stands out: gene `Rv2036` carries two SNVs in the same
codon:

| Positions | Reference Bases | Base Changes | Reference Codon | MNV Codon | AA Changes | Variant Type |
|---|---|---|---|---|---|---|
| `2282376, 2282377` | `T, T` | `C, C` | `GTT` | `GCC` | `Val93Ala` | `SNP/MNV` |

Annotated one SNV at a time, each change reads differently; combined, the codon
`GTT → GCC` is a single Val→Ala substitution. That codon-aware reclassification
is exactly what get_MNV is for. See [Output Formats](output-formats.md) for
every column.

## 4. Add read support (BAM)

Provide an aligned, indexed BAM to count how many reads actually carry the
combined change:

```bash
get_mnv \
  --vcf G35894.var.snp.vcf \
  --fasta MTB_ancestor.fas \
  --genes anot_genes.txt \
  --bam G35894.demo.bam --mnv 1
```

The `Rv2036` row now reports `MNV Reads = 24` (12 forward / 12 reverse) and
`MNV Frequencies = 1.0000`, meaning every read spanning the codon carries both ALT
bases.

!!! tip "BAM requirements"
    The BAM must be coordinate-sorted and indexed (`.bai`), and aligned to the
    same reference. The bundled `G35894.demo.bam` covers only that one locus.

## 5. See the reads in the GUI

The desktop app draws an IGV-style pileup for any row. Load the VCF, FASTA, gene
table and `G35894.demo.bam`, run the analysis, then open the `Rv2036` row to see
all 24 reads with the two ALT bases highlighted. See the
[Desktop GUI](gui.md) guide.

## Next steps

- [Common Recipes](usage.md): ready-to-run commands for the usual jobs.
- [CLI Reference](cli-reference.md): every option, with its default.
- [Input Formats](input-formats.md) and [Output Formats](output-formats.md).
- [Scope and Compatibility](indel-mnv-semantics.md): boundary rules and tuning.
- [Linkage](linkage.md): telling a real haplotype from a coincidence.
