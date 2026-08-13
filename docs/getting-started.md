# Command Line Tutorial

The command line is get_MNV's main interface: everything the tool can do is
reachable from it, and it is what you put in a pipeline. This tutorial runs it
end to end on the bundled *M. tuberculosis* example, and spends most of its time
on the part that matters, reading what comes back.

There is a [Desktop GUI Tutorial](gui-tutorial.md) covering the same ground in
the app.

## 1. Get the example data

The [`example/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/example)
folder of the repository ships everything you need:

| File | Role |
|---|---|
| `MTB_ancestor.fas` | Reference (single contig `MTB_anc`) |
| `anot_genes.txt` | Gene table (`name`, `start`, `end`, `strand`) |
| `G35894.var.snp.vcf` | VarScan-style variant calls |
| `G35894.demo.bam` | Tiny aligned-read subset for read support |

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

get_MNV logs what it decided before it logs what it found. The run prints its
whole configuration, then the lines that matter:

```text
[INFO get_mnv::pipeline] Using genetic code: Bacterial/Archaeal/Plant Plastid (NCBI table 11)
[INFO get_mnv::io::annotation::tsv] TSV gene entries parsed: 4008 | mapped to SNPs: 635 | without SNPs: 3373
[INFO get_mnv::pipeline::processing] Contig 'MTB_anc' -> 950 SNP/variant records in VCF, 635 mapped genes
[INFO get_mnv::pipeline::processing] Contig 'MTB_anc' -> 134 intergenic variant(s) added
[INFO get_mnv::pipeline] Summary 'MTB_anc' -> variants=941 (SNP=797, MNV=0, SNP/MNV=10, INDEL=0, intergenic=134)
[INFO get_mnv::pipeline] Processing complete. Output files generated successfully.
```

It writes `G35894.var.snp.MNV.tsv` in the current directory.

!!! tip "Read the summary line, every time"
    `950` records became `941` variants because a codon carrying two
    substitutions is reported once. `SNP/MNV=10` is the codon-level work: ten
    codons where separate VCF records combine. `MNV=0` means no single record
    already carried more than one substituted base.

    `intergenic=134` is your sanity check on the annotation. get_MNV cannot tell
    a variant outside a gene from a variant whose gene it failed to find, so a
    mismatched annotation does not fail: it succeeds, reports almost everything
    as intergenic, and prints `Processing complete`. If that count is most of
    your variants, your gene coordinates do not line up with your reference.

!!! note "Genetic code"
    get_MNV defaults to NCBI translation table 11 (bacterial), which is correct
    for *M. tuberculosis*. Use `--translation-table` for other organisms.

Use `--dry-run` to get the same summary with nothing written, which is the cheap
way to check inputs before a long run.

## 3. Read the output

Open the TSV. One row stands out: gene `Rv2036` carries two SNVs in the same
codon.

| Column | Value |
|---|---|
| `Positions` | `2282376, 2282377` |
| `Reference Bases` / `Base Changes` | `T, T` → `C, C` |
| `Reference Codon` | `GTT` |
| `SNP Codon` | `GCT, GTC` |
| `MNV Codon` | `GCC` |
| `AA Changes` | `Val93Ala` |
| `SNP AA Changes` | `Val93Ala, Val93Val` |
| `Variant Type` | `SNP/MNV` |
| `SO Term` / `Impact` | `missense_variant` / `MODERATE` |
| `HGVS g.` | `g.[2282376T>C;2282377T>C]` |

Read one base at a time, the two substitutions disagree: `GCT` is Val93Ala and
`GTC` is Val93Val, a silent change. Read together as `GCC`, the codon is a
single Val→Ala substitution. That reclassification is what get_MNV is for, and
it is invisible to a per-SNV annotator. See
[Output Formats](output-formats.md) for every column.

## 4. Add read support

Annotation says what the codon *would* mean. It does not say whether any single
molecule carries both changes. For that, give get_MNV the alignments:

```bash
get_mnv \
  --vcf G35894.var.snp.vcf \
  --fasta MTB_ancestor.fas \
  --genes anot_genes.txt \
  --bam G35894.demo.bam
```

The `Rv2036` row gains its evidence columns:

| Column | Value | Meaning |
|---|---|---|
| `MNV Reads` | `24` | reads carrying the whole haplotype |
| `MNV Forward/Reverse Reads` | `12` / `12` | balanced across strands |
| `Total Reads` | `24` | depth at the codon |
| `MNV Frequencies` | `1.0000` | every read spanning it carries both |
| `SNP Reads` | `0, 0` | reads carrying one change *without* the other |
| `MNV Phasing Support` | `1.0000` | of the reads that could answer, all carried both |
| `MNV Phasing Reads` | `24` | how many reads answered |
| `Haplotype LD` | `-` | see below |

`SNP Reads = 0, 0` next to `MNV Reads = 24` is not a gap: the two columns split
the evidence instead of double-counting it. No read carries one substitution
alone, so nothing lands in the solo counts.

!!! note "Why `Haplotype LD` is empty on a perfect haplotype"
    `D'` measures whether two variants co-occur more than their own frequencies
    predict. Here both are on all 24 molecules, so neither varies and there is
    nothing left to correlate: the honest answer is `-`, not `1`. That is a
    different question from `MNV Phasing Support`, which asks how many of the
    reads that could answer carried the whole haplotype, and answers `1.0000`.
    [Linkage](linkage.md) works through the cases where they diverge.

!!! tip "BAM requirements"
    The BAM must be coordinate-sorted and indexed (`.bai`), and aligned to the
    same reference. The bundled `G35894.demo.bam` covers only the `Rv2036`
    locus, which is why every other row reports no read support.

## 5. Filter on that support

With a BAM you can require evidence. The thing to know before you do is how the
two thresholds combine on a codon-level row:

| Command | Rows | What happened |
|---|---:|---|
| *(no filter)* | 941 | everything |
| `--snp 1` | 10 | every plain SNP row without a supporting read is gone; the 10 codon-level rows stay |
| `--mnv 5` | 941 | **nothing happens** |
| `--snp 1 --mnv 5` | 1 | only `Rv2036`, the one codon with reads |

A `SNP/MNV` row is kept when **either** its individual SNVs clear `--snp` **or**
its haplotype clears `--mnv`. That is deliberate, so a well-supported haplotype
survives even when its individual SNVs do not. The consequence is easy to trip
over: with `--snp` left at its default of `0`, the SNP side passes trivially,
so `--mnv` on its own never removes a row. Raise both, or neither.

## 6. Get a report you can send

The same run writes a single self-contained HTML file:

```bash
get_mnv \
  --vcf G35894.var.snp.vcf \
  --fasta MTB_ancestor.fas \
  --genes anot_genes.txt \
  --bam G35894.demo.bam \
  --report sample.html
```

<div style="text-align: center;" markdown>
![The interactive HTML report get_MNV writes: summary counts, consequence breakdown, search and a variant matrix](assets/cli-01-report.png){ width="840" }
</div>

*The report opens in any browser with no server and travels as one attachment.
It records the exact command that produced it along the top.*

It is searchable and filterable, breaks the calls down by Sequence Ontology term
(here 476 missense, 320 synonymous, 134 intergenic, 10 stop gained, 1 stop
lost), and draws a variant matrix across the contig. For a cohort already
processed one sample per run, build one report from the existing outputs with
`--report-from run1.MNV.tsv run2.MNV.tsv --report cohort.html`.

## 7. The other things it can write

| Flag | Output |
|---|---|
| `--convert` / `--both` | VCF instead of, or alongside, the TSV |
| `--vcf-gz` / `--index-vcf-gz` | BGZF-compressed VCF and its Tabix index |
| `--bcf` | BCF converted from the generated VCF |
| `--summary-json <FILE>` | the summary above, machine-readable |
| `--run-manifest <FILE>` | command, version, inputs, outputs and checksums |
| `--error-json <FILE>` | structured error details when a run fails |

## Next steps

- [Common Recipes](usage.md): ready-to-run commands for the usual jobs.
- [CLI Reference](cli-reference.md): every option, with its default.
- [Input Formats](input-formats.md) and [Output Formats](output-formats.md).
- [Scope and Compatibility](indel-mnv-semantics.md): boundary rules and tuning.
- [Linkage](linkage.md): telling a real haplotype from a coincidence.
- [Desktop GUI Tutorial](gui-tutorial.md): the same run in the app.
