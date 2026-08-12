# Linkage: Telling a Haplotype From a Coincidence

Two variants sitting in one codon are only a multi-nucleotide variant if the
same molecules carry both. get_MNV reports how strongly they travel together, as
`D'` with a p-value, in the `Haplotype LD` and `Haplotype LD p` columns
(`LD` / `LDP` in VCF output). Requires `--bam`.

## Why a co-occurrence ratio is not enough

The obvious measure is the fraction of molecules carrying one variant that also
carry the other. It has a flaw that gets worse the more common the variants are.

Take two substitutions each present on 90% of molecules, with nothing linking
them. By arithmetic alone they are found together on 81% of molecules, and the
ratio reports `0.9000`. Read as linkage, that says "these two nearly always
travel together". What it actually says is "these two are both common".

`MNV Phasing Support` is that ratio and it is still reported, because "how often
do they co-occur" is a fair question with a useful answer. But it cannot answer
the question underneath it, which is whether the co-occurrence means anything.

## What D' measures

`D'` is the excess co-occurrence over what the two frequencies predict,
normalised by the largest excess those frequencies could have produced. That
normalisation is what makes it comparable between sites, and it puts the value
in `[-1, 1]`.

| `D'` | Reading |
|---|---|
| `+1` | The variants travel together as far as their frequencies allow. One haplotype. |
| `~0` | They co-occur exactly as often as chance predicts. Two variants that happen to share a codon, which is not the same thing as an MNV. |
| `-1` | They exclude each other: both present, never on one molecule. |

The same data as above, scored properly:

| Molecules with both | one | the other | neither | Ratio | `D'` | p |
|---|---|---|---|---|---|---|
| 81 | 9 | 9 | 1 | 0.9000 | **0.0000** | 1.0 |
| 90 | 0 | 0 | 10 | 1.0000 | **1.0000** | 5.8e-14 |
| 3 | 0 | 0 | 97 | 1.0000 | **1.0000** | 6.2e-6 |
| 0 | 20 | 20 | 0 | 0.0000 | **-1.0000** | 1.5e-11 |
| 10 | 10 | 10 | 10 | 0.5000 | **0.0000** | 1.0 |

Rows one and five are the point. A ratio of 0.9 and a ratio of 0.5 look like
strong and moderate linkage; both are exactly what independence predicts.

## The negative side

`D'` near `-1` is a finding, not a failure. Both alleles are present in the
sample, at appreciable frequency, and no molecule carries both.

In a haploid organism that is not one variant with two forms. It is **two
lineages**, and at codon level it is the read-level signature of a mixed
infection or of distinct sub-populations. Nothing in the output could express
that before: a co-occurrence ratio of `0.0000` says "they do not co-occur", which
reads the same whether they exclude each other or whether nobody looked.

## The p-value

`D'` of `1.0` from four molecules and `D'` of `1.0` from four hundred are the
same number and not the same evidence. `Haplotype LD p` is the two-tailed Fisher
exact test on the 2x2 table of molecule counts, so it says whether the departure
from independence is more than the depth could produce by chance.

Read them together. A large `|D'|` with a p-value near 1 is a small sample
speaking confidently about nothing.

## When it is absent

`-` means the question could not be answered, and that is deliberately distinct
from an answer of zero:

- No `--bam`, or a single-variant row: nothing to correlate.
- **No molecule observed the variants together.** Common when a codon straddles
  an intron and the fragments are shorter than the intron.
- **One of the alleles is on every observing molecule, or none of them.** With
  no variation at that locus there is no correlation to measure, and `D'`'s
  denominator is zero. Reporting `0.0` there would claim independence was
  measured and found.

That last case is more common than it sounds. If twelve molecules carry three
variants and eight carry only the first two, the first two are on all twenty
molecules that observed the window: the row is reported with its eight
molecules, and its linkage is `-`, because within this window those two never
vary.

## Indel haplotypes get it too

The statistic is not codon-specific. A local indel haplotype row carries it as
well, computed over the variants that row claims travel together.

That row previously had only a read count, which cannot distinguish a real
association from two common variants meeting. An insertion on 40% of molecules
and a substitution on 30%:

| Molecules with both | insertion only | substitution only | neither | Event Reads | `D'` | p |
|---|---|---|---|---|---|---|
| 5 | 15 | 10 | 20 | 5 | -0.1667 | 0.75 |
| 18 | 2 | 2 | 18 | 18 | 0.8000 | 5.3e-7 |

Both rows are real: five and eighteen molecules genuinely carry those two
variants together. The first is what two variants at those frequencies meet at
by chance; the second is a haplotype.

## Reading it beside the read count

They answer different questions and neither replaces the other.

- **`Event Reads` / `MNV Reads`**: how many molecules *are* this combination.
- **`Haplotype LD`**: whether its variants have anything to do with each other.

A haplotype carried by few molecules can be perfectly linked, and one carried by
many can be a coincidence. The first is a real minor species worth reporting;
the second is two variants that happen to be common in the same sample.

## Details of the calculation

- The table is built over **molecules**, not reads: the two mates of a
  paired-end fragment are one molecule. See
  [Indels and Local Haplotypes](indel-haplotypes.md) for what counts as one.
- Only molecules that observed **every** position of the row are included, on
  both sides of the table. A molecule that stopped between two positions saw no
  evidence about the pair.
- With three or more variants, every pair is scored and **the weakest decides**,
  because the row claims one molecule carries all of them, and that claim is
  only as good as its least-linked pair.
- `r²` is computed alongside and is available in the API. `D'` is reported
  because the question here is "are these on the same molecules", which is what
  `D'` answers; `r²` additionally reflects how far apart the two frequencies
  are, which is a different thing to know.

## How this is checked

[`tests/phasing/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/tests/phasing)
builds alignments molecule by molecule, choosing which reads carry which
variants, and compares against `D'` and an exact test derived a second time in
the suite itself. It sweeps 2x2 tables covering independence, both directions of
association, and the degenerate tables where `D'` has no value, alongside the
linkage-support sweeps.

## Limits

- The statistic is **within a window**: a codon, or a local haplotype window.
  get_MNV does not infer long-range phase, and a variant outside the window is
  not part of the table.
- It describes the molecules that were sequenced. Read-level linkage in an
  amplicon scheme is bounded by the amplicon, and a pair separated by more than
  a fragment simply has no molecules in common to count.
- It says nothing about which of two competing lineages is which. `D'` near
  `-1` reports that two exist; identifying them is a job for a phasing or
  strain-deconvolution tool.
