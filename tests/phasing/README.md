# Phasing validation against a known mixture

**English** · [Español](README.es.md)

The [synthetic truth suite](../truth/README.md) validates the codon arithmetic,
but it has no BAM, so it says nothing about the other half of an MNV call:
whether the substitutions sit on the *same molecules*. Two variants in one codon
on two different molecules are not a codon-level MNV at all, and no amount of
codon arithmetic can tell them apart.

This suite supplies that half. It builds alignments molecule by molecule,
choosing which reads carry which variants, and asks whether get_MNV recovers the
mixture it was handed.

```bash
cargo build -p get_mnv --release
python3 tests/phasing/run.py
```

## How the truth is derived

Nothing here is a number copied from a previous run. A case is a list of
molecule classes: how many molecules carry which variants, and whether a read
off that molecule reaches the whole codon. The expectation follows from those
counts alone:

```
denominator = fewest codon-spanning molecules carrying any one constituent SNV
numerator   = codon-spanning molecules carrying all of them
```

If the code and the arithmetic disagree, one of them is wrong.

## What it covers

Five codon geometries, each swept through all 21 mixtures from fully trans to
fully cis:

- a plus-strand codon
- a minus-strand codon, which must give identical answers, because linkage is a
  property of the molecules and not of the strand the gene is read from
- a codon split by an intron, spanned by reads carrying a CIGAR `N`
- all three positions of one codon, which exercises the least-supported rule
  with more than two constituents
- the plus-strand codon again, sequenced as mate pairs whose two mates each
  reach one end of it and not the other, so the answer exists only at the level
  of the molecule

Plus the cases where the honest answer is not a number:

- reads that carry one variant and stop before the other, which must neither
  support nor dilute the ratio
- no read spanning the codon at all, which is unknown (`-`), not zero
- a minority haplotype nested inside a major variant, which is full support on
  few reads, and the reason the read count is reported next to the ratio

120 mixtures, a few seconds.

## Checking that it can still fail

A test that cannot fail is worth nothing. Against the binary from before the
phasing fix it reports 171 disagreements:

```bash
GET_MNV=/path/to/older/get_mnv python3 tests/phasing/run.py
```

Those failures are worth reading, because they are all three ways the old ratio
went wrong at once. It divided by the reads carrying one SNV *alone*, excluding
the reads carrying the whole haplotype, so it overstated linkage everywhere and
pinned to `1.0000` from an even mixture upwards: a coin-flip between cis and
trans was reported as a perfect haplotype. When every read carried the full
haplotype there were no solo reads left to divide by, so flawless linkage came
out as `-`. And when no read spanned the codon it returned `0.0000`, which reads
as proof of trans, for a question nothing had answered.

## What it is worth

Like the truth suite, this is a second implementation by the same hand, so it
catches implementation error rather than a shared misunderstanding. What makes
it stronger than a hand-written expectation is that the mixture is the input:
the answer is fixed before get_MNV runs, and a sweep leaves nowhere for an
off-by-one or a saturating ratio to hide.

What it does not cover: real alignment artefacts. Every read here is a clean
record at MAPQ 60 with Q40 bases, laid out to a chosen geometry. Soft clips,
overlapping mate pairs, and reference bias at the ends of amplicons are outside its reach.
