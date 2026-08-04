# Synthetic ground-truth validation

**English** · [Español](README.es.md)

The [differential suite](../differential/README.md) compares get_MNV with
`bcftools csq` and SnpEff. That is the stronger check, because those tools are
not ours. But it cannot validate the one thing get_MNV exists for: **codon-level
MNVs**. No per-variant annotator combines two SNVs in a codon the way get_MNV
does, so on exactly those calls the oracles disagree by design and cannot say
whether the residue is right.

This suite fills that hole, with a dataset whose answer is known because we
built it.

```bash
cargo build -p get_mnv --release
python3 tests/truth/run.py
```

## How the truth is derived

The expectation must not come from get_MNV's logic, or the test just agrees with
the code. So it takes a different route to the same answer:

```
mutate the genome  ->  re-extract the CDS  ->  translate  ->  compare proteins
```

get_MNV does codon-level surgery in place. This does whole-CDS re-translation
from a codon table written out in `run.py`, using plain string operations on the
reference. An implementation slip shows up as a disagreement.

## What it covers

Four genes, generated deterministically: plus and minus strand, single-exon and
spliced. For every codon of every one of them:

- every single substitution at every position (3 positions x 3 alternates)
- every pair of substitutions in the same codon
- the triple

Around 1560 cases, each run as its own VCF so no two expectations interfere.
Each is checked on gene, amino-acid change, change type, reference codon,
alternate codon, and the HGVS `c.` descriptor.

Plus 324 indels: deletions and insertions of 1, 2 and 3 bases at anchors spread
through every exon, including both edges where the initiator and terminator
codons live. Their expectation comes from the same re-translation, with the
indel applied in CDS coordinate space rather than by cutting the genome, since
the surrounding bases keep their coordinates. Checked on change type, event
class and both codons; the HGVS protein notation for indels is get_MNV's own
formatting, so reproducing it here would test spelling rather than biology.

The whole run takes about six seconds.

## How strong it is

Weaker than a third-party oracle, and worth being clear about why: this is a
second implementation by the same hand. It catches implementation mistakes, not
a shared misunderstanding of the biology. The consequence *labels*
(`Synonymous`, `Stop gained`, ...) follow get_MNV's documented convention, so on
those it validates the code against the documentation rather than against
nature. The amino acid, the codons and the CDS coordinate are the independent
parts.

It is also the reason the differential suite exists alongside it. Neither
replaces the other.

## Checking that it can still fail

A test that cannot fail is worth nothing. `GET_MNV_BIN` points the suite at
another build:

```bash
GET_MNV_BIN=/path/to/older/get_mnv python3 tests/truth/run.py
```

Against the build before the codon-orientation fix this reports 240
disagreements, all on the minus-strand genes. Against the build before the
start-codon fix it reports 11 more, all indels that remove the initiator.

Failures are grouped by category and sampled from each group rather than
printed in order, because one noisy category otherwise buries another: the run
above shows both the 240 ordering failures and the 11 start-codon ones.

## A note on writing the oracle

Its first run reported 665 disagreements. Every one of them was a bug in this
file, not in get_MNV: a leaked loop variable meant the generated VCF carried a
different alternate base from the one the expectation was computed for. The
oracle needs checking as much as the thing it checks.
