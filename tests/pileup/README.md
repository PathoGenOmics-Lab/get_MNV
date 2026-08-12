# Read counting against bcftools mpileup

**English** · [Español](README.es.md)

get_MNV has three suites that check its annotation and one that checks its
phasing, and until this one every check on the *read counting* was written by
the same hand as the code: the scenario expectations, and the arithmetic in the
phasing suite. Both catch implementation error. Neither can catch a shared
misunderstanding of what a read is worth, because both start from the same
understanding.

This suite hands the same BAM to `bcftools mpileup`, which nobody here wrote,
and compares four numbers per site:

| get_MNV | mpileup |
|---|---|
| `Total Reads` | `FORMAT/DP` |
| `SNP Reads` | `FORMAT/AD` for that allele |
| `SNP Forward Reads` | `FORMAT/ADF` |
| `SNP Reverse Reads` | `FORMAT/ADR` |

```bash
cargo build -p get_mnv --release
python3 tests/pileup/run.py
```

## Why those four

They are where this layer's defects have actually been. The strand-attribution
bug fixed on this branch, where a mate hundreds of bases from a variant lent its
strand to it, is a flat disagreement with `ADF`/`ADR`, and it would have been
caught the day it was written instead of surviving until someone thought to look.

## Keeping the two tools comparable

A difference that comes from asking the tools different questions teaches
nothing, so the thresholds are given to both explicitly, and BAQ is disabled
with `-B` because get_MNV does not recalibrate base qualities. Every comparison
site sits in the filler between genes, so each is a single-position row and no
codon grouping stands between the two counts.

The cases cover both strands, an allele that appears on one strand only, mate
pairs whose two mates read different sites, overlapping mates, bases below the
quality floor, reads below the mapping-quality floor, and several sites counted
from one cache.

## Accepted differences

As in the [differential suite](../differential/README.md), the rule is not that
the tools must match. A difference that is understood is written down in the
case that produces it, with the reason; a difference nobody has explained fails
the run.

There is one, on overlapping mates. Both tools count 16 molecules and agree on
depth and support; they split the strands differently. bcftools removes the
overlap by keeping one mate's base and zeroing the other, so each fragment lands
in exactly one arm and `ADF + ADR == AD`. Which arm is effectively arbitrary:
all sixteen fragments are identical and it divides them 5 forward and 11
reverse, a split carrying no information about strand. get_MNV credits both
arms, which is what is actually true of the molecules, and means none of them is
one-strand evidence. The convention differs; neither count is wrong.

With `--count-mates-separately` the disagreement is wholesale, because bcftools
has no equivalent mode. That case is here to pin that the *default* mode is the
one that agrees with an outside tool.

## Checking that it can still fail

Against the build from before the strand fix (`18e426c`):

```bash
GET_MNV=/path/to/older/get_mnv python3 tests/pileup/run.py
```

it reports the bug directly:

```
03_mates_far_apart: reverse at 320 is 15 for get_MNV and 0 for mpileup
03_mates_far_apart: forward at 750 is 15 for get_MNV and 0 for mpileup
```

Fifteen reverse molecules at a position no reverse read ever reached.

## What it does not cover

Indel and complex-event counting. mpileup's indel representation is its own
(`IDV`, `IMF`, and candidate ALTs it chooses itself), so lining it up against
get_MNV's exact-CIGAR event counts would compare two different definitions
rather than one implementation against another. Substitution depth and strand
are the part where the two tools genuinely answer the same question.
