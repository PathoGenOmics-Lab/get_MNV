# Differential annotation suite

get_MNV's own tests can only check that it produces what we wrote down as
expected. When the expectation is wrong the test goes green and the tool stays
wrong.

That is not hypothetical. A ribosomal-slippage CDS join, the shape SARS-CoV-2
ORF1ab is annotated with, was being given `splice_donor_variant` (`HIGH`) on
plain coding bases. The unit test that was supposed to cover splicing built a
"two-exon" gene whose exons abutted, with no intron between them. It passed, and
it confirmed the bug.

This suite takes us out of the loop by comparing get_MNV against annotators we
did not write.

```bash
cargo build -p get_mnv --release
python3 tests/differential/run.py
```

## The two oracles

| Dataset | Oracle | Why |
|---|---|---|
| `slippage` | `bcftools csq`, run live | A slippage join: CDS segments that abut, no intron, nothing spliced out. The exact shape that produced the bug. |
| `spliced` | `bcftools csq`, run live | A genuinely spliced two-exon transcript. The positive control: without it, "never emit splice terms" would also pass. |
| `mtb` | SnpEff 4.1l | The bundled *M. tuberculosis* VCF already carries `ANN=` fields from the original run. Real data, a real annotator, 950 calls. |

For the fixtures, both tools read the **same** reference and the **same** GFF3,
so a disagreement can never be blamed on them seeing different annotation. For
MTB each tool uses its own annotation of the same genome, so gene-model
differences are part of what the baseline records.

The fixtures mirror how these things are annotated in practice. In particular
the slippage fixture has **no `exon` rows**, because the NCBI GFF3 for ORF1ab
has none either: it is one gene and a two-segment CDS with a ribosomal-slippage
exception. Inventing exons there would change what `bcftools csq` reports and
would make the fixture a test of something that does not exist.

## The rule is not "the tools must agree"

get_MNV disagrees on purpose. It reports codon-level MNVs that per-variant
annotators split apart, which is the whole reason it exists. On the MTB dataset
that shows up directly: at position 2282377 SnpEff reports `synonymous_variant`
(`LOW`), because it looks at that variant alone. get_MNV sees that 2282376 sits
in the same codon, so the codon really goes `GTT` to `GCC` and the residue
really changes to Ala. get_MNV is right, and that difference must not fail CI.

So every disagreement is matched against `baseline.tsv`, a versioned list of
accepted differences **with the reason written down**. The suite fails only when
a difference appears that nobody has explained yet.

That file is worth reading on its own: it is the honest, position-by-position
account of where get_MNV differs from the established tools and why.

### When the suite fails

Either get_MNV has a bug, or the behaviour changed on purpose. Look at the
difference before touching anything. If it is intended, add a row to
`baseline.tsv` explaining it. `--update` rewrites the whole file from what is
observed now, with `TODO explain` in every note; use it to bootstrap, then write
the notes by hand. A baseline row with no reason is worse than no row at all.

The suite also reports baseline entries it no longer observes. That usually
means a fix removed a difference, and the row should be deleted.

## Checking that the suite can still fail

A test that cannot fail is worth nothing. `GET_MNV_BIN` points the suite at a
different build, so you can confirm it still catches the bug it was written for:

```bash
git worktree add /tmp/prefix <commit-before-the-fix>
cargo build --release --manifest-path /tmp/prefix/Cargo.toml
GET_MNV_BIN=/tmp/prefix/target/release/get_mnv python3 tests/differential/run.py
```

Against the commit before the slippage fix this exits 1 and points straight at
the boundary positions where get_MNV emitted `splice_donor` and `bcftools csq`
did not.

## Files

| Path | What it is |
|---|---|
| `run.py` | The harness: runs both tools, normalises, compares, checks the baseline |
| `baseline.tsv` | Accepted differences, each with its reason |
| `fixtures/` | The small synthetic genomes, committed so the run is hermetic |
| `make_fixtures.py` | Regenerates `fixtures/`; only run it to change or extend them |

`bcftools` must be on `PATH` for the fixture datasets; they are skipped with a
clear message when it is missing, so the MTB comparison still runs.
