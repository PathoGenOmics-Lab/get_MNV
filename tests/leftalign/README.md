# Indel left-alignment against bcftools norm

**English** · [Español](README.es.md)

`--left-align` moves each insertion and deletion to the leftmost position in its
repeat that describes the same change, and trims the bases REF and ALT share:
the normalization `bcftools norm -f ref.fa` performs. The unit tests pin some
twenty-five records to what bcftools writes for them, and a property test
compares the result with a search over every placement of the indel. Both were
written by the same hand as the code. This suite asks bcftools itself.

```bash
cargo build -p get_mnv --release
python3 tests/leftalign/run.py
```

## What it generates

Every record sits in a repeat of the bundled *M. tuberculosis* genome
(`example/MTB_ancestor.fas`), at least 60 bases from any other, from a fixed
random seed:

- **Indels in short repeats**: homopolymers and tandem repeats of 2 to 6 bases,
  deleting or inserting whole units or any other length, at a random place in
  the repeat or at its right end, where a 3'-shifted HGVS description puts it.
  Most insertions continue the repeat; some are random bases.
- **Whole-unit indels in minisatellites**, units of 7 to 100 bases, the only
  repeats in this genome long enough to move a record hundreds of bases.
- **Complex records**: a deletion and an insertion in one record, which
  bcftools trims but never moves.
- **The two ends of the contig**, where the anchor rule changes.

Each is spelled one of three ways: anchored on the base before the change, on
the base after it, or padded with context on both sides.

## What it compares

The same VCF goes to `bcftools norm -f` and to get_MNV with `--left-align`,
once as a plain `.vcf` and once as a `.vcf.gz`, since the two are read by
different readers. Each record carries an identifier in INFO, which get_MNV
keeps under `--keep-original-info`, so its output VCF is matched to the input
record by record. A record passes when get_MNV writes the POS, REF and ALT that
bcftools writes. The records also go through `--sample all` as a cohort of two
samples, in both containers, since a cohort is read from memory by a path of
its own.

The run also counts how many records get_MNV would write differently without
`--left-align`, and fails if there are none, since a comparison nothing could
fail proves nothing.

## The deliberate difference

bcftools reads an ambiguity code in the reference as its first base (`R` as
`A`), moves an indel through it, and writes that base into REF, which then
disagrees with the FASTA. get_MNV stops the indel at the code instead. The
example genome holds only A, C, G and T, so the suite never meets one; the unit
test `an_ambiguity_code_in_the_reference_stops_the_walk` pins the behaviour.

## Checking that it can still fail

A build with the contig-start rule broken (`start <= 1` for `start == 0` in
`src/io/left_align.rs`) reports the two records that reach the first base:

```
[indels.vcf] edge_start_del (deletion, period 1, right): 2 TG>G -> bcftools 1 TT>T, get_MNV 2 TG>G
[indels.vcf] edge_start_ins (insertion, period 1, left): 2 T>TT -> bcftools 1 T>TT, get_MNV 2 T>TT
```

Point the suite at any other binary with `GET_MNV=/path/to/get_mnv`.
