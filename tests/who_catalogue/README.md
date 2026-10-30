# WHO catalogue converter

**English** · [Español](README.es.md)

`tools/who_catalogue_to_tsv.py` turns the WHO catalogue of mutations in the
*M. tuberculosis* complex (second edition, 2023) into a table for get_MNV's
`--catalogue`, and `example/who_catalogue_2023.tsv` is its output for groups 1
and 2. A converter can be wrong in a way its output does not show: an entry
get_MNV reads but never matches looks like any other. So this suite runs the
entries against WHO's own genomic variants and asks get_MNV whether they match.

```bash
cargo build -p get_mnv
python3 tests/who_catalogue/run.py
```

## What it checks

- **Conversion.** `fixtures/master_rows.tsv` holds real rows of the master
  file, one for each way an entry is rewritten: a protein change, a change
  that needs two substitutions in a codon (rpoB `p.Ser450Phe`, `p.Ser441Gln`),
  a frameshift, `LoF`, a lost start, a lost stop, an in-frame deletion and
  duplication, a same-length deletion-insertion, promoter and rRNA changes, an
  indel outside a coding sequence, a form get_MNV does not write, and a
  whole-gene deletion, which cannot be converted. They must give exactly the
  expected entries, with and without the genomic coordinates, and the
  deletion must be listed with its reason.
- **Matching.** `fixtures/coordinates_rows.tsv` holds WHO's genomic variants
  for those rows. Each is put in a VCF on the example's MTBC ancestor, at most
  one per gene and 300 bp from any other, so none can share a codon or shift
  another's frame, and annotated with the converted table. Each must come out
  flagged with the entry WHO links it to, the two-change codons through the
  codon reading only. One cannot: WHO writes the change behind
  `Rv0565c p.TrpLeu266*` across codons 266 and 267, and get_MNV reports each
  codon on its own row, so no row carries both of the substitutions the
  entry's `g.` names. The suite expects that one miss and fails if it starts
  matching.
- **The shipped table.** `example/who_catalogue_2023.tsv` carries the licence
  notice, agrees with the converter entry by entry, and loads in get_MNV; the
  example sample G35894 carries no graded resistance variant and flags nothing.

## With WHO's files

Download the two files from
[WHO's repository](https://github.com/GTB-tbsequencing/mutation-catalogue-2023)
(`Final Result Files`) into one folder and pass it:

```bash
cargo build -p get_mnv --release
python3 tests/who_catalogue/run.py --who-dir path/to/folder
```

The suite then rebuilds the shipped table and requires it byte for byte, and
annotates every genomic variant WHO links to groups 1 and 2: 26,231 of them.
26,190 come out flagged with their own entry, 198 through the codon reading
only, and 9 more with the gene's `LoF` entry for the same drug. It fails if
more than these 32 are flagged for none of their drugs:

- 28 deletions that take a start codon away along with a whole number of
  codons (in ddn, fbiA, fbiC, fgd1, Rv2983, pepQ, gid, katG, pncA and ethA).
  WHO counts them as losses of function; get_MNV reports them as in-frame
  deletions, and `LoF` does not cover those.
- 2 insertions just before a start codon, tlyA `1917939 C>CAACG` and Rv2983
  `3339117 C>CG`, which WHO counts as losses of function and get_MNV as
  intergenic: the start codon itself is untouched.
- 2 changes whose consequence depends on the base beside them, gid
  `4407589 GC>CA` and tlyA `1917970 CT>TG`. On H37Rv they make a stop codon;
  on the MTBC ancestor, which differs from H37Rv next to them, they make
  `Ala205Trp` and `Leu11Trp`.

## Licence of the fixtures

The fixtures and the shipped table contain information from the WHO mutation
catalogue, which is made available under the
[ODC Attribution License](https://opendatacommons.org/licenses/by/1-0/).
