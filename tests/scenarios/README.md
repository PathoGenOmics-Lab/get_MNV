# End-to-end scenario tests

This directory contains a Python test harness that builds synthetic FASTA,
GFF, VCF (or iVar TSV) and BAM inputs from declarative scenario
definitions, runs `get_mnv` on each one, and checks the output TSV against
pre-defined expected fields.

It is complementary to the Rust unit/integration tests in
[`tests/integration.rs`](../integration.rs). The Rust tests exercise the
library at the function level; these scenario tests exercise the full
binary on realistic inputs, including:

- Read support computed from BAM CIGAR
- Codon grouping across exon junctions (multi-exon CDS models)
- Negative-strand gene annotation
- iVar TSV input with `+SEQ` / `-SEQ` notation
- Multi-contig inputs

## Quick start

Requirements:

- Python 3.10 or newer (only the standard library is used)
- `samtools` on `PATH` (or set the `SAMTOOLS` environment variable)
- A built `get_mnv` binary — by default the suite looks for
  `target/debug/get_mnv`, then `target/release/get_mnv`, then
  `dist/get_mnv`, then `get_mnv` on `PATH`. Override with `GET_MNV=/path/to/get_mnv`.

```bash
# Build get_mnv first
cargo build

# Run all 36 scenarios
python3 tests/scenarios/run.py

# Run a subset by name prefix
python3 tests/scenarios/run.py 18 22 27

# Use a release binary instead of debug
GET_MNV=$(pwd)/target/release/get_mnv python3 tests/scenarios/run.py

# Use samtools from a conda environment
SAMTOOLS=/path/to/conda/envs/bio/bin/samtools python3 tests/scenarios/run.py
```

Intermediate files (FASTA, GFF, VCF, BAM, `get_mnv` output, log) for each
scenario are written under `tests/scenarios/work/<scenario_name>/` and
overwritten on each run. The directory is `.gitignore`d.

## Visual overview (PDF)

[`plot_scenarios.py`](plot_scenarios.py) renders every scenario into a single
illustrative PDF (`scenarios_overview.pdf`). It needs `matplotlib` in addition
to `samtools`.

```bash
SAMTOOLS=/path/to/samtools python3 tests/scenarios/plot_scenarios.py
```

The pages are grouped into thematic **sections** (SNV/MNV, indel
classification, indels combined with SNVs, negative strand, eukaryotic
multi-exon transcripts, input formats, edge cases), preceded by a cover, a
glyph legend and a table of contents. Each scenario page shows:

- the **mapped reads** as a small alignment pileup, marking SNVs, insertions,
  deletions and intron skips (CIGAR `N`), the overlapping genes/CDS, and the
  input variants;
- a **reading-frame / amino-acid track** — alternating shaded codon triplets
  plus the reference amino acid per codon, with the alternate AA in red when a
  variant makes the codon non-synonymous (strand- and phase-aware);
- two tables: **"Without get_mnv"** (the raw caller calls — position, allele,
  length-based frame guess) versus the real **"get_mnv output"**, with a
  one-line **"What get_mnv adds"** summary between them.

The generated PDF and any `plot_*.png` previews are `.gitignore`d; regenerate
them from the script.

## Mini-genome layout

All scenarios share a synthetic reference (`framework.REFERENCE_SEQ`)
composed of two contigs designed to keep codon math explicit:

### `chr_test` (1300 bp)

| Range       | Feature        | Content                                                   |
|-------------|----------------|-----------------------------------------------------------|
| 1–300       | `geneA` (+)    | `ATG` + 98 × `GCT` + `TAA` — Met, 98×Ala, stop           |
| 301–400     | filler         | `A` × 100                                                 |
| 401–700     | `geneB` (−)    | RC of the standard CDS: `TTA` + 98 × `AGC` + `CAT`        |
| 701–800     | filler         | `A` × 100                                                 |
| 801–900     | `geneC` exon 1 | `ATG` + 32 × `GCT` + `G` (codons 1–33 + base 1 of codon 34) |
| 901–1000    | `geneC` intron | `T` × 100                                                 |
| 1001–1200   | `geneC` exon 2 | `CT` + 65 × `GCT` + `TAA` (bases 2–3 of codon 34, codons 35–100) |
| 1201–1300   | filler         | `A` × 100                                                 |

### `chr_test2` (600 bp)

| Range  | Feature     | Content                                |
|--------|-------------|----------------------------------------|
| 1–300  | `geneD` (+) | `ATG` + 98 × `GCT` + `TAA` (same as geneA) |
| 301–600 | filler     | `A` × 300                              |

Codon math summary:

- `geneA` codon N (+): positions `(3N-2, 3N-1, 3N)`
- `geneB` codon N (−): positions `(701-3N, 702-3N, 703-3N)`; mRNA base = RC of genomic
- `geneC` codon 34 spans the exon 1/intron/exon 2 junction: pos 900 (exon 1) + pos 1001-1002 (exon 2)

## Validated scenarios

36 scenarios are currently defined in [`scenarios.py`](scenarios.py).
Each one declares input variants, BAM read groups (with optional
operations: SNV substitution, insertion, deletion, intron skip) and the
expected TSV rows.

### Core SNP / MNV / SNP/MNV

| # | Name | What it validates |
|---|------|-------------------|
| 01 | `snp_simple` | Single SNV `Ala10Thr` (geneA codon 10 GCT→ACT) with 20/20 read support |
| 02a | `snp_mnv_full_phasing` | Two SNVs in the same codon, all reads carry both → `SNP/MNV`, `Ala10Ser`, MNV reads 20, SNP-only reads 0 |
| 02b | `vcf_mnp_decomposed` | A VCF MNP record (`REF=GC ALT=TA`) is decomposed into two SNV components and reported as `SNP/MNV` with `event_class=mnv` |
| 03 | `snp_mnv_mixed` | Partial phasing: 10/10/10 reads (both / only SNV1 / only SNV2) → frequencies 0.3333 / 0.3333 / 0.3333 |

### In-frame and frameshift indels

| # | Name | What it validates |
|---|------|-------------------|
| 04 | `ins_inframe_cds` | In-frame insertion of `GCT` (1 Ala) inside CDS, 20 reads, `INS:30:+GCT` |
| 05 | `del_frameshift_cds` | 1 bp frameshift deletion, `Ala11Leufs`, `DEL:31:G` |
| 24 | `large_inframe_insertion` | 9 bp in-frame insertion (`+GCTGCTGCT` = 3 Ala) |
| 25 | `large_inframe_deletion` | 6 bp in-frame deletion (2 consecutive Ala), `DEL:31-36:GCTGCT` |

### Complex haplotypes (indel + SNV/MNV in cis)

| # | Name | What it validates |
|---|------|-------------------|
| 06 | `indel_plus_snv_haplotype` | Indel + SNV >3 bp apart: NOT merged into `complex_indel` (out of the 3 bp local window) |
| 07 | `fs_del_plus_downstream_snv` | Frameshift propagates to a downstream SNV → `(fs)` suffix in AA Change and `Synonymous (frameshift)` Change Type |
| 08 | `inframe_ins_inside_codon_with_mnv` | In-frame insertion inside a codon + MNV in the same codon: emits 5 rows including 2 `complex_indel`, the MNV row marked `Indel overlap` / `Unknown`, the insertion alone, and `complex_indel` ins+SNV |
| 09 | `fs_del_with_snv_overlap` | Frameshift deletion overlapping the codon + SNV: the SNV row gets `Indel overlap` / `Unknown`; the deletion-alone row has `Event Reads = 0` (no read carries only the del); the `complex_indel` carries the frameshift `Ala10Cysfs` |

### Negative-strand gene (`geneB`)

| # | Name | What it validates |
|---|------|-------------------|
| 10 | `minus_snp_simple` | Single SNV in `geneB` (pos 673 C>T) — mRNA codon 10 `GCT`→`ACT` → `Ala10Thr` |
| 11 | `minus_mnv_same_codon` | MNV in `geneB` codon 10 (pos 671 A>T + pos 673 C>A in cis) — mRNA `GCT`→`TCA` → `Ala10Ser` |
| 12 | `minus_fs_del` | 1 bp frameshift deletion inside `geneB` |

### Multi-exon CDS (`geneC`, with `--gff-features CDS`)

| # | Name | What it validates |
|---|------|-------------------|
| 13 | `multiexon_snp_exon2` | SNV in exon 2 — confirms transcript-aware codon position resolution |
| 14 | `multiexon_junction_snp` | SNV in the junction-spanning codon (base 1 in exon 1 pos 900) — requires the spliced model |
| 15 | `multiexon_junction_mnv` | MNV with one base in exon 1 (pos 900) and another in exon 2 (pos 1002) — adjacent in spliced mRNA but ~100 bp apart in the genome. Validated via N-CIGAR reads spanning the intron |

### Operational edge cases

| # | Name | What it validates |
|---|------|-------------------|
| 16 | `no_bam_coverage` | VCF declares a SNV at a position where the BAM has no overlapping reads — `get_mnv` still emits the row with empty read-support columns |
| 17 | `min_snp_frequency_filter` | SNV present at 10% (2/20 reads) is dropped from the output when `--min-snp-frequency 0.5` is set |

### Amino-acid edge cases

| # | Name | What it validates |
|---|------|-------------------|
| 18 | `stop_gained_via_mnv` | Three-SNV MNV in codon 50 `GCT`→`TAA` → `Ala50Ter`, `Change Type = Stop gained` |
| 19 | `start_codon_altered` | SNV in start codon ATG at pos 2 → `Met1Thr`. NOTE: `get_mnv` does not have a dedicated `Start lost` Change Type; the row is reported as `Non-synonymous` |
| 20 | `stop_lost` | SNV at the stop codon `TAA` (pos 298) → `Change Type = Stop lost` |

### Complex alleles

| # | Name | What it validates |
|---|------|-------------------|
| 21 | `intron_variant` | SNV in the geneC intron (pos 950) → annotated as `intergenic` under `--gff-features CDS` |
| 22 | `multiallelic_split` | Multiallelic VCF `pos 28 REF=G ALT=A,T` with `--split-multiallelic`: each ALT now produces an independent annotation row (regression test for the fix in commit 9ea2aed) |
| 23 | `delins_subst_plus_del` | `REF=GCT ALT=GA`: a substitution + 1 bp deletion compound allele yields a frameshift INDEL row |

### Multi-contig and iVar TSV input

| # | Name | What it validates |
|---|------|-------------------|
| 26 | `multicontig` | Two variants in two different contigs (`chr_test` and `chr_test2`) produce two correctly annotated rows |
| 27 | `ivar_tsv_snv` | iVar TSV input with a simple SNV — same annotated row as the VCF equivalent |
| 28 | `ivar_tsv_insertion` | iVar TSV `+GCT` notation is normalized to the VCF-style anchored allele and annotated as `INS:30:+GCT` |
| 29 | `ivar_tsv_deletion` | iVar TSV `-G` notation is normalized to `DEL:31:G` and produces the same frameshift row |

### Indel refinements (indels branch)

These exercise the indel-handling refinements added on the `indels` branch
(commit 49d2f09): stop detection for in-frame indels and the
`--frameshift-min-freq` downstream-propagation gate.

| # | Name | What it validates |
|---|------|-------------------|
| 30 | `stop_gained_inframe_ins` | In-frame insertion of `TAA` after pos 30 forms a premature stop codon → `Change Type = Stop gained`, `AA Changes = 10_11ins*` (instead of the generic `In-frame Indel`). Driven by `indel_stop_effect`, which compares the number of stop residues in the ref vs alt protein |
| 31 | `stop_lost_inframe_del` | In-frame 3 bp deletion of the terminal stop `TAA` (pos 298-300, `DEL:298-300:TAA`) → `Change Type = Stop lost`, `AA Changes = *100del` |
| 32 | `fs_gate_default_propagates` | Low-frequency upstream frameshift deletion (`AF=0.20`) + a downstream SNV. With the default `--frameshift-min-freq 0.0` the frameshift propagates → the downstream SNV is labelled `Synonymous (frameshift)` / `Ala13Ala (fs)` |
| 33 | `fs_gate_high_freq_suppressed` | Identical inputs to scenario 32 but run with `--frameshift-min-freq 0.5`. The upstream deletion (`AF=0.20`) does not pass the gate, so the frameshift is **not** propagated → the downstream SNV is a plain `Synonymous` / `Ala13Ala` |

Scenarios 32 and 33 are an A/B pair: same VCF (now carrying `AF` in `INFO`)
and same BAM, differing only in the `--frameshift-min-freq` flag. Together
they show the gate affects only downstream frameshift propagation — the
upstream deletion's own row stays `Frameshift Indel` (`Ala11Leufs`) in both.

## Framework architecture

- [`framework.py`](framework.py) defines the building blocks:
  - `REFERENCE_SEQ`, `REFERENCE_SEQ2`, `GFF_GENE_ONLY`, `GFF_CDS_MULTIEXON`
  - Data classes: `Op`, `ReadGroup`, `VcfRecord`, `IvarRecord`, `ExpectedRow`, `Scenario`
  - Builders: `write_fasta`, `write_gff`, `write_vcf`, `write_ivar_tsv`, `write_bam`
  - Driver: `run_get_mnv`, `parse_tsv`, `compare`, `run_scenario`
  - `VcfRecord(af=...)` emits an `AF=` INFO field (default omitted) so
    scenarios can drive frequency-gated logic such as `--frameshift-min-freq`
  - `Scenario.extra_cli_args` passes extra flags through to `get_mnv`
- [`scenarios.py`](scenarios.py) declares the 36 scenarios.
- [`run.py`](run.py) is the CLI driver:
  ```bash
  python3 run.py             # all scenarios
  python3 run.py 18 22       # only scenarios starting with "18" or "22"
  ```

## Adding a new scenario

1. Pick a unique numeric prefix and descriptive name.
2. Choose where the variant falls in the mini-genome (see the codon math
   above). Common positions:
   - `geneA` codon 10: pos 28-30
   - `geneA` codon 50: pos 148-150
   - `geneB` codon 10 (− strand): pos 671-673
   - `geneC` junction codon 34: pos 900 + pos 1001-1002
3. Build the read groups with the operations needed to support the
   variant in BAM. The `Op` kinds are:
   - `Op("snv", pos=P, seq="X")` — substitute base at ref pos `P` with `X`
   - `Op("ins", pos=P, seq="ABC")` — insert `ABC` after ref pos `P`
   - `Op("del", pos=P, length=N)` — delete `N` bases starting at `P`
   - `Op("skip", pos=P, length=N)` — emit a CIGAR `N` (intron skip), for
     spliced reads in multi-exon scenarios
4. Declare the expected TSV rows with the fields you want to verify
   (`positions`, `gene`, `base_changes`, `aa_changes`, `variant_type`,
   `change_type`, `event_class`, `event_components`, `*_reads`,
   `*_frequencies`, etc.). Unset fields are not checked.
5. Optionally set `expected_row_count` to assert the exact total number
   of rows produced — useful when the local haplotype window emits more
   than one row (see scenarios 08 and 09).
6. Append the scenario to `ALL_SCENARIOS` and run `python3 run.py`.

The driver will show `PASS` or `FAIL` for each scenario, and on failure
print the actual rows produced so you can adjust the expectation.

## Known behaviours documented by these tests

- `geneB` (negative strand) — `Reference Codon` and `MNV Codon` columns
  show the **genomic forward** sequence (e.g. `AGC`), not the mRNA
  sequence (`GCT`). The AA effect is still correct.
- Multi-exon CDS without a `Name=` attribute on the CDS feature reports
  the gene column as the first CDS ID (e.g. `cds-geneC-e1`). Standard
  NCBI/Ensembl GFFs include `Name=`, so this only affects hand-rolled
  GFFs.
- A SNV/MNV that alters the initiator Met (protein position 1) is reported
  with Change Type `Start lost` (scenario 19); `Met1` → stop stays
  `Stop gained`, and internal Met substitutions are unaffected.
- `--split-multiallelic`: each ALT at the same codon position now emits
  an independent annotation row.
- In-frame indels that **create or remove a stop codon** are reported as
  `Stop gained` / `Stop lost` rather than the generic `In-frame Indel`.
  Frameshift indels keep the `Frameshift Indel` label (they almost always
  introduce a downstream stop, so the distinction would be uninformative).
- `--frameshift-min-freq F` gates **downstream frameshift propagation**: an
  upstream indel only shifts the frame of downstream codons if its reported
  allele frequency (VCF `AF`/`FREQ`/`AD`, iVar `ALT_FREQ`) is ≥ `F`. The
  default `0.0` reproduces the historical always-propagate behaviour, and
  indels without a known frequency always propagate. The gate never changes
  the indel's own classification.
- When `--gff-features` is not given and the GFF contains `CDS` features, the
  phase/splice-aware CDS model is selected automatically (scenario 34); pass
  `--gff-features gene` to keep the legacy whole-gene behaviour.
- A variant downstream of a frameshift-introduced premature stop is reported
  with Change Type `Downstream of premature stop` instead of carrying an `(fs)`
  annotation as if it were translated (scenario 35). Ordinary frameshift
  propagation, where no early stop is introduced, is unchanged.
- The local haplotype window is 3 bp (`LOCAL_HAPLOTYPE_JOIN_DISTANCE`)
  and accepts up to 8 events
  (`MAX_LOCAL_HAPLOTYPE_VARIANTS`). Events further apart produce
  separate rows instead of being merged into a `complex_indel`.
