# Output Formats

get_MNV can write TSV, VCF, BCF, and JSON metadata files.

## Default TSV Output

Default file name:

```text
<input_name>.MNV.tsv
```

With `--sample all`, one file is written per sample instead, carrying the sample
in its name:

```text
<input_name>.sample_<SAMPLE>.MNV.tsv
```

Characters that cannot appear in a file name are replaced, so a sample called
`sample/1:bad` becomes `sample_1_bad`.

Use this format for spreadsheets, downstream parsing, and quick inspection.

Main columns:

| Column | Meaning |
|---|---|
| `Chromosome` | Contig name |
| `Gene` | Gene or feature name. Intergenic variants are marked as `intergenic`. |
| `Positions` | One position for SNPs, multiple comma-separated positions for MNVs, in ascending genomic order on both strands. The per-SNV columns (`Reference Bases`, `Base Changes`, `SNP AA Changes`, `SNP Codon`, `Event Components`, and the read-support columns) are parallel arrays in that same order. |
| `Reference Bases` | Reference bases at those positions. |
| `Base Changes` | Alternative bases. |
| `AA Changes` | Amino acid change after combining all SNVs in the codon. |
| `SNP AA Changes` | Amino acid change for each SNV considered separately. |
| `Local AA Changes` | Exon-local numbering, which is what get_MNV reported before 1.1.2. Identical to `AA Changes` when the feature has no transcript context (a genes TSV, a prokaryotic or single-exon CDS, and the first exon of a spliced model), and different for the later exons of a spliced transcript, where `AA Changes` counts from the start of the protein and this counts from the start of the exon. |
| `Local SNP AA Changes` | Per-SNP amino acid changes in local numbering. |
| `Variant Type` | `SNP`, `MNV`, `SNP/MNV`, or `INDEL`. |
| `Change Type` | Synonymous, non-synonymous, stop gained/lost, unknown, etc. |
| `Reference Codon` | Original codon, in the transcript's own orientation: on a minus-strand gene these are the coding-strand bases, not the genomic ones, so the codon always translates to the amino acid reported beside it. |
| `SNP Codon` | Codon with individual SNP substitutions, same orientation as `Reference Codon`. |
| `MNV Codon` | Codon with all grouped substitutions, same orientation as `Reference Codon`. |
| `Event Class` | Canonical allele event class: `snp`, `mnv`, `insertion`, `deletion`, `delins`, `complex_indel`, or `symbolic`. |
| `Event Components` | REF/ALT decomposition such as `SNV:10:A>G`, `INS:10:+T`, or `DEL:11-12:TG`. |
| `SO Term` | Sequence Ontology consequence term (`missense_variant`, `synonymous_variant`, `stop_gained`, `start_lost`, `frameshift_variant`, `inframe_deletion`, `intergenic_variant`, …). Variants near an internal exon-exon junction of a spliced transcript also carry a splice term: `splice_donor_variant` / `splice_acceptor_variant` (the two essential intronic bases at each intron end, `HIGH`) or `splice_region_variant` (the exon's first or last 3 bases, or the intron's 3rd-8th bases, `LOW`). An exonic coding change near a junction is combined, e.g. `missense_variant&splice_region_variant`. A junction means a real intron: CDS segments that abut or overlap (a ribosomal-slippage join such as SARS-CoV-2 ORF1ab) are one continuous reading frame and carry no splice terms. A variant inside an intron but away from its splice sites is `intron_variant` (`MODIFIER`), reported against its gene rather than as intergenic, and a variant in a feature declared non-coding is `non_coding_transcript_exon_variant` (`MODIFIER`) with no amino-acid change. |
| `Impact` | Predicted impact following SnpEff/VEP conventions: `HIGH`, `MODERATE`, `LOW` or `MODIFIER`. A combined splice/coding consequence keeps the more severe impact. |
| `Grantham` | Grantham distance and conservation category of a missense change (e.g. `177 (radical)`); `-` for synonymous, nonsense or non-coding changes. |
| `MNV Consequence Shift` | How the combined MNV compares with its individual SNVs: `MNV-gained` (more severe than any single SNV, which is what per-SNV annotators miss), `MNV-masked` (a nonsense SNV rescued by its neighbour) or `Concordant`. `-` for single SNVs. |
| `DBS Class` | COSMIC-style doublet base substitution class for an MNV of two adjacent single-base substitutions, e.g. `CC>TT` (reverse-complement collapsed, so `GG>AA` reports as `CC>TT`). `-` for single SNVs, indels, and non-adjacent or 3-SNV MNVs. |
| `Declared Phase` | The phase the **input VCF** declared for this row's alleles, as `cis:12345` (verdict and `PS` phase set) or just `cis` when the caller phased without a phase set. Read from a `|`-separated `GT`; a `/` genotype is the caller saying it did not resolve the phase, and is reported as `-`. `|contradicted-by-reads` is appended when the BAM leaves the claim no room: a declared cis that not one spanning read carries whole, or a declared trans that every informative read carries. This is the caller's claim, not an observation; the phasing columns beside it are the evidence. `-` for single-position rows and unphased input. |
| `MNV Phasing Support` | BAM-derived phasing (linkage) support: among the reads that observe *every* position of the codon and carry the least-supported constituent SNV, the fraction that also carry the full MNV haplotype. `1.0000` = perfect co-occurrence (a genuine haplotype); low values mean the SNVs largely fall on different molecules (a same-codon coincidence, not a real MNV). `0.0000` is a finding: reads did span the codon and none carried both. `-` means the question could not be answered: no `--bam`, a single SNV, or no read reaching across the codon (common when a codon straddles an intron and the fragments are shorter than the intron). |
| `Haplotype LD` | Read-level linkage disequilibrium (`D'`) between the variants this row claims travel together, over the molecules that observed them. It covers both kinds of multi-variant row: a codon MNV and a local indel haplotype. A co-occurrence ratio cannot separate a haplotype from an accident of frequency: two substitutions on 90% of molecules each are found together on 81% of them by arithmetic alone, and the ratio calls that 0.9. `D'` measures the excess over what the two frequencies predict, normalised by the most it could have been. `+1` = they travel together as far as their frequencies allow, so the MNV is one haplotype. `~0` = they co-occur exactly as chance predicts, so they merely share a codon. `-1` = they exclude each other: both present, never on one molecule, which in a haploid population is two competing lineages rather than one variant. With three or more variants the weakest pair decides, since the row claims one molecule carries all of them. It answers a different question from the read count beside it: the count is how many molecules *are* this combination, while `D'` is whether its variants co-occur more than their own frequencies predict, so a haplotype carried by few molecules can still be perfectly linked and one carried by many can be a coincidence. `-` when no molecule observed the variants together, or when one of them is on every such molecule or none, which leaves nothing to correlate. Only present with `--bam`. |
| `Haplotype LD p` | Two-tailed Fisher exact p-value for that table, so a `D'` of 1.0 from four molecules is not read as one from four hundred. `-` under the same conditions as the column above. Only present with `--bam`. See [Linkage](linkage.md). |
| `MNV Phasing Reads` | How many reads `MNV Phasing Support` was computed from, so `1.0000` from 3 reads is not read as `1.0000` from 300. `-` exactly when that column is `-`, which is not the same condition as the two linkage columns beside it. Only present with `--bam`. |
| `Frameshift Phasing` | What the reads said about this codon sharing molecules with each upstream indel, as `trans:1234:0/18`: the verdict, the indel's position, and the cis reads out of the reads able to answer. Several are joined with ` | `. A codon that is not labelled frameshifted otherwise looks the same whether the reads proved the indel is on other molecules or nobody asked; `-` is that second case. Only present with `--bam`. |
| `NMD Prediction` | Nonsense-mediated decay prediction for a premature stop under the 50-nt rule: `NMD-triggering` when the PTC is more than 50 nt upstream of the last exon-exon junction, `NMD-escaping` when it is in the last exon or within 50 nt of that junction. `-` for variants without a premature stop and for transcripts with no exon-exon junction. Requires a spliced (GFF/GTF transcript) CDS model whose segments are separated by real introns; a single CDS segment, or segments joined by ribosomal slippage, has no junction. |
| `HGVS g.` | HGVS genomic descriptor: `g.100A>G` for an SNV, the allele-bracket `g.[28G>T;30T>A]` for an MNV, and `g.101_102del` / `g.100_101insTG` / `g.101delinsC` for indels. Not 3'-shifted (uses the input allele placement) and carries no reference-accession prefix. |
| `HGVS c.` | HGVS coding descriptor for a coding substitution, numbered from the CDS start with coding-strand bases: `c.30A>G` (SNV) or the allele bracket `c.[28G>A;30T>C]` (MNV). `-` for indels and non-coding variants; the protein change (`p.`) in the AA columns conveys indel consequences. |

Extra columns when `--bam` is used:

| Column | Meaning |
|---|---|
| `SNP Reads` | Reads carrying each individual SNV **without** the full MNV haplotype. On a row where every read carries the whole haplotype this is `0` for each constituent and the count lives in `MNV Reads`, so the two columns partition the support rather than double-count it. |
| `SNP Forward/Reverse Reads` | Strand-specific counts for the reads above. |
| `MNV Reads` | Reads supporting the full MNV haplotype. |
| `MNV Forward/Reverse Reads` | Strand-specific MNV support. |
| `Total Reads` | Depth at the variant positions. |
| `SNP Frequencies` | Per-position SNP frequencies. |
| `MNV Frequencies` | MNV haplotype frequency. |
| `Event Reads` | Exact reads supporting an indel/complex event. |
| `Event Forward/Reverse Reads` | Strand-specific exact event support. |
| `Event Depth` | Reads with an observed allele across the indel/complex event span. |
| `Event Frequency` | Exact event reads divided by event depth. |

Exact event support is CIGAR-aware. A read must reconstruct the same local ALT
sequence and, for complex haplotypes, contain the expected insertion and
deletion components. This prevents net-neutral insertion/deletion combinations
from being counted as support merely because their sequence looks like an MNV.

Frequency columns are calculated from BAM support. `--min-snp-frequency` and
`--min-mnv-frequency` use these same BAM-derived values. The filters are
independent: `--min-snp-frequency` applies to individual SNP observations, and
`--min-mnv-frequency` applies to phased MNV haplotypes. In mixed `SNP/MNV`
calls, a row or VCF record is kept when either component passes its own active
threshold.
Read-count and strand-support filters (`--snp`, `--mnv`, `--min-snp-strand`,
and `--min-mnv-strand`) follow the same independent SNP/MNV behavior.

When a codon-level MNV overlaps an indel, the MNV row is kept as a positional
context row but its amino-acid effect is marked `Unknown` with
`Change Type = Indel overlap`. If BAM reads support the full combined event,
get_MNV emits a separate exact `complex_indel` row with the combined REF/ALT,
event components, and event read support.

Indel overlap follows VCF interbase semantics. Deletions overlap a feature by
their deleted reference span. Insertions overlap a feature only when the inserted
sequence falls between two reference bases inside that feature, so an insertion
anchored at the final feature base is reported outside that feature.

Example:

```text
Chromosome	Gene	Positions	Base Changes	AA Changes	Variant Type	Change Type
MTB_anc	Rv0095c_Rv0095c	104838	T	Asp126Glu	SNP	Non-synonymous
MTB_anc	Rv0095c_Rv0095c	104941, 104942	T, G	Gly92Gln	SNP/MNV	Non-synonymous
```

## VCF Output

Write VCF with:

```bash
--convert
```

or write both TSV and VCF with:

```bash
--both
```

Default file name:

```text
<input_name>.MNV.vcf
```

Use `--vcf-gz` for compressed output:

```text
<input_name>.MNV.vcf.gz
```

Common INFO fields:

| Field | Meaning |
|---|---|
| `GENE` | Gene or feature name |
| `AA` | Amino acid change |
| `CT` | Change type |
| `TYPE` | Variant type |
| `EC` | Canonical allele event class |
| `COMP` | REF/ALT event components |
| `ODP` | Original depth from the input variant file |
| `OFREQ` | Original allele frequency from the input variant file |
| `SR`, `SRF`, `SRR` | SNP reads: total, forward, reverse |
| `MR`, `MRF`, `MRR` | MNV reads: total, forward, reverse |
| `DP` | Depth recalculated from BAM |
| `FREQ` | Frequency recalculated from BAM |
| `ER`, `ERF`, `ERR` | Exact indel/complex event reads: total, forward, reverse |
| `EDP` | Exact event depth for indel/complex alleles |
| `EFREQ` | Exact event frequency for indel/complex alleles |
| `SBP` | SNP strand-bias p-value |
| `MSBP` | MNV strand-bias p-value |
| `SO`, `IMPACT` | Sequence Ontology consequence term and predicted impact |
| `GD` | Grantham distance of a missense change |
| `MNVSHIFT` | Combined MNV consequence vs. its individual SNVs |
| `DBS` | COSMIC-style doublet class for adjacent 2-SNV MNVs (e.g. `CC>TT`) |
| `MNVPS` | MNV phasing support (of the codon-spanning reads carrying the limiting SNV, the fraction carrying the full haplotype) |
| `MNVPR` | Reads that ratio was computed from |
| `FSPH` | Read-level phasing with each upstream indel, as verdict:position:cis/informative |
| `DPHASE` | Phase the input VCF declared for this row, verdict:phase_set |
| `LD` | Linkage disequilibrium D-prime between the codon's substitutions |
| `LDP` | Fisher exact p-value for that linkage table |
| `NMD` | Nonsense-mediated decay prediction for a premature stop (50-nt rule) |
| `HGVSG` | HGVS genomic descriptor (MNV allele-bracket `;` percent-encoded) |
| `HGVSC` | HGVS coding descriptor for a coding substitution (`;` percent-encoded) |

The VCF header records the get_MNV version, command line, and thresholds used.
When `--emit-filtered` is enabled, VCF records below read-support, frequency,
strand-support, or strand-bias thresholds are written with FILTER tags such as
`LowSupport`, `LowFrequency`, `StrandSupport`, or `StrandBias`; otherwise they
are skipped.

## BCF Output

Write BCF with:

```bash
--bcf
```

BCF requires VCF output mode, so use it with `--convert` or `--both`.
This is output conversion only; BCF is not accepted as an input format.

Default file name:

```text
<input_name>.MNV.bcf
```

## JSON Files

### Summary JSON

Write with:

```bash
--summary-json run.summary.json
```

Includes:

- Input file checksums
- Per-contig variant counts
- Global variant counts
- Runtime timings
- Output paths

!!! warning "These counts are what get_MNV produced, before the output filters"
    `produced_variants` and the per-type counts beside it are taken after
    annotation and before the read-support, frequency and strand filters that
    decide what reaches the TSV. On a filtered run the two disagree on purpose:
    the same command can log `produced variants=941` and write a one-row TSV,
    and the HTML report, which counts rows, will say `1`. Read the summary as
    what the annotation found and the TSV as what passed.

### Run Manifest

Write with:

```bash
--run-manifest run.manifest.json
```

Includes the summary plus:

- Command line
- Tool version
- Output file checksums
- Timestamp

### Error JSON

Write errors as JSON with:

```bash
--error-json run.error.json
```

This is useful in automated pipelines.

## Interactive HTML report

`--report <FILE.html>` writes a single self-contained HTML file for exploring the
called variants. It embeds its data, loads no external scripts or fonts, and
therefore opens offline by double-clicking and can be attached to an email or
archived with the results.

[**Open the example report**](assets/example-report.html){ target=_blank }: 941
variants called from the single-sample dataset bundled in `example/`, so it is the
real output of the command below rather than a mock-up. Regenerate it with
`scripts/build_example_report.sh`. The cohort views (the matrix rows, haplotype
recurrence) fill out with `--sample all` or `--report-from` across several samples.

```bash
# Report for one run
get_mnv --vcf sample.vcf --fasta ref.fasta --gff ref.gff --report sample.html

# One report covering every sample of a multi-sample VCF
get_mnv --vcf cohort.vcf --fasta ref.fasta --gff ref.gff --sample all --report cohort.html

# Cohort processed one sample per run: aggregate the TSVs afterwards
get_mnv --report-from results/*.MNV.tsv --report cohort.html
```

With `--report-from` no pipeline runs: the report is built from get_MNV TSV files
that already exist, which is the usual shape in a Nextflow or Snakemake workflow
that calls get_MNV once per sample. Each file becomes one sample, labelled by its
file name with the `.MNV.tsv` suffix removed (`TB-001.MNV.tsv` becomes `TB-001`).

The report contains:

- **A headline count** of the variants currently shown, with supporting tiles for
  samples, genes, MNV rows and high-impact variants, and the **consequence
  distribution** by Sequence Ontology term. All follow the active filters.
- **A variant matrix**, plotted as a small genome browser: samples down the side,
  real genomic coordinates across the top, coloured by the alternate base. Scroll
  to zoom around the cursor, drag to pan, drag on the whole-contig strip to jump,
  and double-click to reset. Dragging across the whole-contig strip selects a
  range; a plain click on it recentres. The region box takes `start-end`,
  `contig:start-end`, a single coordinate, or a gene name, and frames it. Above
  the ruler sit density lanes: one small chart each for all calls, SNP, MNV,
  indel, HIGH impact and the number of distinct samples with a call, binned by
  genomic position. Each lane carries its own scale and its own maximum, so a
  rare class is readable next to a common one instead of being flattened at the
  base of a stack, and the `Tracks` control switches between all lanes, a compact
  set and none. A gene track marks the
  extent of each gene's called sites (which is where calls were made, not the
  annotated gene boundary, which the report does not carry). Zooming in far enough shows
  the base letter inside each cell. Samples can be ordered by shared profile
  (identical patterns land together), by variant count or by name, and each sample
  label carries a bar of how many calls it has in the visible window, which is why
  there is no separate per-sample chart. Colour means one thing only: the
  nucleotide hues belong to the matrix cells, every magnitude is neutral, and the
  reserved status colour marks HIGH impact. Positions
  called together on the same reads (codon MNVs, phased complex indels) are marked
  with a tick above their columns, and a contig selector appears for multi-contig
  data, since a continuous coordinate axis cannot span contigs.

    A cell is only ever "an ALT call" or **"not called"**. get_MNV output cannot
    distinguish a reference base from a position with no coverage, so a blank
    cell is never presented as reference. Read that state as "reference or no
    coverage".

- **Read-backed haplotypes**: the allele combinations get_MNV actually observed
  on the same reads, that is codon-level MNVs and locally phased complex indels,
  ranked by how many samples carry each one and shown with their phasing support
  where a BAM was used. Long-range phasing between distant sites is not inferred
  and is never displayed.
- **A sortable, filterable table** of every variant, virtualised so tens of
  thousands of rows stay responsive. Every column carries its own filter: a
  checkbox list for sample, contig, gene, variant type, consequence and impact
  (several values at once, with a search box when the list is long), a contains
  match for position, base change and amino-acid change, and a min/max range for
  Grantham and frequency. Column filters combine with each other and with the
  free-text search box, and they drive the whole page: the headline counts, the
  charts, the matrix and the haplotype panel all follow the same selection.
- **A detail panel** for the selected variant with its location, consequence,
  HGVS `g.`/`c.`/`p.` descriptors, Grantham distance, DBS class, NMD prediction,
  codons, event components and read support.
- **Export filtered TSV**, which downloads exactly the rows currently shown.
- **Links to the source repository and the documentation** in the masthead. They are
  ordinary hyperlinks: nothing is fetched when the report is opened, so it stays
  self-contained offline.

The report follows the operating-system light or dark theme and has its own
toggle. Because it is built from the TSV, `--report` needs TSV output: it works
with the default output mode and with `--both`, but not with `--convert` alone or
with `--dry-run`. Read-support columns (frequency, depth) are populated only when
the run used `--bam`.

Size scales with the number of variants. Repeated fields are dictionary-encoded,
so a cohort of tens of thousands of variants stays in the low megabytes.

## Notes

- For MNV records, depth and frequency are calculated from reads spanning all
  positions in the grouped haplotype.
- Frequencies are printed with 4 decimal places.
- `--min-snp-frequency` and `--min-mnv-frequency` are values from `0` to `1`
  and require `--bam`.
- SNP and MNV frequency filters are independent, so a strong MNV haplotype is
  not removed by a stricter SNP-frequency threshold.
- SNP and MNV read-support and strand-support filters are also independent.
- `--sample all` writes one output set per VCF sample.
- `--keep-original-info` preserves non-get_MNV INFO fields from the input VCF.
