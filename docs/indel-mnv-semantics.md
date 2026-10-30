# Scope and Compatibility

What get_MNV takes responsibility for when it annotates indels and MNVs, what it
leaves to your variant caller, and where its limits are.

get_MNV is a variant annotator and haplotype summarizer for existing variant
calls. It does not replace a variant caller: the input VCF or iVar TSV is still
responsible for deciding which alleles exist. When a BAM is provided, get_MNV
uses read evidence to count support and to emit exact combined indel/SNV/MNV
haplotypes that are observed on the same reads.

For how an indel is read off the alignments, what each output number counts,
and how local haplotypes are discovered and reported, see
[Indels and Local Haplotypes](indel-haplotypes.md). For whether a haplotype's
variants genuinely travel together, see [Linkage](linkage.md).

## Caller Compatibility

Different callers can represent the same local haplotype in different ways.
For example, FreeBayes is haplotype-based and can emit SNPs, indels, MNPs and
complex events as raw haplotype calls. Its own documentation recommends
decomposing those calls after calling if primitive SNP/indel records are needed.

bcftools `norm` handles the other side of the problem: it can left-align and
normalize indels with a FASTA reference, split multiallelic sites, and atomize
complex variants or MNVs into simpler records. That normalization step is useful
before comparing callsets from different callers.

GATK HaplotypeCaller is a caller, not an annotator. It uses local assembly and
genotyping parameters such as sample ploidy, maximum haplotypes, and MNP merge
distance. get_MNV currently does not perform this de novo haplotype discovery or
genotype likelihood modelling.

## Left-Alignment

An insertion or deletion inside a homopolymer or a tandem repeat can be written
at any of several positions that all describe the same change: deleting one
`GCT` from `GCTGCTGCT` removes the first copy, the second or the third.
Aligners put the gap at the leftmost of them, and get_MNV counts an indel's read
support from reads whose CIGAR carries that exact gap, so an allele written
further right finds no supporting reads at all. A VCF written from a 3'-shifted
HGVS description, or by a caller that places indels its own way, does exactly
that.

`--left-align` moves each indel to its leftmost equivalent position against the
`--fasta` reference before anything else reads it, which is the normalization
`bcftools norm -f ref.fa` performs:

- An insertion or deletion moves left for as long as the result describes the
  same change, and is written with the base before it as its anchor. One that
  reaches the first base of the contig, where there is no base before it, is
  anchored on the base after it, as the VCF specification asks.
- The bases REF and ALT share at either end are trimmed, as
  `--normalize-alleles` trims them, so `--left-align` includes that flag.
- Substitutions (SNVs and MNVs) and complex records (a deletion and an insertion
  in one record, or an indel and a substitution) are trimmed but never moved:
  once trimmed, their two alleles end in different bases, so no other position
  describes them. A complex record is not split into simpler ones, which is
  what `bcftools norm -a` does.
- Left as written: symbolic alleles get_MNV does not write out (`<INV>`,
  `<INS>`, events over 10,000 bases), breakends, `*`, alleles carrying an
  ambiguity code such as `R`, which bcftools refuses as well, and a REF that
  disagrees with the reference, which stops the run as it stops bcftools.
- A `<DEL>` or `<DUP>` that get_MNV writes out from the reference is an explicit
  indel from then on and is left-aligned like any other. bcftools leaves
  symbolic records where they are.
- iVar TSV input is left-aligned too. iVar reads its indels off the alignment,
  so they are usually where the reads have them already.

The protein change and the HGVS descriptions do not depend on where the record
was written, since HGVS places an indel at its 3'-most position whatever the
input says. What moves is `Positions`, `Reference Bases`, `Base Changes`,
`Event Components`, the output VCF record and, with them, the read support.
Left-alignment helps that support only where the aligner left-aligned the gaps,
which most aligners do; when the reads place a gap somewhere else, the warning
about spanning reads without exact support says so.

It is off by default, on the command line and in the desktop app. The one
deliberate difference from bcftools is an ambiguity code in the reference:
bcftools reads `R` as `A`, moves an indel through it and writes that `A` into
REF, which then disagrees with the FASTA, where get_MNV stops the indel at the
`R`. The `tests/leftalign` suite runs both tools on thousands of indels and
complex records placed at random in the repeats of the bundled genome, spelled
in several ways, and requires the same record from each.

## What get_MNV Adds

- Codon-aware grouping of SNVs into SNP, MNV, and mixed `SNP/MNV` rows.
- iVar insertion/deletion normalization into VCF-compatible anchored alleles.
- Left-alignment of indels against the reference with `--left-align`, as
  `bcftools norm -f` does (see [Left-Alignment](#left-alignment)).
- Indel and complex-allele protein effects using the selected feature, strand,
  GFF phase, and transcript protein offset.
- BAM-derived support for individual SNP/MNV observations and exact indel or
  complex-indel events.
- Extra `complex_indel` rows only when the full combined indel/SNV/MNV haplotype
  is observed in BAM reads.
- Local haplotypes read off the molecules. Nearby events are grouped into a
  window by proximity, but proximity only proposes: the reads spanning that
  window decide which combinations are reported, and each is reported with the
  count of molecules carrying it. A combination no read carries is not emitted,
  including one that is a subset of a combination reads do carry, since a
  molecule with A, B and C is not evidence for a molecule with only A and B.
  Two combinations that genuinely coexist are both reported, with their own
  counts. There is no cap on how many variants a window may hold. Each
  candidate is judged on its own reference span, and for an insertion that span
  reaches one base past the anchor, since an insertion lives between two
  reference bases and a read ending on the anchor has seen nothing of the
  junction.

## Boundary Rules

VCF insertions are interbase events anchored on the previous reference base.
An insertion is treated as overlapping a feature only when the inserted sequence
falls between two reference bases inside that feature. An insertion anchored at
the final feature base is outside that feature.

VCF deletions are anchored, but their biological effect is the deleted reference
span. A deletion anchored just before a CDS/gene still affects that feature if
the deleted bases overlap it.

A deletion that runs from the bases before the start codon into the gene is read
at its 3'-most place along the transcript, within the first exon, as HGVS places
it. A deletion in a repeat can be written at several places, and left-aligned,
as `--left-align` and `bcftools norm` write it, deleting the second ATG of
`ATGATG` begins on the base before the gene; read where it is written, it would
keep only two coding bases and come out a frameshift. Placed 3'-most it lies
inside the gene and reads `Met2del`, with the same codons and `c.` as the same
deletion written there. One that still reaches the base before the start codon
in that place loses the start. The row keeps the record as written.

When an MNV row overlaps an indel, get_MNV keeps the MNV row as positional
context and marks its amino-acid effect as `Unknown` with
`Change Type = Indel overlap`. If the BAM supports the full combined event,
get_MNV emits a separate exact `complex_indel` row.

Exact support for complex indels is CIGAR-aware. A read must produce the same
local ALT sequence and contain the expected insertion/deletion components. This
matters for net-neutral combinations, where an insertion plus a deletion can
produce the same sequence as a simple MNV under another alignment.

## Eukaryotic CDS Notes

For eukaryotic GFF/GTF annotations, use `--gff-features CDS` so get_MNV can use
GFF phase and transcript identifiers. When CDS rows carry `transcript_id` or
`Parent`, get_MNV builds a spliced CDS model for each transcript and evaluates
codon grouping, MNV amino-acid effects, and indel frameshift context on that
full coding sequence.

This means codons split across exon junctions can be annotated as one
transcript-level MNV, and downstream SNPs are marked as frameshifted only when
the net upstream coding indel shift remains out of frame. CDS rows without a
usable transcript identifier keep the older per-feature behavior.

## Current Limits

- Indels are left-aligned only when `--left-align` asks for it; without it an
  indel is matched against the reads where the input wrote it. Complex records
  are trimmed, not split into simpler records, and records are not joined
  across lines; `bcftools norm` does both when a comparison needs them.
- Genotypes, ploidy, phase sets, and genotype likelihoods are not re-estimated.
  For diploid/polyploid eukaryotic data, interpret unphased heterozygous
  exon-level MNVs cautiously unless the VCF or BAM evidence confirms that the
  alleles are on the same haplotype.
- Local de novo assembly is out of scope. get_MNV only combines alleles already
  present in the input and, when available, confirmed by BAM read support.
- Local haplotype discovery is bounded to nearby event windows; very large
  haplotype reconstruction should still be handled by a dedicated caller.
- Full local assembly remains out of scope: get_MNV reannotates alleles present
  in the input VCF/iVar file rather than discovering new candidate variants.

## Tuning Knobs

These flags change indel-aware behaviour. Two of them no longer default to the
tool's original behaviour, because that behaviour was wrong more often than it
was right: the frame shift propagates only from a majority upstream indel, and
indel locus depth is counted from the anchor base.

- `--left-align` (default off): move each indel to its leftmost equivalent
  position against the reference before annotating it, as `bcftools norm -f`
  does, so it matches the gap an aligner reports. See
  [Left-Alignment](#left-alignment).

- `--frameshift-min-freq <0.0-1.0>` (default `0.5`): minimum allele frequency an
  *upstream* indel must reach before it shifts the reading frame of downstream
  SNV/MNV codons (the `(fs)` marker and frameshift change types). The default
  propagates only from a majority upstream indel; set `0.0` to propagate from
  every one. For intra-host / mixed populations, this avoids relabelling a
  high-frequency downstream substitution as frameshifted because of a
  low-frequency upstream indel that is almost certainly on a different molecule.

  The frequency compared against the threshold is the one **the reads give**
  when a BAM is provided, which is the same number reported as `EFREQ`. Only
  without a BAM does the caller's declared `AF` stand in for it. This matters
  because many callers write no `AF` at all: while the gate consulted the
  declared value alone, such an input passed every indel at every threshold.
  When a BAM is provided, propagation is additionally **read-phased**: for each
  upstream indel and a downstream SNV within a read's reach, reads spanning both
  loci are inspected, and the frame shift is *not* propagated to that codon when
  the reads carrying the SNV overwhelmingly lack the indel (the two are in trans,
  on different molecules). This is a conservative, suppression-only refinement
  (it never adds propagation the frequency gate would not), and only applies where
  reads span both loci; beyond a read's reach the frequency gate still governs.
- `--legacy-indel-depth` (default off): restrict indel locus depth (the
  `EDP`/`EFREQ` denominator) to reads that fully span the REF allele. By default
  it is counted from every read observing the anchor base, which avoids the
  depth under-counting and `EFREQ` bias a multi-base deletion otherwise suffers,
  where partially-overlapping reads are dropped from the denominator. This flag
  restores the older, narrower denominator.
- `--count-mates-separately` (default off): count the two mates of a paired-end
  fragment as two observations instead of one molecule. By default they are one.
  A fragment is a single DNA molecule sequenced from both ends, so counting the
  mates separately double-counts wherever they overlap, and treats a variant on
  each mate as unrelated when it is in fact proof the two travel together: a
  codon split by an intron can be unanswerable read by read and settled by the
  pair. Where the mates overlap and disagree about a base, one of them is wrong
  and there is no way to tell which, so the molecule is treated as not having
  observed that position. Single-end data is unaffected either way.

  Local haplotype *discovery* is molecule-level too: each candidate variant is
  judged on its own reference span, so a fragment can answer for a variant its
  first mate covers and another its second does. Exact event *counting* is not:
  a read must reconstruct the whole compound allele in one contiguous
  observation, which no single mate of a split haplotype can do. A local
  haplotype whose variants fall on different mates is therefore found but
  cannot be supported, and is dropped by the support threshold rather than
  reported with no evidence.
- `--phased-indel-min-reads <N>` (default `2`) and
  `--phased-indel-min-freq <0.0-1.0>` (default `0.0`): minimum BAM support a
  phased indel/complex haplotype row must have to be emitted. A dense window in
  intra-host data can hold several genuine combinations at low frequency;
  raising these keeps only the well-supported ones.

When an indel that has read coverage produces zero exact-CIGAR support, get_MNV
emits a warning: this usually means the input indel is not left-aligned the same
way as the BAM (common in homopolymers/tandem repeats). Rerun with
`--left-align`, or normalize the input first with `bcftools norm -f ref.fa`, so
the allele matches the read alignment. If the warning stays with
`--left-align`, the reads place the gap elsewhere or do not carry the indel.

## References

- FreeBayes: https://github.com/freebayes/freebayes
- bcftools norm: https://samtools.github.io/bcftools/bcftools.html#norm
- GATK HaplotypeCaller: https://gatk.broadinstitute.org/hc/en-us/articles/30332006386459-HaplotypeCaller
