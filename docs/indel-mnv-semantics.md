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

## What get_MNV Adds

- Codon-aware grouping of SNVs into SNP, MNV, and mixed `SNP/MNV` rows.
- iVar insertion/deletion normalization into VCF-compatible anchored alleles.
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

- Input BCF is not accepted directly. Convert first with `bcftools view`.
- Indel left-alignment and full normalization are not performed automatically
  unless `--normalize-alleles` can trim a shared prefix or suffix. For
  cross-caller comparisons, normalize inputs first with a FASTA-aware tool such
  as `bcftools norm -f ref.fa`.
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
way as the BAM (common in homopolymers/tandem repeats). Normalise the input
first with a FASTA-aware tool (`bcftools norm -f ref.fa`) so the allele matches
the read alignment.

## References

- FreeBayes: https://github.com/freebayes/freebayes
- bcftools norm: https://samtools.github.io/bcftools/bcftools.html#norm
- GATK HaplotypeCaller: https://gatk.broadinstitute.org/hc/en-us/articles/30332006386459-HaplotypeCaller
