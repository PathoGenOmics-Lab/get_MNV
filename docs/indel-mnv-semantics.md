# Indel and MNV Semantics

get_MNV is a variant annotator and haplotype summarizer for existing variant
calls. It does not replace a variant caller: the input VCF or iVar TSV is still
responsible for deciding which alleles exist. When a BAM is provided, get_MNV
uses read evidence to count support and to emit exact combined indel/SNV/MNV
haplotypes that are observed on the same reads.

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
- Bounded local haplotype windows that can combine multiple nearby events, such
  as insertion-plus-deletion haplotypes, when exact read support exists.

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

These optional flags change indel-aware behaviour. All default to the historical
behaviour, so existing pipelines are unaffected unless a flag is set.

- `--frameshift-min-freq <0.0-1.0>` (default `0.0`): minimum allele frequency an
  *upstream* indel must reach before it shifts the reading frame of downstream
  SNV/MNV codons (the `(fs)` marker and frameshift change types). The default
  propagates from every indel regardless of frequency. For intra-host / mixed
  populations, raising this avoids relabelling a high-frequency downstream
  substitution as frameshifted because of a low-frequency upstream indel that is
  almost certainly on a different molecule. Frameshift propagation is positional,
  not read-phased, so this is a coarse but useful guard.
- `--indel-anchor-depth` (default off): count indel locus depth (the `EDP`/`EFREQ`
  denominator) from reads that observe the anchor base, instead of only reads
  that fully span the REF allele. Reduces depth under-counting and `EFREQ` bias
  for multi-base deletions, where partially-overlapping reads are otherwise
  dropped from the denominator.
- `--phased-indel-min-reads <N>` (default `1`) and
  `--phased-indel-min-freq <0.0-1.0>` (default `0.0`): minimum BAM support a
  phased indel/complex haplotype row must have to be emitted. Local windows can
  enumerate many overlapping sub-haplotypes; raising these suppresses
  low-confidence rows from dense variant clusters.

When an indel that has read coverage produces zero exact-CIGAR support, get_MNV
emits a warning: this usually means the input indel is not left-aligned the same
way as the BAM (common in homopolymers/tandem repeats). Normalise the input
first with a FASTA-aware tool (`bcftools norm -f ref.fa`) so the allele matches
the read alignment.

## References

- FreeBayes: https://github.com/freebayes/freebayes
- bcftools norm: https://samtools.github.io/bcftools/bcftools.html#norm
- GATK HaplotypeCaller: https://gatk.broadinstitute.org/hc/en-us/articles/30332006386459-HaplotypeCaller
