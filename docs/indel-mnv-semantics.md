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

## Current Limits

- Input BCF is not accepted directly. Convert first with `bcftools view`.
- Indel left-alignment and full normalization are not performed automatically
  unless `--normalize-alleles` can trim a shared prefix or suffix. For
  cross-caller comparisons, normalize inputs first with a FASTA-aware tool such
  as `bcftools norm -f ref.fa`.
- Genotypes, ploidy, phase sets, and genotype likelihoods are preserved only as
  original INFO/FORMAT context when requested; they are not re-estimated.
- Local de novo assembly is out of scope. get_MNV only combines alleles already
  present in the input and, when available, confirmed by BAM read support.
- Local haplotype discovery is bounded to nearby event windows; very large
  haplotype reconstruction should still be handled by a dedicated caller.

## References

- FreeBayes: https://github.com/freebayes/freebayes
- bcftools norm: https://samtools.github.io/bcftools/bcftools.html#norm
- GATK HaplotypeCaller: https://gatk.broadinstitute.org/hc/en-us/articles/30332006386459-HaplotypeCaller
