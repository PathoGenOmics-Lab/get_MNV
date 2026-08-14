# Glossary

The words this documentation uses, in the sense get_MNV uses them. Terms in
bold appear as columns or INFO fields in the output; see
[Output Formats](output-formats.md) for where each one is written.

## The change itself

**SNV** (single-nucleotide variant)
: One reference base replaced by one other base. The output also calls this a
  **SNP** in the `Variant Type` column.

**MNV** (multi-nucleotide variant)
: Two or more substitutions that fall inside one codon and are read together, so
  the amino acid comes from the whole codon rather than from each base alone.
  This is the thing get_MNV exists to find.

**SNP/MNV**
: A row whose codon holds several substitutions, reported with both readings: the
  amino acid each substitution would give on its own, and the one the codon gives
  with all of them. A row typed plain `MNV` has no separate per-base reading to
  report.

**Codon**
: Three consecutive coding bases, read in the transcript's direction, that
  specify one amino acid. On the minus strand the transcript reads from higher
  coordinates down, so a codon's first base has the highest coordinate of the
  three.

**Indel**
: An insertion or a deletion. A VCF writes both anchored on a base that does not
  change, so `T>TGCT` inserts `GCT` after the anchor and `TGCT>T` deletes the
  three bases after it.

**delins**
: A change that deletes and inserts at once, such as `AT>GGG`. Neither a pure
  substitution nor a pure indel.

**Complex indel**
: An allele whose decomposition holds an indel together with at least one
  substitution, so it cannot be described as one simple event.

**Anchor**
: The base a VCF record names to place an indel: the base the inserted sequence
  follows, or the base before the deleted ones. It is not changed by the indel,
  which is why get_MNV asks where a record *acts* rather than where it begins.

**Event class** and **event components**
: How get_MNV read the `REF`/`ALT` pair: the class is the shape of the whole
  allele (`snp`, `mnv`, `insertion`, `deletion`, `delins`, `complex_indel`), and
  the components are the individual changes it decomposes into, written as
  `SNV:110:A>C` or `INS:200:+GG`.

## Where the change falls

**Intergenic**
: Outside every annotated feature. get_MNV writes the placeholder `intergenic` in
  the `Gene` column for these, which is not the name of a gene, and
  `--exclude-intergenic` removes only these rows.

**Frameshift** and **in-frame**
: An indel whose length is not a multiple of three shifts the reading frame of
  everything downstream in that feature; one whose length is a multiple of three
  leaves it intact and adds or removes whole residues.

**Splice donor**, **splice acceptor**, **splice region**
: The first two bases of an intron, its last two, and the bases immediately
  around a junction. The first two are the essential sites and score HIGH impact;
  the surrounding region scores LOW.

**NMD** (nonsense-mediated decay)
: The cell's route for destroying a transcript that terminates early. get_MNV
  predicts it only where a stop is genuinely premature, using the rule that a
  stop more than 50 nucleotides before the last exon junction triggers decay.

**Translation table**
: Which NCBI genetic code to read codons with, chosen by `--translation-table`.
  It decides both the amino acids and which codons end a protein: `TGA` stops
  under table 11, codes tryptophan under table 2 and glycine under table 25.

## What the evidence says

**Haplotype**
: A set of changes carried by the same molecule. A codon's substitutions form a
  haplotype only if the same DNA molecule carries all of them.

**Phasing**
: Placing changes onto molecules. get_MNV never phases on its own: it reads what
  the caller declared, and counts what the reads show.

**Declared phase**
: What the input VCF claims, taken from a `|`-separated `GT` and its `PS` phase
  set. A `/`-separated genotype is unphased and claims nothing.

**Linkage**, **D'** and **r²**
: How far two changes travel together in the reads, beyond what their individual
  frequencies would produce by chance. Two common variants meeting often is not
  evidence of a haplotype; see [Linkage](linkage.md).

**Read support** and **event support**
: How many reads carry a change. Substitutions are counted per base through a
  window cache; an indel is counted from the CIGAR of each read, and its counts
  are written in the `Event Reads` columns.

**Strand bias**
: Support that arrives almost entirely from one strand, which usually points at
  an artefact rather than a variant. Reported as a Fisher exact p-value with
  `--strand-bias-info` and filtered with `--min-strand-bias-p`.

## How the consequence is labelled

**SO term**
: The Sequence Ontology name for the consequence, such as `missense_variant` or
  `splice_donor_variant`. A row can state more than one at once, joined with an
  ampersand, when a base is both coding and in a splice region.

**Impact**
: The severity bucket the SO term falls in: `HIGH`, `MODERATE`, `LOW` or
  `MODIFIER`. When a row states two consequences it takes the more severe of the
  two.

**Grantham distance**
: A number for how chemically different two amino acids are, reported for a
  genuine missense so that a conservative change can be told from a drastic one.
