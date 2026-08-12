# Indels and Local Haplotypes

How get_MNV reads an indel off the alignments, what it counts, and what each
number in the output means. For the caller-compatibility and boundary rules, see
[Indel & MNV Semantics](indel-mnv-semantics.md); this page is about the reading
and the counting.

## What an indel is, positionally

A VCF indel is **anchored**. `POS=100 REF=A ALT=AGG` does not change position
100: the `A` is padding, and the `GG` is inserted between 100 and 101. A
deletion `POS=100 REF=AGG ALT=A` keeps 100 and removes 101 and 102.

That has a consequence get_MNV has to respect everywhere: an insertion does not
occupy a reference base. It lives *between* two of them.

## What a read has to see before it may speak

A read is entitled to an opinion about an allele only if it observed enough of
the reference to have one. That span is not always the REF allele's own span:

| Allele | REF span | Span a read must observe | Why |
|---|---|---|---|
| Substitution | the base | the base | Nothing lies outside it. |
| Insertion | the anchor | anchor **+ 1** | A read that stops on the anchor has covered the whole REF span and seen nothing of the junction. Without the extra base it reads as "the reference is here", which is evidence the insertion is *absent* that the read never had. |
| Deletion | anchor and deleted bases | that span **+ 1** | Inside its own span, a read deleting more than this allele looks identical to one deleting exactly it. Where the deletion ends can only be seen from outside. |

One function answers this for every path that observes an allele
(`observation_ref_len`). Before it existed, the exact counting, the haplotype
discovery and the indel-versus-substitution phasing each decided it separately,
and disagreed. The worst combination cancelled real support: a mate ending on an
insertion's anchor was read as "reference here" and stripped the support its
partner mate had provided.

## Exact event support

`Event Reads` counts molecules that **reproduce the allele exactly** over its
reference span, CIGAR and all. A read must produce the same local sequence and
carry the expected insertion and deletion components. This matters for
net-neutral combinations, where an insertion plus a deletion can produce the
same sequence as a plain substitution under a different alignment.

Two rules keep that count honest:

- A read whose deletion **runs past** the allele's span carries a longer,
  different deletion, and is not support for this one.
- A molecule whose two mates **read different alleles** over the locus supports
  nothing. One of them is wrong and there is no telling which, which is the same
  view the substitution counting takes of a contradicted overlap. Without it, a
  single mate carrying an alignment artefact was enough to credit the whole
  molecule, and mate overlap is exactly where indel realignment artefacts appear.

`Event Depth` is the locus depth: by default every molecule observing the anchor
base at sufficient quality. `--legacy-indel-depth` restricts it to molecules
that fully span the REF allele, which under-counts depth for a multi-base
deletion and biases `Event Frequency` upward.

## Local haplotypes

When several variants fall close together, the question is no longer what each
one does but **which of them share molecules**. get_MNV answers it by reading
the molecules, not by enumerating possibilities.

### Discovery: what the molecules show

Nearby variants are grouped into a window by proximity, but proximity only
*proposes*. Each candidate variant is judged on its own reference span, so a
paired-end fragment can answer for a variant its first mate covers and another
its second does. A molecule that did not observe every variant of the window is
set aside as partial: it carries what it carries, but placing it in a
combination would assert it lacks the variants nobody looked at.

The molecules that did observe the whole window are grouped by the combination
they show. That is the output of discovery, and it is a **joint distribution**:
how many molecules are each combination.

Two rules follow from reading rather than enumerating:

- A combination no molecule carries is not emitted, **including a subset of one
  that is carried**. A molecule with A, B and C is not evidence for a molecule
  with only A and B.
- Two combinations that genuinely coexist both appear, each with its own count.
  An enumerate-and-test scheme cannot tell that apart from one molecule carrying
  everything.

### Counting: how many molecules are this haplotype

A row's support is the number of molecules whose combination is **exactly** that
one, taken from discovery. The exact-allele count cannot answer this: it matches
over the combination's own reference span, so a molecule carrying that
combination *and something else outside the span* matches too. It keeps the job
it can do, which is confirming the allele is reproducible at all, and supplying
the window depth.

The difference is not academic. With twelve molecules carrying three variants
and eight carrying only the first two:

| Event Components | Event Reads | Event Depth | Event Frequency |
|---|---|---|---|
| `SNV:28:G>T \| INS:29:+GCT \| SNV:30:T>A` | 12 | 20 | 0.6000 |
| `SNV:28:G>T \| INS:29:+GCT` | 8 | 20 | 0.4000 |

Counted by allele match, the second row read `20` at `1.0000`: the whole
species' support, for a haplotype eight molecules carry. The counts now sum to
the depth, which is the coherent picture of a window holding two species.

### What is emitted

- The input's own records, annotated as usual.
- One `complex_indel` row per combination the molecules show, when it holds more
  than one variant and at least one is an indel.
- Nothing else. `--phased-indel-min-reads` (default `2`) and
  `--phased-indel-min-freq` (default `0.0`) set the floor: one read is not
  evidence of a haplotype, since a single sequencing error at a called position
  mints a combination of its own.

A window that shows more than 64 distinct combinations reports the best
supported and says in the log what it set aside, rather than truncating quietly.

## Reading an indel row

For the example above:

```text
Positions          28
Variant Type       INDEL
Change Type        In-frame Indel
AA Changes         Ala10delinsSerLeu
Event Class        complex_indel
Event Components   SNV:28:G>T | INS:29:+GCT | SNV:30:T>A
Event Reads        12
Event Depth        20
Event Frequency    0.6000
SO Term            inframe_insertion
HGVS g.            g.28_30delinsTCGCTA
```

`Positions` is the compound allele's anchor, not one position per component;
`Event Components` is where the parts are listed. The amino-acid change is the
effect of the whole compound allele, so a haplotype of three variants gets one
protein consequence rather than three.

A codon-level MNV row overlapping an indel is kept as positional context with
`Change Type = Indel overlap` and `Unknown` for the residue, because the codon's
protein effect is not defined independently of the indel. The `complex_indel`
row beside it is where the real consequence is.

## Limits worth knowing

- **Exact counting needs one contiguous observation.** Discovery works at the
  level of the molecule, so a haplotype whose variants fall on different mates
  is found; the allele cannot be reproduced by either mate alone, so the row
  cannot be supported and the threshold drops it.
- **Left-alignment is the input's job.** get_MNV matches the allele as given.
  In a homopolymer or tandem repeat, an input indel placed differently from the
  BAM's alignment yields zero exact support; a warning says so, and
  `bcftools norm -f ref.fa` fixes it.
- **An insertion anchored at a feature's last base is outside that feature**,
  because the inserted sequence falls beyond it. On a minus-strand gene that
  same coordinate is the transcript's *first* base, so the insertion lands in
  the 5' UTR. Same rule, opposite reason.
- **No de novo assembly.** get_MNV combines alleles the input already declares.
- **Symbolic alleles** (`<DEL>`, `<DUP>`) have no sequence for a read to
  reproduce, so read evidence cannot speak about them at all.

## How this is checked

- [`tests/scenarios/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/tests/scenarios)
  builds alignments by hand for each rule above, including the ones that only
  bite at an edge: a read ending on an insertion anchor, mates disagreeing about
  an indel, a deletion running past the declared one, and a haplotype nested
  inside another.
- [`tests/truth/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/tests/truth)
  sweeps a deletion and an insertion of 1 to 3 bases at **every position of
  every exon** of four genes, checking the protein consequence against an
  independent re-translation.
- [`tests/pileup/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/main/tests/pileup)
  checks depth and strand against `bcftools mpileup`, which nobody here wrote.
