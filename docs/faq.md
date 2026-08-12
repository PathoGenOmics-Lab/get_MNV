# FAQ

## Is get_MNV a variant caller?

No. get_MNV **annotates and summarizes existing variant calls** (from a VCF or
iVar TSV) against a reference. It does not call variants or assemble reads — it
reinterprets the alleles already present in your input, with codon awareness.

## What inputs do I need?

A reference FASTA, a variant file (`--vcf` or `--tsv`), and a gene annotation
(`--gff` or `--genes`). A BAM (`--bam`) is optional and only needed for read
support and the read-based filters. See [Input Formats](input-formats.md).

## When should I use `--gff` versus `--genes`?

Use `--genes` with a simple four-column TSV (`name`, `start`, `end`, `strand`)
for single-exon / prokaryotic genes. Use `--gff` with a GFF/GFF3 file when you
need spliced transcripts; add `--gff-features CDS` so codons are built from the
joined CDS segments (including codons that span exon junctions).

## Which `--translation-table` should I use?

The default is `11` (bacterial), correct for organisms like *M. tuberculosis*.
For the standard nuclear code use `1`; other supported tables are 2, 3, 4, 5, 6,
12 and 25. See the [CLI Reference](cli-reference.md).

## Why did two of my SNPs become one MNV?

When two or more SNVs fall in the same codon, the amino acid effect depends on
the **combined** codon, not the individual changes. get_MNV groups them and
reports a `MNV` (or `SNP/MNV`) row with the real codon and amino acid. See
[Indel & MNV Semantics](indel-mnv-semantics.md).

## Do I need a BAM?

Only if you want read support. Without `--bam` you still get full codon/MNV and
indel annotation; with it you also get per-event read counts, frequencies and
strand metrics, and you can apply read-based filters.

## Why are the frequencies different from my input VCF?

Read-based filters and the reported support are recalculated from `--bam`, not
taken from the original `OFREQ`/`ODP`. This reflects what the reads actually
show for each SNV, MNV haplotype and indel event.

## Can I use a BCF file as input?

Not directly — convert it to VCF first, for example with
`bcftools view input.bcf -O v -o input.vcf`. get_MNV can *write* BCF output via
`--bcf`. See [Troubleshooting](troubleshooting.md).

## My VCF has multiallelic records — what happens?

By default get_MNV stops so you can decide how to handle them. Pass
`--split-multiallelic` to split each record into independent ALT alleles, or
pre-split with `bcftools norm -m -`.

## Are indels supported?

Yes. get_MNV decomposes `REF/ALT` alleles into SNV, MNV, insertion, deletion,
delins and complex components, reports their protein effect when they overlap a
coding feature, and (with `--bam`) counts exact indel-event support. It does not
left-align or fully normalize indels for you — normalize cross-caller inputs
with `bcftools norm -f ref.fa` (or `--normalize-alleles` for simple trimming).

## What about intergenic variants?

They are included by default and labelled `intergenic`. Use
`--exclude-intergenic` to drop variants outside annotated genes.

## Where does the output go?

By default get_MNV writes `<input_name>.MNV.tsv` into the current working directory, not next to the input. Use
`--convert` for VCF, `--both` for both, and the related output options in the
[CLI Reference](cli-reference.md).

## The macOS app won't open ("unidentified developer")

The desktop app is not signed with an Apple Developer certificate. On first
launch, right-click the app → **Open** → click **Open** in the dialog. See the
[Desktop GUI](gui.md) guide.
