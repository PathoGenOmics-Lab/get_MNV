# Changelog

All notable changes to this project are documented in this file.

## [1.1.5] - unreleased

### Added

- **Interactive HTML report.** New `--report <FILE.html>` writes a single self-contained page for exploring the called variants: it embeds its data, loads nothing external, and opens offline.
  - Contents: a headline count with tiles, the Sequence Ontology consequence distribution, an alignment-style variant matrix, a read-backed haplotype panel, and a virtualised, sortable table with a filter on every column.
  - Per-variant detail panel: HGVS `g.`/`c.`/`p.`, Grantham, DBS class, NMD, codons, event components, read support.
  - The matrix is a genome browser: samples down the side, genomic coordinates across the top, coloured by alternate base, with zoom and pan, density lanes, a gene track and a region box taking `start-end`, `contig:start-end`, a coordinate or a gene name.
  - It orders samples by shared profile, marks positions called together on the same reads, and shows only an ALT call or "not called".
  - The haplotype panel ranks codon-level MNVs and phased complex indels by how many samples share them, with their phasing support; long-range phasing is never inferred.
  - Filters drive the charts as well as the table, the selection exports as TSV, and cohort-sized reports stay in the low megabytes.
  - `--sample all` puts every sample in one report; `--report-from a.MNV.tsv b.MNV.tsv ...` aggregates existing TSV outputs without re-running the pipeline.
- **Sequence Ontology consequence terms and impact levels.** New `SO Term` and `Impact` TSV columns (and `SO` / `IMPACT` VCF INFO fields) annotate every variant with a standard SO term (`missense_variant`, `synonymous_variant`, `stop_gained`, `stop_lost`, `start_lost`, `frameshift_variant`, `inframe_insertion`/`inframe_deletion`, `intergenic_variant`, `coding_sequence_variant`) and a `HIGH`/`MODERATE`/`LOW`/`MODIFIER` impact, following SnpEff/VEP conventions.
- **Grantham distance for missense changes.** A new `Grantham` TSV column (and `GD` VCF INFO field) reports the Grantham (1974) distance of the combined amino-acid substitution together with its conservation category (e.g. `177 (radical)`), for genuine missense changes only.
- **MNV-vs-SNV consequence shift.** A new `MNV Consequence Shift` TSV column (and `MNVSHIFT` VCF INFO field) flags a combined MNV as `MNV-gained` (more severe than any single SNV alone, e.g. two synonymous SNVs producing a non-synonymous residue), `MNV-masked` (a nonsense SNV rescued by its neighbour), or `Concordant`.
- **Doublet base substitution (DBS) class.** A new `DBS Class` TSV column (and `DBS` VCF INFO field) reports the COSMIC-style doublet class for an MNV of two adjacent single-base substitutions, e.g. `CC>TT`, reverse-complement collapsed so both strands map to one class, so mutational signatures can be tallied from codon-level MNV calls.
- **MNV phasing (linkage) support.** With `--bam`, a new `MNV Phasing Support` TSV column (and `MNVPS` VCF INFO field) reports the fraction of reads covering every codon position and carrying the least-supported constituent SNV that also carry the full haplotype, `0.0000` when none do and `-` when nothing could answer the question.
  - A companion `MNV Phasing Reads` column (and `MNVPR` INFO field) reports how many reads the ratio rests on.
- **Phase declared by the input caller is now read and reported.** Both VCF parsers read `GT` and `PS` for the selected sample, and a new `Declared Phase` TSV column (and `DPHASE` VCF INFO field) reports what the input claimed, as `cis:12345`.
  - Only a `|`-separated genotype counts, only records sharing a phase set are compared, and a record is matched to a row by both coordinate and base.
  - The claim sits beside the read evidence, never merged into it; `|contradicted-by-reads` is appended for a declared cis that not one spanning read carries whole, or a declared trans that every informative read carries.
- **Codon-level linkage disequilibrium.** With `--bam`, new `Haplotype LD` and `Haplotype LD p` TSV columns (and `LD` / `LDP` VCF INFO fields) report `D'` between a codon's substitutions and the two-tailed Fisher exact p-value behind it, over the molecules observing every position.
  - `+1` is one haplotype, `~0` is two variants that merely share a codon, and `-1` is mutual exclusion, in a haploid population the read-level signature of a mixed infection.
  - Local indel haplotypes are covered too, not only codon MNVs; with three or more variants the weakest pair decides.
  - The value is absent, not zero, when no molecule spans the codon or when an allele is on every molecule or none.
- **Read counting checked against `bcftools mpileup`** (`tests/pileup/`, new CI job): the same BAM goes to both tools and depth, allele support and the two strand arms are compared against `FORMAT/DP`, `AD`, `ADF` and `ADR`, with an explained difference written down where it arises and an unexplained one failing the run.
- **Read-level phasing validation** (`tests/phasing/`, new CI job): the suite builds alignments molecule by molecule and checks that get_MNV recovers the mixture it was handed across five codon geometries and all 21 mixtures from fully trans to fully cis, including the cases where the honest answer is not a number, checking `D'` and the exact test.
- **Nonsense-mediated decay (NMD) prediction.** A new `NMD Prediction` TSV column (and `NMD` VCF INFO field) applies the 50-nucleotide rule to premature stops in multi-exon transcripts: `NMD-triggering` more than 50 nt upstream of the last exon-exon junction, `NMD-escaping` in the last exon or within 50 nt of that junction; single-exon transcripts are left unannotated.
- **HGVS genomic and coding descriptors.** New `HGVS g.` and `HGVS c.` TSV columns (and `HGVSG` / `HGVSC` VCF INFO fields) join the protein `p.` notation: `g.` covers every variant, including the allele bracket `g.[28G>T;30T>A]`, and `c.` covers substitutions numbered from the CDS start with coding-strand bases on both strands; descriptors are not 3'-shifted and coding indel descriptors are not yet emitted.
- **Splice-site consequence terms.** Variants near an internal exon-exon junction of a multi-exon transcript now get `splice_donor_variant` / `splice_acceptor_variant` (`HIGH`) for the two essential intronic bases at each intron end and `splice_region_variant` (`LOW`) for the exon's last 3 bases or the intron's 3rd-8th bases, in `SO Term` / `Impact` and the `SO` / `IMPACT` VCF INFO fields.
  - Intronic splice variants, previously reported as intergenic, are now associated with their gene.
  - An exonic coding change near a junction is combined, e.g. `missense_variant&splice_region_variant`, keeping the more severe impact; single-exon transcripts and ribosomal-slippage CDS joins are unaffected.
- **Synthetic ground-truth validation** (`tests/truth/`, new CI job): the suite builds a four-gene genome and derives the expected answer by a different route from get_MNV's, covering every substitution, every pair and the triple in a codon plus 1 to 3 base indels at every exon position, each checked on gene, amino-acid change, change type, both codons and the HGVS `c.` descriptor.
- **Differential annotation suite** (`tests/differential/`, new CI job): get_MNV is compared against `bcftools csq` on small committed fixtures and against SnpEff 4.1l through the `ANN=` fields in the bundled *M. tuberculosis* VCF, with every difference matched against `baseline.tsv` and CI failing only on an unexplained one; `GET_MNV_BIN` points the suite at another build.
- Added regression coverage for phased MNV-plus-indel haplotypes, verifying that codon MNV rows overlapping an indel are flagged as `Indel overlap` while BAM-supported combined events are emitted as exact `complex_indel` rows.
- Added an indel/MNV semantics note documenting caller compatibility, boundary rules, current limits, and how exact complex haplotypes are represented.
- Added transcript-aware regression coverage for exon-junction MNV codons and restored-frame indel contexts in multi-exon CDS models.
- **A run now warns when an MNV threshold cannot remove anything.** A codon-level row survives if either its individual SNVs or its haplotype clears its bar, so with the SNP thresholds at their default of `0` a flag like `--mnv 5` removes nothing, and the run now says so, naming the flags to raise; no filtering behaviour changed.
- **Pipeline-level invariants over generated inputs** (`tests/property_pipeline.rs`, run by the existing Rust CI job): properties assert what must hold for any input, padding invariance for substitutions and for indels, record-order independence, one gene spelled as TSV or GFF annotating the same, cohort and multi-gene consistency, iVar-versus-VCF agreement, and the report showing the rows the TSV of the same run holds.
- **Behaviour checks for the JavaScript and TypeScript the project ships** (`tests/js/run.mjs`, new CI job): node exercises the page's filter, search, sort, matrix and canvas drawing plus the desktop form's preset merge, with no new dependencies.
- The scenario suite now covers the splice terms and `--exclude-intergenic`, and every scenario writes the VCF too and checks that its two outputs agree on every base.
- A glossary page in English and Spanish defining the documentation's vocabulary, plus a text version of the home page's worked example and a table of what each step of a run reads and produces.

### Changed

- Mutation testing now covers the read counting: a test constructor feeds observations without an indexed BAM, and fourteen new unit tests pin molecules, per-position strand, the spanning rule, the quality floor and the linkage arithmetic.
- The mutation-testing workflow no longer runs on a schedule. It is `workflow_dispatch` only, like every other manual job here.
- A local haplotype row is now counted by the molecules that are exactly it, not by every molecule matching its allele over its own reference span; discovery supplies the row's support and strand split.
- An audit of the haplotype and read-counting layer fixed 17 defects. Grouped by what went wrong:
  - **Strand was tracked per molecule instead of per position.** It now travels with each observation, so `--min-snp-strand`, `--min-mnv-strand` and the strand-bias test can act on paired-end data.
  - **A whole run aborted on a multi-base substitution outside a gene.** `AC>GT` in an intergenic region with `--bam` failed with "Missing SNP read counts"; such rows now take the placeholder layout in TSV and VCF.
  - **Three copies of "how much must a read see".** One `observation_ref_len` now answers it for exact counting, local-haplotype discovery and indel-versus-SNV phasing, and a read whose deletion runs past the allele is no longer exact support.
  - **Mate disagreement was resolved by OR in two places.** Local-haplotype discovery and the trans test now agree with the substitution counting: a contradicted observation is no observation.
  - **Symbolic SV alleles were always judged trans.** `<DEL>` has no sequence for a read to reproduce, so such pairs are now left out rather than answered wrongly.
  - **Read phasing was keyed on position without the allele.** At a multi-allelic site the last ALT processed overwrote the other's evidence and then decided its codon.
  - **A declared-phase claim spoke for alleles nobody had phased.** Every position must now carry a claim in one phase set, and a no-call (`.`) slot is no longer read as "not this ALT".
  - **Unnamed reads collapsed into one molecule.** With SAM `QNAME` `*` the fallback key was only the segment flags and the start position, so every unnamed read starting there became one molecule.
  - **The exact test was quadratic.** `ln(n!)` is now tabulated once per call instead of recomputed inside the Fisher loop, which also speeds up the strand-bias test.
- Added mutation testing (`cargo mutants`, `.cargo/mutants.toml`, a manually-run advisory CI job) that alters the source and reports whether the suite notices; five new tests now cover the spliced-transcript indel and frameshift-propagation path.
- A paired-end fragment now counts as one molecule instead of two reads in read counting, local haplotype discovery, event counts and the trans test; single-end data is unaffected and `--count-mates-separately` restores the per-record behaviour.
- The read-level phasing behind frameshift suppression is now reported: a new `Frameshift Phasing` TSV column and `FSPH` VCF INFO field carry the verdict, indel position and reads in cis, as `trans:1234:0/18`.
- `--phased-indel-min-reads` now defaults to `2` (was `1`). One read carrying a sequencing error minted a combination of its own and was enough to publish it as a haplotype row; two reads is the least that means more than one molecule. Set `1` to emit every combination any read shows.
- Local haplotypes are now read off the spanning molecules instead of enumerated and queried, each reported with its molecule count; the variant cap is gone and a bound of 64 reported combinations per window remains.
- Where a transcript's introns are is now derived once, by `Gene::introns`, instead of being recomputed by splice-site classification, the intronic-position test and the NMD junction independently.
- The GFF CDS model now rejects and reports annotation it cannot assemble: CDS rows that disagree about the strand stay as per-feature annotations, repeated CDS rows are dropped, both with a warning.
- Added property-based tests (`proptest`): the headline invariant is strand symmetry, a variant on a plus-strand gene must match the mirrored variant on the reverse-complemented genome in consequence, change type and codons.
- Gene fixtures are now built from GFF text through the real parser (`test_support`): `spliced_gene` and `joined_gene` refuse mismatched segments, 40 of 46 `Gene` struct literals are gone, and output is unchanged.
- Consequence classification now lives in one module (`variants::consequence`) instead of two implementations, so indels get the initiator rule that previously reported a start-codon deletion at `MODERATE` impact; output is unchanged.
- `--chrom` now restricts FASTA loading: when a single contig is requested, only that contig is read and IUPAC-validated instead of the whole genome, cutting peak memory for large (e.g. eukaryotic) references. A missing requested contig now fails with a clear "not found in FASTA" error.
- SNVs/MNVs that alter the initiator Met (protein position 1) are now reported with Change Type `Start lost` instead of `Non-synonymous`, matching standard annotators. The reported amino-acid change was already correct; only the classification label changes. `Met1` → stop is still `Stop gained`.
- Added a tuned `[profile.release]` (thin LTO, single codegen unit) and wrapped plain (non-BGZF) VCF output in a buffered writer with an explicit flush, for faster production builds and record emission.
- Pinned `sha2` to the stable `0.10` line (previously a `0.11` pre-release), removing a duplicate hashing dependency stack from the resolved graph.
- `get_mnv_variants_for_gene` and `get_mnv_variants_for_transcript` now build a list of mutually-exclusive codon interpretations per codon start. Multi-allelic positions expand the interpretation set as a Cartesian product, deduplicated by `(position, alt)` keys, so a codon that contains N independent alts emits N output rows.
- Bumped project, GUI, citation, README, and frontend metadata to version 1.1.5.
- Updated the Tauri desktop dependency set to the 2.11 patch line, including `tauri` 2.11.2 and `tauri-plugin-dialog` 2.7.1.
- Updated frontend lockfile dependencies, including `postcss` 8.5.10.
- Added allele-level event decomposition for `snp`, `mnv`, `insertion`, `deletion`, `delins`, `complex_indel`, and symbolic alleles, so length-changing events that also contain SNV/MNV components become a single local haplotype event.
- Phased local haplotype discovery now covers bounded multi-event windows rather than just indel-plus-SNV pairs, so combinations such as insertion-plus-deletion are emitted as exact `complex_indel` rows.
- iVar TSV inputs now keep insertion and deletion rows by converting `+SEQ` and `-SEQ` alleles into VCF-compatible anchored `REF/ALT` alleles using the FASTA reference.
- TSV and VCF outputs now include canonical event class/component annotations plus exact BAM-derived event support metrics for indel/complex alleles.
- The desktop BAM viewer now renders insertion-aware interbase columns, showing inserted bases between reference positions with matching coverage, ruler, reference, and read-pileup alignment instead of hiding them inside the anchor base.
- With a BAM, nearby SNV/MNV rows that phase with an indel on the same reads are now emitted as an additional exact `complex_indel` haplotype row reporting the combined `REF/ALT`, protein effect, and event support, while the original rows are preserved.
- GFF/GTF `CDS` rows with `transcript_id` or `Parent` are now collapsed into spliced transcript CDS models, so codon grouping, MNV amino-acid effects, and indel frameshift context are evaluated against the full coding sequence instead of isolated exon rows.
- `--frameshift-min-freq` now defaults to `0.5` instead of `0.0`: an upstream indel only shifts the frame of downstream SNV/MNV codons when it is the consensus (majority) allele. Indels without a known frequency still propagate. Pass `--frameshift-min-freq 0.0` for the previous behaviour.
- Frameshift propagation is now **read-phased** when a BAM is provided: reads spanning both an upstream indel and a downstream SNV are inspected, and the shift is not propagated to that codon when the SNV-carrying reads overwhelmingly lack the indel. Suppression only, on top of the frequency gate, which still governs beyond a read's reach.
- Indel locus depth (the `EDP`/`EFREQ` denominator) is now counted from reads observing the anchor base by default, removing depth under-counting and `EFREQ` bias for multi-base deletions. Pass `--legacy-indel-depth` to restrict the denominator to reads that fully span the REF allele.
- **An insertion written against the base after it is re-anchored on entry.** Such records are rewritten once when the contig is loaded, the run reports how many it moved, and their POS/REF/ALT is the re-anchored form.
- **The frameshift gate now weighs the frequency the reads give an indel, not the one the caller declared.** With `--bam`, `--frameshift-min-freq` uses the frequency get_MNV counts and publishes as `EFREQ`, falling back to the declared `AF` when reads were not consulted, and an indel with no frequency still propagates.
- **Scientific**: in a multi-sample VCF, each sample is now annotated only for the alleles its own `GT` carries, with `1/1` taking the first ALT, `1/2` both, and an absent or no-call `./.` genotype keeping the allele.
  - Skipped alleles are logged, and Troubleshooting explains why a run can report fewer rows than the VCF has records.
- A sample carrying nothing now writes an empty output file instead of aborting the whole cohort, which previously stopped a `--sample all` run at the first clean sample and left every later sample unwritten. A file with no records at all is still an error.
- A `--sample all` run whose samples would share an output file name now stops and names them both instead of letting one overwrite the other. Each sample writes to a file named after it with file-name-illegal characters replaced, and a caller-supplied output prefix now keeps the per-sample suffix.
- The interactive report's masthead now carries the wordmark and an icon on every action, including `Theme` and `Export filtered TSV`; the inlined mark follows light, dark and print, and every report is about 90 KB smaller.
- The documentation is now a site in English and Spanish at <https://pathogenomics-lab.github.io/get_MNV/>, and the README no longer keeps its own copy of the option table, the output columns and the limitations.

### Fixed

- **Scientific**: a read that skips an intron (CIGAR `N`) no longer counts towards `Event Reads` and `Event Frequency` as carrying a deletion of those bases; `D` still means the base is absent from the molecule.
- **Scientific**: the codon columns (`Reference Codon`, `SNP Codon`, `MNV Codon`) are now reported in the transcript's own orientation on every annotation path. The genomic-coordinate path emitted genomic-strand bases for minus-strand genes, so the same gene and variant gave different codons for a TSV/single-CDS file than for a multi-exon GFF. Amino-acid calls, change types and impacts are unaffected.
- **Scientific**: the TSV annotation format gained an optional sixth column for the feature biotype, so a non-coding feature is reported as `non_coding_transcript_exon_variant` (`MODIFIER`) with no amino-acid change.
  - Four- and five-column files behave as before, an unrecognised biotype is rejected with an error listing the accepted vocabulary, and GFF/GTF input is unaffected.
- **Scientific**: a variant inside a real intron is now reported as `intron_variant` (`MODIFIER`) against the gene it sits in, instead of `intergenic_variant`, which dropped the gene entirely. Splice sites keep their more specific terms, and an unspliced transcript has no intronic positions, so prokaryotic annotation is unchanged.
- **Scientific**: an in-frame indel that destroys the initiator codon is now reported as `Start lost` / `start_lost` (`HIGH`) instead of `In-frame Indel` / `inframe_deletion`-`inframe_insertion` (`MODERATE`), matching substitutions. The amino-acid notation is unchanged (`Met1del`, `Met1delinsArgVal`); frameshifts keep the frameshift label, and a CDS that does not begin at an initiator is left alone.
- **Scientific**: only CDS segment pairs separated by a real gap now count as exon-exon junctions, so a ribosomal-slippage join no longer draws `splice_donor_variant`, `splice_acceptor_variant` or `splice_region_variant`, nor an NMD prediction under the 50-nucleotide rule.
- **Scientific**: `--report-from` no longer merges distinct samples. Labels came from the file name alone, so the per-sample layout `results/<sample>/variants.MNV.tsv` collapsed every input into one row of the report matrix. Colliding inputs are now qualified with as many parent directories as needed, and the chosen label is logged per file.
- **Scientific**: `--split-multiallelic` no longer silently drops alternate ALT alleles when two or more alts share the same codon position. Each alt now produces an independent annotation row with its own AA effect, codon, and BAM-derived read support; true duplicates (same position + same alt) still collapse to one row.
- **Scientific**: a multi-base substitution outside every gene now reaches the VCF whole. The writer read only the first of the row's per-base entries, so the second base of a two-base substitution appeared in the TSV and nowhere in the VCF.
- **Scientific**: an intergenic row is judged in the VCF the way the TSV judges it. Multi-base substitutions were filtered on SNP support, which is zero by construction, and intergenic indels on a hardcoded support of 0 over depth 0 despite writing `ER=20` and `EFREQ=1.0000`. Both now use the haplotype support the TSV uses, and carry it in the INFO.
- **Scientific**: a declared `GT`/`PS` phase is matched to the change it was made about, and is also read for rows the per-gene path does not build, such as an MNP one base outside a gene.
- **Scientific**: a premature stop is compared in the numbering of the rows it judges, undoing the indel's length change and honouring a transcript's `protein_offset`, so translated codons are no longer reported as untranslated `MODIFIER` rows.
- **Scientific**: an NMD verdict is given only where a stop is premature, so in-frame deletions in spliced transcripts no longer get one beside their `inframe_deletion` term, and large in-frame insertions are no longer skipped.
- **Scientific**: the bases a declared phase skips now count towards an indel's length, so in-frame deletions are no longer misreported as frameshifts with `(fs)` downstream, and a row reaching into those bases reports unknown residues.
- **Scientific**: a splice site is now called from the bases a record changes, the most severe consequence winning, so a left-padded record destroying a donor is no longer called intergenic and deleted by `--exclude-intergenic`.
- **Scientific**: a record straddling a gene's edge is annotated once and only where it is: each base is classified on its own and the row built per group, ending the duplicate `Unknown` row that survived `--exclude-intergenic`.
- **Scientific**: a gene keeps its own account of a base another gene annotates, so an overlapping gene's `splice_donor_variant` (`HIGH`) and a gene's `intron_variant` survive; rows are keyed on the consequence, not the `Name`.
- **Scientific**: an indel reaches both codons of a base a ribosomal-slippage join reads twice. Substitutions already used both transcript offsets, indels asked only for the first, so a deletion of the shared base left the second codon believing nothing had happened and its substitution row named a residue for a codon that is no longer whole.
- **Scientific**: an in-frame indel that leaves the made protein unchanged is reported as such. An insertion inside a gene's terminal stop codon that keeps a stop there was reported as a residue inserted past the protein's last one, at `MODERATE` impact, because the comparison ran past the stop; stop gained and stop lost still answer from the whole translation.
- **Scientific**: a strand counts for a haplotype only when a read on it saw the whole of it. Both strands were credited when only the forward mate covered an insertion, so a composite row reported 20 reverse reads for an insertion no reverse read observed and survived `--min-mnv-strand` while the plain insertion row was dropped for having none.
- **Scientific**: the requested translation table reaches the stop rules, and a table-blind translation helper beside the table-aware one is deprecated. It has no callers, and any future call site would silently answer as if the run were table 11 whatever `--translation-table` says.
- **Scientific**: the documented insertion boundary rule now applies in the spliced-transcript path too, and is mirrored on the minus strand. An insertion anchored at a codon's last base leaves both codons intact; it was blanking that codon to `Indel overlap`, demoting a `HIGH` `start_lost` to `MODIFIER` when the only difference between two inputs was a `Parent` attribute.
- **Scientific**: a gene name holding a tab or newline no longer breaks the TSV it is written into: TSV cells now carry the same encoding the VCF writer already applied.
- **Scientific**: the report's haplotype panel no longer merges two haplotypes from different contigs. Grouping by gene, positions and alternate bases left out the contig, so the same codon MNV on a chromosome and on a plasmid came out as one row under one contig, with both samples pooled into its recurrence count.
- **Scientific**: a left-padded insertion's codon now comes from the event's components, matching `Event Components`, `HGVS g.` and the amino-acid columns; `Positions`, `Reference Bases` and `Base Changes` still echo the record, and `--normalize-alleles` still trims them.
- **Scientific**: an insertion that truncates nothing is not given an NMD prediction. Whether the alternate stop lay downstream of an indel was decided against the insertion's anchor offset, which on the plus strand is one too low, so an insertion inside the terminal stop codon was reported as producing a premature termination codon and given an NMD call.
- **Scientific**: an indel in the splice region of a spliced transcript is one row, not two: a row stating a protein consequence now carries both terms at the more severe impact, otherwise the `splice_region_variant` alone.
- **Scientific**: a variant keeps the depth and frequency its own record declared: each `ODP`/`OFREQ` entry now carries its own presence, and a combined row writes the missing one as the VCF's `.`.
- **Scientific**: a symbolic allele no longer claims a frameshift nobody measured: `<INV>`, `<DEL>` and `<DUP>`, whose `SVTYPE`, `END` and `SVLEN` get_MNV does not read, are now reported as an unknown coding change rather than `frameshift_variant` at HIGH.
- **Scientific**: a plain `.vcf` and the same bytes bgzipped annotate the same: the text parser now reads the per-ALT FORMAT `FREQ` element rather than the whole field, so multiallelic records no longer diverge in consequence or `--strict`.
- **Scientific**: a strand-bias p-value is reported at the precision its FILTER was decided on. `SBP` and `MSBP` were written with six fixed decimals, so everything below 5e-7 printed as exactly `0.000000` while the FILTER decision compared the full-precision value. They are written in scientific notation now, as the linkage p-value beside them already was.
- **Scientific**: the report shows the frequency and read count its own TSV holds for a SNP: each row now takes its numbers from the columns for its own type, `--report-from` rebuilds included.
- **Scientific**: a double quote in a GFF3 attribute value no longer swallows the rest of the column, `Parent` included: a quote delimits a value only after a GFF3 `=` or a GTF space.
- **Scientific**: a non-coding gene read from a GFF is not translated: `gene_biotype`, `biotype` or `transcript_biotype` now decides, then the feature type in column 3, so `--gff-features` defaulting to `gene,pseudogene` no longer forces protein-coding.
- **Scientific**: the desktop app no longer pairs a sample with another sample's reads.
  - Pairing is decided across the whole set, strongest match first.
  - A variant file with several equally good candidates is left unpaired.
  - The sample list names the BAM each sample was paired with.
- **Scientific**: the desktop app's read viewer reads a BAM that carries no base qualities: SAM's `QUAL=*` no longer counts as a quality of zero, so those reads are no longer drawn blank and called uninformative.
- **Scientific**: in the desktop app's read viewer, an insertion whose bases fall below the run's base-quality floor now counts as no evidence either way, rather than as support for the reference.
- **Scientific**: the desktop app's read viewer no longer treats CIGAR `N` like CIGAR `D`: a skipped stretch leaves those positions unset, which the viewer draws as empty and the support call reads as no evidence.
- **Scientific**: the desktop app's read viewer no longer answers an unreadable FASTA window with `N`: a file shorter than its `.fai` errors with the rebuild command, and each window stops at the next record header.
- **The report's variant table no longer paints fewer rows than it says it is showing.** Filtering while scrolled down left rows blank and unreachable, so the first visible row is now clamped to the viewport and a new result set starts at its top.
- **The report no longer counts `intergenic` as a gene.** It is the placeholder written in the Gene column for a variant outside every annotated gene, so the "Genes with a call" tile claimed one more gene than the same run's log and summary JSON, and the variant matrix drew a gene bar for it across the whole contig.
- **`HGVS g.` no longer uses the single-nucleotide `>` form for a multi-base allele.** A fallback row (intergenic, intron, splice site, or a position inside a gene without a codon) is now described by the bases that actually differ, not by the record's whole REF and ALT.
- **A multi-base substitution outside a CDS is counted.** An `MNV` row on the intergenic, intron and splice path now carries one entry per base that changes, instead of losing read support and failing with `[E000] Missing SNP read counts`.
- **A base shared by two overlapping CDS rows is annotated in both codons it occupies.** Substitutions now walk every segment of a `join(a..b, b..c)` ribosomal-slippage feature such as ORF1ab, as the indel path already did.
- **A left-padded substitution outside a gene is counted at the base that changes.** The intergenic, intron and splice path now takes positions and alleles from the decomposed components, as the gene path always did, so a record such as `350 AA>AG` is counted at base 351.
- **An insertion re-anchored onto an occupied coordinate no longer disappears.** Whether the gene path had produced a row was decided by coordinate, so an insertion re-anchored onto a position a substitution already held was judged annotated and never emitted: two records in, one row out. The question is now asked about the event a record describes, not its coordinate.
- **A reference call is no longer annotated as a variant.** A record whose ALT repeats its REF describes no change, but it came out as a row attributed to `intergenic` from inside a gene, typed `INDEL`, with an `intergenic_variant` consequence. Since callers do emit such rows, they are skipped and the count is reported rather than failing the run.
- **A record whose POS base does not change no longer produces two rows.** The check for a gene-path row now considers every position the record could have produced, so `28 GC>GA` no longer emits `g.29C>A` alongside `g.28GC>GA`. Introduced by the partial-codon fix earlier in this release.
- **A read marked as failing vendor QC is no longer counted as evidence.** The QC-fail flag (`0x200`) was not excluded where duplicate, secondary and supplementary records were, so such reads inflated support, depth and frequency. `samtools mpileup` excludes it by default: on a BAM where half the reads carry the flag, mpileup reported depth 10 and get_MNV reported 20.
- **A substitution in a trailing partial codon is no longer dropped from every output.** When a feature's length is not a multiple of three, its last position is now reported as `coding_sequence_variant` with an unknown amino-acid effect, instead of vanishing from the TSV, VCF and summary with exit 0.
- **`--exclude-intergenic` no longer deletes variants that are inside a gene.** The flag gated the whole fallback block, so splice-site and intron variants attributed to a named gene disappeared with it, including a HIGH-impact `splice_donor_variant`. Only variants outside every feature are excluded now.
- **The output VCF is written in coordinate order.** A multi-position variant emits one record per constituent position, so a variant between them was written after them and the file came out unsorted. `tabix` refused it, and `--index-vcf-gz` aborted before the BCF, summary JSON, manifest and report were written. Records are now buffered per contig and flushed in order.
- **An iVar indel row is checked against the FASTA like every other row.** The indel branches now read the TSV's own `REF` column instead of taking the anchor base from the reference, so all three row types fail the same way on a mismatched FASTA.
- **A CDS row with a multi-valued `Parent` belongs to every transcript it names.** GFF3 defines `Parent` as a comma-separated list, and each id in it is now honoured, instead of the raw value becoming a transcript called `t1,t2` that removed the exon from both real transcripts.
- **Fisher exact p-values below roughly 1e-12 are computed rather than floored.** The "at least as extreme" test used an absolute tolerance, so a 30/0/0/30 split reported 1.5e-14 where the answer is 1.7e-17, and the ordering the strand-bias filter depends on could invert. The comparison is relative now, checked against scipy.
- **An insertion written against the base after it is no longer placed on the wrong side of that base.** The anchor is now the base the insertion follows, even when it lies before the record's span, so `30 T>AT` reports `INS:29:+A` like `29 C>CA` instead of `INS:30:+A`.
- **An indel outside every annotated feature keeps its read support.** Intergenic indels are now counted one at a time with the same exact-event counter used inside genes and filtered on that support, instead of going out with `Event Reads = 0` and no `ER`/`EDP`/`EFREQ`.
- **`--normalize-alleles` no longer re-anchors an indel and destroys its support.** Suffix trimming now runs only when it leaves the event in place, so `29 CT>CTGCT` trims to `30 T>TGCT` and `5 CG>CTG` still reaches `5 C>CT`.
- **`--run-manifest` records output checksums under `--sample all`.** That mode built its own payload and listed the per-sample output paths with no checksum, contrary to the flag's description. Each sample entry now carries the SHA-256 of the files it wrote.
- **A multi-base substitution inside a gene but outside every codon no longer aborts the run.** The VCF writer now keys on whether read counts exist rather than on whether `--bam` was given, so intron, splice-region and non-coding-transcript alleles no longer fail the run with `[E000] Missing total read depth`.
- **`--report` no longer tells you to do something impossible.** The error, `--help` and the docs said to add `--both` to a `--convert` run, but clap declares the two mutually exclusive, so that advice always failed. They now say to use `--both` instead of `--convert`.
- **Documentation corrections found by running the binary against its own promises**: `--min-strand-bias-p` is VCF-only and never removes a TSV row; `--convert` writes the VCF *instead* of the TSV rather than in addition to it; and `--sample all` writes `<input_name>.sample_<SAMPLE>.MNV.tsv`, a naming pattern that appeared nowhere.
- **A `--bam` run no longer drops VCF records the TSV keeps.** Alleles with zero recomputed support now pass the same filter gate: emitted with zero counts when no threshold is set, skipped under a threshold, or tagged `LowSupport` under `--emit-filtered`.
- `Positions` and the other per-SNV columns of a minus-strand MNV are now emitted in ascending genomic order on both annotation paths, matching the `HGVS g.` descriptor, and the order is documented.
- The report's per-column categorical filters now use a prototype-less lookup, so a value named `constructor` or `toString` can no longer match an inherited property and pass a filter it was not selected in.
- `--report-from` now reports a compressed input as such instead of failing with `invalid UTF-8 in field 0 near byte index 1`; the report reads plain get_MNV TSV files.
- Protein-level insertions now use residue-aware HGVS notation naming both flanking residues (e.g. `Lys2_Phe3insGly`) instead of the bare-position `2_3insGly`. Insertions at the protein N- or C-terminus keep the bare-position form.
- VCF INFO values (`GENE`, `AA`, `CT`, `TYPE`, `EC`, `COMP`) are now percent-encoded for reserved characters (`;`, `=`, `,`, `%`, tab/newline/CR), so GFF gene names containing them can no longer corrupt the INFO column or spawn bogus keys.
- `--keep-original-info` now subsets per-allele (`Number=A/R/G`) INFO fields to the split allele when a multiallelic record is divided, instead of copying the whole array onto each single-ALT record and producing a VCF `bcftools` rejects.
- The BAM is now validated up front (exists, coordinate-sorted/indexed, header readable) before any output file is created, so a missing index fails fast with an actionable message instead of erroring inside a worker thread after partial output was written.
- Output is now transactional: if a run errors after the output files are created, the partial `.MNV.tsv` / `.MNV.vcf` / BCF files are removed on exit so downstream tooling never consumes a truncated file.
- BAM region queries now use the structured noodles `Region` API instead of a `chrom:start-end` string, so contig names containing `:` (e.g. HLA allele contigs) are queried at the correct coordinates.
- The BGZF VCF parser now rejects `POS=0` and a missing `#CHROM` header line, matching the plain-text fast parser.
- When `--gff-features` is not specified and the GFF contains `CDS` features, get_mnv now analyses `CDS` (phase- and splice-aware) automatically instead of whole-gene spans. Passing `--gff-features gene` keeps the previous whole-gene behaviour and still emits the CDS-phase-ignored warning.
- In a spliced transcript (CDS) model, an intronic indel is no longer merged into a phased coding haplotype with nearby exonic SNVs; only exonic variants participate in coding haplotype phasing.
- A variant downstream of a premature stop introduced by an upstream frameshift is no longer labelled `(fs)`; it is reported with Change Type `Downstream of premature stop`. Applies only when a single upstream frameshift indel introduces an early stop; ordinary frameshift propagation is unchanged.
- Resolved the frontend security audit by updating vulnerable transitive packages, including `brace-expansion` 5.0.6.
- Regenerated the Rust lockfile so vulnerable `rand` package entries are no longer present in the resolved dependency graph.
- Corrected CLI, GUI, and documentation wording for BCF input, BAM base-quality filtering, strand-bias INFO tags, and MNV rows that overlap indels.
- VCF records that already encode an MNV as a multi-base `REF/ALT` allele are now decomposed into codon-level haplotypes instead of being treated as generic indels.
- Deletions whose VCF anchor falls just outside a CDS/gene feature now still apply the overlapping deleted bases to the feature sequence, preserving frameshift/in-frame protein effects instead of reporting `Unknown`.
- Insertions anchored at the final base of a CDS/gene feature are no longer treated as inside that feature, and boundary-spanning indels are no longer duplicated as intergenic rows.
- Indel and complex alleles in coding regions now reconstruct the local alternate CDS sequence, respect strand/phase/protein offset, and report in-frame or frameshift protein effects instead of leaving amino-acid changes blank.
- Codons split across neighbouring CDS exons can now produce a single transcript-level MNV annotation when the selected GFF/GTF `CDS` records provide a usable transcript model.
- BAM support for indels is now counted from the CIGAR-derived observed allele across the event span, including inserted sequence and deleted reference bases; exact complex haplotypes also require the expected insertion/deletion components in the read CIGAR, so net-neutral indel complexes are not mistaken for simple MNVs.
- Phased `complex_indel` rows now preserve the original event component coordinates from the input variants, keeping ambiguous repeat-context deletions consistent with the source VCF/iVar event and the original indel row.
- The report's variant matrix now keeps every gene of a site, with bar order and tooltip gene lists sorted independently of the table so overlapping genes each get a bar.
- The desktop app keeps the chosen output formats when a preset is applied. Resetting them while preserving `vcfGz` left a run with compression on and no VCF to compress, and a user who asked for VCF output silently got a TSV.
- The report's matrix keeps its marks inside the plot. Group ticks, which mark the bases one row calls together, lacked the clamp the neighbouring gene bars had, so zooming until the window edge fell between the two bases of an MNV drew the tick over the sample names.
- A matrix cell sits on the base it stands for. The gap separating neighbouring cells was taken off the right side alone, drawing every mark a quarter of a pixel off its own coordinate; it is now taken off both sides.
- The report's region box no longer answers about a contig you did not name. Typing `chr2:150-250` into a report with no `chr2` dropped the name and moved the window to 150-250 on whatever contig was on screen. An unknown name is now reported and the view stays put, as an unknown gene name already was.
- The desktop app now warns before overwriting a cohort's `--sample all` outputs by checking for the files themselves, counting only the formats the run will write, and the previewed pattern includes the output directory.
- The two example scripts run from a clean checkout, and stop writing into the bundled data.
  - `reproduce_example_run.sh` now uses the demo BAM that ships, treats read evidence as optional, and falls back to the plain TSV annotation beside the GFF3.
  - `run_example_with_bam.sh` now computes the repository root correctly and points only at files that ship.
  - Both scripts write to a directory of their own instead of over the committed reference TSV in `example/`.
- The example report shipped with the documentation was regenerated: the page in `docs/assets` predated this release's report work and a new row-schema field, so it was not the page this version produces.
- The synthetic ground-truth suite follows the documented label for an in-frame indel that leaves the made protein unchanged. Its label mapping still followed the pre-change convention; it now derives the judgement from its own translated proteins, as it does everything else.
- The codon columns no longer depend on whether a BAM was supplied. `--bam` also emptied the `MNV Codon` cell of every SNP row, so the report showed a SNP's alternate codon without reads and nothing for the same variant with them, and a `--report-from` cohort mixing the two contradicted itself between samples.
- The desktop app's variant table now judges numeric columns by their values, so `Grantham` sorts, filters and aligns as a number, with `-` before every real number, matching the HTML report's table.
- The report's "Genes with a call" tile counts a gene per contig. It counted distinct gene names, so a report covering a chromosome and its plasmid claimed one gene where the log said two; it now uses the log's measure, in which a gene spread over several exons counts once.
- The `--bcf` and `--index-vcf-gz` helpers now report whether they skipped or converted, so a machine without `bcftools` or `tabix` no longer logs a conversion, records a missing path in the summary JSON, or breaks `--run-manifest`.
- The warning for a variant in the phase-skipped region of a CDS row no longer says the variant was skipped.
  - The variant is still reported against its gene with `-` in every amino-acid column and a `coding_sequence_variant` (`MODIFIER`) call.
  - The troubleshooting page is corrected, and an end-to-end test pins the row in the output file.

## [1.1.4] - 2026-06-09

### Fixed
- Intergenic variants were silently dropped under read or strand thresholds because they were never read-counted. Read counting runs gene by gene, so intergenic positions reached the output filters with a support of 0, and any positive `--snp`, `--mnv`, `--min-snp-strand`, `--min-mnv-strand`, or `--min-snp-frequency` removed them whenever a BAM was supplied (affecting the VCF since 1.1.0 and the TSV since 1.1.3). Intergenic SNPs are now read-counted at their own position and filtered by their real support, exactly like SNPs inside genes, so a threshold applies uniformly to every variant. Use `--exclude-intergenic` to drop intergenic variants on purpose.

## [1.1.3] - 2026-05-11

### Added
- **iVar variant TSV input support** via the new `--tsv <FILE>` option. The new iVar parser reads the standard `REGION`, `POS`, `REF`, `ALT`, `TOTAL_DP`, `REF_DP`, `ALT_DP`, `ALT_FREQ`, and `PASS` columns and maps passing SNVs into the internal variant model used by the VCF pipeline.
- **Automatic iVar TSV detection** for older commands that pass a typical iVar `variants.tsv` file through the existing `--vcf` input argument.
- **Desktop GUI support for iVar TSV variant files.** The GUI now accepts variant `.tsv` files, detects iVar headers separately from annotation TSV files, and sends the selected input format to the Rust backend.
- **BAM-derived SNP and MNV frequency filters** via `--min-snp-frequency` and `--min-mnv-frequency`, including GUI controls and VCF `LowFrequency` FILTER tags when `--emit-filtered` is enabled.

### Changed
- `--vcf` is now documented for VCF/BCF inputs and `--tsv` is the explicit public option for iVar TSV inputs. The legacy `--input-format tsv` and `--input-format ivar` forms remain accepted for backwards compatibility.
- The GUI output-directory logic now defaults to a writable location next to the selected variant input and preserves explicit output directories across reruns.
- Documentation now explains the difference between original input frequency (`OFREQ`) and BAM-derived frequency filtering.
- The macOS GUI build flow now uses `scripts/build_gui_bundle.sh`, which builds the Tauri `.app` first and then creates the `.dmg` with `hdiutil` for reproducible local packaging.

### Fixed
- The GUI no longer treats every `.tsv` file as a gene annotation file; iVar variant TSV headers are recognized as variant input.
- GUI reruns now avoid stale output-path conflicts and read-only working-directory errors.
- `SNP/MNV` TSV rows now apply `--min-snp-frequency` and `--min-mnv-frequency` independently, so a high-frequency MNV haplotype is not dropped only because individual SNP observations are below the SNP-frequency threshold.
- TSV output now also honors read-count and strand-support filters (`--snp`, `--mnv`, `--min-snp-strand`, `--min-mnv-strand`) using independent SNP and MNV components.
- The README now documents the correct mapping-quality flag, `--min-mapq`.
- `scripts/dev.sh build` now points to the current `src-tauri/` application layout instead of the old desktop path.

## [1.1.2] - 2026-04-08

### Fixed
- **GFF phase (column 8) is now honoured for CDS features.** Previously, the parser read all nine GFF fields but silently dropped the phase column, so codon boundaries were computed as if every CDS started in frame 0. For eukaryotic annotations where exons begin at non-zero phase (e.g. GRCh38 Ensembl `GNAQ` CDS with `phase=1`), this caused SNVs to be grouped into the wrong codons on both strands. The GFF parser now reads column 8 (`0`, `1`, `2`, or `.`), stores it in `Gene`, and `codon_bounds_for_position` applies the offset from `gene.start` (plus strand) or `gene.end` (minus strand). TSV gene files and non-CDS features default to `phase = 0`, preserving previous behaviour. Fixes [#12](https://github.com/PathoGenOmics-Lab/get_MNV/issues/12).

### Added
- Regression tests for phase handling on both strands, including the GNAQ case from issue #12.
- TSV gene annotation files now accept an **optional 5th column** with the phase (`0`, `1`, `2` or `.`). When omitted, phase defaults to `0`, preserving the historical 4-column prokaryote-style format. 4- and 5-column rows can be mixed in the same file.
- Documentation: `docs/input-formats.md` describes the GFF phase handling and the new optional TSV phase column.
- **Two new TSV columns: `Local AA Changes` and `Local SNP AA Changes`** which carry the per-feature (exon-local) amino-acid numbering used by versions ≤ 1.1.1, alongside the new full-protein numbering in `AA Changes` / `SNP AA Changes`. For prokaryotes and single-exon features the local columns are identical to the protein columns; for multi-exon eukaryotic CDS the two pairs differ (e.g. GNAQ exon 5: `Met227Ser` vs `Met25Ser`). This lets downstream pipelines that depended on the old numbering keep working without giving up the VEP-compatible reporting.

### Changed
- **Protein amino-acid numbering is now transcript-aware for multi-exon eukaryotic CDS.** Previously, each CDS row in the GFF was treated as an independent feature and the amino-acid position of an MNV was computed locally to that exon (so a variant in GNAQ exon 5 would be reported as `Met25` even though the actual protein residue is `Met224`). The GFF parser now groups CDS rows by `transcript_id` (falling back to `Parent`), sorts the exons in transcript order (ascending for plus strand, descending for minus strand), and assigns each exon a cumulative `protein_offset` equal to the number of complete codons contributed by all prior exons of the same transcript. The final amino-acid position is `protein_offset + local_codon_index + 1`, matching the numbering used by ANNOVAR, SnpEff, Ensembl VEP, and UniProt. Prokaryotic inputs (TSV gene files or single-feature GFFs) are unaffected — `protein_offset` defaults to 0.
- **Warning when `--gff-features` does not include `CDS` but the GFF contains CDS rows with non-zero phase.** In that configuration the phase-aware codon grouping is silently disabled, which was the exact trap reported in issue #12. A single `WARN` log line now tells the user to re-run with `--gff-features CDS` so the phase information is honoured.
- **Gene name extraction from GFF attributes** is smarter for eukaryotic GTF/GFF3 files. `gene_name_from_gff` now prefers `gene_name` (GTF/Ensembl/GENCODE) and `gene_id` in addition to the existing `gene` / `Name` / `locus_tag` / `ID` fallbacks, so a GNAQ CDS row with `gene_name=GNAQ;ID=agat-cds-37838` now shows up as `GNAQ` in the TSV instead of `agat-cds-37838_agat-cds-37838`.
- When no `locus_tag` is present, the returned gene name is no longer duplicated with itself (`primary_primary`) — it is just the primary name.

## [1.1.1] - 2026-03-30

### Added
- `--translation-table <N>` flag to select NCBI genetic code tables for codon translation (default: 11 — Bacterial/Archaeal/Plant Plastid). Supported tables: 1 (Standard), 2 (Vertebrate Mitochondrial), 3 (Yeast Mitochondrial), 4 (Mold/Protozoan Mitochondrial), 5 (Invertebrate Mitochondrial), 6 (Ciliate), 11 (Bacterial), 12 (Alternative Yeast Nuclear), 25 (SR1/Gracilibacteria).
- New `genetic_code` module with `GeneticCode` struct and full NCBI table support.
- 82 new tests (+58%): unit tests for 6 previously untested modules (`io/annotation`, `io/fasta`, `io/validation`, `io/vcf_fast`, `variants/types`, `cli`), plus integration tests for malformed inputs (empty VCF, truncated records, missing headers, error JSON).
- Integration tests for edge-case VCF inputs (empty, truncated, no header).
- **Desktop GUI**: native app (Tauri) with drag-and-drop, parameter presets, genomic track viewer (BAM pileup with codon annotation, IGV-style colors, per-position coverage), and multi-sample batch analysis.
- **Build & Release** GitHub Actions workflow for macOS (ARM + Intel), Linux, and Windows.
- AGPL-3.0 license (previously GPL-3.0).

### Changed
- **Performance: 6.2× faster** (111ms → 18ms on example dataset):
  - Fast text VCF parser bypassing htslib for plain `.vcf` files (9.2× faster parsing).
  - Zero-copy FASTA loading with hand-rolled parser (removed `bio` crate dependency).
  - Lazy SHA-256 checksums: only computed when `--summary-json` or `--run-manifest` is used.
  - O(1) LRU cache (replaced linear `SimpleLruCache` with `lru` crate).
  - `Rc<ReadKey>` shared ownership eliminates per-read qname clones.
- **Architecture: modular split** (max file 901 → 761 LOC):
  - `variants.rs` → `variants/types.rs` (domain types) + `variants/codon.rs` (MNV logic).
  - `output/common.rs` → + `output/stats.rs` (Fisher exact strand bias).
  - `pipeline/mod.rs` → + `pipeline/output_paths.rs` (path resolution).
- CLI migrated from clap builder API (301 LOC) to derive macros (~170 LOC).
- `VcfWriter::new()` 15 positional args replaced with `VcfWriterConfig` builder struct.
- Removed `protein-translate` dependency: inline `translate_codon()` lookup table.
- Overflow-safe arithmetic: `saturating_sub` for amino acid position calculation.
- Bounds check: skip codons exceeding reference sequence length.

### Fixed
- **Scientific**: ambiguous codon (`X`) was classified as "Synonymous" instead of "Unknown".
- **Scientific**: lowercase ALT alleles produced unknown amino acid `X` (now uppercased before translation).
- **Scientific**: duplicate same-position SNPs from `--split-multiallelic` caused incorrect codon grouping.
- Incomplete codon (gene length not multiple of 3) now logs a debug warning instead of silently skipping.
- MNV `original_info` now merges INFO from all SNPs in a codon group (deduplicated, pipe-separated).
- `get_base_name()` strips `.vcf.gz` as compound extension for clean output filenames.
- `sanitized_command_line()` escapes tab/newline/CR to prevent VCF header corruption.
- VCF contig headers: control characters replaced with underscore.
- Clap derive: fixed `required_unless_present = "gff"` → `"gff_file"` (was crashing `--help`).

### Removed
- **`rust-htslib` dependency**: migrated to [noodles](https://github.com/zaeleus/noodles) (pure Rust). Eliminates all C compilation requirements, enabling native Windows builds without POSIX toolchains.
- `bio` crate dependency (replaced with hand-rolled FASTA parser).
- `protein-translate` crate dependency (replaced with inline lookup table).
- BCF as input format (use `bcftools view input.bcf > input.vcf` to convert). BCF output is still supported via external `bcftools`.

## [1.1.0] - 2026-03-12

### Added
- **Intergenic variant support**: SNPs outside annotated gene boundaries are now preserved in the output as intergenic entries (gene = `intergenic`, change type = `Unknown`, no codon/AA annotation). Previously these variants were silently discarded.
- `--exclude-intergenic` flag to exclude intergenic SNPs from the output (default: included).
- `--keep-original-info` flag to preserve all original INFO fields from the input VCF in the output VCF (e.g. SnpEff `ANN`, VEP `CSQ`). Only get_mnv-generated INFO tags are replaced; all others are carried through verbatim.

### Changed
- Desktop app: added "Exclude intergenic" and "Keep original INFO" toggles to the parameter form.
- Desktop app: improved sample switching UI in results view with larger tabs, a "Sample:" label, and better visual hierarchy.

## [1.0.2] - 2026-03-05

### Added
- `--gff-features <FEATURES>` argument for configurable GFF feature type filtering (default: `gene,pseudogene`). Allows analyzing CDS, tRNA, exon, or any other GFF feature type.
- Pseudogene detection in GFF/GFF3 parser (previously only `gene` features were parsed).
- Warning when `--gff-features` is used with TSV annotation format (where it has no effect).
- 11 new unit tests for GFF parsing, percent-decoding, attribute splitting, and feature type filtering.

### Changed
- Refactored GFF parser: extracted shared logic into `GffGeneRecord` struct and `parse_gff_gene_records` helper, eliminating ~40 lines of duplication between `load_genes_from_gff` and `preload_gff_genes`.
- Optimized `iupac_aa` lookup from `HashMap` to direct `match` expression.
- Optimized SNP-gene interval queries with binary search (`partition_point`) instead of linear scan.
- Optimized reference sequence slicing with pre-computed indices.
- GFF genes are now preloaded once and shared across contigs (eliminates redundant file reads).
- Consistent field separators in output: `" ; "` → `"; "` and `","` → `", "`.

### Fixed
- `configure_threads` crash when using `--sample all --threads N` (Rayon global pool was configured multiple times).
- `decode_percent_encoded` produced mojibake for multibyte UTF-8 sequences (e.g. `%C3%A9` now correctly decodes to `é`).

## [1.0.1] - 2026-02-21

### Added
- `--gff` support for gene annotation (GFF/GFF3).
- `--chrom` support for selecting one or multiple contigs.
- `--sample` support for selecting a sample in multi-sample VCFs.
- `--strict` validation mode for requiring original VCF metrics (`ODP`/`OFREQ`).
- `--vcf-gz` support for BGZF-compressed VCF output (`.MNV.vcf.gz`).
- `--index-vcf-gz` support for automatic Tabix index creation (`.tbi`).
- `--split-multiallelic` support for optional in-tool multiallelic expansion.
- `--normalize-alleles` support to trim shared REF/ALT context before processing.
- `--summary-json` structured run summaries with schema version and phase timings.
- `--error-json` structured error output on failures.
- `--run-manifest` reproducibility manifest with command line, summary and output checksums.
- `--bcf` output generation (`.MNV.bcf`) converted from generated VCF.
- `--sample all` mode for per-sample output generation from multi-sample VCFs.
- `--emit-filtered` VCF mode to keep threshold-failing records with FILTER tags.
- Strand-direction thresholds (`--min-snp-strand`, `--min-mnv-strand`) and optional strand-bias INFO p-values (`--strand-bias-info`).
- Strand-bias thresholding with `--min-strand-bias-p`.
- Strand-aware read support counters (`SRF/SRR`, `MRF/MRR`) in TSV and VCF outputs.
- Benchmark threshold guard (`--max-avg-ms`) and benchmark thread control (`--threads`).
- Reproducibility scripts under `/analysis` and benchmark smoke checks in CI.
- Additional E2E tests for:
  - sorted VCF/TSV output
  - strict-mode behavior
  - sample selection errors
  - compressed VCF output
  - caller compatibility (bcftools/freebayes-like/GATK-like/lofreq-like records)

### Changed
- Refactored pipeline orchestration into `src/pipeline.rs` with clearer phase boundaries.
- Improved cache behavior with LRU-style regional BAM observation cache.
- Canonicalized INFO-field output ordering.
- Improved VCF metric parsing fallbacks for incomplete VCF headers (string-typed `DP`/`AF`) and caller-specific FORMAT/INFO combinations (`AD`, `AO/RO`).
- Added input SHA-256 checksums to summary JSON.
- Expanded validation and error messages for contig compatibility and malformed inputs.

### Fixed
- Incorrect handling of original annotation metrics in edge VCF inputs.
- Silent-reporting scenarios by introducing explicit structured summary/error outputs.
- Multi-sample parsing edge cases (including zero-sample headers).
