# Changelog

All notable changes to this project are documented in this file.

## [1.1.5] - 2026-06-09

### Added
- **Sequence Ontology consequence terms and impact levels.** New `SO Term` and `Impact` TSV columns (and `SO` / `IMPACT` VCF INFO fields) annotate every variant with a standard SO term (`missense_variant`, `synonymous_variant`, `stop_gained`, `stop_lost`, `start_lost`, `frameshift_variant`, `inframe_insertion`/`inframe_deletion`, `intergenic_variant`, `coding_sequence_variant`) and a `HIGH`/`MODERATE`/`LOW`/`MODIFIER` impact, following SnpEff/VEP conventions.
- **Grantham distance for missense changes.** A new `Grantham` TSV column (and `GD` VCF INFO field) reports the Grantham (1974) distance of the combined amino-acid substitution together with its conservation category (e.g. `177 (radical)`), for genuine missense changes only.
- **MNV-vs-SNV consequence shift.** A new `MNV Consequence Shift` TSV column (and `MNVSHIFT` VCF INFO field) flags whether a combined MNV is `MNV-gained` (more severe than any single SNV alone — e.g. two individually-synonymous SNVs producing a non-synonymous residue, which per-SNV annotators miss), `MNV-masked` (a nonsense SNV rescued by its neighbour), or `Concordant`.
- **Doublet base substitution (DBS) class.** A new `DBS Class` TSV column (and `DBS` VCF INFO field) reports the COSMIC-style doublet class for an MNV of two adjacent single-base substitutions, e.g. `CC>TT`, reverse-complement collapsed so both strands map to one class. Enables tallying mutational signatures (UV `CC>TT`, APOBEC/ADAR dinucleotide contexts) directly from codon-level MNV calls.
- **MNV phasing (linkage) support.** With `--bam`, a new `MNV Phasing Support` TSV column (and `MNVPS` VCF INFO field) reports the fraction of the least-supported constituent SNV's reads that also carry the full MNV haplotype. Near `1.0` indicates a genuine co-occurring haplotype; low values flag same-codon SNVs that largely fall on different molecules rather than a real MNV, which matters most for intra-host data.
- **Nonsense-mediated decay (NMD) prediction.** A new `NMD Prediction` TSV column (and `NMD` VCF INFO field) applies the 50-nucleotide rule to premature stops in multi-exon transcripts: `NMD-triggering` when the premature stop is more than 50 nt upstream of the last exon-exon junction, `NMD-escaping` when it lies in the last exon or within 50 nt of that junction. Covers nonsense SNVs/MNVs and stop-gaining indels (frameshift and in-frame); single-exon transcripts have no junction and are left unannotated.
- **HGVS genomic and coding descriptors.** New `HGVS g.` and `HGVS c.` TSV columns (and `HGVSG` / `HGVSC` VCF INFO fields) complement the existing protein-level (`p.`) notation. Genomic descriptors cover every variant: `g.100A>G` for SNVs, the allele bracket `g.[28G>T;30T>A]` for MNVs, and `g.101_102del` / `g.100_101insTG` / `g.101delinsC` for indels. Coding descriptors cover substitutions (`c.30A>G`, `c.[28G>A;30T>C]`), numbered from the CDS start with coding-strand bases on both strands. Descriptors are not 3'-shifted and carry no reference-accession prefix; coding indel descriptors are not yet emitted.
- **Splice-site consequence terms.** Variants near an internal exon-exon junction of a multi-exon transcript now receive the standard Sequence Ontology splice terms in the `SO Term` / `Impact` columns (and `SO` / `IMPACT` VCF INFO fields): `splice_donor_variant` / `splice_acceptor_variant` for the two essential intronic bases at each intron end (`HIGH`), and `splice_region_variant` for the exon's last 3 bases or the intron's 3rd-8th bases (`LOW`). Intronic splice variants, previously reported as intergenic, are now associated with their gene; an exonic coding change near a junction is combined, e.g. `missense_variant&splice_region_variant`, keeping the more severe impact. Single-exon transcripts (prokaryotes, most viruses) have no junction and are unaffected.
- Added regression coverage for phased MNV-plus-indel haplotypes, verifying that codon MNV rows overlapping an indel are flagged as `Indel overlap` while BAM-supported combined events are emitted as exact `complex_indel` rows.
- Added an indel/MNV semantics note documenting caller compatibility, boundary rules, current limits, and how exact complex haplotypes are represented.
- Added transcript-aware regression coverage for exon-junction MNV codons and restored-frame indel contexts in multi-exon CDS models.

### Changed
- `--chrom` now restricts FASTA loading: when a single contig is requested, only that contig is read and IUPAC-validated instead of the whole genome, cutting peak memory for large (e.g. eukaryotic) references. A missing requested contig now fails with a clear "not found in FASTA" error.
- SNVs/MNVs that alter the initiator Met (protein position 1) are now reported with Change Type `Start lost` instead of `Non-synonymous`, matching standard annotators. The reported amino-acid change was already correct; only the classification label changes. `Met1` → stop is still `Stop gained`.
- Added a tuned `[profile.release]` (thin LTO, single codegen unit) and wrapped plain (non-BGZF) VCF output in a buffered writer with an explicit flush, for faster production builds and record emission.
- Pinned `sha2` to the stable `0.10` line (previously a `0.11` pre-release), removing a duplicate hashing dependency stack from the resolved graph.
- `get_mnv_variants_for_gene` and `get_mnv_variants_for_transcript` now build a list of mutually-exclusive codon interpretations per codon start. Multi-allelic positions expand the interpretation set as a Cartesian product, deduplicated by `(position, alt)` keys, so a codon that contains N independent alts emits N output rows.
- Bumped project, GUI, citation, README, and frontend metadata to version 1.1.4.
- Updated the Tauri desktop dependency set to the 2.11 patch line, including `tauri` 2.11.2 and `tauri-plugin-dialog` 2.7.1.
- Updated frontend lockfile dependencies, including `postcss` 8.5.10.
- Added allele-level event decomposition for `snp`, `mnv`, `insertion`, `deletion`, `delins`, `complex_indel`, and symbolic alleles so length-changing events that also contain SNV/MNV components are represented as a single local haplotype event.
- Expanded phased local haplotype discovery from indel-plus-SNV pairs to bounded multi-event windows, allowing supported combinations such as insertion-plus-deletion haplotypes to be emitted as exact `complex_indel` rows.
- iVar TSV inputs now keep insertion and deletion rows by converting `+SEQ` and `-SEQ` alleles into VCF-compatible anchored `REF/ALT` alleles using the FASTA reference.
- TSV and VCF outputs now include canonical event class/component annotations plus exact BAM-derived event support metrics for indel/complex alleles.
- The desktop BAM viewer now renders insertion-aware interbase columns, so inserted bases are shown between reference positions with matching coverage, ruler, reference, and read-pileup alignment instead of being hidden inside the anchor base.
- When a BAM is provided, nearby SNV/MNV rows that phase with an indel on the same reads are now emitted as an additional exact `complex_indel` haplotype row, preserving the original rows while reporting the combined `REF/ALT`, protein effect, and event support.
- GFF/GTF `CDS` rows with `transcript_id` or `Parent` are now collapsed into spliced transcript CDS models, so codon grouping, MNV amino-acid effects, and indel frameshift context are evaluated against the full coding sequence instead of isolated exon rows.
- `--frameshift-min-freq` now defaults to `0.5` instead of `0.0`: an upstream indel only shifts the reading frame of downstream SNV/MNV codons when it is the consensus (majority) allele, so a high-frequency downstream substitution is no longer relabelled frameshifted because of a low-frequency upstream indel that is almost certainly on a different molecule (intra-host data). Indels without a known frequency still propagate. Pass `--frameshift-min-freq 0.0` for the previous behaviour.
- Indel locus depth (the `EDP`/`EFREQ` denominator) is now counted from reads observing the anchor base by default, removing the depth under-counting and `EFREQ` bias for multi-base deletions. Pass `--legacy-indel-depth` to restrict the denominator to reads that fully span the REF allele.

### Fixed
- Protein-level insertions now use residue-aware HGVS notation naming both flanking residues (e.g. `Lys2_Phe3insGly`), matching the `del` / `delins` forms, instead of the bare-position `2_3insGly`. Insertions at the protein N- or C-terminus, which have no flanking residue on one side, keep the bare-position form.
- **Scientific**: `--split-multiallelic` no longer silently drops alternate ALT alleles when two or more alts share the same codon position. Each alt now produces an independent annotation row with its own AA effect, codon, and BAM-derived read support; true duplicates (same position + same alt) still collapse to one row.
- VCF output INFO values (`GENE`, `AA`, `CT`, `TYPE`, `EC`, `COMP`) are now percent-encoded for the structurally reserved characters (`;`, `=`, `,`, `%`, and tab/newline/CR), so GFF gene names containing those characters can no longer corrupt the INFO column or spawn bogus keys.
- `--keep-original-info` now subsets per-allele (`Number=A/R/G`) INFO fields to the split allele when a multiallelic record is divided, instead of copying the whole array onto each single-ALT output record (which produced a cardinality-invalid VCF that `bcftools` rejects).
- The BAM is now validated up front (exists, is coordinate-sorted/indexed, header readable) before any output file is created, so a missing index fails fast with an actionable message instead of erroring lazily inside a worker thread after partial output was already written.
- Output is now transactional: if a run errors after the output files are created, the partial `.MNV.tsv` / `.MNV.vcf` / BCF files are removed on exit so downstream tooling never consumes a truncated file.
- BAM region queries are now built with the structured noodles `Region` API instead of a `chrom:start-end` string, so contig names containing `:` (e.g. HLA allele contigs) are queried at the correct coordinates instead of being misparsed.
- The BGZF VCF parser now rejects `POS=0` and a missing `#CHROM` header line, matching the plain-text fast parser so both code paths validate inputs identically.
- When `--gff-features` is not specified and the GFF contains `CDS` features, get_mnv now analyses `CDS` (phase- and splice-aware) automatically instead of whole-gene spans, so eukaryotic/multi-exon annotations are no longer silently mis-numbered over introns. Passing `--gff-features gene` keeps the previous whole-gene behaviour (and still emits the CDS-phase-ignored warning).
- In a spliced transcript (CDS) model, an indel falling in an intron is no longer merged into a phased coding haplotype with nearby exonic SNVs; only exonic variants participate in coding haplotype phasing.
- A variant lying downstream of a premature stop introduced by an upstream frameshift is no longer labelled `(fs)` as if it were translated; it is reported with Change Type `Downstream of premature stop`. Applies only when a single upstream frameshift indel introduces a stop earlier than the natural one; ordinary frameshift propagation (no early stop) is unchanged.
- Resolved the frontend security audit by updating vulnerable transitive packages, including `brace-expansion` 5.0.6.
- Regenerated the Rust lockfile so vulnerable `rand` package entries are no longer present in the resolved dependency graph.
- Corrected CLI, GUI, and documentation wording for BCF input, BAM base-quality filtering, strand-bias INFO tags, and MNV rows that overlap indels.
- VCF records that already encode an MNV as a multi-base `REF/ALT` allele are now decomposed into codon-level haplotypes instead of being treated as generic indels.
- Deletions whose VCF anchor falls just outside a CDS/gene feature now still apply the overlapping deleted bases to the feature sequence, preserving frameshift/in-frame protein effects instead of reporting `Unknown`.
- Insertions anchored at the final base of a CDS/gene feature are no longer treated as if the inserted sequence were inside that feature, and boundary-spanning indels are no longer duplicated as intergenic rows.
- Indel and complex alleles in coding regions now reconstruct the local alternate CDS sequence, respect strand/phase/protein offset, and report in-frame or frameshift protein effects instead of leaving amino-acid changes blank.
- Codons split across neighbouring CDS exons can now produce a single transcript-level MNV annotation when the selected GFF/GTF `CDS` records provide a usable transcript model.
- BAM support for indels is now counted from the CIGAR-derived observed allele across the event span, including inserted sequence and deleted reference bases; exact complex haplotypes now also require the expected insertion/deletion components in the read CIGAR so net-neutral indel complexes are not mistaken for simple MNVs.
- Phased `complex_indel` rows now preserve the original event component coordinates from the input variants, so ambiguous repeat-context deletions remain consistent with the source VCF/iVar event and the original indel row.

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
