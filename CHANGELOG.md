# Changelog

All notable changes to this project are documented in this file.

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
