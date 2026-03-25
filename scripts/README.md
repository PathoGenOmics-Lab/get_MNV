# Analysis and Reproducibility

This folder contains helper scripts to reproduce benchmark and example-run artifacts for development, QA, and publication supplements.

## Scripts

- `run_reproducible_benchmark.sh`
  - Runs `bench_variants` on the example dataset for multiple thread counts.
  - Writes:
    - `benchmark.csv` with one row per benchmark run
    - `environment.txt` with runtime metadata
- `reproduce_example_run.sh`
  - Runs `get_mnv` against bundled example inputs.
  - Produces TSV, VCF, compressed VCF (`.vcf.gz`), Tabix index (`.tbi`), summary JSON, and SHA-256 checksums.

## Usage

```bash
chmod +x analysis/run_reproducible_benchmark.sh analysis/reproduce_example_run.sh
```

```bash
analysis/run_reproducible_benchmark.sh
```

```bash
analysis/reproduce_example_run.sh
```

Optional benchmark environment overrides:
- `THREADS_LIST` (default: `1 2 4`)
- `WARMUP` (default: `3`)
- `ITERS` (default: `10`)
- `MAX_AVG_MS` (optional threshold gate)
- `DATASET_DIR` (default: `example`)
