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
- `build_get_mnv.sh`
  - Builds the release CLI and copies it to `dist/get_mnv`.
- `build_gui_bundle.sh`
  - Builds the React frontend, Tauri `.app`, and macOS `.dmg` bundle.
  - The DMG step uses `hdiutil` directly after Tauri creates the `.app`, avoiding the generated Tauri DMG script when it is invoked without arguments on local builds.

## Usage

```bash
chmod +x scripts/*.sh
```

```bash
bash scripts/run_reproducible_benchmark.sh
```

```bash
bash scripts/reproduce_example_run.sh
```

```bash
bash scripts/build_gui_bundle.sh
```

Optional benchmark environment overrides:
- `THREADS_LIST` (default: `1 2 4`)
- `WARMUP` (default: `3`)
- `ITERS` (default: `10`)
- `MAX_AVG_MS` (optional threshold gate)
- `DATASET_DIR` (default: `example`)
