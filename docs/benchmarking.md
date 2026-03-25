# Benchmarking

get_MNV includes a built-in benchmark binary for measuring variant annotation performance.

## Synthetic benchmark

```bash
cargo run --release --bin bench_variants -- --warmup 5 --iters 30
```

Scale up the synthetic dataset:

```bash
cargo run --release --bin bench_variants -- \
  --warmup 3 --iters 20 --threads 4 --synthetic-scale 4
```

## Dataset benchmark

Run against a real dataset with CSV export:

```bash
cargo run --release --bin bench_variants -- \
  --dataset example \
  --warmup 5 --iters 20 \
  --threads 4 \
  --csv benchmark.csv \
  --max-avg-ms 200
```

This appends one row per run to `benchmark.csv`.

## Benchmark options

| Argument | Default | Description |
|----------|---------|-------------|
| `--warmup <N>` | 5 | Warmup iterations (not measured) |
| `--iters <N>` | 30 | Measured iterations |
| `--threads <N>` | 1 | Worker threads |
| `--dataset <DIR>` | — | Path to dataset directory (FASTA + VCF + annotation) |
| `--contig <NAME>` | all | Restrict to a single contig |
| `--csv <FILE>` | — | Append CSV results |
| `--max-avg-ms <MS>` | — | Fail if average exceeds threshold (regression detection) |
| `--synthetic-scale <N>` | 1 | Scale factor for synthetic gene/SNP count |

## Reproducible scripts

See `scripts/` for batch benchmark runners:
- `scripts/run_reproducible_benchmark.sh`
- `scripts/reproduce_example_run.sh`
