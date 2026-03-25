#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-$ROOT_DIR/analysis/results/bench_$(date -u +%Y%m%dT%H%M%SZ)}"
DATASET_DIR="${DATASET_DIR:-$ROOT_DIR/example}"
WARMUP="${WARMUP:-3}"
ITERS="${ITERS:-10}"
THREADS_LIST="${THREADS_LIST:-1 2 4}"
MAX_AVG_MS="${MAX_AVG_MS:-}"

mkdir -p "$OUT_DIR"
CSV_PATH="$OUT_DIR/benchmark.csv"
ENV_PATH="$OUT_DIR/environment.txt"

{
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "root_dir=$ROOT_DIR"
  echo "dataset_dir=$DATASET_DIR"
  echo "warmup=$WARMUP"
  echo "iters=$ITERS"
  echo "threads_list=$THREADS_LIST"
  echo "rustc=$(rustc --version)"
  echo "cargo=$(cargo --version)"
  echo "uname=$(uname -a)"
} >"$ENV_PATH"

echo "Benchmark output directory: $OUT_DIR"
echo "Benchmark CSV: $CSV_PATH"

for threads in $THREADS_LIST; do
  cmd=(
    cargo run --release --bin bench_variants --
    --dataset "$DATASET_DIR"
    --warmup "$WARMUP"
    --iters "$ITERS"
    --threads "$threads"
    --csv "$CSV_PATH"
  )
  if [[ -n "$MAX_AVG_MS" ]]; then
    cmd+=(--max-avg-ms "$MAX_AVG_MS")
  fi

  echo "Running benchmark with threads=$threads"
  (cd "$ROOT_DIR" && "${cmd[@]}")
done

echo "Done. Results:"
echo "  $CSV_PATH"
echo "  $ENV_PATH"
