#!/usr/bin/env bash
# Regenerate the example HTML report published with the documentation.
#
# The report is a committed static asset (docs/assets/example-report.html) so the
# docs site can link to something clickable. It is generated from the bundled
# example dataset, so anyone can reproduce it byte for byte. Re-run this whenever
# the report template or the annotations change, otherwise the published example
# drifts from what the tool actually emits.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXAMPLE_DIR="$ROOT_DIR/example"
OUT="$ROOT_DIR/docs/assets/example-report.html"

VCF_FILE="$EXAMPLE_DIR/G35894.var.snp.vcf"
FASTA_FILE="$EXAMPLE_DIR/MTB_ancestor.fas"
GENES_FILE="$EXAMPLE_DIR/anot_genes.txt"

if [[ -x "$ROOT_DIR/dist/get_mnv" ]]; then
  BIN="$ROOT_DIR/dist/get_mnv"
else
  cargo build --release --manifest-path "$ROOT_DIR/Cargo.toml" -p get_mnv
  BIN="$ROOT_DIR/target/release/get_mnv"
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Run from a scratch directory so the TSV/VCF side outputs do not land in the repo.
cd "$WORK"
"$BIN" \
  --vcf "$VCF_FILE" \
  --fasta "$FASTA_FILE" \
  --genes "$GENES_FILE" \
  --report "$OUT"

printf 'Wrote %s (%s)\n' "$OUT" "$(du -h "$OUT" | cut -f1)"
