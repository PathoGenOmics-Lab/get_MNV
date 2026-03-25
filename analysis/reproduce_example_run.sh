#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-$ROOT_DIR/analysis/results/example_$(date -u +%Y%m%dT%H%M%SZ)}"
EXAMPLE_DIR="$ROOT_DIR/example"

VCF_FILE="$EXAMPLE_DIR/G35894.var.snp.vcf"
BAM_FILE="$EXAMPLE_DIR/G35894.bam"
FASTA_FILE="$EXAMPLE_DIR/MTB_ancestor.fas"
GFF_FILE="$EXAMPLE_DIR/MTB_ancestor.gff3"

if [[ -x "$ROOT_DIR/dist/get_mnv" ]]; then
  BIN="$ROOT_DIR/dist/get_mnv"
elif [[ -x "$ROOT_DIR/target/release/get_mnv" ]]; then
  BIN="$ROOT_DIR/target/release/get_mnv"
else
  echo "Release binary not found. Building..."
  (cd "$ROOT_DIR" && cargo build --release --locked)
  BIN="$ROOT_DIR/target/release/get_mnv"
fi

for file in "$VCF_FILE" "$BAM_FILE" "$FASTA_FILE" "$GFF_FILE"; do
  if [[ ! -f "$file" ]]; then
    echo "Missing required file: $file"
    exit 1
  fi
done

mkdir -p "$OUT_DIR"
SUMMARY_JSON="$OUT_DIR/summary.json"

echo "Running example with TSV+VCF output"
(cd "$OUT_DIR" && "$BIN" \
  --vcf "$VCF_FILE" \
  --bam "$BAM_FILE" \
  --fasta "$FASTA_FILE" \
  --gff "$GFF_FILE" \
  --both \
  --summary-json "$SUMMARY_JSON")

echo "Running example with compressed VCF output"
(cd "$OUT_DIR" && "$BIN" \
  --vcf "$VCF_FILE" \
  --bam "$BAM_FILE" \
  --fasta "$FASTA_FILE" \
  --gff "$GFF_FILE" \
  --convert \
  --vcf-gz \
  --index-vcf-gz)

{
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "binary=$BIN"
  echo "vcf=$VCF_FILE"
  echo "bam=$BAM_FILE"
  echo "fasta=$FASTA_FILE"
  echo "gff=$GFF_FILE"
  echo "rustc=$(rustc --version)"
  echo "cargo=$(cargo --version)"
  echo "uname=$(uname -a)"
} >"$OUT_DIR/environment.txt"

(cd "$OUT_DIR" && shasum -a 256 ./*.MNV.tsv ./*.MNV.vcf ./*.MNV.vcf.gz ./*.MNV.vcf.gz.tbi summary.json >checksums.sha256)

echo "Done. Outputs in: $OUT_DIR"
