#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-$ROOT_DIR/analysis/results/example_$(date -u +%Y%m%dT%H%M%SZ)}"
EXAMPLE_DIR="$ROOT_DIR/example"

VCF_FILE="$EXAMPLE_DIR/G35894.var.snp.vcf"
BAM_FILE="$EXAMPLE_DIR/G35894.demo.bam"
FASTA_FILE="$EXAMPLE_DIR/MTB_ancestor.fas"
GFF_FILE="$EXAMPLE_DIR/MTB_ancestor.gff3"
TSV_GENES_FILE="$EXAMPLE_DIR/anot_genes.txt"
ANNOTATION_FLAG="--gff"
ANNOTATION_FILE="$GFF_FILE"

if [[ -x "$ROOT_DIR/dist/get_mnv" ]]; then
  BIN="$ROOT_DIR/dist/get_mnv"
elif [[ -x "$ROOT_DIR/target/release/get_mnv" ]]; then
  BIN="$ROOT_DIR/target/release/get_mnv"
else
  echo "Release binary not found. Building..."
  (cd "$ROOT_DIR" && cargo build --release --locked)
  BIN="$ROOT_DIR/target/release/get_mnv"
fi

for file in "$VCF_FILE" "$FASTA_FILE"; do
  if [[ ! -f "$file" ]]; then
    echo "Missing required file: $file"
    exit 1
  fi
done

# The bundle ships a plain TSV annotation; the GFF3 is optional and not part of
# it. Asking for only the GFF3 made this script stop on a clean checkout.
if [[ ! -f "$ANNOTATION_FILE" ]]; then
  if [[ -f "$TSV_GENES_FILE" ]]; then
    ANNOTATION_FILE="$TSV_GENES_FILE"
    ANNOTATION_FLAG="--genes"
  else
    echo "Missing annotation file. Expected either:"
    echo "  $GFF_FILE"
    echo "  $TSV_GENES_FILE"
    exit 1
  fi
fi

# Read evidence is what --bam adds, and the bundle ships a demo BAM covering
# part of the genome. Without one the run still reproduces, with the read
# columns absent rather than empty, so it is reported and not required.
BAM_ARGS=()
if [[ -f "$BAM_FILE" && -f "$BAM_FILE.bai" ]]; then
  BAM_ARGS=(--bam "$BAM_FILE")
else
  echo "No indexed BAM at $BAM_FILE; running without read evidence"
  BAM_FILE="(none)"
fi

mkdir -p "$OUT_DIR"
SUMMARY_JSON="$OUT_DIR/summary.json"

echo "Running example with TSV+VCF output"
echo "Annotation: $ANNOTATION_FILE ($ANNOTATION_FLAG)"
(cd "$OUT_DIR" && "$BIN" \
  --vcf "$VCF_FILE" \
  ${BAM_ARGS[@]+"${BAM_ARGS[@]}"} \
  --fasta "$FASTA_FILE" \
  "$ANNOTATION_FLAG" "$ANNOTATION_FILE" \
  --both \
  --summary-json "$SUMMARY_JSON")

echo "Running example with compressed VCF output"
(cd "$OUT_DIR" && "$BIN" \
  --vcf "$VCF_FILE" \
  ${BAM_ARGS[@]+"${BAM_ARGS[@]}"} \
  --fasta "$FASTA_FILE" \
  "$ANNOTATION_FLAG" "$ANNOTATION_FILE" \
  --convert \
  --vcf-gz \
  --index-vcf-gz)

{
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "binary=$BIN"
  echo "vcf=$VCF_FILE"
  echo "bam=$BAM_FILE"
  echo "fasta=$FASTA_FILE"
  echo "annotation=$ANNOTATION_FILE ($ANNOTATION_FLAG)"
  echo "rustc=$(rustc --version)"
  echo "cargo=$(cargo --version)"
  echo "uname=$(uname -a)"
} >"$OUT_DIR/environment.txt"

(cd "$OUT_DIR" && shasum -a 256 ./*.MNV.tsv ./*.MNV.vcf ./*.MNV.vcf.gz ./*.MNV.vcf.gz.tbi summary.json >checksums.sha256)

echo "Done. Outputs in: $OUT_DIR"
