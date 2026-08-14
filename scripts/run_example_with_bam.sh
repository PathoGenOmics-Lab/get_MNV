#!/usr/bin/env bash
set -euo pipefail

# A run over the bundled example with read evidence, which is what --bam adds:
# read counts, frequencies, the strand arms, phasing and linkage.
#
# Outputs go to a directory of their own. Writing them into example/ overwrote
# the committed reference output that the tests compare against.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXAMPLE_DIR="$ROOT_DIR/example"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/analysis/results/example_bam_$(date -u +%Y%m%dT%H%M%SZ)}"

VCF_FILE="$EXAMPLE_DIR/G35894.var.snp.vcf"
BAM_FILE="$EXAMPLE_DIR/G35894.demo.bam"
FASTA_FILE="$EXAMPLE_DIR/MTB_ancestor.fas"
GFF3_FILE="$EXAMPLE_DIR/MTB_ancestor.gff3"
TSV_GENES_FILE="$EXAMPLE_DIR/anot_genes.txt"
GENES_FILE="$GFF3_FILE"
ANNOTATION_FLAG="--gff"

if [[ -x "$ROOT_DIR/dist/get_mnv" ]]; then
  BIN="$ROOT_DIR/dist/get_mnv"
elif [[ -x "$ROOT_DIR/target/release/get_mnv" ]]; then
  BIN="$ROOT_DIR/target/release/get_mnv"
elif [[ -x "$ROOT_DIR/target/debug/get_mnv" ]]; then
  BIN="$ROOT_DIR/target/debug/get_mnv"
else
  echo "Error: get_mnv binary not found."
  echo "Build it first with: $ROOT_DIR/scripts/build_get_mnv.sh"
  exit 1
fi

for file in "$VCF_FILE" "$BAM_FILE" "$FASTA_FILE" "$BAM_FILE.bai"; do
  if [[ ! -f "$file" ]]; then
    echo "Error: missing required file: $file"
    exit 1
  fi
done

# The GFF3 is not part of the bundle; the plain TSV annotation beside it is.
if [[ ! -f "$GENES_FILE" ]]; then
  if [[ -f "$TSV_GENES_FILE" ]]; then
    GENES_FILE="$TSV_GENES_FILE"
    ANNOTATION_FLAG="--genes"
  else
    echo "Error: missing annotation file. Expected either:"
    echo "  $GFF3_FILE"
    echo "  $TSV_GENES_FILE"
    exit 1
  fi
fi

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

extra_args=()
for arg in "$@"; do
  if [[ "$arg" != "--convert" && "$arg" != "--both" ]]; then
    extra_args+=("$arg")
  fi
done

echo "Running get_mnv over the bundled example with its demo BAM (TSV + VCF)"
echo "Binary:     $BIN"
echo "Annotation: $GENES_FILE ($ANNOTATION_FLAG)"
echo "Reads:      $BAM_FILE"

"$BIN" \
  --vcf "$VCF_FILE" \
  --bam "$BAM_FILE" \
  --fasta "$FASTA_FILE" \
  "$ANNOTATION_FLAG" "$GENES_FILE" \
  --both \
  ${extra_args[@]+"${extra_args[@]}"}

echo "Done. Outputs in: $OUT_DIR"
