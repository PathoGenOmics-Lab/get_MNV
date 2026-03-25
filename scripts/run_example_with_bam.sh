#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXAMPLE_DIR="$ROOT_DIR/example"

VCF_FILE="$EXAMPLE_DIR/MIP00022.MTB_anc.ann.vcf"
BAM_FILE="$EXAMPLE_DIR/MIP00022.MTB_anc.final.bam"
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
  echo "Build it first with: $ROOT_DIR/build_get_mnv.sh"
  exit 1
fi

for file in "$VCF_FILE" "$BAM_FILE" "$FASTA_FILE" "$BAM_FILE.bai"; do
  if [[ ! -f "$file" ]]; then
    echo "Error: missing required file: $file"
    exit 1
  fi
done

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

cd "$EXAMPLE_DIR"
extra_args=()
for arg in "$@"; do
  if [[ "$arg" != "--convert" && "$arg" != "--both" ]]; then
    extra_args+=("$arg")
  fi
done

echo "Running get_mnv with example files (TSV + VCF)..."
echo "Binary: $BIN"
echo "Annotation file: $GENES_FILE ($ANNOTATION_FLAG)"

if "$BIN" --help 2>/dev/null | grep -q -- "--both"; then
  if [[ ${#extra_args[@]} -gt 0 ]]; then
    "$BIN" \
      --vcf "$VCF_FILE" \
      --bam "$BAM_FILE" \
      --fasta "$FASTA_FILE" \
      "$ANNOTATION_FLAG" "$GENES_FILE" \
      --both \
      "${extra_args[@]}"
  else
    "$BIN" \
      --vcf "$VCF_FILE" \
      --bam "$BAM_FILE" \
      --fasta "$FASTA_FILE" \
      "$ANNOTATION_FLAG" "$GENES_FILE" \
      --both
  fi
else
  echo "Warning: binary does not support --both. Falling back to two runs."
  if [[ ${#extra_args[@]} -gt 0 ]]; then
    "$BIN" \
      --vcf "$VCF_FILE" \
      --bam "$BAM_FILE" \
      --fasta "$FASTA_FILE" \
      "$ANNOTATION_FLAG" "$GENES_FILE" \
      "${extra_args[@]}"
    "$BIN" \
      --vcf "$VCF_FILE" \
      --bam "$BAM_FILE" \
      --fasta "$FASTA_FILE" \
      "$ANNOTATION_FLAG" "$GENES_FILE" \
      "${extra_args[@]}" \
      --convert
  else
    "$BIN" \
      --vcf "$VCF_FILE" \
      --bam "$BAM_FILE" \
      --fasta "$FASTA_FILE" \
      "$ANNOTATION_FLAG" "$GENES_FILE"
    "$BIN" \
      --vcf "$VCF_FILE" \
      --bam "$BAM_FILE" \
      --fasta "$FASTA_FILE" \
      "$ANNOTATION_FLAG" "$GENES_FILE" \
      --convert
  fi
fi

echo "Done. Output written in: $EXAMPLE_DIR"
