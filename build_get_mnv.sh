#!/usr/bin/env bash
set -euo pipefail

BIN_NAME="get_mnv"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${1:-$ROOT_DIR/dist}"

if ! command -v cargo >/dev/null 2>&1; then
  echo "Error: cargo is not installed or not in PATH."
  exit 1
fi

cd "$ROOT_DIR"

echo "Building release binary..."
cargo build --release --locked

mkdir -p "$OUT_DIR"
cp "target/release/$BIN_NAME" "$OUT_DIR/$BIN_NAME"
chmod +x "$OUT_DIR/$BIN_NAME"

if command -v strip >/dev/null 2>&1; then
  strip "$OUT_DIR/$BIN_NAME" 2>/dev/null || true
fi

echo "Build complete:"
echo "  Binary: $OUT_DIR/$BIN_NAME"
shasum -a 256 "$OUT_DIR/$BIN_NAME"
