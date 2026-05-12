#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
FRONTEND_DIR="$ROOT_DIR/frontend"
TAURI_DIR="$ROOT_DIR/src-tauri"
BUNDLE_DIR="$ROOT_DIR/target/release/bundle"
MACOS_DIR="$BUNDLE_DIR/macos"
DMG_DIR="$BUNDLE_DIR/dmg"
PRODUCT_NAME="get_MNV"
APP_NAME="$PRODUCT_NAME.app"

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Error: required command '$1' is not installed or not in PATH." >&2
    exit 1
  fi
}

require_command cargo
require_command npm
require_command node
require_command hdiutil
require_command ditto

VERSION="$(
  cd "$ROOT_DIR"
  node -e 'const fs = require("fs"); const cfg = JSON.parse(fs.readFileSync("src-tauri/tauri.conf.json", "utf8")); process.stdout.write(cfg.version);'
)"

case "$(uname -m)" in
  arm64) ARCH_LABEL="aarch64" ;;
  x86_64) ARCH_LABEL="x64" ;;
  *) ARCH_LABEL="$(uname -m)" ;;
esac

APP_PATH="$MACOS_DIR/$APP_NAME"
DMG_PATH="$DMG_DIR/${PRODUCT_NAME}_${VERSION}_${ARCH_LABEL}.dmg"
STAGE_DIR="${TMPDIR:-/tmp}/get_mnv_gui_dmg_stage_$$"

cleanup() {
  rm -rf "$STAGE_DIR"
}
trap cleanup EXIT

if [ ! -d "$FRONTEND_DIR/node_modules" ]; then
  echo "==> Installing frontend dependencies..."
  npm install --prefix "$FRONTEND_DIR"
fi

echo "==> Building frontend..."
npm run build --prefix "$FRONTEND_DIR"

echo "==> Building Tauri app bundle..."
(
  cd "$TAURI_DIR"
  cargo tauri build --bundles app
)

if [ ! -d "$APP_PATH" ]; then
  echo "Error: expected app bundle was not created: $APP_PATH" >&2
  exit 1
fi

echo "==> Creating DMG..."
mkdir -p "$DMG_DIR"
rm -rf "$STAGE_DIR"
mkdir -p "$STAGE_DIR"
ditto "$APP_PATH" "$STAGE_DIR/$APP_NAME"
ln -s /Applications "$STAGE_DIR/Applications"
rm -f "$DMG_PATH"
hdiutil create -volname "$PRODUCT_NAME" -srcfolder "$STAGE_DIR" -ov -format UDZO "$DMG_PATH"
hdiutil verify "$DMG_PATH"

echo "Build complete:"
echo "  App: $APP_PATH"
echo "  DMG: $DMG_PATH"
