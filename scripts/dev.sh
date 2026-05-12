#!/usr/bin/env bash
# dev.sh — Launch get_MNV desktop (Tauri + React) in development mode.
# Usage: ./dev.sh         (dev mode with hot-reload)
#        ./dev.sh build   (production build)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FRONTEND="$ROOT/frontend"
DESKTOP="$ROOT/src-tauri"

# Ensure frontend dependencies are installed
if [ ! -d "$FRONTEND/node_modules" ]; then
  echo "==> Installing frontend dependencies..."
  (cd "$FRONTEND" && npm install)
fi

if [ "${1:-}" = "build" ]; then
  "$SCRIPT_DIR/build_gui_bundle.sh"
else
  # Start Vite dev server in background
  echo "==> Starting Vite dev server..."
  (cd "$FRONTEND" && npm run dev) &
  VITE_PID=$!

  # Give Vite a moment to start
  sleep 2

  echo "==> Starting Tauri..."
  (cd "$DESKTOP" && cargo tauri dev)

  # Clean up Vite when Tauri exits
  kill "$VITE_PID" 2>/dev/null || true
fi
