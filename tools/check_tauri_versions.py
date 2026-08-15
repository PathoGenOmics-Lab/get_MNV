#!/usr/bin/env python3
"""The Tauri crate and the Tauri npm package are on the same minor.

Tauri refuses to bundle when they are not:

    Found version mismatched Tauri packages. Make sure the NPM package and Rust
    crate versions are on the same major/minor releases:
    tauri (v2.11.2) : @tauri-apps/api (v2.10.1)

That is one dependency spelled in two ecosystems, and nothing bumps them
together: Dependabot treats cargo and npm separately, so it can move one and
leave the other. Nothing in the repository noticed, because the check lives
inside `tauri build`, which only runs when a release is being published. The
1.1.5 release failed on all four platforms for exactly this, after the tag was
already pushed.

Reads the resolved versions, not the requirements: a caret range says what is
allowed, and the lockfiles say what would actually be built.

    python3 tools/check_tauri_versions.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Each pair is (Rust crate, npm package). Tauri checks the first one itself; the
# plugins are held to the same rule because they are the same coupling.
PAIRS = (
    ("tauri", "@tauri-apps/api"),
    ("tauri-plugin-dialog", "@tauri-apps/plugin-dialog"),
    ("tauri-plugin-shell", "@tauri-apps/plugin-shell"),
)


def crate_versions() -> dict[str, str]:
    """Every crate in Cargo.lock, by name."""
    text = (REPO / "Cargo.lock").read_text()
    out: dict[str, str] = {}
    for block in text.split("[[package]]"):
        name = re.search(r'^name = "([^"]+)"', block, re.M)
        version = re.search(r'^version = "([^"]+)"', block, re.M)
        if name and version:
            out[name.group(1)] = version.group(1)
    return out


def npm_versions() -> dict[str, str]:
    """Every package in the frontend lockfile, by name."""
    data = json.loads((REPO / "frontend/package-lock.json").read_text())
    out: dict[str, str] = {}
    for path, entry in data.get("packages", {}).items():
        if not path.startswith("node_modules/"):
            continue
        version = entry.get("version")
        if version:
            out[path[len("node_modules/") :]] = version
    return out


def minor(version: str) -> tuple[str, str]:
    parts = version.split(".")
    return (parts[0], parts[1] if len(parts) > 1 else "0")


def main() -> int:
    crates, packages = crate_versions(), npm_versions()
    problems: list[str] = []
    checked = 0

    for crate, package in PAIRS:
        rust, node = crates.get(crate), packages.get(package)
        if rust is None or node is None:
            # A pair that is not in both lockfiles is not a mismatch; the shell
            # plugin, for one, may be present on only one side.
            continue
        checked += 1
        print(f"  {crate} {rust}  <->  {package} {node}")
        if minor(rust) != minor(node):
            problems.append(
                f"  {crate} is on {rust} and {package} on {node}: Tauri will "
                f"refuse to bundle, and only a release build would say so"
            )

    if not checked:
        print("no Tauri pair found in both lockfiles; nothing to compare", file=sys.stderr)
        return 1
    if problems:
        print("\nFAIL: the two halves of a dependency disagree.\n")
        print("\n".join(problems))
        return 1
    print(f"\nOK: {checked} Tauri pair(s), each on the same minor.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
