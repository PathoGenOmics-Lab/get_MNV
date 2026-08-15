#!/usr/bin/env python3
"""Packages that constrain each other end up in the same Dependabot group.

Two npm packages can be joined by a `peerDependencies` range, and when they are,
they have to move together. Dependabot does not know that. It groups by whatever
`.github/dependabot.yml` says, and grouping the frontend by update type alone
split one such pair down the middle:

    eslint 9.39.4 -> 10.8.1                     went to npm-major   (a major)
    eslint-plugin-react-hooks 7.0.1 -> 7.1.1    went to npm-routine (a patch)

The patch was the half that widened the plugin's peer range to accept eslint 10.
So the major pull request died at `npm ci` with ERESOLVE against a range only the
other pull request fixed, and neither could be judged without the other.

The fix was to give each coupled family its own group, listed before the
catch-alls. What this guards is that the fix stays fixed, because the failure
mode is quiet: the groups are matched in file order, so moving a family group
below `npm-routine` silently restores the original bug, and adding a new
dependency that peers on an existing one silently creates a new split. Neither
shows up until a pull request fails a month later.

The families are read from `frontend/package-lock.json` rather than listed here.
A list would be the same fact written in two places, and would go stale the first
time a dependency gained a peer range.

    python3 tools/check_dependency_groups.py
"""

from __future__ import annotations

import json
import sys
from fnmatch import fnmatch
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
FRONTEND = REPO / "frontend"
CONFIG = REPO / ".github" / "dependabot.yml"

# The update types Dependabot can put on a group. A family has to land in one
# group for all three, since the constraint is between the packages and has
# nothing to do with how big either version bump happens to be.
UPDATE_TYPES = ("patch", "minor", "major")


def direct_dependencies() -> set[str]:
    """The packages Dependabot will actually open pull requests for."""
    manifest = json.loads((FRONTEND / "package.json").read_text())
    return set(manifest.get("dependencies", {})) | set(manifest.get("devDependencies", {}))


def peer_edges(direct: set[str]) -> list[tuple[str, str]]:
    """Every `a peers on b` where both are direct dependencies.

    Read from the lockfile, so this needs no network and describes the versions
    that would actually be installed.
    """
    lock = json.loads((FRONTEND / "package-lock.json").read_text())
    edges: list[tuple[str, str]] = []
    for path, entry in lock.get("packages", {}).items():
        if not path.startswith("node_modules/"):
            continue
        name = path[len("node_modules/") :]
        if name not in direct:
            continue
        for peer in entry.get("peerDependencies", {}):
            if peer in direct:
                edges.append((name, peer))
    return edges


def families(direct: set[str], edges: list[tuple[str, str]]) -> list[set[str]]:
    """Connected components of the peer graph, largest first.

    Components rather than pairs, because the couplings chain: typescript-eslint
    peers on both eslint and typescript, which puts all three in one family even
    though eslint and typescript constrain each other not at all.
    """
    parent = {name: name for name in direct}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    groups: dict[str, set[str]] = {}
    for name in direct:
        groups.setdefault(find(name), set()).add(name)
    return sorted((f for f in groups.values() if len(f) > 1), key=len, reverse=True)


def npm_groups() -> dict[str, dict]:
    """The npm groups, in the order the file lists them.

    The order is the load-bearing part: Dependabot puts a dependency in the first
    group it matches, so a family group below a `patterns: ["*"]` catch-all never
    gets to claim anything.
    """
    config = yaml.safe_load(CONFIG.read_text())
    npm = [u for u in config["updates"] if u["package-ecosystem"] == "npm"]
    if len(npm) != 1:
        sys.exit(f"expected exactly one npm entry in {CONFIG.name}, found {len(npm)}")
    return npm[0].get("groups") or {}


def group_for(name: str, update_type: str, groups: dict[str, dict]) -> str | None:
    """The group this package joins for this kind of bump, or None if ungrouped."""
    for group, spec in groups.items():
        if not any(fnmatch(name, pattern) for pattern in spec.get("patterns", [])):
            continue
        allowed = spec.get("update-types")
        if allowed and update_type not in allowed:
            continue
        return group
    return None


def main() -> int:
    direct = direct_dependencies()
    found = families(direct, peer_edges(direct))
    groups = npm_groups()

    print(f"{len(direct)} direct npm dependencies, {len(found)} coupled families")
    print(f"groups, in the order they are matched: {', '.join(groups) or '(none)'}\n")

    problems: list[str] = []
    for family in found:
        landings = {
            member: {t: group_for(member, t, groups) for t in UPDATE_TYPES}
            for member in sorted(family)
        }
        distinct = {g for per_type in landings.values() for g in per_type.values()}

        if distinct == {None}:
            problems.append(
                f"  {', '.join(sorted(family))}\n"
                f"    constrain each other but match no group at all, so each one\n"
                f"    would arrive in a pull request of its own"
            )
        elif len(distinct) > 1:
            detail = "\n".join(
                f"      {member:<32} "
                + ", ".join(f"{t}: {g or 'ungrouped'}" for t, g in per_type.items())
                for member, per_type in landings.items()
            )
            problems.append(
                f"  {', '.join(sorted(family))}\n"
                f"    constrain each other but can land in {len(distinct)} different\n"
                f"    groups, so one half can be updated without the other:\n{detail}"
            )
        else:
            print(f"  OK  {distinct.pop()}: {', '.join(sorted(family))}")

    if problems:
        print("\nFAIL: a package can move without the package that constrains it.\n")
        print("\n\n".join(problems))
        print(
            "\nGive the family one group in .github/dependabot.yml, listed before\n"
            "any group with a catch-all pattern, and with no update-types filter."
        )
        return 1

    if not found:
        print("no peer couplings between direct dependencies; nothing to keep together")
    print("\nOK: every coupled family travels as one pull request.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
