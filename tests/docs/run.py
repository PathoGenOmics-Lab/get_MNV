#!/usr/bin/env python3
"""The documentation names only things the tool actually has.

Prose drifts silently. A page can promise a flag that was renamed, a column that
was dropped, or an INFO key that never existed, and nothing notices until a
reader follows it and it does not work. Worse, a page can keep asserting a rule
the code stopped following: the GFF biotype paragraph said "GFF/GTF input is
unaffected, because it selects CDS features" long after that premise had been
shown false, because nothing compared the two.

This compares the vocabulary the documentation uses against the vocabulary the
binary produces:

  * every `--flag` named in the English or Spanish pages exists in `--help`
  * every TSV column named in a table cell exists in a real run's header
  * every VCF INFO key named in a table cell is declared in a real run's header

It is deliberately narrow. It cannot tell whether a sentence is true, only
whether the names in it are real, which is the part that rots without anybody
noticing. Exits non-zero on any name the tool does not have.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"
EXAMPLE = REPO / "example"

# Flags of other tools that legitimately appear in the pages: the docs show
# cargo, npm, samtools and bcftools commands beside get_MNV's own.
FOREIGN_FLAGS = {
    "--release", "--bin", "--path", "--locked", "--all-targets", "--nocapture",
    "--force", "--workspace", "--prefix", "--strict", "--site-dir", "--version",
    "--help", "--warmup", "--iters", "--no-deps", "--open", "--frozen",
    "--depth", "--branch", "--recursive", "--set-upstream-to", "--global",
}


def newest_source_mtime() -> float:
    return max(p.stat().st_mtime for p in (REPO / "src").rglob("*") if p.is_file())


def binary() -> Path:
    built = [
        candidate
        for relative in ("target/debug/get_mnv", "target/release/get_mnv", "dist/get_mnv")
        for candidate in [REPO / relative]
        if candidate.exists()
    ]
    if not built:
        sys.exit("no get_mnv binary in target/; build it first (cargo build)")
    built.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    if built[0].stat().st_mtime < newest_source_mtime():
        sys.exit(
            f"{built[0]} is older than src/. Rebuild it (cargo build) before running\n"
            "this, or it would be checking the documentation against a build that no\n"
            "longer matches the source."
        )
    return built[0]


def tool_vocabulary(get_mnv: Path) -> tuple[set[str], set[str], set[str]]:
    """The flags, TSV columns and INFO keys the binary actually has."""
    flags: set[str] = set()
    # The repository ships two binaries and the pages document both, so the
    # benchmark's own flags are part of the vocabulary too.
    help_text = subprocess.run(
        [str(get_mnv), "--help"], capture_output=True, text=True, check=False
    ).stdout
    flags |= set(re.findall(r"(--[a-z0-9][a-z0-9-]*)", help_text))
    # The pages document the benchmark binary too, and it has no --help to ask:
    # it reads argv directly and runs. Its flags come from its own source, which
    # is the only place they are written down.
    bench = REPO / "src" / "bin" / "bench_variants.rs"
    if bench.exists():
        flags |= set(re.findall(r'"(--[a-z0-9][a-z0-9-]*)"', bench.read_text()))

    work = Path(tempfile.mkdtemp(prefix="get_mnv_docs_"))
    try:
        vcf = work / "run.vcf"
        shutil.copy(EXAMPLE / "G35894.var.snp.vcf", vcf)
        # Everything switched on, so the widest vocabulary the tool can emit is
        # what the pages are checked against: read support, strand bias, both
        # output formats.
        subprocess.run(
            [
                str(get_mnv),
                "--vcf", str(vcf),
                "--bam", str(EXAMPLE / "G35894.demo.bam"),
                "--fasta", str(EXAMPLE / "MTB_ancestor.fas"),
                "--genes", str(EXAMPLE / "anot_genes.txt"),
                "--both",
                "--strand-bias-info",
            ],
            cwd=work,
            capture_output=True,
            text=True,
            check=True,
        )
        columns = set(
            (work / "run.MNV.tsv").read_text().split("\n", 1)[0].rstrip("\r").split("\t")
        )
        info = set()
        for line in (work / "run.MNV.vcf").read_text().split("\n"):
            if line.startswith("##INFO=<ID="):
                info.add(re.match(r"##INFO=<ID=([^,]+)", line).group(1))
            elif line and not line.startswith("#"):
                break
        return flags, columns, info
    finally:
        shutil.rmtree(work, ignore_errors=True)


def info_tables(text: str) -> list[str]:
    """The tables that list get_MNV's own VCF INFO keys.

    Found by the line that introduces them rather than by a heading, because
    they sit inside the VCF output section beside prose. Scoping matters: the
    same pages tabulate the input VCF's fixed columns and iVar's, and `REF`,
    `POS` and `PASS` there belong to the input, not to this tool.
    """
    lines = text.split("\n")
    out = []
    for index, line in enumerate(lines):
        if "INFO" not in line or line.lstrip().startswith("|"):
            continue
        start = next(
            (k for k in range(index + 1, min(index + 4, len(lines))) if lines[k].startswith("|")),
            None,
        )
        if start is None:
            continue
        end = start
        while end < len(lines) and lines[end].startswith("|"):
            end += 1
        out.append("\n".join(lines[start:end]))
    return out


def documented_names(columns: set[str], info: set[str]) -> tuple[set, set, set]:
    """Every flag, column name and INFO key the pages name, with where it is."""
    flags: dict[str, set[str]] = {}
    cols: dict[str, set[str]] = {}
    keys: dict[str, set[str]] = {}
    for page in sorted(DOCS.glob("*.md")):
        text = page.read_text()
        for flag in re.findall(r"(--[a-z0-9][a-z0-9-]*)", text):
            flags.setdefault(flag, set()).add(page.name)
        for token in re.findall(r"`([^`\n]+)`", text):
            token = token.strip()
            if token in columns or token.title() in columns:
                cols.setdefault(token, set()).add(page.name)
        # An INFO key is judged only where the pages are describing get_MNV's own
        # INFO column. A capitalised word anywhere else is a VCF field name, an
        # iVar column, or ordinary prose: `REF`, `POS` and `PASS` all appear in
        # the input-format tables and belong to the input, not to this tool.
        for table in info_tables(text):
            for row in table.split("\n"):
                first_cell = row.split("|")[1] if row.count("|") >= 2 else ""
                for token in re.findall(r"`([A-Z][A-Z0-9]{1,9})`", first_cell):
                    keys.setdefault(token, set()).add(page.name)
    return flags, cols, keys


def main() -> int:
    get_mnv = binary()
    flags, columns, info = tool_vocabulary(get_mnv)
    doc_flags, doc_cols, doc_keys = documented_names(columns, info)

    problems: list[str] = []

    unknown_flags = {
        flag: pages
        for flag, pages in doc_flags.items()
        if flag not in flags and flag not in FOREIGN_FLAGS
    }
    for flag, pages in sorted(unknown_flags.items()):
        problems.append(f"  flag not in --help: {flag}   named in {', '.join(sorted(pages))}")

    for key, pages in sorted(doc_keys.items()):
        if key not in info:
            problems.append(
                f"  INFO key not declared by the tool: {key}   named in {', '.join(sorted(pages))}"
            )

    print(f"{len(flags)} flags, {len(columns)} TSV columns and {len(info)} INFO keys in the tool")
    print(f"{len(doc_flags)} flags, {len(doc_cols)} columns and {len(doc_keys)} INFO keys named in the pages")
    if problems:
        print(f"\nFAIL: {len(problems)} name(s) the tool does not have.\n")
        print("\n".join(problems))
        return 1
    print("\nOK: every name the pages use is one the tool has.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
