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
  * every JSON key tabulated under "JSON files" is in a real run's JSON, for
    both the single-sample and the `--sample all` shapes
  * every page is reachable from the nav and has a twin in the other language
  * the light and dark wordmarks are the same drawing, differing only in ink
  * the light and dark diagrams cover exactly the same pixels
  * every screenshot the pages show has a dark twin of the same size

It is deliberately narrow. It cannot tell whether a sentence is true, only
whether the names in it are real, which is the part that rots without anybody
noticing. Exits non-zero on any name the tool does not have.
"""

from __future__ import annotations

import json
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
# The lines that introduce a table of TSV column names, in both languages.
# The headings that open the JSON reference, in both languages. Everything
# tabulated under them names a key of a payload the tool writes.
JSON_SECTION_HEADINGS = ("## JSON Files", "## Archivos JSON")

COLUMN_TABLE_MARKERS = (
    "Main columns",
    "Extra columns when",
    "Columnas principales",
    "Columnas extra",
    "Columnas adicionales cuando",
)

# A flag begins a token: the two dashes are not preceded by a word character or
# by another dash. Without that anchor `.md-button--primary`, a CSS class in a
# Markdown attribute list, is read as a flag named `--primary`, and the pages
# get blamed for naming something the tool does not have.
FLAG_PATTERN = r"(?<![\w-])(--[a-z0-9][a-z0-9-]*)"

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
    flags |= set(re.findall(FLAG_PATTERN, help_text))
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


def json_key_universe(get_mnv: Path) -> set[str]:
    """Every key name the tool puts in a JSON payload.

    Runs the tool four times because the payloads differ: a single-sample run
    writes one shape, `--sample all` writes another with `aggregate` and
    `samples` in place of the run-wide keys, and a failed run writes the error
    payload. Documenting the second shape is what this is for: the pages tell a
    reader that `data["global"]` becomes `data["aggregate"]["global"]`, and a
    reader who follows that gets a KeyError if it ever stops being true.

    Compares key names rather than paths. A key that moves between nesting
    levels is a real change this will not catch, but a key that is renamed or
    dropped is the way these payloads actually rot.
    """
    names: set[str] = set()

    def collect(value) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                names.add(key)
                collect(child)
        elif isinstance(value, list):
            for child in value:
                collect(child)

    work = Path(tempfile.mkdtemp(prefix="get_mnv_docs_json_"))
    try:
        vcf = work / "run.vcf"
        shutil.copy(EXAMPLE / "G35894.var.snp.vcf", vcf)
        common = [
            "--fasta", str(EXAMPLE / "MTB_ancestor.fas"),
            "--genes", str(EXAMPLE / "anot_genes.txt"),
        ]
        subprocess.run(
            [str(get_mnv), "--vcf", str(vcf), *common,
             "--summary-json", "s.json", "--run-manifest", "m.json"],
            cwd=work, capture_output=True, text=True, check=True,
        )
        # A two-sample VCF built from the example, so the `--sample all` shape
        # is read off a real run rather than trusted from the prose.
        header, records = [], []
        for line in vcf.read_text().split("\n"):
            if line.startswith("##"):
                header.append(line)
            elif line.startswith("#CHROM"):
                header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
                header.append("\t".join(line.split("\t")[:8] + ["FORMAT", "SAMPLE_A", "SAMPLE_B"]))
            elif line.strip() and len(records) < 8:
                records.append("\t".join(line.split("\t")[:8] + ["GT", "1/1", "1/1"]))
        multi = work / "multi.vcf"
        multi.write_text("\n".join(header + records) + "\n")
        subprocess.run(
            [str(get_mnv), "--vcf", str(multi), *common, "--sample", "all",
             "--summary-json", "sall.json", "--run-manifest", "mall.json"],
            cwd=work, capture_output=True, text=True, check=True,
        )
        # And a run that fails, because the error payload is only written then.
        subprocess.run(
            [str(get_mnv), "--vcf", str(work / "absent.vcf"), *common,
             "--error-json", "e.json"],
            cwd=work, capture_output=True, text=True, check=False,
        )
        for name in ("s.json", "m.json", "sall.json", "mall.json", "e.json"):
            path = work / name
            if path.exists():
                collect(json.loads(path.read_text()))
        return names
    finally:
        shutil.rmtree(work, ignore_errors=True)


def tables_introduced_by(text: str, needles: tuple[str, ...]) -> list[str]:
    """The tables introduced by a line naming one of `needles`.

    Scoping by the line above the table rather than by a heading, because these
    tables sit inside a section beside prose, and because the same pages also
    tabulate the input formats' own columns, which are not this tool's to
    produce.
    """
    lines = text.split("\n")
    out = []
    for index, line in enumerate(lines):
        if line.lstrip().startswith("|") or not any(n in line for n in needles):
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


def json_section(text: str) -> str:
    """The part of a page that documents the JSON payloads, or nothing."""
    for heading in JSON_SECTION_HEADINGS:
        if heading in text:
            rest = text.split(heading, 1)[1]
            end = re.search(r"^## ", rest, re.M)
            return rest[: end.start()] if end else rest
    return ""


def documented_names(columns: set[str], info: set[str]) -> tuple[set, set, set, set]:
    """Every flag, column name, INFO key and JSON key the pages name."""
    flags: dict[str, set[str]] = {}
    cols: dict[str, set[str]] = {}
    keys: dict[str, set[str]] = {}
    json_keys: dict[str, set[str]] = {}
    for page in sorted(DOCS.glob("*.md")):
        text = page.read_text()
        for flag in re.findall(FLAG_PATTERN, text):
            flags.setdefault(flag, set()).add(page.name)
        # The tables that write down the TSV schema. Every first cell in them is
        # a column name, so a name that is not in the real header is a promise
        # the tool does not keep: a reader indexing that column gets a KeyError
        # from the one page meant to prevent exactly that. Collecting only names
        # that were already columns, which is what this did before, made the
        # check unable to fail.
        for table in tables_introduced_by(text, COLUMN_TABLE_MARKERS):
            for row in table.split("\n"):
                first_cell = row.split("|")[1] if row.count("|") >= 2 else ""
                for token in re.findall(r"`([^`]+)`", first_cell):
                    cols.setdefault(token.strip(), set()).add(page.name)
        # An INFO key is judged only where the pages are describing get_MNV's own
        # INFO column. A capitalised word anywhere else is a VCF field name, an
        # iVar column, or ordinary prose: `REF`, `POS` and `PASS` all appear in
        # the input-format tables and belong to the input, not to this tool.
        for table in info_tables(text):
            for row in table.split("\n"):
                first_cell = row.split("|")[1] if row.count("|") >= 2 else ""
                for token in re.findall(r"`([A-Z][A-Z0-9]{1,9})`", first_cell):
                    keys.setdefault(token, set()).add(page.name)
        # The JSON reference tabulates one key per row, sometimes several in a
        # row when they share a description. Only first cells, and only inside
        # that section: `null`, `sample_all` and the exit codes live in the
        # description column and in prose, and none of them is a key.
        for row in json_section(text).split("\n"):
            if not row.startswith("|") or row.count("|") < 2:
                continue
            for token in re.findall(r"`([a-z][a-z0-9_]{2,})`", row.split("|")[1]):
                json_keys.setdefault(token, set()).add(page.name)
    return flags, cols, keys, json_keys


def nav_pages() -> set[str]:
    """The pages the nav actually links to.

    Reads the `nav:` tree rather than scanning the file for `.md`, because the
    file mentions pages outside the nav: a comment naming `linkage.md` was
    enough to make a scan of the whole file report every page as reachable,
    which is a check that cannot fail. Loaded with a loader that tolerates
    mkdocs' `!!python/name:` tags, which `safe_load` refuses.
    """
    import yaml

    class Loader(yaml.SafeLoader):
        pass

    Loader.add_multi_constructor("tag:yaml.org,2002:python/name:", lambda *_: None)
    Loader.add_multi_constructor("!", lambda *_: None)

    config = yaml.load((REPO / "mkdocs.yml").read_text(), Loader=Loader)
    pages: set[str] = set()

    def walk(node) -> None:
        if isinstance(node, str):
            pages.add(node)
        elif isinstance(node, list):
            for child in node:
                walk(child)
        elif isinstance(node, dict):
            for child in node.values():
                walk(child)

    walk(config.get("nav") or [])
    return {page for page in pages if page.endswith(".md")}


def wordmark_problems() -> list[str]:
    """The two wordmarks are one drawing in two inks, and stay that way.

    `logo.svg` and `logo-dark.svg` exist as separate files because a mark that
    decides its own colour from `prefers-color-scheme` follows the reader's
    system rather than the background it was put on, and paints itself white on
    white the moment a site with its own toggle disagrees. The cost of two files
    is that they can drift: someone re-exports one from the drawing program and
    the other keeps yesterday's artwork. This compares them with their
    stylesheets removed, so the only difference allowed is the ink.
    """
    problems: list[str] = []
    light = DOCS / "assets" / "logo.svg"
    dark = DOCS / "assets" / "logo-dark.svg"
    for path in (light, dark):
        if not path.exists():
            problems.append(f"  missing wordmark: {path.name}")
    if problems:
        return problems

    def drawing(path: Path) -> str:
        """Everything but the stylesheet: the shapes themselves."""
        return re.sub(r"<style[^>]*>.*?</style>", "", path.read_text(), flags=re.S)

    if drawing(light) != drawing(dark):
        problems.append(
            "  logo.svg and logo-dark.svg are no longer the same drawing; "
            "re-cut the dark one from the light one rather than editing it apart"
        )

    # Looking for "fill:" anywhere in the stylesheet would pass on the drawing's
    # own twenty fills and could never fail. The rule has to be the one that
    # repaints the five near-black shapes, and it has to repaint them light.
    ink_rule = re.search(
        r"([^{}]*#text2924[^{}]*)\{([^}]*)\}", dark.read_text()
    )
    if not ink_rule:
        problems.append(
            "  logo-dark.svg has no rule for the ink shapes, so it is the light mark under another name"
        )
    else:
        selector, body = ink_rule.groups()
        missing = [i for i in ("text2924", "path3", "rect2936", "rect2938", "rect2940")
                   if i not in selector]
        if missing:
            problems.append(
                f"  logo-dark.svg leaves {', '.join(missing)} in the light ink, so part of the mark "
                "disappears on a dark ground"
            )
        colour = re.search(r"fill:\s*#([0-9a-fA-F]{6})", body)
        if not colour:
            problems.append("  logo-dark.svg's ink rule sets no fill")
        else:
            r, g, b = (int(colour.group(1)[i:i+2], 16) / 255 for i in (0, 2, 4))
            lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
            if lum < 0.5:
                problems.append(
                    f"  logo-dark.svg paints its ink #{colour.group(1)}, which is dark: "
                    "on a dark ground the mark loses its lettering"
                )
    return problems


def artwork_problems() -> list[str]:
    """The diagram's two variants are one drawing in two inks.

    `get_mnv_aa.png` is line art on transparent paper; `get_mnv_aa-dark.png` is
    the same art with the near-black ink inverted and the two saturated inks
    replaced, because pure blue on a dark ground is 1.9:1 and inverting does not
    fix that. Their alpha channels are the shape of the drawing, so comparing
    those catches the case where one is regenerated and the other is not.
    """
    from PIL import Image

    problems: list[str] = []
    light = DOCS / "assets" / "get_mnv_aa.png"
    dark = DOCS / "assets" / "get_mnv_aa-dark.png"
    for path in (light, dark):
        if not path.exists():
            problems.append(f"  missing diagram: {path.name}")
    if problems:
        return problems

    a, b = (Image.open(p).convert("RGBA") for p in (light, dark))
    if a.size != b.size:
        problems.append(f"  the diagrams are different sizes: {a.size} and {b.size}")
        return problems
    if a.getchannel("A").tobytes() != b.getchannel("A").tobytes():
        problems.append(
            "  the light and dark diagrams no longer cover the same pixels; "
            "re-cut the dark one from the light one"
        )

    # Screenshots are photographs of the application, so their two versions are
    # different renderings rather than one drawing in two inks: only their size
    # can be compared, and a missing twin is the failure that actually happens.
    for shot in sorted(DOCS.glob("assets/*.png")):
        if shot.stem.endswith("-dark") or not shot.stem.startswith(("gui-", "cli-")):
            continue
        twin = shot.with_name(f"{shot.stem}-dark.png")
        if not twin.exists():
            problems.append(f"  {shot.name} has no dark twin, so the dark pages show a light screenshot")
            continue
        one, other = Image.open(shot).size, Image.open(twin).size
        if one != other:
            problems.append(
                f"  {shot.name} and its dark twin are different sizes, {one} and {other}, "
                "so the page reflows when the theme changes"
            )

    # And the dark one has to be light: the whole point is that its ink shows on
    # a dark page. Measured on the pixels that are actually drawn.
    for path, image, wanted in ((light, a, "dark"), (dark, b, "light")):
        pixels = [p for p in image.getdata() if p[3] > 200]
        if not pixels:
            problems.append(f"  {path.name} has no opaque ink at all")
            continue
        mid = sorted(sum(p[:3]) / 3 for p in pixels)[len(pixels) // 2]
        if wanted == "light" and mid < 128:
            problems.append(f"  {path.name} is inked dark, so it vanishes on the dark page it is for")
        if wanted == "dark" and mid > 128:
            problems.append(f"  {path.name} is inked light, so it vanishes on the light page it is for")
    return problems


def page_structure_problems() -> list[str]:
    """No page is unreachable, and no page exists in one language only.

    A page that nothing links to is a page nobody reads, and a page with no twin
    is a reader in the other language hitting a wall. Both are invisible to a
    build: mkdocs is happy to ship an orphan.
    """
    problems: list[str] = []
    nav = nav_pages()
    for page in sorted(DOCS.glob("*.md")):
        if page.name.endswith(".es.md"):
            english = DOCS / page.name.replace(".es.md", ".md")
            if not english.exists():
                problems.append(f"  Spanish page with no English original: {page.name}")
            continue
        if page.name not in nav:
            problems.append(f"  page missing from the nav in mkdocs.yml: {page.name}")
        spanish = DOCS / page.name.replace(".md", ".es.md")
        if not spanish.exists():
            problems.append(f"  page with no Spanish translation: {page.name}")
    return problems


def main() -> int:
    get_mnv = binary()
    flags, columns, info = tool_vocabulary(get_mnv)
    json_keys = json_key_universe(get_mnv)
    doc_flags, doc_cols, doc_keys, doc_json = documented_names(columns, info)

    problems: list[str] = []

    unknown_flags = {
        flag: pages
        for flag, pages in doc_flags.items()
        if flag not in flags and flag not in FOREIGN_FLAGS
    }
    for flag, pages in sorted(unknown_flags.items()):
        problems.append(f"  flag not in --help: {flag}   named in {', '.join(sorted(pages))}")

    for column, pages in sorted(doc_cols.items()):
        if column not in columns:
            problems.append(
                f"  TSV column not in the header: {column}   named in {', '.join(sorted(pages))}"
            )

    for key, pages in sorted(doc_keys.items()):
        if key not in info:
            problems.append(
                f"  INFO key not declared by the tool: {key}   named in {', '.join(sorted(pages))}"
            )

    for key, pages in sorted(doc_json.items()):
        if key not in json_keys:
            problems.append(
                f"  JSON key not in any payload: {key}   named in {', '.join(sorted(pages))}"
            )

    problems.extend(page_structure_problems())
    problems.extend(wordmark_problems())
    problems.extend(artwork_problems())

    print(
        f"{len(flags)} flags, {len(columns)} TSV columns, {len(info)} INFO keys "
        f"and {len(json_keys)} JSON keys in the tool"
    )
    print(
        f"{len(doc_flags)} flags, {len(doc_cols)} columns, {len(doc_keys)} INFO keys "
        f"and {len(doc_json)} JSON keys named in the pages"
    )
    if problems:
        print(f"\nFAIL: {len(problems)} name(s) the tool does not have.\n")
        print("\n".join(problems))
        return 1
    print("\nOK: every name the pages use is one the tool has.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
