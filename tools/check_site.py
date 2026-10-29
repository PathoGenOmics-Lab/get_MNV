#!/usr/bin/env python3
"""Check the documentation site the way `mkdocs build --strict` used to.

The pages are written by tools/build_docs.py and committed. What a strict build
caught, this catches, over the pages as a reader receives them:

  * every link and image inside the site points at a file that exists, and
    every `#anchor` at a heading or element that the target page has;
  * every markdown page in docs/ is in the sidebar, in both languages, and
    nothing in the sidebar is missing its page;
  * no MkDocs Material syntax is left in the markdown, where GitHub and this
    builder would both show it as text (`!!!`, `=== "`, `{ .class }`,
    `:material-icon:`, `markdown` attributes);
  * the two languages have the same headings structure page by page, so a
    section added in one is not forgotten in the other.

    python3 tools/check_site.py

Standard library only.
"""

from __future__ import annotations

import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
sys.path.insert(0, str(ROOT / "tools"))
import build_docs  # noqa: E402

LEFTOVERS = [
    (re.compile(r"^\s*(!!!|\?\?\?)\s"), "an MkDocs admonition (use a GitHub alert: > [!NOTE])"),
    (re.compile(r'^\s*=== "'), "an MkDocs content tab"),
    (re.compile(r"\{\s*[.#][\w-]+[^}]*\}"), "an attribute list ({ .class })"),
    (re.compile(r":(material|octicons|fontawesome)-[\w-]+:"), "an MkDocs icon shortcode"),
    (re.compile(r"\smarkdown(=\"1\")?>"), "a markdown attribute on raw HTML"),
    (re.compile(r"#only-(light|dark)\)"), "a Material theme fragment (use #gh-light-mode-only / #gh-dark-mode-only)"),
]


class Page(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.ids: set[str] = set()
        self.refs: list[tuple[str, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if values.get("id"):
            self.ids.add(values["id"])
        for key in ("href", "src", "srcset"):
            if values.get(key) and not (tag == "link" and values.get("rel") in ("canonical", "preconnect", "alternate")):
                self.refs.append((tag, values[key]))


def parse(path: Path) -> Page:
    page = Page()
    page.feed(path.read_text(encoding="utf-8"))
    return page


def link_problems(pages: dict[Path, Page]) -> list[str]:
    problems = []
    for path, page in pages.items():
        here = path.relative_to(DOCS)
        for tag, ref in page.refs:
            parts = urlsplit(ref)
            if parts.scheme or ref.startswith(("//", "mailto:", "data:")):
                continue
            target = (path.parent / unquote(parts.path)).resolve() if parts.path else path.resolve()
            if not target.exists():
                problems.append(f"  {here}: <{tag}> points at {ref}, which does not exist")
                continue
            if parts.fragment and target.suffix == ".html":
                target_page = pages.get(target) or parse(target)
                if unquote(parts.fragment) not in target_page.ids:
                    problems.append(f"  {here}: {ref} names an anchor its page does not have")
    return problems


def sidebar_problems() -> list[str]:
    problems = []
    listed = {source for _, source, _ in build_docs.PAGES if source.endswith(".md")}
    for md in sorted(DOCS.glob("*.md")):
        english = md.name[: -len(".es.md")] + ".md" if md.name.endswith(".es.md") else md.name
        if english not in listed:
            problems.append(f"  docs/{md.name} is not in the sidebar (PAGES in tools/build_docs.py)")
    for source in sorted(listed):
        for lang in ("en", "es"):
            if not (DOCS / build_docs.localized(source, lang)).exists():
                problems.append(f"  the sidebar lists {source} but docs/{build_docs.localized(source, lang)} does not exist")
    for _, source, _ in build_docs.PAGES:
        if not source.endswith(".md") and not source.startswith("http") and not (DOCS / source).exists():
            problems.append(f"  the sidebar links docs/{source}, which does not exist")
    return problems


def syntax_problems() -> list[str]:
    problems = []
    for md in sorted(DOCS.glob("*.md")):
        fence = False
        for number, line in enumerate(md.read_text(encoding="utf-8").splitlines(), 1):
            if re.match(r"^\s*(```|~~~)", line):
                fence = not fence
                continue
            if fence:
                continue
            for pattern, what in LEFTOVERS:
                if pattern.search(line):
                    problems.append(f"  docs/{md.name}:{number}: {what}")
    return problems


def twin_problems() -> list[str]:
    problems = []
    for _, source, _ in build_docs.PAGES:
        if not source.endswith(".md"):
            continue
        levels = {}
        for lang in ("en", "es"):
            text = (DOCS / build_docs.localized(source, lang)).read_text(encoding="utf-8")
            text = re.sub(r"(?ms)^(```|~~~).*?^\1", "", text)
            levels[lang] = [len(m.group(1)) for m in re.finditer(r"(?m)^(#{1,4})\s", text)]
        if levels["en"] != levels["es"]:
            problems.append(
                f"  {source}: the English page has headings {levels['en']} and the Spanish one {levels['es']}"
            )
    return problems


def main() -> int:
    # 404.html links absolutely, to the published site, so it is not followed here.
    generated = sorted(p for p in DOCS.rglob("*.html") if p.name == "index.html")
    pages = {p.resolve(): parse(p) for p in generated}
    checks = [
        ("links and anchors", link_problems(pages)),
        ("sidebar", sidebar_problems()),
        ("leftover MkDocs syntax", syntax_problems()),
        ("the two languages", twin_problems()),
    ]
    failed = False
    for name, problems in checks:
        if problems:
            failed = True
            print(f"{name}:")
            print("\n".join(problems))
    if failed:
        return 1
    refs = sum(len(p.refs) for p in pages.values())
    print(f"OK: {len(pages)} pages, {refs} links and images, every one resolves.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
