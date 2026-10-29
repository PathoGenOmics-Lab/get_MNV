#!/usr/bin/env python3
"""Build the documentation site from the markdown in docs/.

The markdown stays the source: GitHub renders it, and this writes one HTML page
per markdown file with the site's header, theme, sidebar, page contents and
search, plus the search index. The home page is the exception: its figure and
codon playground do not fit in markdown, so it is written from
tools/site/home.en.html and home.es.html, inside the same header and footer. Pages land where the old MkDocs site served
them, `docs/<page>/index.html` and `docs/es/<page>/index.html`, so every link
already out in the world still resolves. The pages are committed, so the site
opens from a clone with no server and no build step; --check fails when they
no longer match the markdown, which is what CI runs.

    python3 tools/build_docs.py           # write the pages
    python3 tools/build_docs.py --check   # exit 1 if any page is stale

Only the markdown these pages use is supported, read the way GitHub reads it:
headings, paragraphs, lists (nested, loose or tight, with code inside), tables,
block quotes and GitHub's alerts (> [!NOTE], titled by a first line that is
bold and nothing else), fenced code, rules, images (a `#gh-light-mode-only` /
`#gh-dark-mode-only` fragment shows the image in that theme only, as GitHub
does), inline code, backslash escapes, emphasis, links, URLs written out, a few
inline tags and raw HTML blocks. Reference links, footnotes and task lists are
not, and tools/check_site.py stops on them; a table row with more cells than
its header stops the build.

Standard library only.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import re
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
REPO = "https://github.com/PathoGenOmics-Lab/get_MNV"
SITE_URL = "https://pathogenomics-lab.github.io/get_MNV/"
ASK = f"{REPO}/discussions/categories/q-a"

# The sidebar, in reading order: (section, source, title). A source ending in
# .md is built into a page; anything else is linked as it is. Each language
# has its own titles; the Spanish page of `x.md` is `x.es.md`.
PAGES = [
    # The home page is written from tools/site/home.<lang>.html, not markdown.
    ("Home", "index.html", "Home"),
    ("Get started", "getting-started.md", "Command Line Tutorial"),
    ("Get started", "usage.md", "Common Recipes"),
    ("Get started", "gui.md", "Desktop GUI"),
    ("Get started", "gui-tutorial.md", "Desktop GUI Tutorial"),
    ("Reference", "cli-reference.md", "CLI Reference"),
    ("Reference", "input-formats.md", "Input Formats"),
    ("Reference", "output-formats.md", "Output Formats"),
    ("Reference", "assets/example-report.html", "Example report"),
    ("How it works", "indel-mnv-semantics.md", "Scope and Compatibility"),
    ("How it works", "indel-haplotypes.md", "Indels and Local Haplotypes"),
    ("How it works", "linkage.md", "Linkage"),
    ("How it works", "benchmarking.md", "Benchmarking"),
    ("Help", "glossary.md", "Glossary"),
    ("Help", "troubleshooting.md", "Troubleshooting"),
    ("Help", "faq.md", "FAQ"),
    # The only entry that leaves the site, and it is deliberate: someone reading
    # Help has already failed to find their answer, and the place to ask should
    # not need finding.
    ("Help", ASK, "Ask a question"),
]

SPANISH = {
    "Home": "Inicio",
    "Get started": "Empezar",
    "Reference": "Referencia",
    "How it works": "Cómo funciona",
    "Help": "Ayuda",
    "Command Line Tutorial": "Tutorial de línea de comandos",
    "Common Recipes": "Recetas habituales",
    "Desktop GUI": "GUI de escritorio",
    "Desktop GUI Tutorial": "Tutorial de la GUI",
    "CLI Reference": "Referencia de CLI",
    "Input Formats": "Formatos de entrada",
    "Output Formats": "Formatos de salida",
    "Example report": "Informe de ejemplo",
    "Scope and Compatibility": "Alcance y compatibilidad",
    "Indels and Local Haplotypes": "Indels y haplotipos locales",
    "Linkage": "Ligamiento",
    "Benchmarking": "Rendimiento",
    "Glossary": "Glosario",
    "Troubleshooting": "Resolución de problemas",
    "FAQ": "Preguntas frecuentes",
    "Ask a question": "Pregunta algo",
}

# Interface strings, per language.
UI = {
    "en": {
        "lang_name": "English",
        "other": "es",
        "other_label": "Español",
        "search": "Search",
        "search_label": "Search the documentation",
        "site_nav": "Main",
        "start": "Start here",
        "try": "Try it",
        "reference": "Reference",
        "report": "Example report",
        "home_title": "get_MNV: codon-level variant annotation",
        "on_page": "On this page",
        "prev": "Previous",
        "next": "Next",
        "edit": "Edit this page on GitHub",
        "theme": "Switch between light and dark",
        "skip": "Skip to the content",
        "menu": "Pages",
        "docs_label": "Documentation",
        "description": "Codon-level Multi-Nucleotide Variant detection and indel annotation from VCF or iVar TSV. Pure Rust, no C dependencies.",
        "cite": "Cite",
        "license": "AGPL-3.0 · PathoGenOmics Lab",
        "project": "Project",
        "issues": "Report a bug",
        "ask": "Ask a question",
        "releases": "Releases",
        "notes": {"NOTE": "Note", "TIP": "Tip", "IMPORTANT": "Important", "WARNING": "Warning", "CAUTION": "Caution"},
    },
    "es": {
        "lang_name": "Español",
        "other": "en",
        "other_label": "English",
        "search": "Buscar",
        "search_label": "Buscar en la documentación",
        "site_nav": "Principal",
        "start": "Empieza aquí",
        "try": "Pruébalo",
        "reference": "Referencia",
        "report": "Informe de ejemplo",
        "home_title": "get_MNV: anotación de variantes a nivel de codón",
        "on_page": "En esta página",
        "prev": "Anterior",
        "next": "Siguiente",
        "edit": "Editar esta página en GitHub",
        "theme": "Cambiar entre claro y oscuro",
        "skip": "Saltar al contenido",
        "menu": "Páginas",
        "docs_label": "Documentación",
        "description": "Detección de variantes multinucleotídicas a nivel de codón y anotación de indels desde VCF o TSV de iVar. Rust puro, sin dependencias de C.",
        "cite": "Cita",
        "license": "AGPL-3.0 · PathoGenOmics Lab",
        "project": "Proyecto",
        "issues": "Informar de un fallo",
        "ask": "Pregunta algo",
        "releases": "Versiones",
        "notes": {"NOTE": "Nota", "TIP": "Consejo", "IMPORTANT": "Importante", "WARNING": "Aviso", "CAUTION": "Cuidado"},
    },
}

SEARCH_INDEX = DOCS / "assets" / "search-index.js"

GH_PATH = (
    "M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94"
    "-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07"
    "-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0"
    " 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25"
    ".54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"
)
GH_SVG = f'<svg viewBox="0 0 16 16" fill="currentColor" aria-hidden="true"><path d="{GH_PATH}"/></svg>'
SEARCH_SVG = (
    '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6" aria-hidden="true">'
    '<circle cx="7" cy="7" r="4.6"/><path d="m10.4 10.4 3.6 3.6" stroke-linecap="round"/></svg>'
)
THEME_SVG = (
    '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.4" aria-hidden="true">'
    '<circle cx="8" cy="8" r="5.9"/><path d="M8 2.1a5.9 5.9 0 0 0 0 11.8z" fill="currentColor" stroke="none"/></svg>'
)
# Before the stylesheet paints anything: the reader's own choice if they made
# one, otherwise dark, as the lab's other sites open.
THEME_BOOT = (
    "<script>(function(){var t;try{t=localStorage.getItem(\"get_mnv.theme\")}catch(e){}"
    "if(t!==\"light\"&&t!==\"dark\"){t=\"dark\"}"
    "document.documentElement.dataset.theme=t})();</script>"
)
FONTS = (
    '<link rel="preconnect" href="https://fonts.googleapis.com">\n'
    '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n'
    '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&'
    'family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&display=swap">\n'
)
INLINE_TAGS = re.compile(r"</?(kbd|br|sup|sub|strong|em|b|i|code)\s*/?>", re.I)


# ── Paths ──


def localized(source: str, lang: str) -> str:
    """The markdown file for `source` in `lang`."""
    return source if lang == "en" else source[: -len(".md")] + ".es.md"


def output_of(source: str, lang: str) -> str:
    """Where the page for an English `source` is written, relative to docs/."""
    stem = source.rsplit(".", 1)[0]
    prefix = "" if lang == "en" else "es/"
    return f"{prefix}index.html" if stem == "index" else f"{prefix}{stem}/index.html"


def relative(from_page: str, to_path: str) -> str:
    """A link from the page at `from_page` to `to_path`, both relative to docs/."""
    return os.path.relpath(to_path, os.path.dirname(from_page) or ".").replace(os.sep, "/")


def title_in(title: str, lang: str) -> str:
    return title if lang == "en" else SPANISH.get(title, title)


# ── Markdown ──


class DocsError(Exception):
    """Markdown the site cannot render the way GitHub does."""


ASCII_PUNCTUATION = set("!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~")

FENCE = re.compile(r"^(\s*)(```+|~~~+)\s*(.*)$")
# A closing run of #s only counts after a space, so `## C#` keeps its #.
HEADING = re.compile(r"^ {0,3}(#{1,6})[ \t]+(.*?)(?:[ \t]+#+)?[ \t]*$")
SETEXT = re.compile(r"^ {0,3}(=+|-+)[ \t]*$")
RULE = re.compile(r"^ {0,3}((\*[ \t]*){3,}|(-[ \t]*){3,}|(_[ \t]*){3,})$")
LIST_ITEM = re.compile(r"^(\s*)([-*+]|\d+\.)(\s+)")
# The tags that open or close an HTML block, which runs to the next blank line.
HTML_BLOCK = re.compile(
    r"^</?(address|article|aside|blockquote|caption|center|dd|details|dialog|div|dl|dt|fieldset"
    r"|figcaption|figure|footer|form|h[1-6]|header|hr|li|main|nav|ol|p|picture|section|summary"
    r"|table|tbody|td|tfoot|th|thead|tr|ul)(\s|/?>|$)",
    re.I,
)
TABLE_DELIMITER = re.compile(r"^\|?\s*:?-+:?\s*(\|\s*:?-+:?\s*)*\|?$")
# A link's destination may hold one level of parentheses, as Wikipedia's do
# (Codon_(biology)), and may be followed by a "title".
DESTINATION = r"((?:[^()\s]|\([^()\s]*\))+)"
TITLE = r'(?:\s+"([^"]*)")?'
IMAGE = re.compile(r"!\[([^\]]*)\]\(\s*" + DESTINATION + TITLE + r"\s*\)")
LINK = re.compile(r"\[([^\]]+)\]\(\s*" + DESTINATION + TITLE + r"\s*\)")
# A URL written out in the text, which GitHub links.
BARE_URL = re.compile(r"(?<![\w/@])https?://[^\s<\x00-\x03]+")


def code_and_escapes(text: str, code, escaped) -> str:
    """Reads code spans and backslash escapes left to right, as CommonMark
    does: a backslash makes any ASCII punctuation literal, a run of backticks
    opens a code span only where a run of the same length closes it, and
    nothing inside a code span is markup. `code` and `escaped` say what each
    becomes."""
    out: list[str] = []
    i, n = 0, len(text)
    while i < n:
        char = text[i]
        if char == "\\" and i + 1 < n and text[i + 1] in ASCII_PUNCTUATION:
            out.append(escaped(text[i + 1]))
            i += 2
        elif char == "`":
            opening = re.match(r"`+", text[i:]).group(0)
            closing = re.compile(rf"(?<!`){opening}(?!`)").search(text, i + len(opening))
            if closing is None:
                out.append(opening)
                i += len(opening)
                continue
            content = text[i + len(opening) : closing.start()].replace("\n", " ")
            # One space comes off each end when both ends have one, which is
            # how a span can start with a backtick; any other space is the
            # code's own, as in `; `.
            if content.startswith(" ") and content.endswith(" ") and content.strip(" "):
                content = content[1:-1]
            out.append(code(content))
            i = closing.end()
        else:
            out.append(char)
            i += 1
    return "".join(out)


def table_cells(row: str) -> list[str]:
    """A table row's cells. A | inside a cell, in code too, is written \\|, as
    GitHub asks."""
    row = row.strip()
    row = row[1:] if row.startswith("|") else row
    row = row[:-1] if row.endswith("|") and not row.endswith("\\|") else row
    return [cell.strip().replace("\\|", "|") for cell in re.split(r"(?<!\\)\|", row)]


def plain(markdown_inline: str) -> str:
    """The text of an inline span, without its markup: what a heading's anchor,
    the page contents and the search index are made of."""
    keep: list[str] = []

    def stash(fragment: str) -> str:
        keep.append(fragment)
        return f"\x00{len(keep) - 1}\x00"

    text = code_and_escapes(markdown_inline, stash, stash)
    text = re.sub(r"<(https?://[^>\s]+)>", r"\1", text)
    text = IMAGE.sub(lambda m: m.group(1), text)
    text = LINK.sub(lambda m: m.group(1), text)
    text = INLINE_TAGS.sub("", text)
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"(?<![\w_])__(?!\s)(.+?)(?<!\s)__(?![\w_])", r"\1", text)
    text = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\1", text)
    text = re.sub(r"(?<!\w)_(?!\s)(.+?)(?<!\s)_(?!\w)", r"\1", text)
    text = html.unescape(text)
    return re.sub(r"\x00(\d+)\x00", lambda m: keep[int(m.group(1))], text)


def github_slug(text: str) -> str:
    """GitHub's anchor for a heading, so links written for GitHub land here too."""
    return re.sub(r"[^\w\- ]", "", text.strip().lower()).replace(" ", "-")


def unique_slug(slug: str, seen: dict[str, int]) -> str:
    """`slug`, numbered the way GitHub numbers a repeated heading. A numbered
    anchor is taken too, so "Foo", "Foo" and "Foo 1" get foo, foo-1 and foo-1-1
    rather than two foo-1."""
    result = slug
    while result in seen:
        seen[slug] += 1
        result = f"{slug}-{seen[slug]}"
    seen[result] = 0
    return result


def mkdocs_slug(text: str, taken: set[str]) -> str:
    """The anchor the MkDocs site gave the same heading: accents folded to
    ASCII, runs of spaces and hyphens joined, and a repeat numbered _1, _2.
    Links to the old site's sections are out in the world, and each page keeps
    these as well, so they still land on their section."""
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    slug = re.sub(r"[-\s]+", "-", re.sub(r"[^\w\s-]", "", text).strip().lower())
    while slug in taken or not slug:
        numbered = re.match(r"^(.*)_([0-9]+)$", slug)
        slug = f"{numbered.group(1)}_{int(numbered.group(2)) + 1}" if numbered else f"{slug}_1"
    taken.add(slug)
    return slug


class Renderer:
    def __init__(self, source: str, lang: str):
        self.source = source  # the English source, e.g. linkage.md
        self.lang = lang
        self.page = output_of(source, lang)
        self.seen: dict[str, int] = {}
        self.old_anchors: set[str] = set()
        self.headings: list[tuple[int, str, str]] = []
        self.sections: list[dict] = []

    @property
    def markdown(self) -> str:
        return f"docs/{localized(self.source, self.lang)}"

    # Links

    def link(self, target: str) -> str:
        if re.match(r"^(https?:|mailto:|#)", target):
            return target
        path, _, anchor = target.partition("#")
        suffix = f"#{anchor}" if anchor else ""
        if anchor in ("gh-light-mode-only", "gh-dark-mode-only"):
            suffix = ""
        # A path from / starts at the repository's root, as GitHub reads it.
        resolved = (ROOT / path.lstrip("/") if path.startswith("/") else DOCS / path).resolve()
        if resolved.is_relative_to(DOCS):
            inside = resolved.relative_to(DOCS).as_posix()
            if inside.endswith(".md"):
                spanish = inside.endswith(".es.md")
                english = inside[: -len(".es.md")] + ".md" if spanish else inside
                if english in {p[1] for p in PAGES}:
                    return relative(self.page, output_of(english, "es" if spanish else "en")) + suffix
            return relative(self.page, inside) + suffix
        if resolved.is_relative_to(ROOT):
            return f"{REPO}/blob/main/{resolved.relative_to(ROOT).as_posix()}" + suffix
        raise DocsError(f"{self.markdown} links to {target}, which is outside the repository")

    # Inline

    def inline(self, text: str) -> str:
        keep: list[str] = []

        def stash(fragment: str) -> str:
            keep.append(fragment)
            return f"\x00{len(keep) - 1}\x00"

        text = code_and_escapes(
            text, lambda code: stash(f"<code>{html.escape(code)}</code>"), lambda char: stash(html.escape(char))
        )
        text = re.sub(r"<(https?://[^>\s]+)>", lambda m: stash(f'<a href="{html.escape(m.group(1))}">{html.escape(m.group(1))}</a>'), text)
        text = INLINE_TAGS.sub(lambda m: stash(m.group(0)), text)

        def image(m: re.Match) -> str:
            alt, target, title = m.groups()
            theme = ""
            if target.endswith("#gh-light-mode-only"):
                theme = ' class="only-light"'
            elif target.endswith("#gh-dark-mode-only"):
                theme = ' class="only-dark"'
            titled = f' title="{html.escape(title)}"' if title else ""
            src = self.link(target)
            return stash(f'<img{theme} src="{html.escape(src)}" alt="{html.escape(plain(alt))}"{titled} loading="lazy">')

        text = IMAGE.sub(image, text)

        def link(m: re.Match) -> str:
            label, target, title = m.groups()
            titled = f' title="{html.escape(title)}"' if title else ""
            # The label is fenced off, so a URL written in it is not linked twice.
            return stash(f'<a href="{html.escape(self.link(target))}"{titled}>') + f"\x02{label}\x03" + stash("</a>")

        def bare(m: re.Match) -> str:
            # As on GitHub, trailing punctuation and an unmatched ) are the
            # sentence's, not the URL's.
            url, trail = m.group(0), ""
            while url and (url[-1] in "?!.,:*_~" or (url[-1] == ")" and url.count(")") > url.count("("))):
                url, trail = url[:-1], url[-1] + trail
            return stash(f'<a href="{html.escape(url)}">{html.escape(url)}</a>') + trail

        text = LINK.sub(link, text)
        text = "".join(
            part if part.startswith("\x02") else BARE_URL.sub(bare, part)
            for part in re.split(r"(\x02.*?\x03)", text, flags=re.S)
        )
        text = re.sub(r"&(?!#?\w+;)", "&amp;", text).replace("<", "&lt;").replace(">", "&gt;")
        text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
        text = re.sub(r"(?<![\w_])__(?!\s)(.+?)(?<!\s)__(?![\w_])", r"<strong>\1</strong>", text)
        text = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", text)
        text = re.sub(r"(?<![\w])_(?!\s)(.+?)(?<!\s)_(?![\w])", r"<em>\1</em>", text)
        text = text.replace("\x02", "").replace("\x03", "")
        while "\x00" in text:
            text = re.sub(r"\x00(\d+)\x00", lambda m: keep[int(m.group(1))], text)
        return text

    # Blocks

    def blocks(self, lines: list[str]) -> str:
        out: list[str] = []
        i = 0
        while i < len(lines):
            line = lines[i]
            stripped = line.strip()
            if not stripped:
                i += 1
                continue
            if stripped.startswith("<!--"):
                while i < len(lines) and "-->" not in lines[i]:
                    i += 1
                i += 1
                continue
            fence = FENCE.match(line)
            if fence:
                indent, marker, info = len(fence.group(1)), fence.group(2), fence.group(3).strip()
                body = []
                i += 1
                while i < len(lines) and not lines[i].strip().startswith(marker):
                    body.append(lines[i][indent:] if lines[i][:indent].strip() == "" else lines[i])
                    i += 1
                i += 1
                out.append(self.code(info, body))
                continue
            heading = HEADING.match(line)
            if heading:
                out.append(self.heading(len(heading.group(1)), heading.group(2)))
                i += 1
                continue
            if RULE.match(line):
                out.append("<hr>")
                i += 1
                continue
            if self.starts_table(lines, i):
                rows = [lines[i], lines[i + 1]]
                i += 2
                while i < len(lines) and lines[i].strip() and not self.starts_block(lines, i):
                    rows.append(lines[i])
                    i += 1
                out.append(self.table(rows))
                continue
            if stripped.startswith(">"):
                quoted = []
                while i < len(lines) and lines[i].strip().startswith(">"):
                    quoted.append(re.sub(r"^\s*> ?", "", lines[i]))
                    i += 1
                out.append(self.quote(quoted))
                continue
            # An HTML block runs to the next blank line, and its closing tag can
            # open one of its own, so markdown can sit between the two, as in
            # GitHub's <details> with a blank line after its <summary>.
            if HTML_BLOCK.match(stripped):
                chunk = []
                while i < len(lines) and lines[i].strip():
                    chunk.append(lines[i])
                    i += 1
                out.append(self.raw_html("\n".join(chunk)))
                continue
            if LIST_ITEM.match(line):
                i, rendered = self.list_block(lines, i)
                out.append(rendered)
                continue
            para: list[str] = []
            level = 0
            while i < len(lines) and lines[i].strip():
                if para and SETEXT.match(lines[i]):
                    level = 1 if lines[i].strip().startswith("=") else 2
                    i += 1
                    break
                if para and self.starts_block(lines, i):
                    break
                para.append(lines[i].strip())
                i += 1
            text = " ".join(para)
            if level:
                out.append(self.heading(level, text))
                continue
            self.note(text)
            rendered = self.inline(text)
            # A paragraph of nothing but images is a figure.
            if re.fullmatch(r"(\s*<img [^>]+>\s*)+", rendered):
                out.append(f'<p class="figure">{rendered}</p>')
            else:
                out.append(f"<p>{rendered}</p>")
        return "\n".join(out)

    @staticmethod
    def starts_table(lines: list[str], i: int) -> bool:
        if i + 1 >= len(lines) or "|" not in lines[i] or "|" not in lines[i + 1]:
            return False
        delimiter = lines[i + 1].strip()
        return bool(TABLE_DELIMITER.match(delimiter)) and len(table_cells(lines[i])) == len(table_cells(delimiter))

    def starts_block(self, lines: list[str], i: int) -> bool:
        """Whether line `i` ends the paragraph or list item above it. GitHub
        starts a new block at a heading, rule, fence, quote, comment, HTML
        block or table, and at a list item with text, an ordered one only when
        it counts from 1, so a wrapped line that begins `2024.` stays text."""
        line = lines[i]
        stripped = line.strip()
        if (
            FENCE.match(line)
            or HEADING.match(line)
            or RULE.match(line)
            or stripped.startswith((">", "<!--"))
            or HTML_BLOCK.match(stripped)
            or self.starts_table(lines, i)
        ):
            return True
        item = LIST_ITEM.match(line)
        if not item:
            return False
        marker = item.group(2)
        return bool(line[item.end() :].strip()) and (not marker[0].isdigit() or marker == "1.")

    def heading(self, level: int, text: str) -> str:
        name = plain(text)
        slug = unique_slug(github_slug(name), self.seen)
        old = mkdocs_slug(name, self.old_anchors)
        self.headings.append((level, name, slug))
        self.sections.append({"heading": name, "slug": slug, "level": level, "text": []})
        alias = f'<span class="alias" id="{old}"></span>' if old != slug else ""
        anchor = f'<a class="anchor" href="#{slug}" aria-hidden="true" tabindex="-1">#</a>' if level > 1 else ""
        return f'<h{level} id="{slug}">{alias}{self.inline(text)}{anchor}</h{level}>'

    def note(self, text: str) -> None:
        if self.sections:
            self.sections[-1]["text"].append(plain(text))

    def quote(self, quoted: list[str]) -> str:
        """A block quote, or one of GitHub's alerts (> [!WARNING]) as a callout
        whose title is its first line when that line is bold and nothing else."""
        first = next((q.strip() for q in quoted if q.strip()), "")
        alert = re.match(r"^\[!(\w+)\]$", first)
        notes = UI[self.lang]["notes"]
        if not alert or alert.group(1).upper() not in notes:
            return f"<blockquote>{self.blocks(quoted)}</blockquote>"
        kind = alert.group(1).upper()
        body = quoted[quoted.index(next(q for q in quoted if q.strip())) + 1 :]
        while body and not body[0].strip():
            body = body[1:]
        title = notes[kind]
        bold = re.fullmatch(r"\*\*((?:(?!\*\*).)+)\*\*", body[0].strip()) if body else None
        if bold:
            title = self.inline(bold.group(1))
            self.note(bold.group(1))
            body = body[1:]
        return (
            f'<div class="callout callout-{kind.lower()}"><p class="callout-title">{title}</p>'
            f"{self.blocks(body)}</div>"
        )

    def list_block(self, lines: list[str], i: int) -> tuple[int, str]:
        first = LIST_ITEM.match(lines[i])
        base = len(first.group(1))
        ordered = first.group(2)[0].isdigit()
        start = int(first.group(2)[:-1]) if ordered else 1
        items: list[list[str]] = []
        loose = False
        while i < len(lines):
            line = lines[i]
            m = LIST_ITEM.match(line)
            if m and len(m.group(1)) == base and (m.group(2)[0].isdigit()) == ordered:
                content_indent = len(m.group(1)) + len(m.group(2)) + len(m.group(3))
                items.append([line[content_indent:]])
                i += 1
                while i < len(lines):
                    nxt = lines[i]
                    if not nxt.strip():
                        j = i + 1
                        while j < len(lines) and not lines[j].strip():
                            j += 1
                        if j < len(lines) and (len(lines[j]) - len(lines[j].lstrip())) > base:
                            items[-1].append("")
                            loose = loose or not re.match(r"^\s*([-*+]|\d+\.)\s", lines[j])
                            i += 1
                            continue
                        if j < len(lines) and re.match(rf"^\s{{{base}}}([-*+]|\d+\.)\s", lines[j]):
                            loose = True
                            i = j
                        break
                    indent = len(nxt) - len(nxt.lstrip())
                    if indent > base:
                        items[-1].append(nxt[min(indent, content_indent):])
                        i += 1
                        continue
                    # A line as far left as the marker continues the item's
                    # text, unless it starts an item or a block of its own.
                    if LIST_ITEM.match(nxt) or self.starts_block(lines, i):
                        break
                    items[-1].append(nxt.strip())
                    i += 1
                continue
            break
        tag = "ol" if ordered else "ul"
        attrs = f' start="{start}"' if ordered and start != 1 else ""
        rendered = []
        for content in items:
            body = self.blocks(content)
            if not loose and body.startswith("<p>") and body.count("<p>") == 1:
                body = body[3:].replace("</p>", "", 1)
            rendered.append(f"<li>{body}</li>")
        return i, f"<{tag}{attrs}>" + "".join(rendered) + f"</{tag}>"

    def table(self, rows: list[str]) -> str:
        head = table_cells(rows[0])
        aligns = table_cells(rows[1])
        body = []
        for row in rows[2:]:
            cells = table_cells(row)
            # GitHub drops the cells past the header's; here they would break
            # the table, so a row like that is an error to fix in the markdown.
            if len(cells) > len(head):
                raise DocsError(
                    f"{self.markdown}: a table row has {len(cells)} cells under a header of {len(head)}. "
                    f"Write a | inside a cell, in code too, as \\|:\n  {row.strip()[:120]}"
                )
            body.append(cells + [""] * (len(head) - len(cells)))
        for row in body:
            self.note(" ".join(row))

        def cls(k: int) -> str:
            return ' class="num"' if k < len(aligns) and aligns[k].endswith(":") and not aligns[k].startswith(":") else ""

        out = ['<div class="table-wrap"><table><thead><tr>']
        out += [f"<th{cls(k)}>{self.inline(c)}</th>" for k, c in enumerate(head)]
        out.append("</tr></thead><tbody>")
        for row in body:
            out.append("<tr>" + "".join(f"<td{cls(k)}>{self.inline(c)}</td>" for k, c in enumerate(row)) + "</tr>")
        out.append("</tbody></table></div>")
        return "".join(out)

    def code(self, info: str, body: list[str]) -> str:
        lang = info.split()[0] if info else ""
        data_lang = f' data-lang="{html.escape(lang)}"' if lang else ""
        text = html.escape("\n".join(body))
        self.note(" ".join(body))
        return f'<div class="code"{data_lang}><pre><code>{text}</code></pre></div>'

    def raw_html(self, chunk: str) -> str:
        """Raw HTML passes through, with its asset paths made relative to the
        page. A <picture> that swaps on the system colour scheme follows the
        site's own theme instead, which is the one the reader chose."""
        picture = re.search(
            r'<source media="\(prefers-color-scheme: dark\)" srcset="([^"]+)">\s*<img src="([^"]+)" alt="([^"]*)">', chunk
        )
        if picture:
            dark, light, alt = picture.groups()
            cls = re.search(r'^<p class="([^"]+)"', chunk.strip())
            extra = f" {cls.group(1)}" if cls else ""
            return (
                f'<p class="figure{extra}"><img class="only-dark" src="{self.link(dark)}" alt="{alt}">'
                f'<img class="only-light" src="{self.link(light)}" alt="{alt}"></p>'
            )
        return re.sub(r'(src|href)="(?!https?:|#|mailto:)([^"]+)"', lambda m: f'{m.group(1)}="{self.link(m.group(2))}"', chunk)


# ── Page ──


def header(page: str, lang: str, twin: str, current: str = "") -> str:
    ui = UI[lang]
    home = relative(page, output_of("index.html", lang))
    emblem = relative(page, "assets/emblem.svg")
    other = ui["other"]
    links = [
        ("getting-started.md", relative(page, output_of("getting-started.md", lang)), ui["start"]),
        ("#try", f"{home}#try", ui["try"]),
        ("cli-reference.md", relative(page, output_of("cli-reference.md", lang)), ui["reference"]),
        ("example-report", relative(page, "assets/example-report.html"), ui["report"]),
    ]
    here = ' aria-current="page"'
    nav = "".join(
        f'      <a href="{href}"{here if key == current else ""}>{label}</a>\n' for key, href, label in links
    )
    return (
        '<header class="top">\n  <div class="top-in">\n'
        f'    <a class="brand" href="{home}"><img src="{emblem}" alt="" width="28" height="28"><span>get_MNV</span></a>\n'
        f'    <nav class="top-links" aria-label="{ui["site_nav"]}">\n{nav}    </nav>\n'
        '    <div class="top-tools">\n'
        '      <div class="search" role="search">\n'
        f'        {SEARCH_SVG}<input type="search" placeholder="{ui["search"]}" aria-label="{ui["search_label"]}" data-search>'
        '<kbd aria-hidden="true">/</kbd>\n'
        f'        <button class="icon-btn search-open" type="button" aria-label="{ui["search_label"]}" data-search-open>{SEARCH_SVG}</button>\n'
        '        <div class="results" data-results hidden></div>\n      </div>\n'
        f'      <a class="lang" href="{twin}" hreflang="{other}" lang="{other}">{ui["other_label"]}</a>\n'
        f'      <a class="gh" href="{REPO}" aria-label="get_MNV on GitHub">{GH_SVG}<span>GitHub</span></a>\n'
        f'      <button class="icon-btn" type="button" data-theme-toggle aria-label="{ui["theme"]}">{THEME_SVG}</button>\n'
        "    </div>\n  </div>\n</header>\n"
    )


def sidebar(page: str, lang: str, current: str) -> str:
    ui = UI[lang]
    out = [
        f'<aside class="side" aria-label="{ui["docs_label"]}">\n',
        f'  <details class="menu" open>\n    <summary>{ui["menu"]}</summary>\n',
    ]
    section = None
    for sec, source, title in PAGES:
        if sec != section:
            if section is not None:
                out.append("    </ul>\n")
            out.append(f'    <p class="menu-section">{title_in(sec, lang)}</p>\n    <ul>\n')
            section = sec
        if source.endswith(".md") or source == "index.html":
            target = relative(page, output_of(source, lang))
        elif source.startswith("http"):
            target = source
        else:
            target = relative(page, source)
        here = ' aria-current="page"' if source == current else ""
        out.append(f'      <li><a href="{target}"{here}>{title_in(title, lang)}</a></li>\n')
    out.append("    </ul>\n  </details>\n</aside>\n")
    return "".join(out)


def footer(page: str, lang: str) -> str:
    ui = UI[lang]
    return (
        '<footer class="foot">\n  <div class="foot-in">\n'
        f'    <div><p class="foot-h">{ui["cite"]}</p>\n'
        "      <p>Ruiz-Rodriguez P, Coscolla M. get_MNV. "
        '<a href="https://doi.org/10.5281/zenodo.13907422">doi:10.5281/zenodo.13907422</a></p></div>\n'
        f'    <div><p class="foot-h">{ui["project"]}</p>\n      <ul>\n'
        f'        <li><a href="{REPO}">GitHub</a></li>\n'
        f'        <li><a href="{REPO}/releases">{ui["releases"]}</a></li>\n'
        f'        <li><a href="{ASK}">{ui["ask"]}</a></li>\n'
        f'        <li><a href="{REPO}/issues">{ui["issues"]}</a></li>\n'
        "      </ul></div>\n"
        f'    <p class="foot-note">Copyright &copy; Paula Ruiz-Rodriguez &amp; Mireia Coscolla, I&sup2;SysBio (University of Valencia-CSIC). {ui["license"]}</p>\n'
        "  </div>\n</footer>\n"
    )


def page_html(
    page: str, lang: str, title: str, description: str, body: str, twin: str, scripts: tuple[str, ...] = ()
) -> str:
    css = relative(page, "assets/site.css")
    js = relative(page, "assets/site.js")
    index = relative(page, "assets/search-index.js")
    icon = relative(page, "assets/emblem.svg")
    root = relative(page, ".") + "/"
    canonical = SITE_URL + (page[: -len("index.html")] if page.endswith("index.html") else page)
    other_lang = UI[lang]["other"]
    return (
        "<!doctype html>\n"
        f'<html lang="{lang}" data-root="{root}">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{html.escape(title)}</title>\n"
        f'<meta name="description" content="{html.escape(description)}">\n'
        f'<link rel="canonical" href="{canonical}">\n'
        f'<link rel="alternate" hreflang="{other_lang}" href="{twin}">\n'
        f'<link rel="icon" href="{icon}" type="image/svg+xml">\n'
        '<meta property="og:type" content="website">\n<meta property="og:site_name" content="get_MNV">\n'
        f'<meta property="og:title" content="{html.escape(title)}">\n'
        f'<meta property="og:description" content="{html.escape(description)}">\n'
        f'<meta property="og:url" content="{canonical}">\n'
        f"{THEME_BOOT}\n{FONTS}"
        f'<link rel="stylesheet" href="{css}">\n</head>\n<body>\n'
        f'<a class="skip" href="#content">{UI[lang]["skip"]}</a>\n'
        "<!-- Generated by tools/build_docs.py from the markdown in docs/. Edit that, then run it. -->\n"
        + body
        + f'<script src="{index}"></script>\n<script src="{js}"></script>\n'
        + "".join(f'<script src="{relative(page, s)}"></script>\n' for s in scripts)
        + "</body>\n</html>\n"
    )


def built_pages() -> list[tuple[int, str]]:
    return [(n, p[1]) for n, p in enumerate(PAGES) if p[1].endswith(".md")]


def doc_page(position: int, lang: str) -> tuple[Path, str, list[dict]]:
    section, source, title = PAGES[position]
    page = output_of(source, lang)
    renderer = Renderer(source, lang)
    lines = (DOCS / localized(source, lang)).read_text(encoding="utf-8").splitlines()
    content = renderer.blocks(lines)
    # An old anchor that is now a heading's own stays that heading's.
    taken = {slug for _, _, slug in renderer.headings}
    content = re.sub(
        r'<span class="alias" id="([^"]+)"></span>', lambda m: "" if m.group(1) in taken else m.group(0), content
    )
    toc = "".join(
        f'<li class="toc-{level}"><a href="#{slug}">{html.escape(text)}</a></li>'
        for level, text, slug in renderer.headings
        if level in (2, 3)
    )
    ui = UI[lang]
    order = built_pages()
    at = [n for n, _ in order].index(position)
    pager = '<nav class="pager" aria-label="Previous and next page">'
    if at > 0:
        prev_n, prev_src = order[at - 1]
        pager += f'<a class="prev" href="{relative(page, output_of(prev_src, lang))}"><span>{ui["prev"]}</span>{title_in(PAGES[prev_n][2], lang)}</a>'
    if at + 1 < len(order):
        next_n, next_src = order[at + 1]
        pager += f'<a class="next" href="{relative(page, output_of(next_src, lang))}"><span>{ui["next"]}</span>{title_in(PAGES[next_n][2], lang)}</a>'
    pager += "</nav>"
    twin = relative(page, output_of(source, ui["other"]))
    body = (
        header(page, lang, twin, source)
        + '<div class="layout">\n'
        + sidebar(page, lang, source)
        + f'<main class="prose" id="content">\n'
        + f'<p class="crumb">{title_in(section, lang)}</p>\n'
        + content
        + f'\n<p class="edit"><a href="{REPO}/blob/main/docs/{localized(source, lang)}">{GH_SVG}{ui["edit"]}</a></p>\n'
        + pager
        + "\n</main>\n"
        + (f'<aside class="toc" aria-label="{ui["on_page"]}"><p>{ui["on_page"]}</p><ul>{toc}</ul></aside>\n' if toc else "")
        + "</div>\n"
        + footer(page, lang)
    )
    shown = title_in(title, lang)
    full_title = f"{shown} · get_MNV"
    entries = [
        {
            "lang": lang,
            "page": shown,
            "section": s["heading"] if s["level"] > 1 else "",
            "url": page + (f"#{s['slug']}" if s["level"] > 1 else ""),
            "text": " ".join(s["text"]),
        }
        for s in renderer.sections
    ]
    return DOCS / page, page_html(page, lang, full_title, ui["description"], body, twin), entries


HOME_SOURCES = ROOT / "tools" / "site"
HOME_LINK = re.compile(r"\{\{(page|asset|repo|home)(?::([^}]+))?\}\}")


def home_page(lang: str) -> tuple[Path, str, list[dict]]:
    """The front page, from tools/site/home.<lang>.html: its own layout, the
    site's header and footer. Links in the template are written as
    {{page:getting-started.md}}, {{asset:example-report.html}}, {{home}} or
    {{repo}}, so each language's page gets paths that resolve from where it is
    written, and tools/check_site.py follows them like any other."""
    ui = UI[lang]
    page = output_of("index.html", lang)
    twin = relative(page, output_of("index.html", ui["other"]))

    def link(match: re.Match) -> str:
        kind, target = match.groups()
        if kind == "page":
            return relative(page, output_of(target, lang))
        if kind == "asset":
            return relative(page, f"assets/{target}")
        if kind == "home":
            return relative(page, page)
        return REPO

    template = (HOME_SOURCES / f"home.{lang}.html").read_text(encoding="utf-8")
    body = header(page, lang, twin, "home") + HOME_LINK.sub(link, template) + footer(page, lang)

    # Its sections go into the search index like any page's.
    def text_of(chunk: str) -> str:
        return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", HOME_LINK.sub("", chunk))).split())

    shown = title_in("Home", lang)
    hero = re.search(r'<p class="lede">(.*?)</p>', template, re.S)
    entries = [{"lang": lang, "page": shown, "section": "", "url": page, "text": text_of(hero.group(1))}]
    for slug, heading, rest in re.findall(
        r'<section class="sec" id="([\w-]+)".*?<h2[^>]*>(.*?)</h2>(.*?)</section>', template, re.S
    ):
        entries.append(
            {"lang": lang, "page": shown, "section": text_of(heading), "url": f"{page}#{slug}", "text": text_of(rest)}
        )
    html_text = page_html(page, lang, ui["home_title"], ui["description"], body, twin, ("assets/home.js",))
    html_text = html_text.replace(
        "<!-- Generated by tools/build_docs.py from the markdown in docs/. Edit that, then run it. -->",
        f"<!-- Generated by tools/build_docs.py from tools/site/home.{lang}.html. Edit that, then run it. -->",
    )
    return DOCS / page, html_text, entries


def not_found_page() -> str:
    """The page Pages serves for an address it does not have. It can be served
    at any depth, so every link on it is absolute."""
    body = (
        '<header class="top">\n  <div class="top-in">\n'
        f'    <a class="brand" href="{SITE_URL}"><img src="{SITE_URL}assets/emblem.svg" alt="" width="28" height="28">'
        "<span>get_MNV</span></a>\n  </div>\n</header>\n"
        '<main class="prose notfound" id="content">\n<h1>Page not found</h1>\n'
        "<p>Nothing lives at this address. The documentation starts at the "
        f'<a href="{SITE_URL}">home page</a>; '
        f'<a href="{SITE_URL}es/">la documentación en español</a> también.</p>\n'
        f'<p>If a link brought you here, <a href="{REPO}/issues">say where it was</a> and it will be fixed.</p>\n'
        "</main>\n"
    )
    return (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        "<title>Page not found · get_MNV</title>\n"
        f'<link rel="icon" href="{SITE_URL}assets/emblem.svg" type="image/svg+xml">\n'
        f"{THEME_BOOT}\n{FONTS}"
        f'<link rel="stylesheet" href="{SITE_URL}assets/site.css">\n</head>\n<body>\n'
        "<!-- Generated by tools/build_docs.py. -->\n" + body + "</body>\n</html>\n"
    )


def build() -> dict[Path, str]:
    files: dict[Path, str] = {DOCS / "404.html": not_found_page()}
    entries: list[dict] = []
    for lang in ("en", "es"):
        path, text, home_entries = home_page(lang)
        files[path] = text
        entries += home_entries
        for position, _ in built_pages():
            path, text, page_entries = doc_page(position, lang)
            files[path] = text
            entries += page_entries
    files[SEARCH_INDEX] = (
        "/* Generated by tools/build_docs.py: the documentation, section by section, for the search box. */\n"
        "window.GET_MNV_DOCS_INDEX = " + json.dumps(entries, ensure_ascii=False, indent=0) + ";\n"
    )
    return files


GENERATED = "Generated by tools/build_docs.py"


def leftovers(files: dict[Path, str]) -> list[Path]:
    """Pages this script wrote that nothing writes now, such as the old page
    of a markdown file that was renamed. Committed, they would still be served
    and found, with nothing left to update them."""
    return sorted(
        p
        for p in DOCS.rglob("*.html")
        if p not in files and GENERATED in p.read_text(encoding="utf-8", errors="replace")
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="exit 1 if a generated page differs from what is committed")
    args = parser.parse_args()
    try:
        files = build()
    except DocsError as error:
        print(error, file=sys.stderr)
        return 1
    stale = [p for p, text in files.items() if not p.exists() or p.read_text(encoding="utf-8") != text]
    orphans = leftovers(files)
    if args.check:
        if stale:
            print("These pages no longer match their markdown; run python3 tools/build_docs.py:", file=sys.stderr)
            for p in stale:
                print(f"  {p.relative_to(ROOT)}", file=sys.stderr)
        if orphans:
            print("These pages have no markdown any more; run python3 tools/build_docs.py to remove them:", file=sys.stderr)
            for p in orphans:
                print(f"  {p.relative_to(ROOT)}", file=sys.stderr)
        if stale or orphans:
            return 1
        print(f"docs/: {len(files)} generated files match their markdown")
        return 0
    for p in stale:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(files[p], encoding="utf-8")
    for p in orphans:
        p.unlink()
        folder = p.parent
        while folder != DOCS and not any(folder.iterdir()):
            folder.rmdir()
            folder = folder.parent
    removed = f", removed {len(orphans)} with no markdown" if orphans else ""
    print(f"docs/: wrote {len(stale)} of {len(files)} generated files{removed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
