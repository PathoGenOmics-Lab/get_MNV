#!/usr/bin/env python3
"""Generate an illustrative PDF with, for every scenario:

  - The mapped reads (a small alignment pileup) against the reference,
    marking SNVs / insertions / deletions / introns (CIGAR N), the genes
    that overlap, and the input variants (VCF or iVar).
  - Two output tables:
      * "Without get_mnv": each input variant annotated INDEPENDENTLY
        (no MNV grouping) — computed by running get_mnv on one variant at a
        time, which is what a plain per-call annotation would yield.
      * "get_mnv output": the real, codon-aware grouped annotation.

Designed to avoid label collisions:
  - Variant / insertion labels are spread across "levels" (de-collision).
  - Per-read-group labels live in the left margin, outside the data area.
  - Output is shown in grid tables.

Usage:
  SAMTOOLS=.../samtools BGZIP=.../bgzip TABIX=.../tabix \\
      python3 tests/scenarios/plot_scenarios.py [--png NAME ...]
"""

from __future__ import annotations

import os
import shutil
import sys
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# --- Make samtools/bgzip/tabix locatable BEFORE importing framework ---
_BIOCONDA = Path.home() / "miniconda3" / "envs" / "bioconda" / "bin"
for _tool, _env in (("samtools", "SAMTOOLS"), ("bgzip", "BGZIP"), ("tabix", "TABIX")):
    if not os.environ.get(_env) and shutil.which(_tool) is None:
        cand = _BIOCONDA / _tool
        if cand.exists():
            os.environ[_env] = str(cand)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Polygon
from matplotlib.transforms import blended_transform_factory

import framework as F
from framework import run_scenario, run_get_mnv, parse_tsv, contig_len_for
from scenarios import ALL_SCENARIOS

# ----------------------------------------------------------------------------
# Style
# ----------------------------------------------------------------------------
BASE_COLOR = {"A": "#3FA45B", "C": "#3C6FE0", "G": "#E8A33D", "T": "#D6453D", "N": "#999999"}
GENE_COLORS = ["#6F398D", "#149389", "#C2185B", "#E8843C", "#2D7DD2"]
ACCENT = "#2D7DD2"
NAIVE_HDR = "#8A8F98"
DEL_COLOR = "#D6453D"
INS_COLOR = "#7B1FA2"
READ_FILL = "#C9D6E8"
READ_EDGE = "#4A6079"

READ_PITCH = 0.78           # vertical spacing between stacked reads (lanes)
READ_H = 0.52               # read bar height
READ_BUDGET = 12            # max reads drawn per view (shared across groups)

GETMNV_COLS = ["Positions", "Gene", "Variant Type", "Change Type",
               "AA Changes", "Event Class", "Event Components"]
GETMNV_W = [0.08, 0.085, 0.09, 0.135, 0.15, 0.10, 0.36]

RAW_COLS = ["Position", "Allele", "Class (caller)", "Frame (by length)"]
RAW_W = [0.12, 0.26, 0.30, 0.24]

# English one-line descriptions (the PDF is fully English; scenarios.py keeps
# its Spanish docstrings for the test harness).
DESC = {
    "01_snp_simple": "Simple non-synonymous SNP: pos 28 G>A in codon 10 (GCT→ACT, Ala10Thr); 20/20 alt reads.",
    "02a_snp_mnv_full_phasing": "SNP/MNV: pos 28 G>T + pos 30 T>A together on all reads (full MNV haplotype).",
    "02b_vcf_mnp_decomposed": "VCF MNP REF=GC ALT=TA at pos 28-29 decomposed into a SNP/MNV (Ala10Tyr); event_class=mnv.",
    "03_snp_mnv_mixed": "Mixed SNP/MNV: 10 reads carry pos28+pos30, 10 only pos28, 10 only pos30 (30 total).",
    "04_ins_inframe_cds": "In-frame insertion of GCT (1 Ala) between pos 30 and 31; 20 reads.",
    "05_del_frameshift_cds": "Frameshift 1 bp deletion at pos 31; 20 reads (Ala11Leufs).",
    "06_indel_plus_snv_haplotype": "Combined haplotype: in-frame insertion at pos 30 + SNV pos 36 G>A (Ala12Thr) in cis.",
    "07_fs_del_plus_downstream_snv": "Frameshift deletion at pos 31 + downstream SNV pos 39 (>3 bp away): the SNV gets the (fs) marker.",
    "08_inframe_ins_inside_codon_with_mnv": "In-frame insertion INSIDE codon 10 + MNV pos28+pos30, all in cis: get_mnv emits 5 rows.",
    "09_fs_del_with_snv_overlap": "Frameshift 2 bp deletion over codon 10 + SNV pos28 in cis: SNP-overlap, lone INDEL (no reads), complex_indel.",
    "10_minus_snp_simple": "Non-synonymous SNP on geneB (− strand): pos 673 C>T → mRNA codon 10 GCT→ACT (Ala10Thr). Codons shown genomic-forward.",
    "11_minus_mnv_same_codon": "MNV on geneB codon 10 (− strand): pos 671 A>T + pos 673 C>A in cis (mRNA GCT→TCA, Ala10Ser); genomic-forward codons.",
    "12_minus_fs_del": "Frameshift 1 bp deletion at pos 671 (geneB, − strand): codon 10 frameshifted.",
    "13_multiexon_snp_exon2": "Non-synonymous SNP in exon 2 of geneC: pos 1048 G>A in codon 50 (GCT→ACT, Ala50Thr).",
    "14_multiexon_junction_snp": "SNP at the geneC junction codon: pos 900 G>A (base 1 of codon 34 GCT→ACT, Ala34Thr); spliced N-CIGAR reads.",
    "15_multiexon_junction_mnv": "MNV crossing the junction: pos 900 G>A + pos 1002 T>C in cis (mRNA codon 34 GCT→ACC, Ala34Thr).",
    "16_no_bam_coverage": "VCF SNV at pos 28 (geneA) with 0 reads covering it in the BAM.",
    "17_min_snp_frequency_filter": "SNV at 10% frequency (2/20 reads) filtered out by --min-snp-frequency 0.5 → 0 rows.",
    "18_stop_gained_via_mnv": "Stop-gained via a 3-SNV MNV: codon 50 GCT→TAA (Ala50Ter).",
    "19_start_codon_altered": "SNV in the start codon ATG: pos 2 T>C → ACG (Met1Thr). NOTE: reported as Non-synonymous, not 'Start lost'.",
    "20_stop_lost": "Stop-lost: codon 100 TAA→CAA (pos 298 T>C, Ter100Gln).",
    "21_intron_variant": "SNV in the geneC intron (pos 950 T>A): reported as intergenic under --gff-features CDS.",
    "22_multiallelic_split": "Multiallelic SNV pos 28 G>A,T with --split-multiallelic: 2 independent SNV rows (G>A Ala10Thr + G>T Ala10Ser).",
    "23_delins_subst_plus_del": "Delins REF=GCT ALT=GA at pos 28 (substitution pos 29 C>A + deletion pos 30): frameshift.",
    "24_large_inframe_insertion": "In-frame 9 bp insertion (GCTGCTGCT = 3 Ala) between pos 30 and 31.",
    "25_large_inframe_deletion": "In-frame 6 bp deletion (2 consecutive Ala) at pos 31-36.",
    "26_multicontig": "Multi-contig: SNV on chr_test2 geneD codon 10 (GCT→ACT, Ala10Thr) + SNV on chr_test geneA.",
    "27_ivar_tsv_snv": "iVar TSV input: SNV pos 28 G>A (Ala10Thr).",
    "28_ivar_tsv_insertion": "iVar TSV: in-frame insertion +GCT at pos 30 (+SEQ notation).",
    "29_ivar_tsv_deletion": "iVar TSV: frameshift deletion -G at pos 31 (-SEQ notation).",
    "30_stop_gained_inframe_ins": "In-frame insertion of TAA after pos 30 creates a premature stop → Change Type 'Stop gained' (not 'In-frame Indel').",
    "31_stop_lost_inframe_del": "In-frame 3 bp deletion removing the stop TAA (298-300) → Change Type 'Stop lost'.",
    "32_fs_gate_default_propagates": "Upstream frameshift deletion AF=0.20 + downstream SNV: by default the frameshift propagates → SNV marked (fs).",
    "33_fs_gate_high_freq_suppressed": "Same inputs as scenario 32 but with --frameshift-min-freq 0.5: the AF=0.20 deletion fails the gate → downstream SNV WITHOUT (fs).",
}

# Thematic sections (the PDF is grouped by case, not by raw number).
# Each entry: (full title, short tag, blurb, [scenario codes]).
CAT_COLORS = ["#2D7DD2", "#149389", "#E8843C", "#C2185B", "#6F398D", "#3FA45B", "#8A6D3B"]
CATEGORIES = [
    ("SNVs and MNVs (single-codon substitutions)", "SNV / MNV",
     "Single-nucleotide and multi-nucleotide variants within one codon. get_mnv's headline feature: "
     "grouping co-occurring SNVs into the correct combined amino-acid change, plus stop/start effects.",
     ["01", "02a", "02b", "03", "18", "19", "20"]),
    ("Indels — classification and amino-acid effect", "Indels: classification",
     "Single insertions/deletions: in-frame vs frameshift, large in-frame indels, delins, and in-frame "
     "indels that create or remove a stop codon.",
     ["04", "05", "23", "24", "25", "30", "31"]),
    ("Indels combined with SNVs (complex events, frameshift propagation)", "Indels: complex + frameshift",
     "Where get_mnv's indel handling is unique: joining an indel with in-cis SNVs into one complex_indel, "
     "and propagating a frameshift to downstream variants (with an optional frequency gate).",
     ["06", "07", "08", "09", "32", "33"]),
    ("Negative-strand gene", "Negative strand",
     "Genes on the minus strand (reverse-complement codon maths). Reported codons are shown in "
     "genomic-forward orientation — a documented behaviour.",
     ["10", "11", "12"]),
    ("Eukaryotic multi-exon transcripts (phase and splicing)", "Eukaryotic / multi-exon",
     "Spliced CDS models (gene → mRNA → CDS with Parent), non-zero phase, intron skips, and codons / MNVs "
     "that span an exon–exon junction.",
     ["13", "14", "15", "21"]),
    ("Input representation (multiallelic, multi-contig, iVar TSV)", "Input formats",
     "How get_mnv ingests different inputs: multiallelic VCF records, multiple contigs, and iVar TSV "
     "+/-SEQ notation.",
     ["22", "26", "27", "28", "29"]),
    ("Operational edge cases (coverage and filtering)", "Edge cases",
     "BAM coverage gaps and frequency-based filtering.",
     ["16", "17"]),
]


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def parse_gff(content: str) -> list[dict]:
    feats = []
    for line in content.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        feats.append({
            "contig": f[0], "type": f[2],
            "start": int(f[3]), "end": int(f[4]), "strand": f[6],
            "phase": int(f[7]) if f[7] in ("0", "1", "2") else 0,
            "name": attrs.get("Name", attrs.get("ID", "")),
            "parent": attrs.get("Parent", ""),
        })
    return feats


# Standard genetic code (table 1; identical to table 11 for the internal codons
# used by the synthetic genes). Used only to label the reference reading frame.
_GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_AA3 = {
    "F": "Phe", "L": "Leu", "I": "Ile", "M": "Met", "V": "Val", "S": "Ser",
    "P": "Pro", "T": "Thr", "A": "Ala", "Y": "Tyr", "*": "Ter", "H": "His",
    "Q": "Gln", "N": "Asn", "K": "Lys", "D": "Asp", "E": "Glu", "C": "Cys",
    "W": "Trp", "R": "Arg", "G": "Gly", "?": "?",
}
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _revcomp(s):
    return "".join(_COMP.get(b, "N") for b in reversed(s.upper()))


def _aa3(codon):
    return _AA3.get(_GENETIC_CODE.get(codon.upper(), "?"), "?")


def codon_track(view, features, snv_map):
    """Full codons overlapping the window, per CDS/gene feature (strand + phase
    aware). Returns dicts {lo, hi, ref_aa, alt_aa, hit} where ref_aa/alt_aa are
    3-letter amino acids and `hit` is True if an input SNV falls in the codon.
    Codons split across an exon junction are skipped (no clean triplet)."""
    x0, x1 = view["x0"], view["x1"]
    contig = view["contig"]
    ref = F.reference_for(contig)
    out = []
    for f in features:
        if f["contig"] != contig or f["end"] < x0 or f["start"] > x1:
            continue
        p = f.get("phase", 0) or 0
        if f["strand"] == "+":
            s = f["start"] + (p % 3)
            while s + 2 <= f["end"]:
                if s + 2 >= x0 and s <= x1:
                    out.append(_codon_entry(ref, s, s + 2, [s, s + 1, s + 2], snv_map, contig, "+"))
                s += 3
        else:
            h = f["end"] - (p % 3)
            while h - 2 >= f["start"]:
                if h >= x0 and h - 2 <= x1:
                    out.append(_codon_entry(ref, h - 2, h, [h - 2, h - 1, h], snv_map, contig, "-"))
                h -= 3
    return out


def _codon_entry(ref, lo, hi, positions, snv_map, contig, strand):
    tri = [ref[g - 1] for g in positions]
    alt = list(tri)
    hit = False
    for i, g in enumerate(positions):
        b = snv_map.get((contig, g))
        if b:
            alt[i] = b
            hit = True
    ref_g, alt_g = "".join(tri), "".join(alt)
    if strand == "-":
        ref_g, alt_g = _revcomp(ref_g), _revcomp(alt_g)
    return {"lo": lo, "hi": hi, "ref_aa": _aa3(ref_g), "alt_aa": _aa3(alt_g), "hit": hit}


def assign_levels(items: list[tuple[float, float]]) -> list[int]:
    """items = [(left, right)] in data coords; return a level per item so that
    items on the same level never overlap horizontally."""
    levels: list[float] = []
    out = []
    for left, right in items:
        placed = None
        for lvl, occupied in enumerate(levels):
            if left > occupied:
                placed = lvl
                levels[lvl] = right
                break
        if placed is None:
            placed = len(levels)
            levels.append(right)
        out.append(placed)
    return out


def read_op_positions(rg) -> list[int]:
    pos = []
    for op in rg.ops:
        if op.kind in ("snv", "ins"):
            pos.append(op.pos)
        elif op.kind in ("del", "skip"):
            pos.append(op.pos)
            pos.append(op.pos + max(1, op.length) - 1)
    return pos


def scenario_inputs(scn):
    out = []
    if scn.ivar_records:
        for r in scn.ivar_records:
            out.append({"pos": r.pos, "chrom": r.chrom, "label": r.alt})
    else:
        for v in scn.variants:
            out.append({"pos": v.pos, "chrom": v.chrom, "label": f"{v.ref}→{v.alt}"})
    return out


def build_views(scn):
    inputs = scenario_inputs(scn)
    contigs = []
    for it in inputs:
        if it["chrom"] not in contigs:
            contigs.append(it["chrom"])
    if not contigs:
        for rg in scn.reads:
            c = getattr(rg, "chrom", F.CONTIG)
            if c not in contigs:
                contigs.append(c)

    views = []
    for c in contigs:
        reads = [rg for rg in scn.reads if getattr(rg, "chrom", F.CONTIG) == c]
        vins = [it for it in inputs if it["chrom"] == c]
        poi = [it["pos"] for it in vins]
        for rg in reads:
            poi += read_op_positions(rg)
        if not poi:
            poi = [rg.start for rg in reads] + [rg.start + rg.length - 1 for rg in reads]
        lo, hi = min(poi), max(poi)
        pad = max(7, int((hi - lo) * 0.25) + 3)
        views.append({
            "contig": c,
            "x0": max(1, lo - pad),
            "x1": min(contig_len_for(c), hi + pad),
            "show_bases": (min(contig_len_for(c), hi + pad) - max(1, lo - pad)) <= 58,
            "reads": reads, "inputs": vins,
        })
    return views


def _frame_by_length(delta: int) -> str:
    return "frameshift" if delta % 3 else "in-frame"


def classify_vcf(pos, ref, alt):
    """Naive, get_mnv-free classification of a single VCF allele, using only
    REF/ALT lengths (what a plain caller record tells you)."""
    dl = len(alt) - len(ref)
    if len(ref) == 1 and len(alt) == 1:
        return (pos, f"{ref}→{alt}", "substitution", "—")
    if len(ref) == len(alt):
        return (pos, f"{ref}→{alt}", f"MNP ({len(ref)} bp)", "—")
    if dl > 0 and alt.startswith(ref):
        return (pos, f"+{alt[len(ref):]}", f"insertion (+{dl} bp)", _frame_by_length(dl))
    if dl < 0 and ref.startswith(alt):
        return (pos, f"−{ref[len(alt):]}", f"deletion (−{-dl} bp)", _frame_by_length(-dl))
    return (pos, f"{ref}→{alt}", "delins", _frame_by_length(abs(dl)))


def raw_calls(scn):
    """The raw variant calls as the caller emits them — no annotation at all.
    This is what you have *before* get_mnv."""
    out = []
    if scn.ivar_records:
        for r in scn.ivar_records:
            a = r.alt
            if a.startswith("+"):
                out.append((r.pos, f"+{a[1:]}", f"insertion (+{len(a) - 1} bp)",
                            _frame_by_length(len(a) - 1)))
            elif a.startswith("-"):
                out.append((r.pos, f"−{a[1:]}", f"deletion (−{len(a) - 1} bp)",
                            _frame_by_length(len(a) - 1)))
            else:
                out.append((r.pos, f"{r.ref}→{a}", "substitution", "—"))
    else:
        for v in scn.variants:
            for alt in v.alt.split(","):
                out.append(classify_vcf(v.pos, v.ref, alt))
    return out


def getmnv_adds(rows):
    """One honest line summarising what get_mnv contributes beyond the raw call."""
    if not rows:
        return "the variant is filtered out / intergenic — get_mnv emits no row here"
    cls = {r.get("Event Class", "") for r in rows}
    ct = {r.get("Change Type", "") for r in rows}
    adds = []
    if "complex_indel" in cls:
        adds.append("joins the indel with the in-cis SNV(s) into one complex_indel event")
    if "mnv" in cls or any(r.get("Variant Type") == "SNP/MNV" for r in rows):
        adds.append("groups co-occurring SNVs in the same codon into a single MNV codon")
    if any("(fs)" in (r.get("AA Changes", "") or "") for r in rows):
        adds.append("propagates the frameshift downstream — marks later variants (fs)")
    if ("Stop gained" in ct or "Stop lost" in ct) and any(r.get("Variant Type") == "INDEL" for r in rows):
        which = "Stop gained" if "Stop gained" in ct else "Stop lost"
        adds.append(f"detects the in-frame indel as {which} (a length-based guess would say 'in-frame')")
    if not adds:
        aa = next((r.get("AA Changes", "") for r in rows if r.get("AA Changes")), "")
        cty = next((r.get("Change Type", "") for r in rows if r.get("Change Type")), "")
        if any(r.get("Variant Type") == "INDEL" for r in rows):
            adds.append(f"classifies the indel ({cty}) and resolves its amino-acid effect"
                        + (f" ({aa})" if aa else ""))
        else:
            adds.append("resolves the codon-level amino-acid change" + (f" ({aa})" if aa else ""))
    return "; ".join(adds[:3])


# ----------------------------------------------------------------------------
# Read-map drawing (with a small alignment pileup)
# ----------------------------------------------------------------------------
def draw_readmap(ax, view, features, gene_color, snv_map):
    x0, x1 = view["x0"], view["x1"]
    contig = view["contig"]
    ref = F.reference_for(contig)
    reads = view["reads"]

    n_groups = max(1, len(reads))
    per_cap = max(3, READ_BUDGET // n_groups)
    display_n = [min(rg.count, per_cap) for rg in reads]
    n_inst = sum(display_n) if reads else 1

    base_top = 5.7 + max(0, n_inst - 1) * READ_PITCH
    y_aa = base_top              # codon / amino-acid + reading-frame track
    y_ref = base_top - 1
    y_gene = base_top - 2
    y_var = base_top - 4
    read_base = base_top - 5

    ax.set_xlim(x0 - 0.5, x1 + 0.5)
    ax.set_ylim(0.2, base_top + 1.1)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", labelsize=9, length=3)
    ax.set_xlabel(f"genomic position · {contig}", fontsize=10)

    trans = blended_transform_factory(ax.transAxes, ax.transData)
    span = (x1 - x0) or 1
    charw = span * 0.016

    for it in view["inputs"]:
        ax.axvline(it["pos"], color=ACCENT, lw=0.6, ls=(0, (2, 2)), alpha=0.30, zorder=0)

    # ---- reading frame (codon triplets) + amino-acid track ----
    codons = codon_track(view, features, snv_map) if view["show_bases"] else []
    codons.sort(key=lambda c: c["lo"])
    ax.text(-0.012, y_aa, "codon / aa", transform=trans, ha="right", va="center",
            fontsize=9, color="#444", fontweight="bold")
    for i, c in enumerate(codons):
        # triplet shading (the reading frame), behind everything
        if i % 2 == 0:
            ax.add_patch(Rectangle((c["lo"] - 0.5, 0.2), (c["hi"] - c["lo"] + 1),
                                   y_ref + 0.45 - 0.2, facecolor="#eef0f4",
                                   edgecolor="none", zorder=0))
        mid = (c["lo"] + c["hi"]) / 2
        if not (c["lo"] >= x0 and c["hi"] <= x1):
            continue  # partial codon at the window edge → shade only, no AA label
        nonsyn = c["hit"] and c["alt_aa"] != c["ref_aa"]
        if nonsyn:
            ax.text(mid, y_aa + 0.20, c["ref_aa"], ha="center", va="center", fontsize=8.5,
                    fontweight="bold", color="#333")
            ax.text(mid, y_aa - 0.28, c["alt_aa"], ha="center", va="center", fontsize=8.5,
                    fontweight="bold", color=DEL_COLOR)
            ax.annotate("", xy=(mid, y_aa - 0.13), xytext=(mid, y_aa + 0.05),
                        arrowprops=dict(arrowstyle="-|>", color=DEL_COLOR, lw=0.8))
        else:
            col = "#C9820B" if c["hit"] else "#333"   # synonymous hit = amber
            ax.text(mid, y_aa, c["ref_aa"], ha="center", va="center", fontsize=8.5,
                    fontweight="bold", color=col)

    # reference bases
    ax.text(-0.012, y_ref, "ref", transform=trans, ha="right", va="center",
            fontsize=9, color="#444", fontweight="bold")
    if view["show_bases"]:
        fs = 9 if span <= 40 else 7
        for p in range(x0, x1 + 1):
            b = ref[p - 1]
            ax.text(p, y_ref, b, ha="center", va="center", fontsize=fs,
                    family="monospace", color=BASE_COLOR.get(b, "#333"))
    else:
        ax.plot([x0, x1], [y_ref, y_ref], color="#bbb", lw=1.2)

    # genes / CDS
    ax.text(-0.012, y_gene, "genes", transform=trans, ha="right", va="center",
            fontsize=9, color="#444", fontweight="bold")
    drawn = [f for f in features if f["contig"] == contig and f["end"] >= x0 and f["start"] <= x1]
    by_parent: dict[str, list[dict]] = {}
    for f in drawn:
        by_parent.setdefault(f["parent"] or f["name"], []).append(f)
    for blocks in by_parent.values():
        blocks.sort(key=lambda b: b["start"])
        for a, b in zip(blocks, blocks[1:]):
            ax.plot([min(max(a["end"], x0), x1), max(min(b["start"], x1), x0)],
                    [y_gene, y_gene], color=gene_color, lw=0.8, ls=":", zorder=1)
    for f in drawn:
        gs, ge = max(f["start"], x0), min(f["end"], x1)
        ax.add_patch(Rectangle((gs - 0.5, y_gene - 0.28), (ge - gs + 1), 0.56,
                               facecolor=gene_color, edgecolor="none", alpha=0.45, zorder=2))
        if f["strand"] == "+":
            ax.annotate("", xy=(min(ge, gs + (ge - gs) * 0.6) + 0.5, y_gene), xytext=(gs, y_gene),
                        arrowprops=dict(arrowstyle="-|>", color=gene_color, lw=1.1))
        else:
            ax.annotate("", xy=(gs - 0.5, y_gene), xytext=(min(ge, gs + (ge - gs) * 0.4), y_gene),
                        arrowprops=dict(arrowstyle="-|>", color=gene_color, lw=1.1))
        ax.text((gs + ge) / 2, y_gene + 0.46, f["name"] or f["parent"], ha="center",
                va="bottom", fontsize=8.5, fontweight="bold", color="#222", clip_on=True)

    # input variants
    ax.text(-0.012, y_var, "input", transform=trans, ha="right", va="center",
            fontsize=9, color="#444", fontweight="bold")
    inputs = sorted(view["inputs"], key=lambda it: it["pos"])
    boxes = [(it["pos"] - max(len(it["label"]), 3) * charw / 2,
              it["pos"] + max(len(it["label"]), 3) * charw / 2) for it in inputs]
    for it, lvl in zip(inputs, assign_levels(boxes)):
        ax.add_patch(Polygon([(it["pos"], y_var - 0.18), (it["pos"] - 0.34, y_var - 0.54),
                              (it["pos"] + 0.34, y_var - 0.54)], closed=True,
                             facecolor=ACCENT, edgecolor="none", zorder=4))
        ly = y_var + 0.34 + lvl * 0.46
        ax.plot([it["pos"], it["pos"]], [y_var - 0.05, ly - 0.12], color=ACCENT, lw=0.6,
                alpha=0.7, zorder=3)
        ax.text(it["pos"], ly, it["label"], ha="center", va="bottom", fontsize=8.5,
                color="white", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.22", fc=ACCENT, ec="none"), zorder=5)

    # ---- read pileup ----
    inst = 0
    for rg, dn in zip(reads, display_n):
        snv = {o.pos: o.seq for o in rg.ops if o.kind == "snv"}
        ins = {o.pos: o.seq for o in rg.ops if o.kind == "ins"}
        dels = {o.pos: o.length for o in rg.ops if o.kind == "del"}
        skips = {o.pos: o.length for o in rg.ops if o.kind == "skip"}
        rs = max(rg.start, x0)
        re = min(rg.start + rg.length - 1, x1)
        ys = [read_base - (inst + k) * READ_PITCH for k in range(dn)]
        inst += dn
        top_y = ys[0]

        for k, y in enumerate(ys):
            is_top = (k == 0)
            seg_start = None
            p = rs
            while p <= re:
                if p in skips:
                    if seg_start is not None:
                        ax.add_patch(Rectangle((seg_start - 0.5, y - READ_H / 2),
                                               (p - seg_start), READ_H, facecolor=READ_FILL,
                                               edgecolor=READ_EDGE, lw=0.5, zorder=3))
                        seg_start = None
                    ln = skips[p]
                    ax.plot([p - 0.5, p + ln - 0.5], [y, y], color="#9aa6b2", lw=0.9,
                            ls=(0, (4, 2)), zorder=2)
                    if is_top:
                        ax.text((2 * p + ln - 1) / 2, top_y + 0.34, "intron (N)", ha="center",
                                va="bottom", fontsize=7, color="#6b7884", style="italic")
                    p += ln
                    continue
                if p in dels:
                    if seg_start is not None:
                        ax.add_patch(Rectangle((seg_start - 0.5, y - READ_H / 2),
                                               (p - seg_start), READ_H, facecolor=READ_FILL,
                                               edgecolor=READ_EDGE, lw=0.5, zorder=3))
                        seg_start = None
                    ln = dels[p]
                    ax.add_patch(Rectangle((p - 0.5, y - READ_H / 2), ln, READ_H, facecolor="white",
                                           edgecolor=DEL_COLOR, lw=0.7, ls="--", hatch="////", zorder=3))
                    if is_top:
                        ax.text((2 * p + ln - 1) / 2, top_y - 0.38, f"del {ln}bp", ha="center",
                                va="top", fontsize=7, color=DEL_COLOR, fontweight="bold")
                    p += ln
                    continue
                if seg_start is None:
                    seg_start = p
                p += 1
            if seg_start is not None:
                ax.add_patch(Rectangle((seg_start - 0.5, y - READ_H / 2), (re - seg_start + 1),
                                       READ_H, facecolor=READ_FILL, edgecolor=READ_EDGE,
                                       lw=0.5, zorder=3))

            for sp, alt in snv.items():
                if x0 <= sp <= x1:
                    ax.add_patch(Rectangle((sp - 0.5, y - READ_H / 2), 1, READ_H,
                                           facecolor=BASE_COLOR.get(alt, "#333"),
                                           edgecolor="white", lw=0.4, zorder=4))
                    if is_top and span <= 48:
                        ax.text(sp, top_y, alt, ha="center", va="center", fontsize=7,
                                family="monospace", color="white", fontweight="bold", zorder=5)
            for sp, seq in ins.items():
                if x0 <= sp <= x1:
                    ax.plot([sp + 0.5], [y + READ_H / 2], marker="v", ms=4, color=INS_COLOR, zorder=5)

        ins_items = sorted(ins.items())
        ins_boxes = [(sp + 0.5 - max(len(seq) + 1, 3) * charw / 2,
                      sp + 0.5 + max(len(seq) + 1, 3) * charw / 2) for sp, seq in ins_items]
        for (sp, seq), lvl in zip(ins_items, assign_levels(ins_boxes)):
            if x0 <= sp <= x1:
                ax.text(sp + 0.5, top_y + 0.32 + lvl * 0.34, f"+{seq}", ha="center", va="bottom",
                        fontsize=7, color="white", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.18", fc=INS_COLOR, ec="none"), zorder=6)

        cy = sum(ys) / len(ys)
        extra = f"\n{dn}/{rg.count} shown" if dn < rg.count else ""
        ax.text(-0.012, cy, f"{rg.name_prefix}\n×{rg.count} ({rg.strand}){extra}",
                transform=trans, ha="right", va="center", fontsize=7.5, color="#333")


# ----------------------------------------------------------------------------
# Tables
# ----------------------------------------------------------------------------
def _style_table(tab, header_color, row_lines, header_lines=1.5):
    total = sum(row_lines) + header_lines
    tab.auto_set_font_size(False)
    tab.set_fontsize(7.6)
    for (r, c), cell in tab.get_celld().items():
        cell.set_linewidth(0.4)
        cell.get_text().set_verticalalignment("center")
        cell.PAD = 0.03
        if r == 0:
            cell.set_facecolor(header_color)
            cell.set_height(header_lines / total)
            t = cell.get_text(); t.set_color("white"); t.set_fontweight("bold")
        else:
            cell.set_height(row_lines[r - 1] / total)
            cell.set_facecolor("#f3f6fb" if r % 2 else "white")
            cell.set_edgecolor("#d4dde8")


def draw_raw_table(ax, calls):
    ax.axis("off")
    ax.set_title("Without get_mnv — raw variant calls (caller output: position + allele, no annotation)",
                 fontsize=10.5, fontweight="bold", loc="left", color="#555", pad=4)
    cell_text = [[str(p), allele, cls, frame] for (p, allele, cls, frame) in calls]
    row_lines = [1] * len(cell_text)
    tab = ax.table(cellText=cell_text, colLabels=RAW_COLS, colWidths=RAW_W,
                   cellLoc="left", loc="upper center", bbox=[0.0, 0.0, 1.0, 0.90])
    _style_table(tab, NAIVE_HDR, row_lines)


def draw_adds_note(ax, note):
    ax.axis("off")
    ax.text(0.0, 0.5, "▶  What get_mnv adds:  " + textwrap.fill(note, 100),
            transform=ax.transAxes, ha="left", va="center", fontsize=10,
            color="#1d5e2a", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.5", fc="#e9f6ec", ec="#3FA45B", lw=0.8))


def draw_getmnv_table(ax, rows):
    ax.axis("off")
    ax.set_title("get_mnv output — codon-aware grouping", fontsize=10.5, fontweight="bold",
                 loc="left", color=ACCENT, pad=4)
    if not rows:
        ax.text(0.5, 0.4, "get_mnv produced no rows for this scenario (e.g. filtered or intergenic)",
                ha="center", va="center", fontsize=9.5, style="italic", color="#666")
        return
    cell_text, row_lines = [], []
    for r in rows:
        cells, maxl = [], 1
        for key in GETMNV_COLS:
            v = (r.get(key, "") or "").strip()
            if key == "Event Components":
                v = "\n".join(v.split(" | "))
            elif key == "AA Changes":
                v = textwrap.fill(v, 15) if v else ""
            elif key == "Change Type":
                v = textwrap.fill(v, 11) if v else ""
            elif key == "Variant Type":
                v = textwrap.fill(v, 8) if v else ""
            cells.append(v)
            maxl = max(maxl, v.count("\n") + 1)
        cell_text.append(cells)
        row_lines.append(maxl)
    tab = ax.table(cellText=cell_text, colLabels=GETMNV_COLS, colWidths=GETMNV_W,
                   cellLoc="left", loc="upper center", bbox=[0.0, 0.0, 1.0, 0.90])
    _style_table(tab, ACCENT, row_lines)


# ----------------------------------------------------------------------------
# Per-scenario page
# ----------------------------------------------------------------------------
def make_figure(scn, rows, cat_tag=None, cat_color=ACCENT):
    gff = scn.gff_content or F.GFF_GENE_ONLY
    want = "CDS" if (scn.gff_features and "CDS" in scn.gff_features) else "gene"
    features = [f for f in parse_gff(gff) if f["type"] == want]

    views = build_views(scn)
    names = []
    for f in features:
        nm = f["name"] or f["parent"]
        if nm not in names:
            names.append(nm)
    gcol = {nm: GENE_COLORS[i % len(GENE_COLORS)] for i, nm in enumerate(names)}

    def color_for_view(v):
        for f in features:
            if f["contig"] == v["contig"] and f["end"] >= v["x0"] and f["start"] <= v["x1"]:
                return gcol.get(f["name"] or f["parent"], GENE_COLORS[0])
        return GENE_COLORS[0]

    n_groups_per_view = [max(1, len(v["reads"])) for v in views]
    inst_per_view = []
    for v in views:
        ng = max(1, len(v["reads"]))
        cap = max(3, READ_BUDGET // ng)
        inst_per_view.append(sum(min(rg.count, cap) for rg in v["reads"]) or 1)

    # single-base substitutions → alt base, used to label the alt amino acid
    snv_map = {}
    if scn.ivar_records:
        for r in scn.ivar_records:
            if not r.alt.startswith(("+", "-")) and len(r.ref) == 1 and len(r.alt) == 1:
                snv_map[(r.chrom, r.pos)] = r.alt
    else:
        for v in scn.variants:
            if len(v.ref) == 1:
                for a in v.alt.split(","):
                    if len(a) == 1:
                        snv_map[(v.chrom, v.pos)] = a

    raw = raw_calls(scn)
    adds = getmnv_adds(rows)
    view_ratios = [4.6 + 0.78 * n for n in inst_per_view]
    raw_ratio = max(1.5, 0.9 + 0.5 * max(1, len(raw)))
    note_ratio = 0.9
    getmnv_ratio = max(1.8, 0.9 + 0.62 * max(1, len(rows)))

    total_ratio = sum(view_ratios) + raw_ratio + note_ratio + getmnv_ratio
    fig_h = min(17.0, max(7.8, 1.9 + 0.40 * total_ratio))
    fig = plt.figure(figsize=(12.4, fig_h))
    gs = fig.add_gridspec(len(views) + 3, 1,
                          height_ratios=view_ratios + [raw_ratio, note_ratio, getmnv_ratio],
                          left=0.135, right=0.985, top=0.885, bottom=0.05, hspace=0.55)

    for i, v in enumerate(views):
        draw_readmap(fig.add_subplot(gs[i, 0]), v, features, color_for_view(v), snv_map)
    draw_raw_table(fig.add_subplot(gs[len(views), 0]), raw)
    draw_adds_note(fig.add_subplot(gs[len(views) + 1, 0]), adds)
    draw_getmnv_table(fig.add_subplot(gs[len(views) + 2, 0]), rows)

    fig.suptitle(scn.name, fontsize=16, fontweight="bold", x=0.135, ha="left", y=0.975)
    if cat_tag:
        fig.text(0.985, 0.975, cat_tag, ha="right", va="center", fontsize=9, color="white",
                 fontweight="bold", bbox=dict(boxstyle="round,pad=0.35", fc=cat_color, ec="none"))
    flags = ""
    if scn.gff_features:
        flags += f"   ·  --gff-features {scn.gff_features}"
    if scn.extra_cli_args:
        flags += "   ·  " + " ".join(scn.extra_cli_args)
    if scn.ivar_records:
        flags += "   ·  input: iVar TSV"
    desc = textwrap.fill(DESC.get(scn.name, scn.description), 125)
    fig.text(0.135, 0.918, desc + flags, fontsize=9.5, color="#444", va="top")
    return fig


def make_legend_page():
    fig = plt.figure(figsize=(11.7, 8.3))
    ax = fig.add_axes([0.06, 0.06, 0.88, 0.86])
    ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    fig.suptitle("Diagram legend", fontsize=18, fontweight="bold", x=0.06, ha="left")
    y = 9.4

    def row(draw, text):
        nonlocal y
        draw(y)
        ax.text(2.45, y, text, va="center", fontsize=11)
        y -= 0.98

    row(lambda y: ax.add_patch(Rectangle((0.3, y - 0.18), 1.4, 0.36, facecolor=READ_FILL,
        edgecolor=READ_EDGE, lw=0.7)),
        "Mapped read (aligned segment, CIGAR M). Reads are drawn as a small pileup; the left label is "
        "prefix  ×total (strand).")
    row(lambda y: ax.add_patch(Rectangle((0.6, y - 0.18), 0.4, 0.36, facecolor=BASE_COLOR["A"],
        edgecolor="white")),
        "SNV on a read: block coloured by the alternate base (A green, C blue, G amber, T red).")
    row(lambda y: ax.add_patch(Rectangle((0.6, y - 0.18), 0.8, 0.36, facecolor="white",
        edgecolor=DEL_COLOR, lw=0.7, ls="--", hatch="////")),
        "Deletion: red hatched gap over the read (labelled 'del Nbp').")
    row(lambda y: (ax.plot([0.6], [y + 0.22], marker="v", ms=7, color=INS_COLOR),
                   ax.text(0.62, y + 0.5, "+SEQ", ha="center", fontsize=7, color="white",
                           bbox=dict(boxstyle="round,pad=0.15", fc=INS_COLOR, ec="none"))),
        "Insertion: purple triangle between bases with the inserted sequence (+SEQ).")
    row(lambda y: ax.plot([0.3, 1.7], [y, y], color="#9aa6b2", lw=1.0, ls=(0, (4, 2))),
        "Intron (CIGAR N): dashed grey line; the read skips that region (multi-exon genes).")
    row(lambda y: ax.add_patch(Polygon([(1.0, y + 0.16), (0.7, y - 0.16), (1.3, y - 0.16)],
        closed=True, facecolor=ACCENT)),
        "Input variant (VCF REF→ALT, or iVar +/-SEQ): blue triangle on the 'input' track.")
    row(lambda y: (ax.add_patch(Rectangle((0.3, y - 0.16), 1.4, 0.32, facecolor=GENE_COLORS[0],
        alpha=0.45)),
        ax.annotate("", xy=(1.5, y), xytext=(0.4, y),
                    arrowprops=dict(arrowstyle="-|>", color=GENE_COLORS[0], lw=1.1))),
        "Gene/CDS overlapping the window; the arrow shows the strand (+ right, − left).")
    row(lambda y: (ax.add_patch(Rectangle((0.30, y - 0.18), 0.45, 0.36, facecolor="#eef0f4", ec="#ccc", lw=0.5)),
                   ax.add_patch(Rectangle((0.75, y - 0.18), 0.45, 0.36, facecolor="white", ec="#ccc", lw=0.5)),
                   ax.add_patch(Rectangle((1.20, y - 0.18), 0.45, 0.36, facecolor="#eef0f4", ec="#ccc", lw=0.5)),
                   ax.text(0.97, y + 0.46, "Ala", ha="center", fontsize=8.5, fontweight="bold")),
        "Reading frame & amino acid: alternating shaded triplets are codons; the 'codon/aa' track shows the "
        "reference amino acid, with a red AA below (↓) when a variant makes the codon non-synonymous.")

    ax.text(0.06, y - 0.05,
            "Tables:  'Without get_mnv' shows the raw variant calls as a caller emits them — position, "
            "allele and a length-based frame guess, with no annotation.\n'get_mnv output' is the real "
            "codon-aware result; the green 'What get_mnv adds' line summarises the difference (MNV / "
            "complex-indel grouping, frameshift propagation, in-frame stop detection, amino-acid effect).\n\n"
            "Notes:  the window is centred on the region of interest (variants ± margin); if it spans "
            ">58 bp the reference letters are hidden.\nFor genes on the '−' strand the codons reported "
            "by get_mnv are the genomic-forward sequence (documented behaviour).",
            fontsize=10, color="#555", va="top")
    return fig


def make_title_page(n, n_sections):
    fig = plt.figure(figsize=(11.7, 8.3))
    fig.text(0.5, 0.62, "get_MNV — validation scenarios", ha="center", fontsize=22, fontweight="bold")
    fig.text(0.5, 0.54, f"Read mapping and annotated output · {n} scenarios in {n_sections} sections",
             ha="center", fontsize=15, color="#444")
    fig.text(0.5, 0.47, "indels branch · generated from tests/scenarios/", ha="center",
             fontsize=10, color="#777")
    return fig


def make_toc_page(categories):
    fig = plt.figure(figsize=(11.7, 8.3))
    fig.text(0.08, 0.92, "Contents", fontsize=20, fontweight="bold")
    fig.text(0.08, 0.873, "Scenarios are grouped by case, not by raw number.", fontsize=11.5, color="#666")
    y = 0.80
    for i, (title, short, blurb, codes) in enumerate(categories, 1):
        color = CAT_COLORS[(i - 1) % len(CAT_COLORS)]
        fig.text(0.08, y, f"{i}.", fontsize=14, fontweight="bold", color=color)
        fig.text(0.12, y, title, fontsize=14, fontweight="bold")
        fig.text(0.12, y - 0.032, "scenarios:  " + ",  ".join(codes), fontsize=11, color="#666")
        y -= 0.099
    return fig


def make_section_page(num, total, title, blurb, scns, color):
    fig = plt.figure(figsize=(11.7, 8.3))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 0.035, 1.0, color=color, transform=ax.transAxes))
    fig.text(0.085, 0.85, f"Section {num} of {total}", fontsize=13, color=color, fontweight="bold")
    fig.text(0.085, 0.80, title, fontsize=21, fontweight="bold", va="top")
    fig.text(0.085, 0.71, textwrap.fill(blurb, 84), fontsize=12, color="#444", va="top")
    y = 0.58
    for scn in scns:
        fig.text(0.10, y, "●", fontsize=9, color=color, va="center")
        fig.text(0.125, y, scn.name, fontsize=12, fontweight="bold", va="center")
        fig.text(0.125, y - 0.034, textwrap.fill(DESC.get(scn.name, scn.description), 90),
                 fontsize=9.5, color="#555", va="top")
        y -= 0.086
    return fig


def main(argv):
    png_filters, args = [], list(argv)
    while "--png" in args:
        i = args.index("--png")
        png_filters.append(args[i + 1])
        del args[i:i + 2]

    base_work = HERE / "work_plots"
    if base_work.exists():
        shutil.rmtree(base_work)
    base_work.mkdir(parents=True)

    # Group scenarios by section; warn about any not assigned to a category.
    code_to_scn = {scn.name.split("_", 1)[0]: scn for scn in ALL_SCENARIOS}
    assigned = {c for (_, _, _, codes) in CATEGORIES for c in codes}
    sections = [(title, short, blurb, [code_to_scn[c] for c in codes if c in code_to_scn])
                for (title, short, blurb, codes) in CATEGORIES]
    leftover = [scn for code, scn in code_to_scn.items() if code not in assigned]
    if leftover:
        sections.append(("Other scenarios", "Other", "Scenarios not assigned to a section.", leftover))

    out_pdf = HERE / "scenarios_overview.pdf"
    ok, failed = 0, []

    def render(scn, cat_short, color):
        nonlocal ok
        _, _, out = run_scenario(scn, base_work)
        rows = parse_tsv(out)[1]
        fig = make_figure(scn, rows, cat_tag=f"§{cat_short}", cat_color=color)
        pdf.savefig(fig)
        if any(scn.name.startswith(f) for f in png_filters):
            fig.savefig(HERE / f"plot_{scn.name}.png", dpi=130)
        plt.close(fig)
        ok += 1
        print(f"  ok  {scn.name}  (get_mnv {len(rows)} row/s)")

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(make_title_page(len(ALL_SCENARIOS), len(sections))); plt.close()
        pdf.savefig(make_legend_page()); plt.close()
        pdf.savefig(make_toc_page(CATEGORIES)); plt.close()
        for si, (title, short, blurb, scns) in enumerate(sections, 1):
            color = CAT_COLORS[(si - 1) % len(CAT_COLORS)]
            pdf.savefig(make_section_page(si, len(sections), title, blurb, scns, color)); plt.close()
            print(f"— Section {si}: {title}")
            for scn in scns:
                try:
                    render(scn, short, color)
                except Exception as e:  # noqa: BLE001
                    failed.append((scn.name, str(e)))
                    print(f"  FAIL  {scn.name}: {e}")

    shutil.rmtree(base_work, ignore_errors=True)
    print(f"\nPDF: {out_pdf}  ({ok} scenario pages, {len(failed)} failures)")
    for name, err in failed:
        print(f"  - {name}: {err}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
