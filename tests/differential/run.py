#!/usr/bin/env python3
"""Differential annotation suite: get_MNV against independent oracles.

get_MNV's own tests can only check that it produces what we wrote down as
expected. When the expectation is wrong the test goes green and the tool stays
wrong, which is exactly how a ribosomal-slippage CDS join came to be annotated
with splice sites: the fixture that was supposed to prove splicing worked built
two exons with no intron between them. This suite takes us out of the loop by
comparing against annotators we did not write.

Two oracles, chosen for what each is good at:

* `bcftools csq`, run live on the small synthetic fixtures. Those are where the
  structural edge cases live (slippage joins, real introns) and where we control
  the gene model exactly. Both tools read the SAME reference and the SAME GFF3,
  so a disagreement can never be blamed on them seeing different annotation.
* SnpEff, already present in the bundled M. tuberculosis VCF as `ANN=` fields
  from the original 4.1l run. Real data, a real annotator and a real database,
  with no GFF conversion in between. Here the two tools do use their own
  annotation of the same genome, so gene-model differences are part of what the
  baseline records.

The rule is NOT "the tools must agree". get_MNV disagrees on purpose: it reports
codon-level MNVs that per-variant annotators split apart, and it annotates
intergenic variants that some annotators omit. Every disagreement is matched
against a versioned baseline of accepted differences, and the suite fails only
when a difference appears that nobody has explained yet.

Usage:
    python3 tests/differential/run.py            # compare against the baseline
    python3 tests/differential/run.py --update   # rewrite the baseline

Needs `bcftools` on PATH for the fixture datasets; those are skipped with a
clear message when it is missing.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FIXTURES = HERE / "fixtures"
WORK = HERE / "work"
BASELINE = HERE / "baseline.tsv"

NONE = "-"  # "this tool reported nothing here"

# Genuine naming differences between vocabularies, applied after the shared
# `_variant` suffix has been stripped. Kept deliberately tiny: anything longer
# is a synonym table that becomes its own source of mistakes.
ALIASES = {"intergenic_region": "intergenic"}


# --------------------------------------------------------------------------
# datasets
# --------------------------------------------------------------------------
@dataclass
class Dataset:
    name: str
    fasta: Path
    vcf: Path
    oracle: str  # "bcftools" or "snpeff"
    gff: Path | None  # shared with the oracle; None means get_MNV uses `genes`
    genes: Path | None
    note: str


DATASETS = [
    Dataset(
        "slippage",
        FIXTURES / "slippage.fa",
        FIXTURES / "slippage.vcf",
        "bcftools",
        FIXTURES / "slippage.gff",
        None,
        "ribosomal-slippage CDS join, the SARS-CoV-2 ORF1ab shape: no intron, so no splice sites",
    ),
    Dataset(
        "spliced",
        FIXTURES / "spliced.fa",
        FIXTURES / "spliced.vcf",
        "bcftools",
        FIXTURES / "spliced.gff",
        None,
        "a genuinely spliced two-exon transcript: the positive control, splice terms must fire here",
    ),
    Dataset(
        "mtb",
        ROOT / "example" / "MTB_ancestor.fas",
        ROOT / "example" / "G35894.var.snp.vcf",
        "snpeff",
        None,
        ROOT / "example" / "anot_genes.txt",
        "the bundled M. tuberculosis dataset: 950 real calls already annotated by SnpEff 4.1l",
    ),
]


# --------------------------------------------------------------------------
# normalisation
# --------------------------------------------------------------------------
def stem(term: str) -> str:
    """Reduce a consequence to a vocabulary the tools share.

    get_MNV emits full Sequence Ontology names (`missense_variant`) while
    `bcftools csq` emits the short form (`missense`), so dropping a trailing
    `_variant` lines them up without a hand-maintained synonym table.
    """
    term = term.strip().lstrip("*&")
    if term.endswith("_variant"):
        term = term[: -len("_variant")]
    return ALIASES.get(term, term)


@dataclass
class Calls:
    """Per-position consequence terms and impact, as one tool sees them."""

    terms: dict[int, set[str]]
    impact: dict[int, str]


def get_mnv_calls(tsv: Path) -> Calls:
    terms: dict[int, set[str]] = defaultdict(set)
    impact: dict[int, str] = {}
    with tsv.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row_terms = {stem(t) for t in row["SO Term"].split("&") if t.strip()}
            for token in row["Positions"].split(","):
                token = token.strip()
                if not token.isdigit():
                    continue
                position = int(token)
                terms[position] |= row_terms
                impact[position] = row.get("Impact", "").strip() or NONE
    return Calls(terms, impact)


BCSQ_POS = re.compile(r"(\d+)[ACGTNacgtn*]+>")


def bcftools_calls(annotated: Path) -> Calls:
    """Parse a `bcftools csq` annotated VCF.

    A BCSQ entry is `consequence|gene|transcript|biotype[|strand|aa|dna]`. The
    seven-field form carries the DNA change of the whole haplotype
    (`401C>A+402C>A`), so its consequence is attributed to every position it
    names, which is what makes a get_MNV MNV row comparable to a haplotype call.
    Shorter entries apply to the record itself. `@POS` back-references are
    skipped: the record they point at already names the positions it covers.

    bcftools csq has no impact column, so impact is not compared here.
    """
    terms: dict[int, set[str]] = defaultdict(set)
    with annotated.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            position, info = int(fields[1]), fields[7]
            terms.setdefault(position, set())
            for chunk in info.split(";"):
                if not chunk.startswith("BCSQ="):
                    continue
                for entry in chunk[len("BCSQ=") :].split(","):
                    if not entry or entry.startswith("@"):
                        continue
                    parts = entry.split("|")
                    entry_terms = {stem(t) for t in parts[0].split("&") if stem(t)}
                    if len(parts) >= 7 and parts[6]:
                        for touched in BCSQ_POS.findall(parts[6]) or [position]:
                            terms[int(touched)] |= entry_terms
                    else:
                        terms[position] |= entry_terms
    return Calls(terms, {})


def snpeff_calls(vcf: Path) -> Calls:
    """Parse the `ANN=` fields SnpEff already wrote into the bundled VCF.

    SnpEff orders its annotations by severity, so the first one is the call it
    stands behind. Entries from a custom interval file (`custom|MODIFIER|`,
    carrying no gene) are not consequence predictions and are dropped.
    """
    terms: dict[int, set[str]] = defaultdict(set)
    impact: dict[int, str] = {}
    with vcf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            position, info = int(fields[1]), fields[7]
            terms.setdefault(position, set())
            impact.setdefault(position, NONE)
            for chunk in info.split(";"):
                if not chunk.startswith("ANN="):
                    continue
                for entry in chunk[len("ANN=") :].split(","):
                    parts = entry.split("|")
                    if len(parts) < 3 or parts[1] == "custom":
                        continue
                    terms[position] |= {stem(t) for t in parts[1].split("&") if stem(t)}
                    impact[position] = parts[2].strip() or NONE
                    break  # the first real annotation is SnpEff's primary call
    return Calls(terms, impact)


def key_of(terms: set[str]) -> str:
    return ",".join(sorted(terms)) if terms else NONE


# --------------------------------------------------------------------------
# running the tools
# --------------------------------------------------------------------------
def binary() -> Path:
    """The get_mnv build under test.

    `GET_MNV_BIN` overrides the search, which is how you point the suite at an
    older build to check that it still catches a bug it is supposed to catch.
    """
    override = os.environ.get("GET_MNV_BIN")
    if override:
        return Path(override)
    for candidate in (
        ROOT / "target" / "release" / "get_mnv",
        ROOT / "target" / "debug" / "get_mnv",
    ):
        if candidate.exists():
            return candidate
    sys.stderr.write("get_mnv binary not found; run `cargo build -p get_mnv` first\n")
    raise SystemExit(1)


def run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(
            f"\ncommand failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stdout}\n{proc.stderr}\n"
        )
        raise SystemExit(1)
    return proc


def annotate(dataset: Dataset, work: Path) -> tuple[Calls, Calls]:
    case = work / dataset.name
    case.mkdir(parents=True, exist_ok=True)
    local_vcf = case / "input.vcf"
    shutil.copy(dataset.vcf, local_vcf)

    cmd = [str(binary()), "--vcf", "input.vcf", "--fasta", str(dataset.fasta.resolve())]
    if dataset.gff is not None:
        cmd += ["--gff", str(dataset.gff.resolve())]
    else:
        cmd += ["--genes", str(dataset.genes.resolve())]
    run(cmd, cwd=case)
    mine = get_mnv_calls(case / "input.MNV.tsv")

    if dataset.oracle == "bcftools":
        proc = run(
            [
                "bcftools", "csq",
                "-f", str(dataset.fasta.resolve()),
                "-g", str(dataset.gff.resolve()),
                str(local_vcf),
            ]
        )
        annotated = case / "csq.vcf"
        annotated.write_text(proc.stdout)
        (case / "csq.stderr").write_text(proc.stderr)
        return mine, bcftools_calls(annotated)

    return mine, snpeff_calls(local_vcf)


# --------------------------------------------------------------------------
# baseline
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Difference:
    dataset: str
    field: str  # "so" or "impact"
    position: int
    mine: str
    theirs: str

    def rule(self) -> tuple[str, str, str, str]:
        return (self.dataset, self.field, self.mine, self.theirs)


def load_baseline() -> tuple[dict[tuple[str, str, str, str], str], set[tuple[str, str, int]]]:
    rules: dict[tuple[str, str, str, str], str] = {}
    sites: set[tuple[str, str, int]] = set()
    if not BASELINE.exists():
        return rules, sites
    with BASELINE.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["kind"] == "rule":
                rules[(row["dataset"], row["field"], row["get_mnv"], row["oracle"])] = row["note"]
            elif row["kind"] == "site":
                sites.add((row["dataset"], row["field"], int(row["position"])))
    return rules, sites


def write_baseline(differences: list[Difference]) -> None:
    counts: dict[tuple[str, str, str, str], int] = defaultdict(int)
    for diff in differences:
        counts[diff.rule()] += 1
    with BASELINE.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["kind", "dataset", "field", "position", "get_mnv", "oracle", "count", "note"])
        for (dataset, field, mine, theirs), count in sorted(counts.items()):
            writer.writerow(
                ["rule", dataset, field, "*", mine, theirs, count, "TODO explain"]
            )


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def compare(dataset: Dataset, mine: Calls, theirs: Calls) -> tuple[list[Difference], int]:
    differences: list[Difference] = []
    positions = sorted(set(mine.terms) | set(theirs.terms))
    for position in positions:
        mine_key = key_of(mine.terms.get(position, set()))
        their_key = key_of(theirs.terms.get(position, set()))
        if mine_key != their_key:
            differences.append(Difference(dataset.name, "so", position, mine_key, their_key))
        # Impact is only comparable when the oracle reports one.
        if theirs.impact:
            mine_impact = mine.impact.get(position, NONE)
            their_impact = theirs.impact.get(position, NONE)
            if mine_impact != their_impact:
                differences.append(
                    Difference(dataset.name, "impact", position, mine_impact, their_impact)
                )
    return differences, len(positions)


def main() -> int:
    parser = argparse.ArgumentParser(description="Differential annotation suite")
    parser.add_argument("--update", action="store_true", help="rewrite the baseline from what is observed now")
    args = parser.parse_args()

    have_bcftools = shutil.which("bcftools") is not None
    if have_bcftools:
        # Recorded because a different bcftools can legitimately shift the
        # comparison, and a failure here should be diagnosable at a glance.
        version = subprocess.run(
            ["bcftools", "--version"], capture_output=True, text=True
        ).stdout.splitlines()
        print(f"oracle: {version[0] if version else 'bcftools (version unknown)'}")
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)

    differences: list[Difference] = []
    skipped: list[str] = []
    print(f"{'dataset':<10} {'oracle':<9} {'positions':>9} {'agree':>7} {'differ':>7}")
    for dataset in DATASETS:
        if dataset.oracle == "bcftools" and not have_bcftools:
            skipped.append(dataset.name)
            continue
        mine, theirs = annotate(dataset, WORK)
        found, positions = compare(dataset, mine, theirs)
        differences += found
        so_diffs = len({d.position for d in found if d.field == "so"})
        print(f"{dataset.name:<10} {dataset.oracle:<9} {positions:>9} {positions - so_diffs:>7} {len(found):>7}")

    if skipped:
        print(f"\nSKIPPED (bcftools not on PATH): {', '.join(skipped)}")

    if args.update:
        write_baseline(differences)
        print(f"\nbaseline rewritten from {len(differences)} differences; explain every note before committing")
        return 0

    rules, sites = load_baseline()
    unexplained = [
        d for d in differences
        if d.rule() not in rules and (d.dataset, d.field, d.position) not in sites
    ]
    observed = {d.rule() for d in differences}
    stale = [r for r in rules if r not in observed and r[0] not in skipped]

    print(f"\n{len(differences)} difference(s), {len(differences) - len(unexplained)} explained by the baseline")

    if stale:
        print("\nNOTE: baseline entries no longer observed. If a fix removed them, delete them from baseline.tsv:")
        for dataset, field, mine_key, their_key in sorted(stale):
            print(f"  {dataset:<10} {field:<7} get_MNV={mine_key}  oracle={their_key}")

    if unexplained:
        print(f"\nFAIL: {len(unexplained)} unexplained difference(s).")
        print("Either get_MNV has a bug, or the behaviour changed on purpose and the baseline needs an entry.\n")
        print(f"  {'dataset':<10} {'field':<7} {'pos':>9}  {'get_MNV':<40} oracle")
        for diff in unexplained[:40]:
            print(f"  {diff.dataset:<10} {diff.field:<7} {diff.position:>9}  {diff.mine:<40} {diff.theirs}")
        if len(unexplained) > 40:
            print(f"  ... and {len(unexplained) - 40} more")
        return 1

    print("\nOK: every difference is a known, documented one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
