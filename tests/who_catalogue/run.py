#!/usr/bin/env python3
"""The WHO catalogue converter, checked on real rows of the catalogue.

tools/who_catalogue_to_tsv.py turns the WHO catalogue of M. tuberculosis
resistance mutations (2023) into a table for get_MNV's `--catalogue`. A
converter can be wrong in two ways that a look at its output does not show: it
can write an entry get_MNV reads but never matches, and it can drift away from
the table shipped in example/. This checks both, without the WHO files:

  * conversion: real rows of the master file (fixtures/master_rows.tsv), one
    for each way an entry is rewritten, give exactly the expected entries, and
    the whole-gene deletion is listed as not converted, with its reason;
  * matching: WHO's own genomic variants for those rows
    (fixtures/coordinates_rows.tsv) are annotated on the example's MTBC
    ancestor with the converted table, and each comes out flagged with the
    entry WHO links it to; rpoB p.Ser450Phe and p.Ser441Gln, which need two
    changes in one codon, through the codon reading only;
  * the shipped table: example/who_catalogue_2023.tsv carries the licence
    notice, agrees with the converter entry by entry, and loads in get_MNV; the
    example sample, which carries no graded resistance variant, flags nothing.

With --who-dir naming a folder that holds the two WHO files, it also rebuilds
the shipped table and requires it byte for byte, and annotates every genomic
variant WHO links to groups 1 and 2, failing if more of them go unflagged than
the ones README.md explains.

Usage:
    cargo build -p get_mnv
    python3 tests/who_catalogue/run.py
    python3 tests/who_catalogue/run.py --who-dir ~/Downloads/who-2023
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FIXTURES = HERE / "fixtures"
EXAMPLE = ROOT / "example"
SHIPPED = EXAMPLE / "who_catalogue_2023.tsv"
FASTA = EXAMPLE / "MTB_ancestor.fas"
GENES = EXAMPLE / "anot_genes.txt"
MASTER_FILE = "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
COORDINATES_FILE = "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"

# WHO genomic variants of groups 1 and 2 that get_MNV, run with the shipped
# table on the MTBC ancestor, flags for none of their drugs. Every one is
# explained in README.md; more than this many is a regression.
UNFLAGGED_IN_GROUPS_1_2 = 4

spec = importlib.util.spec_from_file_location("converter", ROOT / "tools" / "who_catalogue_to_tsv.py")
converter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(converter)

# What the fixture's rows become, all five groups kept: gene, variant, label.
EXPECTED = [
    ("rrs", "n.1401A>G", "AMK: Assoc w R"),
    ("rrl", "g.1475926_1475927insT", "CAP: Uncertain significance, WHO n.2269_2270insT"),
    ("rrs", "n.1401A>G", "CAP: Assoc w R"),
    ("aftB", "g.4266954dup", "EMB: Uncertain significance, WHO p.Ter628ext*?"),
    ("inhA", "c.-777C>T", "ETO: Assoc w R"),
    ("Rv0565c", "g.[656672G>C;656674C>T]", "ETO: Uncertain significance, WHO p.TrpLeu266*"),
    ("inhA", "c.-777C>T", "INH: Assoc w R"),
    ("katG", "c.[421T>C;423G>C]", "INH: Not assoc w R - Interim, WHO c.421_423delTTGinsCTC"),
    ("katG", "LoF", "INH: Assoc w R"),
    ("katG", "p.Ala16fs", "INH: Assoc w R - Interim"),
    ("katG", "p.Ser315Thr", "INH: Assoc w R"),
    ("rrs", "n.1401A>G", "KAN: Assoc w R"),
    ("gyrA", "p.Asp94Gly", "LFX: Assoc w R"),
    ("gyrA", "p.Glu21Gln", "LFX: Not assoc w R"),
    ("rrl", "g.1475926_1475927insT", "LZD: Assoc w R - Interim, WHO n.2269_2270insT"),
    ("gyrA", "p.Asp94Gly", "MFX: Assoc w R"),
    ("gyrA", "p.Glu21Gln", "MFX: Not assoc w R"),
    ("pncA", "g.2289247del", "PZA: Assoc w R, WHO c.-5delG"),
    ("pncA", "p.Met1?", "PZA: Assoc w R"),
    ("pncA", "p.Ter187Arg", "PZA: Assoc w R - Interim, WHO p.Ter187Argext*?"),
    ("rpoB", "p.Asn438del", "RIF: Assoc w R - Interim"),
    ("rpoB", "p.Met434Ile", "RIF: Assoc w R - Interim"),
    ("rpoB", "p.Phe433_Met434insPhe", "RIF: Assoc w R, WHO p.Phe433dup"),
    ("rpoB", "p.Ser441Gln", "RIF: Assoc w R"),
    ("rpoB", "p.Ser450Leu", "RIF: Assoc w R"),
    ("rpoB", "p.Ser450Phe", "RIF: Assoc w R"),
    ("rrs", "n.1401A>G", "STM: Uncertain significance"),
]

# Genomic variants of the fixture that get_MNV cannot flag with their entry,
# and why. WHO writes this change across codons 266 and 267 as one stop; get_MNV
# reads each codon on its own row (Trp266*, Leu267Val), so no row carries both
# substitutions that the entry's `g.` names.
EXPECTED_MISSES = {
    ("Rv0565c_p.TrpLeu266*", 656672, "GCC", "CCT"),
}

# Two changes in one codon: matched through the codon, never through one base.
CODON_ONLY = {"rpoB_p.Ser450Phe", "rpoB_p.Ser441Gln"}

failures: list[str] = []


def check(condition: bool, message: str) -> None:
    if not condition:
        failures.append(message)


def newest_source_mtime() -> float:
    return max(p.stat().st_mtime for p in (ROOT / "src").rglob("*") if p.is_file())


def binary() -> Path:
    built = [
        ROOT / relative
        for relative in ("target/debug/get_mnv", "target/release/get_mnv")
        if (ROOT / relative).exists()
    ]
    if not built:
        sys.exit("no get_mnv binary in target/; build it first (cargo build -p get_mnv)")
    built.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    if built[0].stat().st_mtime < newest_source_mtime():
        sys.exit(f"{built[0]} is older than src/; rebuild it (cargo build -p get_mnv)")
    return built[0]


def read_table(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def convert_cli(*args: str) -> tuple[list[dict], list[dict], str]:
    """Run the converter as a user would, and read what it wrote."""
    with tempfile.TemporaryDirectory() as work:
        out, unconverted = Path(work) / "out.tsv", Path(work) / "unconverted.tsv"
        run = subprocess.run(
            [sys.executable, str(ROOT / "tools" / "who_catalogue_to_tsv.py"), *args,
             "-o", str(out), "--unconverted", str(unconverted)],
            capture_output=True, text=True, check=True,
        )
        return read_table(out), read_table(unconverted), out.read_text()


def test_helpers() -> None:
    check(converter.grade_of("2) Assoc w R - Interim") == (2, "Assoc w R - Interim"), "grade_of")
    check(converter.parse_grades("1,2") == {1, 2}, "parse_grades 1,2")
    check(converter.parse_grades("1-3") == {1, 2, 3}, "parse_grades 1-3")
    check(converter.parse_grades("all") == {1, 2, 3, 4, 5}, "parse_grades all")
    for bad in ("0", "6", "2-1", "x", ""):
        try:
            converter.parse_grades(bad)
            failures.append(f"parse_grades accepted {bad!r}")
        except argparse.ArgumentTypeError:
            pass
    rewrite = converter.CODING_DELINS.fullmatch("c.421_423delTTGinsCTC")
    check(converter.coding_brackets(rewrite) == "c.[421T>C;423G>C]", "coding_brackets")
    single = converter.CODING_DELINS.fullmatch("c.10_12delAGCinsAGT")
    check(converter.coding_brackets(single) == "c.12C>T", "coding_brackets, one base")

    residues = {433: "Phe", 434: "Met", 68: "His", 69: "Gln", 70: "Leu", 71: "Ala"}
    check(converter.convert("p.Phe433dup", residues) == ["p.Phe433_Met434insPhe"], "dup of one")
    check(
        converter.convert("p.His68_Leu70dup", residues) == ["p.Leu70_Ala71insHisGlnLeu"],
        "dup of a run",
    )
    try:
        converter.convert("p.Phe433dup", {433: "Phe"})
        failures.append("a duplication beside an unnamed residue was converted")
    except converter.Unconverted as why:
        check("--coordinates" in str(why), f"dup without residues: {why}")
    check(converter.convert("p.Ter187Argext*?") == ["p.Ter187Arg"], "stop loss")
    for kept in ("p.Ser450Leu", "p.Trp122*", "p.Asp379fs", "p.Met1?", "p.Asn438del",
                 "p.Gln432_Asp435delinsHis", "p.Leu430_Ser431insArg", "c.-15C>T",
                 "n.1401A>G", "c.102G>A", "LoF"):
        check(converter.convert(kept) == [kept], f"{kept} is kept as written")
    try:
        converter.convert("deletion")
        failures.append("a whole-gene deletion was converted")
    except converter.NeedsGenome:
        failures.append("a whole-gene deletion was sent to the genome")
    except converter.Unconverted:
        pass

    # get_MNV's own `HGVS g.` on the example reference: the pncA promoter
    # deletion sits in a run of two C and is written at the 3' one, the aftB
    # insertion repeats the base before it, the rrl one does not.
    sequence = converter.read_fasta(FASTA)["MTB_anc"]
    for (position, ref, alt), expected in {
        (2289245, "TC", "T"): "g.2289247del",
        (4266953, "T", "TC"): "g.4266954dup",
        (1475926, "A", "AT"): "g.1475926_1475927insT",
        (761155, "CG", "TT"): "g.[761155C>T;761156G>T]",
        (761155, "C", "T"): "g.761155C>T",
    }.items():
        got = converter.genomic_descriptor(position, ref, alt, sequence)
        check(got == expected, f"genomic_descriptor {position} {ref}>{alt}: {got}, want {expected}")
    toy = "AAACTGCTGGGG"
    check(converter.genomic_descriptor(3, "ACTG", "A", toy) == "g.7_9del", "deletion in a repeat")
    check(converter.genomic_descriptor(3, "ACT", "AG", toy) == "g.4_5delinsG", "delins stays")
    check(converter.genomic_descriptor(3, "ACT", "AGG", toy) == "g.[4C>G;5T>G]", "same length")


def test_conversion() -> list[dict]:
    master = str(FIXTURES / "master_rows.tsv")
    coordinates = ["--coordinates", str(FIXTURES / "coordinates_rows.tsv"), "--fasta", str(FASTA)]
    table, unconverted, text = convert_cli(master, "--grades", "all", *coordinates)
    got = [(row["gene"], row["variant"], row["label"]) for row in table]
    check(got == EXPECTED, "fixture conversion:\n  " + "\n  ".join(
        f"{'+' if entry not in EXPECTED else '-' if entry not in got else ' '} {entry}"
        for entry in sorted(set(got) | set(EXPECTED))
        if entry not in got or entry not in EXPECTED
    ))
    check(
        [(row["who_variant"], row["reason"].split(":")[0]) for row in unconverted]
        == [("pncA_deletion", "whole-gene deletion")],
        f"unconverted: {unconverted}",
    )
    check("ODC Attribution License" in text and converter.LICENCE_URL in text, "licence notice")
    for row in table:
        check(row["who_variant"].startswith(row["gene"] + "_"), f"who_variant of {row}")

    # The default keeps groups 1 and 2 only.
    default, _, _ = convert_cli(master, *coordinates)
    check(
        {row["grade"][0] for row in default} == {"1", "2"} and len(default) == 20,
        f"default grades: {len(default)} entries, groups {sorted({r['grade'] for r in default})}",
    )
    # Without the coordinates, what needs the genome is listed, not dropped.
    plain, missing, _ = convert_cli(master, "--grades", "all")
    check(
        sorted(row["who_variant"] for row in missing)
        == sorted(["rrl_n.2269_2270insT"] * 2 + ["aftB_p.Ter628ext*?", "Rv0565c_p.TrpLeu266*",
                   "pncA_c.-5delG", "pncA_deletion"]),
        f"unconverted without coordinates: {[row['who_variant'] for row in missing]}",
    )
    check(len(plain) == len(EXPECTED) - 5, f"converted without coordinates: {len(plain)}")
    return table


def batches_of(records: list[tuple]) -> list[list[tuple]]:
    """Genomic variants in batches that cannot interact: at most one per gene,
    and every two at least 300 bp apart, so none share a codon or shift the
    frame of another."""
    batches: list[list[tuple]] = []
    # The first batch each gene may still join; it only moves forward, so a
    # gene never lands twice in one batch.
    next_free: dict[str, int] = defaultdict(int)
    for record in sorted(records):
        position, ref = record[0], record[1]
        gene = record[3].split("_")[0]
        index = next_free[gene]
        while True:
            if index == len(batches):
                batches.append([])
            batch = batches[index]
            if all(position + len(ref) + 300 < other[0] or other[0] + len(other[1]) + 300 < position
                   for other in batch):
                batch.append(record)
                next_free[gene] = index + 1
                break
            index += 1
    return batches


def annotate(get_mnv: Path, catalogue: Path, records: list[tuple], vcf: Path | None = None) -> list[dict]:
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        if vcf is None:
            vcf = work / "calls.vcf"
            lines = [
                "##fileformat=VCFv4.2",
                "##contig=<ID=MTB_anc,length=4411532>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
            ]
            lines += [f"MTB_anc\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t1/1" for p, r, a, *_ in sorted(records)]
            vcf.write_text("\n".join(lines) + "\n")
        run = subprocess.run(
            [str(get_mnv), "--vcf", str(vcf), "--fasta", str(FASTA), "--genes", str(GENES),
             "--catalogue", str(catalogue), "--output-dir", str(work), "--output-prefix", "run"],
            capture_output=True, text=True,
        )
        if run.returncode != 0:
            raise RuntimeError(run.stderr.strip().splitlines()[-1])
        return read_table(work / "run.MNV.tsv")


def flagged(rows: list[dict], position: int, ref: str, needles: set[str]) -> tuple[bool, str]:
    """Whether the rows of one genomic variant carry one of its entries, and
    the note get_MNV gave the match."""
    near = [row for row in rows
            if any(position - 150 <= int(p) <= position + len(ref) + 150
                   for p in row["Positions"].split(", "))]
    cells = " ; ".join(row["Catalogue"] for row in near)
    for needle in needles:
        match = re.search(re.escape(needle) + r"(?:, ([^)]*))?\)", cells)
        if match:
            return True, match.group(1) or ""
    return False, cells


def match_records(get_mnv: Path, catalogue: Path, coordinates: Path, wanted=None):
    """Every WHO genomic variant of the table's entries, annotated with the
    table: (who_variant, position, ref, alt) -> (flagged, note or cells)."""
    entries = read_table(catalogue)
    needles: dict[str, set[str]] = defaultdict(set)
    for row in entries:
        needles[row["who_variant"]].add(f"({row['gene']} {row['variant']}")
    sequence = converter.read_fasta(FASTA)["MTB_anc"]
    linked: dict[tuple, set[str]] = defaultdict(set)
    for who, records in converter.read_coordinates(coordinates, set(needles)).items():
        if wanted is not None and who not in wanted:
            continue
        for position, ref, alt in records:
            if sequence[position - 1 : position - 1 + len(ref)] == ref:
                linked[(position, ref, alt)].add(who)
    records = [(p, r, a, sorted(whos)[0], whos) for (p, r, a), whos in linked.items()]
    results: dict[tuple, tuple[bool, str]] = {}

    def run(batch):
        return batch, annotate(get_mnv, catalogue, batch)

    with ThreadPoolExecutor(max(1, min(8, os.cpu_count() or 1))) as pool:
        for batch, rows in pool.map(run, batches_of(records)):
            for position, ref, alt, _, whos in batch:
                for who in whos:
                    results[(who, position, ref, alt)] = flagged(rows, position, ref, needles[who])
    return results


def test_matching(get_mnv: Path) -> None:
    with tempfile.TemporaryDirectory() as work:
        table = Path(work) / "fixture.tsv"
        subprocess.run(
            [sys.executable, str(ROOT / "tools" / "who_catalogue_to_tsv.py"),
             str(FIXTURES / "master_rows.tsv"), "--grades", "all",
             "--coordinates", str(FIXTURES / "coordinates_rows.tsv"), "--fasta", str(FASTA),
             "-o", str(table)],
            capture_output=True, text=True, check=True,
        )
        results = match_records(get_mnv, table, FIXTURES / "coordinates_rows.tsv")
    check(len(results) == 52, f"{len(results)} fixture genomic variants annotated, want 52")
    for (who, position, ref, alt), (hit, note) in sorted(results.items()):
        expected_miss = (who, position, ref, alt) in EXPECTED_MISSES
        if expected_miss:
            check(not hit, f"{who} {position} {ref}>{alt} is now flagged; update EXPECTED_MISSES")
        elif not hit:
            failures.append(f"{who} {position} {ref}>{alt} not flagged; its rows say: {note}")
        elif who in CODON_ONLY and len(ref) == len(alt) and sum(a != b for a, b in zip(ref, alt)) == 2:
            check(note == "codon reading only", f"{who} {position} {ref}>{alt}: note {note!r}")


def test_shipped_table(get_mnv: Path) -> None:
    text = SHIPPED.read_text(encoding="utf-8")
    check(converter.LICENCE_URL in text and converter.SOURCE_URL in text, "shipped table: notice")
    entries = read_table(SHIPPED)
    check(len(entries) == 1379, f"shipped table: {len(entries)} entries, want 1379")
    check({row["grade"][0] for row in entries} == {"1", "2"}, "shipped table: groups 1 and 2")
    residues: dict[str, dict[int, str]] = {}
    for row in entries:
        mutation = row["who_variant"].split("_", 1)[1]
        for residue, position in converter.NAMED_RESIDUE.findall(mutation):
            residues.setdefault(row["gene"], {}).setdefault(int(position), residue)
    drifted = []
    for row in entries:
        mutation = row["who_variant"].split("_", 1)[1]
        try:
            ours = converter.name_level(mutation, residues.get(row["gene"], {}))
        except converter.NeedsGenome:
            ours = None
        if ours is None:
            # A duplication may take its neighbour from a group this table
            # leaves out; anything else with no form of its own is a `g.`.
            if not (row["variant"].startswith("g.") or converter.DUPLICATION.fullmatch(mutation)):
                drifted.append(f"{row['who_variant']}: {row['variant']}, want a g. entry")
        elif row["variant"] not in ours:
            drifted.append(f"{row['who_variant']}: {row['variant']}, the converter writes {ours}")
        label = f"{converter.DRUGS[row['drug']]}: {row['grade'].split(') ', 1)[1]}"
        if not row["label"].startswith(label):
            drifted.append(f"{row['who_variant']}: label {row['label']!r}")
    check(not drifted, "shipped table drifted from the converter:\n  " + "\n  ".join(drifted[:20]))

    rows = annotate(get_mnv, SHIPPED, [], vcf=EXAMPLE / "G35894.var.snp.vcf")
    check(len(rows) == 941, f"example run with the shipped table: {len(rows)} rows, want 941")
    hits = [row for row in rows if row["Catalogue"] != "-"]
    check(not hits, f"the example sample flags {len(hits)} rows: {[r['Catalogue'] for r in hits[:3]]}")


def test_with_who_files(get_mnv: Path, folder: Path) -> None:
    master, coordinates = folder / MASTER_FILE, folder / COORDINATES_FILE
    for path in (master, coordinates):
        if not path.exists():
            sys.exit(f"--who-dir: {path} is missing")
    with tempfile.TemporaryDirectory() as work:
        rebuilt = Path(work) / SHIPPED.name
        subprocess.run(
            [sys.executable, str(ROOT / "tools" / "who_catalogue_to_tsv.py"), str(master),
             "--coordinates", str(coordinates), "--fasta", str(FASTA), "-o", str(rebuilt)],
            capture_output=True, text=True, check=True,
        )
        check(rebuilt.read_bytes() == SHIPPED.read_bytes(),
              f"{SHIPPED.relative_to(ROOT)} differs from what the converter writes from {MASTER_FILE}")
    results = match_records(get_mnv, SHIPPED, coordinates)
    drugs: dict[str, set[str]] = defaultdict(set)
    for row in read_table(SHIPPED):
        drugs[row["who_variant"]].add(row["label"].split(":")[0])
    by_variant: dict[tuple, list[tuple[str, bool, str]]] = defaultdict(list)
    for (who, position, ref, alt), (hit, note) in results.items():
        by_variant[(position, ref, alt)].append((who, hit, note))
    own = sum(hit for hit, _ in results.values())
    codon_only = sum(hit and note == "codon reading only" for hit, note in results.values())
    unflagged = []
    for (position, ref, alt), outcomes in sorted(by_variant.items()):
        for who, hit, note in outcomes:
            if not hit and not any(f"{drug}: " in note for drug in drugs[who]):
                unflagged.append(f"{who} {position} {ref[:20]}>{alt[:20]}: {note[:100]}")
    print(f"  {len(results)} WHO genomic variants of groups 1 and 2 annotated: {own} flagged with "
          f"their own entry ({codon_only} through the codon reading only), "
          f"{len(unflagged)} flagged for none of their drugs")
    check(len(unflagged) <= UNFLAGGED_IN_GROUPS_1_2,
          f"{len(unflagged)} unflagged, more than the {UNFLAGGED_IN_GROUPS_1_2} explained:\n  "
          + "\n  ".join(unflagged))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--who-dir", type=Path,
                        help=f"folder holding {MASTER_FILE} and {COORDINATES_FILE}")
    args = parser.parse_args()
    get_mnv = binary()
    steps = [("converter helpers", test_helpers),
             ("conversion of real rows", test_conversion),
             ("matching WHO's genomic variants", lambda: test_matching(get_mnv)),
             ("the shipped table", lambda: test_shipped_table(get_mnv))]
    if args.who_dir:
        steps.append(("the WHO files", lambda: test_with_who_files(get_mnv, args.who_dir)))
    for name, step in steps:
        before = len(failures)
        step()
        print(f"{'PASS' if len(failures) == before else 'FAIL'}  {name}")
    for failure in failures:
        print(f"  - {failure}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
