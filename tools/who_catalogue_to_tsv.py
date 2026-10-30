#!/usr/bin/env python3
"""Turn the WHO catalogue of M. tuberculosis resistance mutations into a table
for get_MNV's `--catalogue`.

The catalogue (WHO, second edition, 2023) grades about 48,000 variant and drug
pairs into five groups, from "1) Assoc w R" to "5) Not assoc w R". It names a
variant the way get_MNV mostly does already: `rpoB_p.Ser450Leu`,
`inhA_c.-777C>T`, `rrs_n.1401A>G`, `katG_LoF`. Most entries therefore carry
over as written, split at the first underscore into a gene and a variant, and
the rest are rewritten into the form get_MNV writes for the same change:

  * a same-length change in a coding sequence (`c.421_423delTTGinsCTC`) is the
    substitutions it is made of, `c.[421T>C;423G>C]`;
  * an in-frame duplication (`p.Phe433dup`) is the insertion get_MNV names
    between the residues around it, `p.Phe433_Met434insPhe`, and a lost stop
    (`p.Ter187Argext*?`) is the residue that replaces it, `p.Ter187Arg`. The
    residues come from the catalogue itself, whose entries name the reference
    residue at nearly every position of its genes;
  * an insertion or deletion outside a coding sequence (`c.-5delG` upstream of
    pncA, `n.2269_2270insT` in rrl) has no `c.` descriptor in get_MNV, so it is
    written as the `g.` one get_MNV gives its row: placed from WHO's genomic
    coordinates and moved 3'-most along the reference, as HGVS asks. That needs
    the coordinates file and the FASTA of the run (`--coordinates`, `--fasta`),
    which also place the few entries the rewrites above cannot name.

Matching protein changes rather than genomic records is what lets a codon with
two substitutions be found: rpoB p.Ser450Phe needs two changes in codon 450,
and get_MNV reads them as one codon where a base-by-base annotator reads
p.Ser450Leu and a synonymous change.

What cannot be converted is listed with the reason, never dropped silently:
the whole-gene deletions (`pncA_deletion`), which no VCF record of a variant
inside the gene describes, and anything written in a form get_MNV does not
produce. `--unconverted` writes that list as a table.

The label of each entry is the drug's WHO abbreviation and the grade, `RIF:
Assoc w R`, which is what get_MNV prints in its `Catalogue` column; the other
columns keep the drug, the grade and WHO's own variant name, which get_MNV
ignores and a reader can join on.

The input is the master file as WHO publishes it in text form,
`WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt`, from
https://github.com/GTB-tbsequencing/mutation-catalogue-2023 (folder `Final
Result Files`). Those data are made available under the Open Data Commons
Attribution License, and the table written here carries that notice.

Usage:
    python3 tools/who_catalogue_to_tsv.py MASTER_FILE -o who_catalogue.tsv
    python3 tools/who_catalogue_to_tsv.py MASTER_FILE --grades 1,2,5 \\
        --coordinates WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt \\
        --fasta example/MTB_ancestor.fas -o who_catalogue.tsv

Standard library only, Python 3.7 or later.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import re
import sys
from collections import Counter
from pathlib import Path

# The drugs as the catalogue names them, with the abbreviations its report uses.
DRUGS = {
    "Amikacin": "AMK",
    "Bedaquiline": "BDQ",
    "Capreomycin": "CAP",
    "Clofazimine": "CFZ",
    "Delamanid": "DLM",
    "Ethambutol": "EMB",
    "Ethionamide": "ETO",
    "Isoniazid": "INH",
    "Kanamycin": "KAN",
    "Levofloxacin": "LFX",
    "Linezolid": "LZD",
    "Moxifloxacin": "MFX",
    "Pyrazinamide": "PZA",
    "Rifampicin": "RIF",
    "Streptomycin": "STM",
}

# The columns read from the master file. Its other hundred-odd columns are the
# statistics behind each grade.
DRUG, GENE, MUTATION, VARIANT, GRADE = (
    "drug",
    "gene",
    "mutation",
    "variant",
    "FINAL CONFIDENCE GRADING",
)
NEEDED = (DRUG, GENE, MUTATION, VARIANT, GRADE)

LICENCE_URL = "https://opendatacommons.org/licenses/by/1-0/"
SOURCE_URL = "https://github.com/GTB-tbsequencing/mutation-catalogue-2023"

RESIDUE = r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)"

# The protein forms get_MNV writes itself, compared as written once it has
# read them: a substitution or a premature stop, a frameshift (matched by its
# first residue whatever follows `fs`), a lost start, and an in-frame deletion,
# insertion or deletion-insertion named by the residues around it.
PROTEIN_FORMS = (
    re.compile(rf"p\.{RESIDUE}\d+(?:{RESIDUE}|\*|Ter|=)"),
    re.compile(rf"p\.{RESIDUE}\d+(?:{RESIDUE})?fs(?:\*\d*|Ter\d*|\*\?)?"),
    re.compile(r"p\.Met1\?"),
    re.compile(rf"p\.{RESIDUE}\d+(?:_{RESIDUE}\d+)?del"),
    re.compile(rf"p\.{RESIDUE}\d+(?:_{RESIDUE}\d+)?delins(?:{RESIDUE})+"),
    re.compile(rf"p\.{RESIDUE}\d+_{RESIDUE}\d+ins(?:{RESIDUE})+"),
)

# An in-frame duplication, of one residue or of a run of them.
DUPLICATION = re.compile(rf"p\.({RESIDUE})(\d+)(?:_({RESIDUE})(\d+))?dup")

# A lost stop that names the residue taking its place. WHO also writes losses
# that name none (`p.Ter628ext*?`), which leave nothing to compare with.
STOP_LOSS = re.compile(rf"p\.Ter(\d+)({RESIDUE})ext\*\?")

# Where a protein change names a residue: `p.Ser450Leu`, both ends of a range.
NAMED_RESIDUE = re.compile(rf"(?:p\.|_)({RESIDUE})(\d+)")

# A single-base substitution, numbered along a gene: `c.` for one that is
# translated, `n.` for an rRNA, a minus sign upstream of the first base.
NUCLEOTIDE_SUBSTITUTION = re.compile(r"[cn]\.-?\d+[ACGT]>[ACGT]")

# A same-length change spelt as a deletion-insertion, inside a coding sequence.
CODING_DELINS = re.compile(r"c\.(\d+)_(\d+)del([ACGT]+)ins([ACGT]+)")

# An insertion, deletion or duplication numbered along a gene. WHO writes those
# inside a coding sequence as protein changes, so these are upstream
# (`c.-30_-29insGCG`) or in an rRNA (`n.1077dupT`).
NUCLEOTIDE_INDEL = re.compile(r"[cn]\.-?\d+(?:_-?\d+)?(?:del|ins|dup)[ACGT]*(?:ins[ACGT]+)?")


class Unconverted(Exception):
    """An entry that has no get_MNV form, and why."""


class NeedsGenome(Unconverted):
    """An entry with no get_MNV form of its own that WHO's genomic variants can
    stand in for: each is written as the `g.` get_MNV gives its row."""


def grade_of(text: str) -> tuple[int, str]:
    """`1) Assoc w R` as (1, "Assoc w R")."""
    match = re.match(r"\s*([1-5])\)\s*(.+?)\s*$", text)
    if not match:
        raise ValueError(f"unrecognised confidence grading: {text!r}")
    return int(match.group(1)), match.group(2)


def parse_grades(text: str) -> set[int]:
    """`1,2`, `1-3` or `all` as the set of groups to keep."""
    if text.strip().lower() == "all":
        return {1, 2, 3, 4, 5}
    grades: set[int] = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        low, _, high = part.partition("-")
        try:
            span = range(int(low), int(high or low) + 1)
        except ValueError:
            raise argparse.ArgumentTypeError(f"not a grade or range of grades: {part!r}")
        if not span or min(span) < 1 or max(span) > 5:
            raise argparse.ArgumentTypeError(f"grades run from 1 to 5, got {part!r}")
        grades.update(span)
    if not grades:
        raise argparse.ArgumentTypeError("no grade given")
    return grades


def open_text(path: Path):
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path), encoding="utf-8", newline="")
    return open(path, encoding="utf-8", newline="")


def read_master(path: Path) -> list[dict]:
    """The master file's rows, from its header line on.

    The text WHO publishes starts with the header; a sheet saved from the Excel
    workbook has two title rows above it, and a fixture may have `#` comments.
    Whatever comes before the line naming the needed columns is skipped.
    """
    if path.suffix.lower() in (".xlsx", ".xls"):
        sys.exit(
            f"{path}: this reads the master file as text. WHO publishes it next to "
            "the workbook as WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
        )
    with open_text(path) as handle:
        lines = handle.read().splitlines()
    for at, line in enumerate(lines):
        columns = line.split("\t")
        if all(name in columns for name in NEEDED):
            break
    else:
        sys.exit(f"{path}: no header line naming the columns {', '.join(NEEDED)}")
    return list(csv.DictReader(lines[at:], delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current = None
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                sequences[current] = []
            elif current is not None:
                sequences[current].append(line.upper())
    return {name: "".join(parts) for name, parts in sequences.items()}


def read_coordinates(path: Path, wanted: set[str]) -> dict[str, list[tuple[int, str, str]]]:
    """WHO's genomic variants (VCF-style, on NC_000962.3) for the variants named."""
    found: dict[str, list[tuple[int, str, str]]] = {}
    with open_text(path) as handle:
        lines = (line for line in handle if not line.startswith("#"))
        for row in csv.DictReader(lines, delimiter="\t"):
            name = row.get("variant", "")
            if name in wanted:
                found.setdefault(name, []).append(
                    (
                        int(row["position"]),
                        row["reference_nucleotide"].upper(),
                        row["alternative_nucleotide"].upper(),
                    )
                )
    return found


def trim_common(ref: str, alt: str) -> tuple[int, int]:
    """The shared prefix, then the shared suffix of what is left, as get_MNV
    trims a REF/ALT pair."""
    most = min(len(ref), len(alt))
    prefix = 0
    while prefix < most and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    while suffix < most - prefix and ref[-1 - suffix] == alt[-1 - suffix]:
        suffix += 1
    return prefix, suffix


def genomic_descriptor(position: int, ref: str, alt: str, sequence: str) -> str:
    """The HGVS `g.` get_MNV writes for a VCF-style change on `sequence`.

    A deletion or insertion moves to its 3'-most equivalent place, and an
    insertion that repeats the bases before it is a `dup`, as HGVS asks and as
    get_MNV's own `HGVS g.` column does (src/variants/hgvs.rs). A same-length
    change is its substitutions, in brackets when there are several.
    """
    prefix, suffix = trim_common(ref, alt)
    ref_core = ref[prefix : len(ref) - suffix]
    alt_core = alt[prefix : len(alt) - suffix]
    start = position + prefix
    if len(ref_core) == len(alt_core):
        changed = [
            f"{start + i}{r}>{a}" for i, (r, a) in enumerate(zip(ref_core, alt_core)) if r != a
        ]
        if not changed:
            raise Unconverted("a genomic variant that changes nothing")
        return f"g.{changed[0]}" if len(changed) == 1 else f"g.[{';'.join(changed)}]"
    if not ref_core:
        after, bases = start - 1, alt_core
        while after < len(sequence) and sequence[after] == bases[0]:
            bases = bases[1:] + bases[0]
            after += 1
        length = len(bases)
        if after >= length and sequence[after - length : after] == bases:
            first = after - length + 1
            return f"g.{after}dup" if first == after else f"g.{first}_{after}dup"
        return f"g.{after}_{after + 1}ins{bases}"
    if not alt_core:
        first, last = start, start + len(ref_core) - 1
        while last < len(sequence) and sequence[first - 1] == sequence[last]:
            first += 1
            last += 1
        return f"g.{first}del" if first == last else f"g.{first}_{last}del"
    end = start + len(ref_core) - 1
    span = f"{start}" if start == end else f"{start}_{end}"
    return f"g.{span}delins{alt_core}"


def coding_brackets(match: re.Match) -> str:
    """`c.421_423delTTGinsCTC` as the substitutions it makes, `c.[421T>C;423G>C]`."""
    first, last, deleted, inserted = match.groups()
    first, last = int(first), int(last)
    if len(deleted) != len(inserted) or last - first + 1 != len(deleted):
        raise NeedsGenome("a coding deletion-insertion that changes the length")
    changed = [
        f"{first + i}{r}>{a}" for i, (r, a) in enumerate(zip(deleted, inserted)) if r != a
    ]
    if not changed:
        raise Unconverted("a deletion-insertion that changes nothing")
    return f"c.{changed[0]}" if len(changed) == 1 else f"c.[{';'.join(changed)}]"


def residue_map(rows: list[dict]) -> dict[str, dict[int, str]]:
    """Each gene's reference residues, as the catalogue's protein entries name
    them: `rpoB_p.Met434Ile` says residue 434 of rpoB is Met."""
    residues: dict[str, dict[int, str]] = {}
    for row in rows:
        mutation = row[MUTATION]
        if mutation.startswith("p."):
            for residue, position in NAMED_RESIDUE.findall(mutation):
                residues.setdefault(row[GENE], {}).setdefault(int(position), residue)
    return residues


def duplication_as_insertion(match: re.Match, residues: dict[int, str]) -> str:
    """`p.Phe433dup` as get_MNV writes it, `p.Phe433_Met434insPhe`.

    get_MNV names an in-frame insertion by the residues on either side of it
    and the residues inserted, never as a duplication, so the copied residues
    and the one after them have to be known.
    """
    first_residue, first, last_residue, last = match.groups()
    first = int(first)
    last = int(last) if last else first
    copied = [residues.get(position) for position in range(first, last + 1)]
    after = residues.get(last + 1)
    if None in copied or after is None:
        raise NeedsGenome("in-frame duplication beside a residue the catalogue names nowhere")
    if copied[0] != first_residue or copied[-1] != (last_residue or first_residue):
        raise NeedsGenome("in-frame duplication whose residues disagree with the catalogue's")
    return f"p.{copied[-1]}{last}_{after}{last + 1}ins{''.join(copied)}"


def convert(mutation: str, residues=None, place=None) -> list[str]:
    """The get_MNV variant(s) for one WHO mutation, or `Unconverted`.

    `residues` are the reference residues of the entry's gene, by position.
    `place` returns the `g.` descriptors of the entry's genomic variants, for
    one with no get_MNV form of its own; it is None when no coordinates and
    FASTA were given.
    """
    try:
        return name_level(mutation, residues or {})
    except NeedsGenome as why:
        if place is None:
            raise Unconverted(
                f"{why}: placed from WHO's genomic coordinates, which needs "
                "--coordinates and --fasta"
            )
        return place()


def name_level(mutation: str, residues: dict[int, str]) -> list[str]:
    """The entry rewritten from its name alone, or the reason it cannot be."""
    if mutation == "LoF":
        return ["LoF"]
    if mutation == "deletion":
        raise Unconverted(
            "whole-gene deletion: WHO gives it no genomic variant, and get_MNV reports "
            "changes inside a gene, not a gene missing from the sample"
        )
    if any(form.fullmatch(mutation) for form in PROTEIN_FORMS):
        return [mutation]
    duplication = DUPLICATION.fullmatch(mutation)
    if duplication:
        return [duplication_as_insertion(duplication, residues)]
    stop_loss = STOP_LOSS.fullmatch(mutation)
    if stop_loss:
        # get_MNV writes the lost stop as the residue replacing it (`*187Arg`).
        position, residue = stop_loss.groups()
        return [f"p.Ter{position}{residue}"]
    if mutation.startswith("p."):
        if "ext" in mutation or "???" in mutation:
            raise NeedsGenome("stop loss that names no residue in its place")
        raise NeedsGenome("protein change in a form get_MNV does not write")
    if NUCLEOTIDE_SUBSTITUTION.fullmatch(mutation):
        return [mutation]
    coding = CODING_DELINS.fullmatch(mutation)
    if coding:
        return [coding_brackets(coding)]
    if NUCLEOTIDE_INDEL.fullmatch(mutation):
        # get_MNV writes no `c.` for a change outside a coding sequence.
        raise NeedsGenome("indel outside a coding sequence")
    raise NeedsGenome("nucleotide change in a form get_MNV does not read")


class Placement:
    """WHO's genomic variants for the entries with no get_MNV form of their own,
    placed on the run's reference."""

    def __init__(self, coordinates: dict, sequence: str) -> None:
        self.coordinates = coordinates
        self.sequence = sequence

    def descriptors(self, who_variant: str) -> list[str]:
        records = self.coordinates.get(who_variant, [])
        if not records:
            raise Unconverted("no genomic variant in WHO's coordinates to place it by")
        descriptors: list[str] = []
        for position, ref, alt in records:
            # A record whose bases this FASTA does not have was written against
            # H37Rv where the two differ, and would describe another change here.
            if self.sequence[position - 1 : position - 1 + len(ref)] != ref:
                continue
            descriptor = genomic_descriptor(position, ref, alt, self.sequence)
            if descriptor not in descriptors:
                descriptors.append(descriptor)
        if not descriptors:
            raise Unconverted("its genomic variants' reference bases differ from the FASTA's")
        return descriptors


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the WHO catalogue of M. tuberculosis resistance mutations "
            "(2023) into a table for get_MNV --catalogue."
        )
    )
    parser.add_argument(
        "master",
        type=Path,
        help="WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt (or .txt.gz)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="where to write the table (default: standard output)",
    )
    parser.add_argument(
        "--grades",
        type=parse_grades,
        default=parse_grades("1,2"),
        help="confidence groups to keep: a list such as 1,2,5, a range such as 1-3, "
        "or all (default: 1,2, the groups associated with resistance)",
    )
    parser.add_argument(
        "--coordinates",
        type=Path,
        help="WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt, to place indels outside "
        "coding sequences (needs --fasta)",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        help="the reference of the run, on H37Rv coordinates (NC_000962.3 or the "
        "MTBC ancestor), to place those indels as get_MNV writes them",
    )
    parser.add_argument(
        "--unconverted",
        type=Path,
        help="also write the entries that could not be converted, with the reason",
    )
    args = parser.parse_args(argv)
    if (args.coordinates is None) != (args.fasta is None):
        parser.error("--coordinates and --fasta go together")

    rows = read_master(args.master)
    residues = residue_map(rows)
    kept = []
    for row in rows:
        number, name = grade_of(row[GRADE])
        if number in args.grades:
            kept.append((row, number, name))

    placement = None
    if args.coordinates is not None:
        sequences = read_fasta(args.fasta)
        if len(sequences) == 1:
            sequence = next(iter(sequences.values()))
        elif "NC_000962.3" in sequences:
            sequence = sequences["NC_000962.3"]
        else:
            parser.error(
                f"{args.fasta} holds {len(sequences)} sequences; give one on H37Rv "
                "coordinates, or one named NC_000962.3 among them"
            )
        wanted = {row[VARIANT] for row, _, _ in kept}
        placement = Placement(read_coordinates(args.coordinates, wanted), sequence)

    entries: list[tuple[str, ...]] = []
    seen: set[tuple[str, ...]] = set()
    unconverted: list[tuple[str, ...]] = []
    reasons: Counter = Counter()
    for row, number, name in kept:
        gene, mutation, drug = row[GENE], row[MUTATION], row[DRUG]
        abbreviation = DRUGS.get(drug, drug)
        grade = f"{number}) {name}"
        place = None
        if placement is not None:
            place = lambda who_variant=row[VARIANT]: placement.descriptors(who_variant)
        try:
            variants = convert(mutation, residues.get(gene), place)
        except Unconverted as why:
            reasons[str(why)] += 1
            unconverted.append((gene, mutation, drug, grade, row[VARIANT], str(why)))
            continue
        for variant in variants:
            label = f"{abbreviation}: {name}"
            if variant != mutation:
                # Written in get_MNV's form, so say what WHO called it.
                label += f", WHO {mutation}"
            entry = (gene, variant, label, drug, grade, row[VARIANT])
            if entry not in seen:
                seen.add(entry)
                entries.append(entry)

    groups = ", ".join(str(g) for g in sorted(args.grades))
    header = [
        "# WHO catalogue of mutations in Mycobacterium tuberculosis complex and their",
        "# association with drug resistance, second edition (2023), for get_MNV --catalogue.",
        f"# Contains information from the WHO mutation catalogue ({SOURCE_URL}),",
        f"# which is made available under the ODC Attribution License ({LICENCE_URL}).",
        f"# Converted by tools/who_catalogue_to_tsv.py from {args.master.name}",
        f"# (sha256 {sha256(args.master)}); groups {groups}.",
    ]
    if placement is not None:
        header.append(
            f"# Entries with no get_MNV form of their own placed on {args.fasta.name} "
            f"from {args.coordinates.name}."
        )
    out = open(args.output, "w", encoding="utf-8", newline="") if args.output else sys.stdout
    try:
        out.write("\n".join(header) + "\n")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "variant", "label", "drug", "grade", "who_variant"])
        writer.writerows(entries)
    finally:
        if args.output:
            out.close()

    if args.unconverted:
        with open(args.unconverted, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene", "mutation", "drug", "grade", "who_variant", "reason"])
            writer.writerows(unconverted)

    print(
        f"WHO catalogue: {len(kept)} entries in groups {groups}; "
        f"{len(kept) - len(unconverted)} converted into {len(entries)} get_MNV entries, "
        f"{len(unconverted)} not converted",
        file=sys.stderr,
    )
    for reason, count in reasons.most_common():
        print(f"  {count:6d}  {reason}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
