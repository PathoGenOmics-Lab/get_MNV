#!/usr/bin/env python3
"""Read-level phasing validated against a BAM whose phase we chose ourselves.

The synthetic truth suite validates the codon arithmetic but has no BAM, so it
cannot say anything about linkage. This one supplies the missing half: it builds
alignments molecule by molecule, deciding which reads carry which variants, and
then asks whether get_MNV recovers the mixture it was handed.

Every expectation here is arithmetic on the molecule counts, never a number read
off a previous run.

    cargo build -p get_mnv --release
    python3 tests/phasing/run.py

Point it at another binary with GET_MNV to confirm it can still fail:

    GET_MNV=/path/to/older/get_mnv python3 tests/phasing/run.py
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "tests" / "scenarios"))

from framework import (  # noqa: E402
    Op,
    ReadGroup,
    Scenario,
    VcfRecord,
    parse_tsv,
    run_scenario,
)
from scenarios import GFF_CDS_MULTIEXON  # noqa: E402

# Absent values in the TSV.
MISSING = "-"


@dataclass(frozen=True)
class Molecule:
    """One class of identical molecules in the sample.

    `carries[i]` says whether this molecule has the ALT at codon position `i`.
    `spans` says whether a read off this molecule reaches every position of the
    codon; a molecule sequenced only in part testifies about no pair at all.
    """

    carries: tuple[bool, ...]
    count: int
    spans: bool = True


@dataclass(frozen=True)
class Codon:
    """A codon under test, and how to lay a read across it."""

    label: str
    gene: str
    positions: tuple[int, ...]
    refs: tuple[str, ...]
    alts: tuple[str, ...]
    # Geometry of a read that reaches every position.
    span_start: int
    span_length: int
    # When set, a spanning molecule is sequenced as a mate pair instead of one
    # read: the first mate covers the leading positions and the second the rest,
    # so neither alone reaches across the codon. Value is the read length used
    # for each mate.
    paired_mate_length: int | None = None
    intron: tuple[int, int] | None = None  # (start, length) of a CIGAR N block
    gff_content: str | None = None
    gff_features: str | None = None


# geneA, plus strand, single exon. Codon 10 = pos 28-30 = GCT.
PLUS = Codon(
    label="plus_strand",
    gene="geneA",
    positions=(28, 30),
    refs=("G", "T"),
    alts=("A", "A"),
    span_start=11,
    span_length=40,
)

# geneB, minus strand, single exon. Codon 10 = pos 671-673, GCT in transcript
# orientation. Linkage is a property of the molecules, so the answer must not
# depend on which strand the gene is read from.
MINUS = Codon(
    label="minus_strand",
    gene="geneB",
    positions=(671, 673),
    refs=("A", "C"),
    alts=("T", "A"),
    span_start=655,
    span_length=40,
)

# geneC, spliced. Codon 34 straddles the intron: pos 900 sits at the end of
# exon 1 and pos 1002 near the start of exon 2, 102 bp away in the genome.
SPLICED = Codon(
    label="spliced_junction",
    gene="geneC",
    positions=(900, 1002),
    refs=("G", "T"),
    alts=("A", "C"),
    span_start=851,
    span_length=200,
    intron=(901, 100),
    gff_content=GFF_CDS_MULTIEXON,
    gff_features="CDS",
)

# All three positions of one codon, to exercise the least-supported rule with
# more than two constituents.
TRIPLE = Codon(
    label="triple_plus_strand",
    gene="geneA",
    positions=(28, 29, 30),
    refs=("G", "C", "T"),
    alts=("A", "G", "C"),
    span_start=11,
    span_length=40,
)


# The same plus-strand codon, sequenced as mate pairs. Each mate reaches one
# end of the codon and not the other, so the answer exists only at the level of
# the molecule. The linkage arithmetic must come out identical to PLUS: a
# fragment is one molecule however it was read.
PAIRED = Codon(
    label="paired_plus_strand",
    gene="geneA",
    positions=(28, 30),
    refs=("G", "T"),
    alts=("A", "A"),
    span_start=11,
    span_length=40,
    paired_mate_length=20,
)


def read_groups(codon: Codon, molecules: list[Molecule]) -> list[ReadGroup]:
    """Turn molecule classes into alignments.

    A spanning molecule becomes one read over the whole codon (skipping the
    intron where there is one). A non-spanning molecule becomes a read that
    stops at the first position it carries, so it never observes the rest.
    """
    groups: list[ReadGroup] = []
    for index, molecule in enumerate(molecules):
        if molecule.count == 0:
            continue
        edits = [
            Op(kind="snv", pos=position, seq=alt)
            for position, alt, carried in zip(codon.positions, codon.alts, molecule.carries)
            if carried
        ]
        if molecule.spans and codon.paired_mate_length is not None:
            # One molecule, two mates. The first ends before the last position
            # of the codon and the second starts after the first, so no single
            # read sees the pair and only the fragment can answer.
            first, last = codon.positions[0], codon.positions[-1]
            mate_length = codon.paired_mate_length
            groups.append(
                ReadGroup(
                    name_prefix=f"m{index}",
                    start=first - mate_length + 1,
                    length=mate_length,
                    ops=[op for op in edits if op.pos <= first],
                    count=molecule.count,
                    mate_start=last,
                    mate_length=mate_length,
                    mate_ops=[op for op in edits if op.pos >= last],
                )
            )
            continue

        if molecule.spans:
            ops = list(edits)
            if codon.intron is not None:
                intron_start, intron_length = codon.intron
                ops.append(Op(kind="skip", pos=intron_start, length=intron_length))
            groups.append(
                ReadGroup(
                    name_prefix=f"m{index}",
                    start=codon.span_start,
                    length=codon.span_length,
                    ops=ops,
                    count=molecule.count,
                )
            )
            continue

        # Not spanning: end the read on the first position this molecule
        # carries, so the partner positions are never observed.
        carried_positions = [
            position
            for position, carried in zip(codon.positions, molecule.carries)
            if carried
        ]
        if not carried_positions:
            continue
        stop = carried_positions[0]
        groups.append(
            ReadGroup(
                name_prefix=f"m{index}_short",
                start=codon.span_start,
                length=stop - codon.span_start + 1,
                ops=[op for op in edits if op.pos <= stop],
                count=molecule.count,
            )
        )
    return groups


def expected_phasing(molecules: list[Molecule], width: int) -> tuple[str, str]:
    """The answer implied by the mixture, computed from the counts alone.

    Only molecules read end to end may testify. Among those, the least-supported
    constituent SNV fixes the denominator; the molecules carrying the whole
    haplotype fix the numerator. With no witnesses at all the answer is unknown,
    which is not the same as no linkage.
    """
    informative = [m for m in molecules if m.spans]
    carrying = [
        sum(m.count for m in informative if m.carries[i]) for i in range(width)
    ]
    least_supported = min(carrying)
    if least_supported == 0:
        return MISSING, MISSING
    full_haplotype = sum(m.count for m in informative if all(m.carries))
    return f"{full_haplotype / least_supported:.4f}", str(least_supported)


def fisher_exact_two_tailed(a: int, b: int, c: int, d: int) -> float:
    """Two-tailed Fisher exact p for the 2x2 table, from first principles.

    Sum the hypergeometric probability of every table with the same margins that
    is no more likely than the observed one. Written out here rather than
    imported so it is a genuinely separate derivation from the Rust.
    """
    row1, row2, col1 = a + b, c + d, a + c
    total = row1 + row2
    if total == 0:
        return 1.0

    def pmf(k: int) -> float:
        if k < 0 or k > row1 or col1 - k < 0 or col1 - k > row2:
            return 0.0
        return math.comb(row1, k) * math.comb(row2, col1 - k) / math.comb(total, col1)

    observed = pmf(a)
    return min(1.0, sum(
        pmf(k) for k in range(max(0, col1 - row2), min(row1, col1) + 1)
        if pmf(k) <= observed + 1e-12
    ))


def expected_linkage(both: int, first: int, second: int, neither: int):
    """D-prime and its p-value for a 2x2 table, derived from the counts alone.

    D-prime is the excess co-occurrence over what the two frequencies predict,
    divided by the most that excess could have been. Undefined when either
    allele is on every molecule or on none, since nothing varies to correlate.
    """
    total = both + first + second + neither
    if total == 0:
        return MISSING, None
    n = float(total)
    p1 = (both + first) / n
    p2 = (both + second) / n
    if p1 in (0.0, 1.0) or p2 in (0.0, 1.0):
        return MISSING, None
    disequilibrium = both / n - p1 * p2
    if disequilibrium >= 0:
        maximum = min(p1 * (1 - p2), (1 - p1) * p2)
    else:
        maximum = min(p1 * p2, (1 - p1) * (1 - p2))
    if maximum <= 0:
        return MISSING, None
    d_prime = max(-1.0, min(1.0, disequilibrium / maximum))
    return f"{d_prime:.4f}", fisher_exact_two_tailed(both, first, second, neither)


def run_linkage_case(both: int, first: int, second: int, neither: int, work: Path) -> list[str]:
    """One 2x2 mixture of molecules, checked against the arithmetic above."""
    groups = []
    edits = {
        "both": [Op(kind="snv", pos=p, seq=a) for p, a in zip(PLUS.positions, PLUS.alts)],
        "first": [Op(kind="snv", pos=PLUS.positions[0], seq=PLUS.alts[0])],
        "second": [Op(kind="snv", pos=PLUS.positions[1], seq=PLUS.alts[1])],
        "neither": [],
    }
    for name, count in [("both", both), ("first", first), ("second", second),
                        ("neither", neither)]:
        if count:
            groups.append(ReadGroup(name_prefix=f"r_{name}", start=PLUS.span_start,
                                    length=PLUS.span_length, ops=edits[name], count=count))

    name = f"linkage_{both}_{first}_{second}_{neither}"
    scenario = Scenario(
        name=name,
        description="linkage from a 2x2 mixture chosen by construction",
        variants=[
            VcfRecord(pos=position, ref=ref, alt=alt)
            for position, ref, alt in zip(PLUS.positions, PLUS.refs, PLUS.alts)
        ],
        reads=groups,
        expected=[],
    )
    _, _, out = run_scenario(scenario, work)
    _, rows = parse_tsv(out)
    row = next((r for r in rows if r.get("Positions") == "28, 30"), None)
    if row is None:
        return [f"{name}: no row for the codon"]

    want_d, want_p = expected_linkage(both, first, second, neither)
    problems = []
    got_d = row.get("Codon LD", "<absent>")
    if got_d != want_d:
        problems.append(f"{name}: Codon LD expected {want_d}, got {got_d}")
    got_p = row.get("Codon LD p", "<absent>")
    if want_p is None:
        if got_p != MISSING:
            problems.append(f"{name}: Codon LD p expected {MISSING}, got {got_p}")
    else:
        try:
            parsed = float(got_p)
        except ValueError:
            problems.append(f"{name}: Codon LD p expected ~{want_p:.3e}, got {got_p}")
        else:
            if abs(parsed - want_p) > 1e-3 * max(want_p, 1e-12):
                problems.append(f"{name}: Codon LD p expected {want_p:.3e}, got {got_p}")
    return problems


def run_case(
    codon: Codon,
    molecules: list[Molecule],
    name: str,
    work: Path,
    extra_args: list[str] | None = None,
) -> list[str]:
    scenario = Scenario(
        name=name,
        description=f"{codon.label}: molecule-level phase chosen by construction",
        variants=[
            VcfRecord(pos=position, ref=ref, alt=alt)
            for position, ref, alt in zip(codon.positions, codon.refs, codon.alts)
        ],
        reads=read_groups(codon, molecules),
        expected=[],
        gff_content=codon.gff_content,
        gff_features=codon.gff_features,
        extra_cli_args=extra_args,
    )
    _, _, out = run_scenario(scenario, work)
    _, rows = parse_tsv(out)

    wanted = ", ".join(str(position) for position in codon.positions)
    row = next((r for r in rows if r.get("Positions") == wanted), None)
    if row is None:
        return [f"{name}: no row for positions {wanted}"]

    if extra_args and "--count-mates-separately" in extra_args:
        # Neither mate reaches across the codon, so read by read nobody can
        # answer. This is what the pair-aware sweep would degrade to, and it is
        # why those cases have teeth.
        want_support, want_reads = MISSING, MISSING
    else:
        want_support, want_reads = expected_phasing(molecules, len(codon.positions))
    problems = []
    got_support = row.get("MNV Phasing Support", "<absent>")
    got_reads = row.get("MNV Phasing Reads", "<absent>")
    if got_support != want_support:
        problems.append(
            f"{name}: MNV Phasing Support expected {want_support}, got {got_support}"
        )
    if got_reads != want_reads:
        problems.append(
            f"{name}: MNV Phasing Reads expected {want_reads}, got {got_reads}"
        )
    return problems


def linkage_sweep(codon: Codon) -> list[tuple[str, list[Molecule]]]:
    """Every mixture from fully trans to fully cis, at a fixed total of 20
    molecules carrying each constituent SNV."""
    width = len(codon.positions)
    both = tuple([True] * width)
    cases = []
    for full in range(0, 21):
        solo_count = 20 - full
        molecules = [Molecule(carries=both, count=full)]
        for i in range(width):
            carries = tuple(j == i for j in range(width))
            molecules.append(Molecule(carries=carries, count=solo_count))
        cases.append((f"{codon.label}_cis{full:02d}", molecules))
    return cases


def all_cases() -> list[tuple[Codon, str, list[Molecule]]]:
    cases: list[tuple[Codon, str, list[Molecule]]] = []

    for codon in (PLUS, MINUS, SPLICED, TRIPLE, PAIRED):
        for name, molecules in linkage_sweep(codon):
            cases.append((codon, name, molecules))

    # Reads that carry one variant and stop before the other are not evidence
    # against linkage: every molecule read end to end here carries both, so the
    # answer is 1.0 out of 20, undiluted by the 20 partial reads.
    cases.append(
        (
            PLUS,
            "plus_strand_partial_reads_do_not_dilute",
            [
                Molecule(carries=(True, True), count=20),
                Molecule(carries=(True, False), count=10, spans=False),
                Molecule(carries=(False, True), count=10, spans=False),
            ],
        )
    )

    # Nothing spans the codon: both variants are seen, never together. Unknown,
    # and reported as such rather than as zero linkage.
    for codon in (PLUS, SPLICED):
        cases.append(
            (
                codon,
                f"{codon.label}_nothing_spans_the_codon",
                [
                    Molecule(carries=(True,) + (False,) * (len(codon.positions) - 1),
                             count=20, spans=False),
                    Molecule(carries=(False,) * (len(codon.positions) - 1) + (True,),
                             count=20, spans=False),
                ],
            )
        )

    # A minority haplotype riding inside a majority variant: every molecule with
    # the rare SNV also has the common one. The documented rule keys on the
    # least-supported constituent, so this is full support at a low count.
    cases.append(
        (
            PLUS,
            "plus_strand_minor_haplotype_nested_in_major",
            [
                Molecule(carries=(True, True), count=3),
                Molecule(carries=(True, False), count=97),
            ],
        )
    )

    return cases


def main() -> int:
    work = HERE / "work"
    work.mkdir(exist_ok=True)

    cases = all_cases()
    problems: list[str] = []
    for codon, name, molecules in cases:
        problems.extend(run_case(codon, molecules, name, work))

    # The contrast that gives the paired sweep its teeth: read by read, no mate
    # reaches across the codon, so the same alignments answer nothing.
    contrast = [
        Molecule(carries=(True, True), count=10),
        Molecule(carries=(True, False), count=10),
        Molecule(carries=(False, True), count=10),
    ]
    problems.extend(
        run_case(
            PAIRED,
            contrast,
            "paired_plus_strand_without_pairing",
            work,
            extra_args=["--count-mates-separately"],
        )
    )

    # Linkage: 2x2 mixtures spanning independence, both directions of
    # association, and the degenerate tables where D-prime has no value.
    linkage_cases = [
        (81, 9, 9, 1),      # both common, exactly what chance predicts
        (90, 0, 0, 10),     # both common and genuinely together
        (3, 0, 0, 97),      # a rare pair, perfectly linked
        (0, 20, 20, 0),     # never together: competing haplotypes
        (10, 10, 10, 10),   # independent at intermediate frequency
        (18, 2, 2, 18),     # strong but imperfect association
        (2, 18, 18, 2),     # strong repulsion, imperfect
        (5, 0, 0, 1),       # linked but on few molecules: p should be weak
        (20, 0, 20, 0),     # first allele on every molecule: undefined
        (20, 20, 0, 0),     # second allele on every molecule: undefined
    ]
    for mixture in linkage_cases:
        problems.extend(run_linkage_case(*mixture, work))

    total_cases = len(cases) + 1 + len(linkage_cases)
    print(f"{total_cases} phased mixtures built, {total_cases} rows compared")
    print("  codons: " + ", ".join(
        f"{c.label}({c.gene}, {len(c.positions)} SNV(s))"
        for c in (PLUS, MINUS, SPLICED, TRIPLE, PAIRED)
    ))
    print()

    if not problems:
        print("OK: every linkage call matches the mixture the BAM was built with.")
        return 0

    print(f"{len(problems)} disagreement(s):")
    for problem in problems[:40]:
        print(f"  {problem}")
    if len(problems) > 40:
        print(f"  ... and {len(problems) - 40} more")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
