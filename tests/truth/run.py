#!/usr/bin/env python3
"""Validation against a synthetic dataset whose answer is known by construction.

The differential suite compares get_MNV with `bcftools csq` and SnpEff, which is
the stronger check because those tools are not ours. But it cannot validate the
one thing get_MNV exists for: **codon-level MNVs**. No per-variant annotator
combines two SNVs in a codon the way get_MNV does, so on exactly those calls the
oracles disagree by design and cannot say whether the residue is right.

This suite fills that hole. It builds a genome, places genes on both strands
(single-exon and spliced), introduces variants, and derives the expected answer
independently:

    mutate the genome -> re-extract the CDS -> translate -> compare proteins

That is a different route from get_MNV's, which does codon-level surgery in
place. The amino acids and codons here come from a codon table written out in
this file and from string operations on the reference, not from any get_MNV
logic, so an implementation slip shows up as a disagreement.

Being honest about its strength: this is a second implementation by the same
hand, which is weaker than a third-party oracle. It catches implementation
mistakes, not a shared misunderstanding of the biology. The consequence *labels*
(`Synonymous`, `Stop gained`, ...) follow get_MNV's documented convention, so on
those this validates the code against the documentation rather than against
nature. The amino acid itself is the independent part.

Usage:
    python3 tests/truth/run.py            # build, run, compare
    python3 tests/truth/run.py --keep     # leave the generated files in place
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent / "work"

# The standard genetic code, written out rather than imported, so this file does
# not share a table with the implementation under test. Internal codons are
# identical in NCBI tables 1 and 11; every gene here starts with ATG so the
# alternative bacterial initiation codons never come up.
BASES = "TCAG"
AMINO_ACIDS = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
CODON_TABLE = {
    a + b + c: AMINO_ACIDS[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in BASES for b in BASES for c in BASES
    )
}

THREE_LETTER = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly",
    "H": "His", "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn",
    "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val",
    "W": "Trp", "Y": "Tyr", "*": "*",
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def revcomp(sequence: str) -> str:
    return "".join(COMPLEMENT[base] for base in reversed(sequence))


def translate(sequence: str) -> str:
    return "".join(
        CODON_TABLE.get(sequence[i : i + 3], "X")
        for i in range(0, len(sequence) - len(sequence) % 3, 3)
    )


@dataclass
class Gene:
    name: str
    strand: str
    # CDS segments in ascending genomic order, 1-based inclusive.
    segments: list[tuple[int, int]]

    def cds(self, genome: str) -> str:
        """The coding sequence, spliced and in transcript orientation."""
        pieces = [genome[start - 1 : end] for start, end in self.segments]
        if self.strand == "+":
            return "".join(pieces)
        # Transcript order on the minus strand runs from the highest segment
        # down, and each piece is read on the other strand.
        return "".join(revcomp(piece) for piece in reversed(pieces))

    def cds_positions(self) -> list[int]:
        """Genomic position of every CDS base, in transcript order."""
        positions: list[int] = []
        if self.strand == "+":
            for start, end in self.segments:
                positions.extend(range(start, end + 1))
        else:
            for start, end in reversed(self.segments):
                positions.extend(range(end, start - 1, -1))
        return positions

    def cds_with_indel(
        self, genome: str, deleted: set[int], inserted_after: dict[str, str]
    ) -> str:
        """The CDS after removing `deleted` genomic bases and inserting after others.

        Walking the CDS position by position keeps the coordinates of the
        untouched bases fixed, which is what an indel actually does; shortening
        the genome and re-slicing by the old coordinates does not.

        A VCF insertion places its bases between genomic p and p+1. On the minus
        strand the transcript walks p+1 before p, so the inserted sequence is
        emitted before the base at p, reverse-complemented.
        """
        out: list[str] = []
        for position in self.cds_positions():
            insertion = inserted_after.get(str(position))
            if insertion and self.strand == "-":
                out.append(revcomp(insertion))
            if position not in deleted:
                base = genome[position - 1]
                out.append(base if self.strand == "+" else COMPLEMENT[base])
            if insertion and self.strand == "+":
                out.append(insertion)
        return "".join(out)

    def codon_of(self, position: int) -> int | None:
        """1-based codon number containing a genomic position, or None."""
        offset = self.transcript_offset(position)
        return None if offset is None else offset // 3 + 1

    def transcript_offset(self, position: int) -> int | None:
        """0-based offset of a genomic position in the spliced CDS."""
        if self.strand == "+":
            walked = 0
            for start, end in self.segments:
                if start <= position <= end:
                    return walked + position - start
                walked += end - start + 1
            return None
        walked = 0
        for start, end in reversed(self.segments):
            if start <= position <= end:
                return walked + end - position
            walked += end - start + 1
        return None


@dataclass
class IndelCase:
    """An indel with the answer derived by re-translating the mutated CDS."""

    position: int
    ref_allele: str
    alt_allele: str
    gene: str
    change_type: str
    event_class: str
    ref_codon: str | None
    alt_codon: str | None
    note: str = ""


@dataclass
class Case:
    """One variant with the answer we expect, derived independently below."""

    positions: list[int]
    ref_bases: list[str]
    alt_bases: list[str]
    gene: str
    aa_change: str
    change_type: str
    ref_codon: str
    alt_codon: str
    hgvs_c: str
    note: str = ""
    records: list[tuple[int, str, str]] = field(default_factory=list)


def build_genome(rng: random.Random) -> tuple[str, list[Gene]]:
    """A genome with four genes: both strands, single-exon and spliced."""

    def sense_codons(count: int) -> str:
        out = []
        while len(out) < count:
            codon = "".join(rng.choice("ACGT") for _ in range(3))
            if codon not in ("TAA", "TAG", "TGA"):
                out.append(codon)
        return "".join(out)

    def orf(codons: int) -> str:
        return "ATG" + sense_codons(codons - 2) + "TAA"

    filler = lambda n: "".join(rng.choice("ACGT") for _ in range(n))

    parts: list[str] = []
    genes: list[Gene] = []
    cursor = 1

    def place(name: str, strand: str, layout: list[int]) -> None:
        """`layout` alternates exon lengths and intron lengths."""
        nonlocal cursor
        total_coding = sum(layout[::2])
        coding = orf(total_coding // 3)
        if strand == "-":
            coding = revcomp(coding)
        segments: list[tuple[int, int]] = []
        taken = 0
        # On the minus strand the transcript starts at the highest coordinate, so
        # the coding sequence laid on the plus strand is already reversed above
        # and is written left to right here.
        for index, length in enumerate(layout):
            if index % 2 == 0:
                parts.append(coding[taken : taken + length])
                segments.append((cursor, cursor + length - 1))
                taken += length
                cursor += length
            else:
                parts.append(filler(length))
                cursor += length
        genes.append(Gene(name, strand, segments))

    parts.append(filler(30))
    cursor += 30
    place("geneP", "+", [90])
    parts.append(filler(40))
    cursor += 40
    place("geneM", "-", [90])
    parts.append(filler(40))
    cursor += 40
    place("geneS", "+", [45, 50, 45])
    parts.append(filler(40))
    cursor += 40
    place("geneMS", "-", [45, 50, 45])
    parts.append(filler(30))

    return "".join(parts), genes


def expected_for(genome: str, gene: Gene, changes: list[tuple[int, str]]) -> Case | None:
    """Derive the expected annotation by mutating and re-translating.

    `changes` are (genomic position, alternate base). Everything below works on
    strings and the codon table above; nothing here knows how get_MNV computes.
    """
    offsets = [gene.transcript_offset(position) for position, _ in changes]
    if any(offset is None for offset in offsets):
        return None
    codons = {offset // 3 for offset in offsets}
    if len(codons) != 1:
        return None  # only same-codon groups are compared as one row
    codon_index = codons.pop()

    mutated = list(genome)
    for position, alternate in changes:
        mutated[position - 1] = alternate
    mutated = "".join(mutated)

    reference_protein = translate(gene.cds(genome))
    alternate_protein = translate(gene.cds(mutated))
    if codon_index >= len(reference_protein):
        return None

    reference_aa = reference_protein[codon_index]
    alternate_aa = alternate_protein[codon_index]
    position_1based = codon_index + 1

    if reference_aa == alternate_aa:
        change_type = "Synonymous"
    elif alternate_aa == "*":
        change_type = "Stop gained"
    elif reference_aa == "*":
        change_type = "Stop lost"
    elif position_1based == 1 and reference_aa == "M":
        change_type = "Start lost"
    else:
        change_type = "Non-synonymous"

    ordered = sorted(changes, key=lambda item: item[0])

    # HGVS c. is derived from the CDS coordinate and the coding-strand bases,
    # both taken from the spliced CDS directly. This is where a minus-strand
    # coordinate slip would show up.
    reference_cds = gene.cds(genome)
    alternate_cds = gene.cds(mutated)
    entries = sorted(
        (gene.transcript_offset(position) for position, _ in ordered)
    )
    parts = [
        f"{offset + 1}{reference_cds[offset]}>{alternate_cds[offset]}"
        for offset in entries
    ]
    hgvs_c = f"c.{parts[0]}" if len(parts) == 1 else "c.[" + ";".join(parts) + "]"

    return Case(
        positions=[position for position, _ in ordered],
        ref_bases=[genome[position - 1] for position, _ in ordered],
        alt_bases=[alternate for _, alternate in ordered],
        gene=gene.name,
        aa_change=(
            f"{THREE_LETTER[reference_aa]}{position_1based}{THREE_LETTER[alternate_aa]}"
        ),
        change_type=change_type,
        ref_codon=reference_cds[codon_index * 3 : codon_index * 3 + 3],
        alt_codon=alternate_cds[codon_index * 3 : codon_index * 3 + 3],
        hgvs_c=hgvs_c,
        records=[
            (position, genome[position - 1], alternate) for position, alternate in ordered
        ],
    )


def expected_indel(
    genome: str, gene: Gene, position: int, ref_allele: str, alt_allele: str
) -> IndelCase | None:
    """Derive an indel's expected annotation by re-translating the mutated CDS.

    The proteins and the frameshift determination come from string operations on
    the reference; only the label mapping follows get_MNV's documented
    convention.
    """
    # A VCF indel is anchored: the shared first base stays, the rest of REF is
    # deleted and the rest of ALT is inserted after the anchor.
    deleted = set(range(position + 1, position + len(ref_allele)))
    inserted_after = {str(position): alt_allele[1:]} if len(alt_allele) > 1 else {}

    if inserted_after:
        # An insertion goes into the gap between genomic p and p+1, and that gap
        # is only inside the coding sequence when both bases are. At the CDS
        # terminus it is not: on the plus strand the highest coordinate is the
        # last transcript base and the gap above it is 3' UTR, while on the minus
        # strand that same coordinate is the *first* transcript base and the gap
        # is 5' UTR. Either way the inserted bases are not translated, and there
        # is no coding consequence for this suite to assert.
        coding = set(gene.cds_positions())
        # Any other exon end is an internal junction: in the spliced mRNA the
        # inserted bases land between two exons and are transcribed. Only the
        # gene's highest coordinate is outside, on either strand, because an
        # insertion anchored there always inserts above the gene.
        highest = gene.segments[-1][1]
        internal_junction = any(
            position == end and end != highest for _, end in gene.segments
        )
        if not ((position in coding and position + 1 in coding) or internal_junction):
            return None
    reference_cds = gene.cds(genome)
    alternate_cds = gene.cds_with_indel(genome, deleted, inserted_after)

    delta = len(alternate_cds) - len(reference_cds)
    if delta == 0:
        return None  # not length-changing inside this gene: nothing to assert
    frameshift = delta % 3 != 0

    reference_protein = translate(reference_cds)
    alternate_protein = translate(alternate_cds)

    # What the ribosome makes stops at the first stop codon; everything past it
    # is translated by nobody, so two alternates that differ only there produce
    # the same protein.
    def made(protein: str) -> str:
        stop = protein.find("*")
        return protein if stop < 0 else protein[: stop + 1]

    if frameshift:
        change_type = "Frameshift Indel"
    elif reference_protein[:1] == "M" and alternate_protein[:1] not in ("M", "*"):
        change_type = "Start lost"
    elif alternate_protein.count("*") > reference_protein.count("*"):
        change_type = "Stop gained"
    elif alternate_protein.count("*") < reference_protein.count("*"):
        change_type = "Stop lost"
    elif made(reference_protein) == made(alternate_protein):
        # An in-frame indel that leaves the made protein byte for byte the same
        # changed nothing. An insertion inside the terminal stop codon that keeps
        # a stop there is the case, and it used to be reported as a residue
        # inserted one position past the protein's last one.
        change_type = "Synonymous"
    else:
        change_type = "In-frame Indel"

    event_class = "deletion" if len(alt_allele) < len(ref_allele) else "insertion"

    # The codon shown is the one containing the first CDS base that actually
    # changes. In a VCF indel the anchor base is padding and is unchanged, so
    # for a deletion that is the first deleted base, and for an insertion it is
    # the anchor the new bases hang off.
    changed_positions = sorted(deleted) if deleted else [position]
    touched = [gene.transcript_offset(pos) for pos in changed_positions]
    touched = sorted(offset for offset in touched if offset is not None)
    if not touched:
        return None
    codon_start = (touched[0] // 3) * 3
    ref_codon = reference_cds[codon_start : codon_start + 3]
    alt_codon = alternate_cds[codon_start : codon_start + 3]
    if len(ref_codon) < 3:
        return None

    return IndelCase(
        position=position,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
        gene=gene.name,
        change_type=change_type,
        event_class=event_class,
        ref_codon=ref_codon,
        alt_codon=alt_codon if len(alt_codon) == 3 else None,
    )


def build_indel_cases(genome: str, genes: list[Gene]) -> list[IndelCase]:
    """Deletions and insertions of 1, 2 and 3 bases across each gene."""
    cases: list[IndelCase] = []
    for gene in genes:
        for segment_start, segment_end in gene.segments:
            # Anchors spread through the exon, kept clear of its edges so the
            # deleted span stays inside the CDS.
            # Every position of every exon. Where an insertion goes is the
            # subtlest decision in the whole indel path: the bases must be
            # reverse-complemented on the minus strand, they must land between
            # the right two residues, and the gap above the gene's highest
            # coordinate is outside the CDS on *either* strand, because there it
            # is the 3' end of a plus-strand transcript and the 5' end of a minus
            # one. Sampling anchors would leave those edges untested.
            for anchor in range(segment_start, segment_end + 1):
                for length in (1, 2, 3):
                    reference = genome[anchor - 1 : anchor + length]
                    case = expected_indel(
                        genome, gene, anchor, reference, genome[anchor - 1]
                    )
                    if case:
                        case.note = f"del{length} {gene.name} at {anchor}"
                        cases.append(case)
                    inserted = "ACG"[:length]
                    case = expected_indel(
                        genome, gene, anchor, genome[anchor - 1], genome[anchor - 1] + inserted
                    )
                    if case:
                        case.note = f"ins{length} {gene.name} at {anchor}"
                        cases.append(case)
    return cases


def build_cases(genome: str, genes: list[Gene]) -> list[Case]:
    """Every single-base change in a stretch of codons, plus same-codon groups."""
    cases: list[Case] = []
    for gene in genes:
        coding_positions = [
            position
            for start, end in gene.segments
            for position in range(start, end + 1)
        ]
        by_codon: dict[int, list[int]] = {}
        for position in coding_positions:
            by_codon.setdefault(gene.codon_of(position), []).append(position)

        for codon_number, positions in sorted(by_codon.items()):
            if len(positions) != 3:
                continue
            positions.sort()
            # Every single substitution in this codon.
            for position in positions:
                for alternate in "ACGT":
                    if alternate == genome[position - 1]:
                        continue
                    case = expected_for(genome, gene, [(position, alternate)])
                    if case:
                        case.note = f"SNV {gene.name} codon {codon_number}"
                        cases.append(case)
            # Every pair, and the triple: the codon-level MNVs no other tool
            # combines the way get_MNV does.
            pairs = [(0, 1), (0, 2), (1, 2)]
            for left, right in pairs:
                a, b = positions[left], positions[right]
                alt_a = "ACGT".replace(genome[a - 1], "")[0]
                alt_b = "ACGT".replace(genome[b - 1], "")[1]
                case = expected_for(genome, gene, [(a, alt_a), (b, alt_b)])
                if case:
                    case.note = f"MNV {gene.name} codon {codon_number}"
                    cases.append(case)
            triple = [
                (position, "ACGT".replace(genome[position - 1], "")[index % 3])
                for index, position in enumerate(positions)
            ]
            case = expected_for(genome, gene, triple)
            if case:
                case.note = f"MNV3 {gene.name} codon {codon_number}"
                cases.append(case)
    return cases


def write_inputs(work: Path, genome: str, genes: list[Gene], cases: list[Case]) -> None:
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / "truth.fa"
    fasta.write_text(
        ">chrT\n" + "\n".join(genome[i : i + 60] for i in range(0, len(genome), 60)) + "\n"
    )

    rows = ["##gff-version 3"]
    for gene in genes:
        for start, end in gene.segments:
            rows.append(
                f"chrT\ttruth\tCDS\t{start}\t{end}\t.\t{gene.strand}\t0\t"
                f"ID=cds:{gene.name};Parent={gene.name};"
                f"transcript_id={gene.name};gene={gene.name}"
            )
    (work / "truth.gff").write_text("\n".join(rows) + "\n")


def write_vcf(path: Path, genome: str, records: list[tuple[int, str, str]]) -> None:
    header = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=chrT,length={len(genome)}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    body = [
        f"chrT\t{position}\t.\t{reference}\t{alternate}\t.\tPASS\t."
        for position, reference, alternate in sorted(records)
    ]
    path.write_text("\n".join(header + body) + "\n")


def binary() -> Path:
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Synthetic ground-truth validation")
    parser.add_argument("--keep", action="store_true", help="leave the generated files behind")
    args = parser.parse_args()

    rng = random.Random(20260804)
    genome, genes = build_genome(rng)
    cases = build_cases(genome, genes)
    indel_cases = build_indel_cases(genome, genes)

    if WORK.exists():
        shutil.rmtree(WORK)
    write_inputs(WORK, genome, genes, cases)

    # One VCF per case keeps each expectation independent: variants in different
    # codons would otherwise be grouped or shift each other's frame.
    failures: list[str] = []
    checked = 0
    checked_indels = 0
    for index, case in enumerate(cases):
        case_dir = WORK / f"case_{index:04d}"
        case_dir.mkdir(parents=True, exist_ok=True)
        vcf = case_dir / "in.vcf"
        write_vcf(vcf, genome, case.records)
        proc = subprocess.run(
            [
                str(binary()),
                "--vcf", "in.vcf",
                "--fasta", str((WORK / "truth.fa").resolve()),
                "--gff", str((WORK / "truth.gff").resolve()),
            ],
            cwd=case_dir,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            failures.append(f"{case.note}: get_mnv failed\n{proc.stderr.strip()[:400]}")
            continue

        with (case_dir / "in.MNV.tsv").open() as handle:
            rows = [row for row in csv.DictReader(handle, delimiter="\t")]
        wanted = ", ".join(str(position) for position in case.positions)
        row = next((r for r in rows if r["Positions"] == wanted), None)
        if row is None:
            failures.append(
                f"{case.note}: no row for positions {wanted}; got "
                + ", ".join(repr(r["Positions"]) for r in rows)
            )
            continue

        checked += 1
        for label, expected, actual in (
            ("gene", case.gene, row["Gene"]),
            ("AA change", case.aa_change, row["AA Changes"]),
            ("change type", case.change_type, row["Change Type"]),
            ("reference codon", case.ref_codon, row["Reference Codon"]),
            ("alternate codon", case.alt_codon, row["MNV Codon"]),
            ("HGVS c.", case.hgvs_c, row["HGVS c."]),
        ):
            if expected != actual:
                failures.append(
                    f"{case.note} at {wanted}: {label} expected {expected!r}, got {actual!r}"
                )

    # Indels, checked on the labels and codons that can be derived without
    # reimplementing HGVS protein notation.
    for index, indel in enumerate(indel_cases):
        case_dir = WORK / f"indel_{index:04d}"
        case_dir.mkdir(parents=True, exist_ok=True)
        vcf = case_dir / "in.vcf"
        write_vcf(vcf, genome, [(indel.position, indel.ref_allele, indel.alt_allele)])
        proc = subprocess.run(
            [
                str(binary()),
                "--vcf", "in.vcf",
                "--fasta", str((WORK / "truth.fa").resolve()),
                "--gff", str((WORK / "truth.gff").resolve()),
            ],
            cwd=case_dir,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            failures.append(f"{indel.note}: get_mnv failed\n{proc.stderr.strip()[:400]}")
            continue
        with (case_dir / "in.MNV.tsv").open() as handle:
            rows = [row for row in csv.DictReader(handle, delimiter="\t")]
        row = next((r for r in rows if r["Gene"] == indel.gene), None)
        if row is None:
            failures.append(f"{indel.note}: no row for gene {indel.gene}")
            continue
        checked_indels += 1
        checks = [
            ("change type", indel.change_type, row["Change Type"]),
            ("event class", indel.event_class, row["Event Class"]),
            ("reference codon", indel.ref_codon, row["Reference Codon"]),
        ]
        if indel.alt_codon is not None:
            checks.append(("alternate codon", indel.alt_codon, row["MNV Codon"]))
        for label, expected, actual in checks:
            if expected != actual:
                failures.append(
                    f"{indel.note}: {label} expected {expected!r}, got {actual!r}"
                )

    print(f"{len(cases)} substitution cases built, {checked} rows compared")
    print(f"{len(indel_cases)} indel cases built, {checked_indels} rows compared")
    print(f"  genes: {', '.join(f'{g.name}({g.strand}, {len(g.segments)} exon(s))' for g in genes)}")

    if not args.keep:
        shutil.rmtree(WORK, ignore_errors=True)

    if failures:
        print(f"\nFAIL: {len(failures)} disagreement(s) with the constructed truth.\n")
        # Grouped, and sampled per group. Printing the first N instead would let
        # one noisy category bury another entirely: a run that broke both the
        # substitution ordering and the start-codon label showed 40 of the
        # former and none of the latter.
        buckets: dict[str, list[str]] = {}
        for failure in failures:
            kind = failure.split()[0].rstrip(":")
            kind = "".join(ch for ch in kind if not ch.isdigit()) or kind
            buckets.setdefault(kind, []).append(failure)
        for kind, group in sorted(buckets.items(), key=lambda item: -len(item[1])):
            print(f"  {len(group):5d}  {kind}")
        print()
        for kind, group in sorted(buckets.items(), key=lambda item: -len(item[1])):
            for failure in group[:4]:
                print(f"  {failure}")
            if len(group) > 4:
                print(f"  ... and {len(group) - 4} more {kind}")
        return 1

    print("\nOK: every call matches the answer the dataset was built with.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
