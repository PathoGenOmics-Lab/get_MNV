#!/usr/bin/env python3
"""Read counting checked against bcftools mpileup.

Everything that validates get_MNV's read counting has so far been written by
the same hand as the code: the scenario expectations, and the phasing suite's
arithmetic. Both catch implementation error and neither can catch a shared
misunderstanding of what a read is worth. This suite hands the same BAM to
`bcftools mpileup`, which nobody here wrote, and compares four numbers per
site: depth, allele support, and support split by strand.

Those four are exactly where the layer's bugs have been living. The
strand-attribution defect fixed on this branch, where a mate hundreds of bases
away lent its strand to a variant it never read, is a straight disagreement with
`ADF`/`ADR` and would have been caught here the day it was written.

    cargo build -p get_mnv --release
    python3 tests/pileup/run.py

Point it at another binary with GET_MNV to confirm it can still fail.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "tests" / "scenarios"))

from framework import (  # noqa: E402
    CONTIG,
    GET_MNV,
    Op,
    ReadGroup,
    VcfRecord,
    reference_for,
    write_bam,
    write_fasta,
    write_gff,
    write_vcf,
)

BCFTOOLS = shutil.which("bcftools") or "bcftools"
SAMTOOLS = shutil.which("samtools") or "samtools"

# Thresholds handed to both tools, so neither is filtering on a rule the other
# does not know about.
MIN_BASE_QUALITY = 20
MIN_MAPPING_QUALITY = 30


@dataclass
class Site:
    """A position to compare, and the alternate allele to compare it on."""

    position: int
    alt: str

    @property
    def ref(self) -> str:
        return reference_for(CONTIG)[self.position - 1]


@dataclass
class Case:
    name: str
    description: str
    sites: list[Site]
    reads: list[ReadGroup]
    extra_args: list[str] = field(default_factory=list)
    # Differences we expect and accept, keyed "position:field".
    known_differences: dict[str, str] = field(default_factory=dict)


def run_mpileup(work: Path, bam: Path, fasta: Path, sites: list[Site]) -> dict:
    """Per-position depth and allele counts, straight from bcftools."""
    regions = work / "regions.txt"
    regions.write_text("".join(f"{CONTIG}\t{site.position}\n" for site in sites))
    out = work / "pileup.vcf"
    subprocess.run(
        [
            BCFTOOLS, "mpileup",
            # No BAQ: get_MNV does not recalibrate base qualities, so leaving it
            # on would have bcftools filtering bases get_MNV keeps.
            "-B",
            "-a", "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP",
            "-Q", str(MIN_BASE_QUALITY),
            "-q", str(MIN_MAPPING_QUALITY),
            "-d", "100000",
            "-R", str(regions),
            "-f", str(fasta),
            "-o", str(out),
            str(bam),
        ],
        check=True, capture_output=True,
    )
    query = subprocess.run(
        [BCFTOOLS, "query", "-f", "%POS\t%REF\t%ALT\t[%DP\t%AD\t%ADF\t%ADR]\n", str(out)],
        check=True, capture_output=True, text=True,
    ).stdout

    parsed = {}
    for line in query.splitlines():
        pos, ref, alts, depth, ad, adf, adr = line.split("\t")
        # Allele order is REF then the ALTs bcftools chose to report.
        alleles = [ref] + [a for a in alts.split(",") if a not in (".", "<*>")]
        counts = {}
        for values, key in ((ad, "ad"), (adf, "adf"), (adr, "adr")):
            numbers = [int(v) for v in values.split(",")]
            counts[key] = {
                allele: numbers[i] for i, allele in enumerate(alleles) if i < len(numbers)
            }
        parsed[int(pos)] = {"depth": int(depth), **counts}
    return parsed


def run_get_mnv(work: Path, bam: Path, fasta: Path, gff: Path, vcf: Path,
                extra: list[str]) -> dict:
    cmd = [
        GET_MNV, "--vcf", str(vcf), "--fasta", str(fasta), "--gff", str(gff),
        "--bam", str(bam),
        "--quality", str(MIN_BASE_QUALITY),
        "--min-mapq", str(MIN_MAPPING_QUALITY),
        *extra,
    ]
    proc = subprocess.run(cmd, cwd=str(work), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"get_mnv failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")

    rows = {}
    with (work / f"{vcf.stem}.MNV.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if "," in row["Positions"]:
                continue  # a grouped codon row is not a per-position count
            rows[int(row["Positions"])] = row
    return rows


def compare(case: Case, pileup: dict, rows: dict) -> list[str]:
    problems = []
    for site in case.sites:
        row = rows.get(site.position)
        if row is None:
            problems.append(f"{case.name}: get_MNV reported no row at {site.position}")
            continue
        seen = pileup.get(site.position)
        if seen is None:
            problems.append(f"{case.name}: mpileup reported nothing at {site.position}")
            continue

        checks = [
            ("depth", seen["depth"], int(row["Total Reads"])),
            ("support", seen["ad"].get(site.alt, 0), int(row["SNP Reads"])),
            ("forward", seen["adf"].get(site.alt, 0), int(row["SNP Forward Reads"])),
            ("reverse", seen["adr"].get(site.alt, 0), int(row["SNP Reverse Reads"])),
        ]
        for field_name, expected, got in checks:
            key = f"{site.position}:{field_name}"
            if expected == got:
                continue
            if key in case.known_differences:
                continue
            problems.append(
                f"{case.name}: {field_name} at {site.position} "
                f"is {got} for get_MNV and {expected} for mpileup"
            )
    return problems


def run_case(case: Case, base: Path) -> list[str]:
    work = base / case.name
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    fasta, gff, bam, vcf = (work / n for n in
                            ("ref.fasta", "genes.gff3", "reads.bam", "variants.vcf"))
    write_fasta(fasta)
    write_gff(gff, None)
    write_bam(bam, work, case.reads)
    write_vcf(vcf, [VcfRecord(pos=s.position, ref=s.ref, alt=s.alt) for s in case.sites])
    subprocess.run([SAMTOOLS, "faidx", str(fasta)], check=True, capture_output=True)

    return compare(case, run_mpileup(work, bam, fasta, case.sites),
                   run_get_mnv(work, bam, fasta, gff, vcf, case.extra_args))


# The comparison positions sit in the filler between genes, so every row is a
# single-position count and no codon grouping stands between the two tools.
FILLER_A, FILLER_B, FILLER_C = 320, 350, 750


def all_cases() -> list[Case]:
    def snv(pos: int, alt: str) -> Op:
        return Op(kind="snv", pos=pos, seq=alt)

    return [
        Case(
            name="01_mixed_strands",
            description="alt and reference reads on both strands",
            sites=[Site(FILLER_A, "G")],
            reads=[
                ReadGroup(name_prefix="alt_f", start=300, length=60,
                          ops=[snv(FILLER_A, "G")], count=12, strand="+"),
                ReadGroup(name_prefix="alt_r", start=300, length=60,
                          ops=[snv(FILLER_A, "G")], count=8, strand="-"),
                ReadGroup(name_prefix="ref_f", start=300, length=60,
                          ops=[], count=5, strand="+"),
                ReadGroup(name_prefix="ref_r", start=300, length=60,
                          ops=[], count=7, strand="-"),
            ],
        ),
        Case(
            name="02_single_strand_only",
            description="the alt appears on reverse reads only",
            sites=[Site(FILLER_A, "T")],
            reads=[
                ReadGroup(name_prefix="alt_r", start=300, length=60,
                          ops=[snv(FILLER_A, "T")], count=20, strand="-"),
                ReadGroup(name_prefix="ref_f", start=300, length=60,
                          ops=[], count=20, strand="+"),
            ],
        ),
        Case(
            name="03_mates_far_apart",
            description="each mate reads a different site; neither lends the other its strand",
            sites=[Site(FILLER_A, "G"), Site(FILLER_C, "G")],
            reads=[
                ReadGroup(
                    name_prefix="pair", start=300, length=60,
                    ops=[snv(FILLER_A, "G")], count=15,
                    mate_start=720, mate_length=60, mate_ops=[snv(FILLER_C, "G")],
                ),
            ],
        ),
        Case(
            name="04_overlapping_mates",
            description="both mates cover the site: one molecule, not two observations",
            sites=[Site(FILLER_A, "G")],
            reads=[
                ReadGroup(
                    name_prefix="pair", start=300, length=60,
                    ops=[snv(FILLER_A, "G")], count=16,
                    mate_start=300, mate_length=60, mate_ops=[snv(FILLER_A, "G")],
                ),
            ],
            # Depth and support agree: both tools see 16 molecules. The strand
            # split does not, and the two conventions are both defensible.
            # bcftools removes the overlap by keeping one mate's base and
            # zeroing the other, so each fragment lands in exactly one arm and
            # ADF + ADR == AD. Which arm is effectively arbitrary here: all 16
            # fragments are identical and it splits them 5 forward, 11 reverse,
            # a division carrying no information about strand at all. get_MNV
            # credits both arms, saying what is actually true, that all 16
            # molecules were read from both directions and none of them is
            # one-strand evidence. Keeping the difference visible and explained
            # rather than matching for the sake of matching.
            known_differences={
                f"{FILLER_A}:forward": "bcftools keeps one mate of an overlap, get_MNV keeps both",
                f"{FILLER_A}:reverse": "bcftools keeps one mate of an overlap, get_MNV keeps both",
            },
        ),
        Case(
            name="05_low_base_quality",
            description="bases below the floor are dropped by both tools",
            sites=[Site(FILLER_A, "G")],
            reads=[
                ReadGroup(name_prefix="good", start=300, length=60,
                          ops=[snv(FILLER_A, "G")], count=10),
                # '#' is Q2, far below the floor both tools are given.
                ReadGroup(name_prefix="poor", start=300, length=60,
                          ops=[snv(FILLER_A, "G")], count=10, base_quality="#"),
            ],
        ),
        Case(
            name="06_low_mapping_quality",
            description="reads below the mapping-quality floor are dropped by both tools",
            sites=[Site(FILLER_A, "G")],
            reads=[
                ReadGroup(name_prefix="good", start=300, length=60,
                          ops=[snv(FILLER_A, "G")], count=10, mapq=60),
                ReadGroup(name_prefix="poor", start=300, length=60,
                          ops=[snv(FILLER_A, "G")], count=10, mapq=5),
            ],
        ),
        Case(
            name="07_several_sites_one_gene_region",
            description="three sites counted from one cache, each on its own reads",
            sites=[Site(FILLER_A, "G"), Site(FILLER_B, "C"), Site(FILLER_C, "T")],
            reads=[
                ReadGroup(name_prefix="a", start=300, length=40,
                          ops=[snv(FILLER_A, "G")], count=9, strand="+"),
                ReadGroup(name_prefix="b", start=330, length=40,
                          ops=[snv(FILLER_B, "C")], count=11, strand="-"),
                ReadGroup(name_prefix="c", start=730, length=40,
                          ops=[snv(FILLER_C, "T")], count=13, strand="+"),
                ReadGroup(name_prefix="cover", start=300, length=60, ops=[], count=4),
            ],
        ),
        Case(
            name="08_mates_counted_separately",
            description="with --count-mates-separately the two mates are two observations",
            sites=[Site(FILLER_A, "G")],
            reads=[
                ReadGroup(
                    name_prefix="pair", start=300, length=60,
                    ops=[snv(FILLER_A, "G")], count=16,
                    mate_start=300, mate_length=60, mate_ops=[snv(FILLER_A, "G")],
                ),
            ],
            extra_args=["--count-mates-separately"],
            # bcftools always removes the overlap of a mate pair, so it has no
            # equivalent of this mode: it reports 16 molecules where get_MNV was
            # asked for 32 observations. The disagreement is the flag doing its
            # job, and it is here to pin that the default mode is the one that
            # agrees.
            known_differences={
                f"{FILLER_A}:depth": "bcftools has no per-record mode",
                f"{FILLER_A}:support": "bcftools has no per-record mode",
                f"{FILLER_A}:forward": "bcftools has no per-record mode",
                f"{FILLER_A}:reverse": "bcftools has no per-record mode",
            },
        ),
    ]


def main() -> int:
    for tool, name in ((BCFTOOLS, "bcftools"), (SAMTOOLS, "samtools")):
        if shutil.which(tool) is None:
            print(f"{name} not found on PATH; skipping the pileup comparison.")
            return 0

    work = HERE / "work"
    work.mkdir(exist_ok=True)

    cases = all_cases()
    problems: list[str] = []
    accepted = 0
    for case in cases:
        problems.extend(run_case(case, work))
        accepted += len(case.known_differences)

    checks = sum(len(case.sites) for case in cases) * 4
    print(f"{len(cases)} pile-ups compared, {checks} counts checked against bcftools mpileup")
    print(f"  {accepted} difference(s) accepted as documented")
    print()

    if not problems:
        print("OK: get_MNV and bcftools mpileup agree on every count.")
        return 0

    print(f"{len(problems)} disagreement(s):")
    for problem in problems[:40]:
        print(f"  {problem}")
    if len(problems) > 40:
        print(f"  ... and {len(problems) - 40} more")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
