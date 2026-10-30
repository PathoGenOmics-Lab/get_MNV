#!/usr/bin/env python3
"""Indel left-alignment checked against bcftools norm.

`--left-align` moves each insertion and deletion to the leftmost place in its
repeat that describes the same change, which is what `bcftools norm -f ref.fa`
does, and trims the bases both alleles share. The unit tests pin a handful of
records to what bcftools writes; this suite asks bcftools itself, on thousands.

It writes indels into the homopolymers and tandem repeats of the bundled
M. tuberculosis genome, at random places in each repeat (its right end, where
HGVS puts them, more often than chance), with random lengths (whole repeat
units and not), and in random spellings: anchored on the base before, on the
base after, or padded with context on both sides. Some are complex records, a
deletion and an insertion in one, which bcftools trims but never moves. Both
tools read the same VCF and the same FASTA, and every record must come out of
get_MNV with the POS, REF and ALT bcftools gives it, from the plain `.vcf`
reader and from the `.vcf.gz` one.

    cargo build -p get_mnv --release
    python3 tests/leftalign/run.py

get_MNV is also run without `--left-align`, and the records it then leaves
where bcftools would not are counted: a suite that cannot fail proves nothing,
and this number says it can. Point it at another binary with GET_MNV to see it
fail outright. Skips, with a message, when bcftools is not on PATH.
"""

from __future__ import annotations

import random
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "tests" / "scenarios"))

from framework import GET_MNV  # noqa: E402

BCFTOOLS = shutil.which("bcftools") or "bcftools"
FASTA = REPO_ROOT / "example" / "MTB_ancestor.fas"
GENES = REPO_ROOT / "example" / "anot_genes.txt"

SEED = 20261030
INDELS = 3000
COMPLEX = 400
# Whole-unit indels in the minisatellites (units of 7 to 100 bases), the only
# repeats in this genome long enough to move a record hundreds of bases.
LONG_UNITS = 150
# Records sit at least this far apart, so no two can reach each other however
# far either moves.
SPACING = 60
TAG = "LAID"
# `--sample all` reads a cohort through a path of its own, the file held in
# memory and each sample's column picked out of it, so the records go through
# that path too, as a cohort of two samples that both carry every one.
COHORT = ("S1", "S2")


@dataclass
class Record:
    ident: str
    kind: str  # deletion, insertion or complex
    period: int
    spelling: str  # left, right or padded
    pos: int  # 1-based
    ref: str
    alt: str


def read_genome(path: Path) -> tuple[str, str]:
    name, parts = "", []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                raise SystemExit(f"{path} holds more than one contig; this suite expects one")
            name = line[1:].split()[0]
        else:
            parts.append(line.strip().upper())
    return name, "".join(parts)


def find_repeats(seq: str) -> list[tuple[int, int, int]]:
    """Homopolymers of 4 or more and tandem repeats of 2 to 6 bases, as
    0-based half-open spans extended over any partial unit at either end."""
    found = []
    for period, copies in ((1, 4), (2, 3), (3, 3), (4, 2), (5, 2), (6, 2)):
        pattern = re.compile(r"([ACGT]{%d})\1{%d,}" % (period, copies - 1))
        for match in pattern.finditer(seq):
            start, end = match.start(), match.end()
            while start > 0 and seq[start - 1] == seq[start - 1 + period]:
                start -= 1
            while end < len(seq) and seq[end] == seq[end - period]:
                end += 1
            found.append((start, end, period))
    return found


def find_long_repeats(seq: str) -> list[tuple[int, int, int]]:
    """Stretches holding two or more copies of a unit of 7 to 100 bases.

    `seq[i] == seq[i + k]` along a whole stretch is what a repeat of period k
    is, and XOR-ing the sequence with itself shifted by k turns every such
    stretch into a run of zero bytes, which a regex finds at C speed. A loop
    over the bases in Python takes a minute for the same answer."""
    data = seq.encode()
    found = []
    for period in range(7, 101):
        width = len(data) - period
        shifted = int.from_bytes(data[:width], "big") ^ int.from_bytes(data[period:], "big")
        for match in re.finditer(rb"\x00{%d,}" % period, shifted.to_bytes(width, "big")):
            found.append((match.start(), match.end() + period, period))
    return found


def spaced(repeats: list[tuple[int, int, int]], rng: random.Random, wanted: int,
           length: int, taken: list[tuple[int, int]]) -> list[tuple[int, int, int]]:
    """A random choice of repeats, none within SPACING of another, of one
    already `taken`, or of either end of the contig."""
    rng.shuffle(repeats)
    chosen = []
    for start, end, period in repeats:
        if start < SPACING or end > length - SPACING:
            continue
        if any(start < e + SPACING and end > s - SPACING for s, e in taken):
            continue
        taken.append((start, end))
        chosen.append((start, end, period))
        if len(chosen) == wanted:
            break
    return chosen


def spell(seq: str, rng: random.Random, first: int, last: int, gap: int, removed: int,
          inserted: str) -> tuple[str, int, str, str]:
    """The change `seq[gap:gap+removed]` -> `inserted` as a VCF record, spelled
    one of three ways, as (spelling, pos, ref, alt)."""
    spelling = rng.choices(["left", "right", "padded"], weights=[5, 2, 3])[0]
    if spelling == "left":
        lo, hi = gap - 1, gap + removed
    elif spelling == "right":
        lo, hi = gap, gap + removed + 1
    else:
        lo, hi = gap - rng.randint(1, 3), gap + removed + rng.randint(0, 3)
    lo, hi = max(lo, first), min(hi, last)
    ref = seq[lo:hi]
    alt = seq[lo:gap] + inserted + seq[gap + removed:hi]
    return spelling, lo + 1, ref, alt


def make_records(seq: str, rng: random.Random) -> list[Record]:
    taken: list[tuple[int, int]] = []
    long_units = spaced(find_long_repeats(seq), rng, LONG_UNITS, len(seq), taken)
    short = spaced(find_repeats(seq), rng, INDELS + COMPLEX, len(seq), taken)
    if len(long_units) < LONG_UNITS or len(short) < INDELS + COMPLEX:
        raise SystemExit(f"only {len(long_units)} + {len(short)} spaced repeats found")
    records = []
    for number, (start, end, period) in enumerate(short + long_units):
        unit = seq[start:start + period]
        # Where in the repeat: anywhere, or its right end more often than
        # chance, which is where a 3'-shifted (HGVS) description puts it.
        if number < INDELS or number >= INDELS + COMPLEX:
            kind = rng.choice(["deletion", "insertion"])
            if number >= INDELS + COMPLEX:
                length = period
            elif rng.random() < 0.7:
                length = period * rng.randint(1, 3)
            else:
                length = rng.randint(1, 2 * period + 1)
            if kind == "deletion":
                length = min(length, end - start)
                rightmost = end - length
                gap = rightmost if rng.random() < 0.3 else rng.randint(start, rightmost)
                removed, inserted = length, ""
            else:
                gap = end if rng.random() < 0.3 else rng.randint(start, end)
                if rng.random() < 0.7:
                    # The repeat grows: the bases it would carry on with here.
                    phase = (gap - start) % period
                    inserted = (unit * (length // period + 2))[phase:phase + length]
                else:
                    inserted = "".join(rng.choice("ACGT") for _ in range(length))
                removed = 0
        else:
            kind = "complex"
            gap = rng.randint(start, end - 1)
            removed = rng.randint(1, min(6, end - gap))
            inserted = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))
            if len(inserted) == removed:
                inserted += rng.choice("ACGT")
        spelling, pos, ref, alt = spell(seq, rng, start - SPACING // 2, end + SPACING // 2,
                                        gap, removed, inserted)
        if ref == alt or len(ref) == len(alt) or not ref or not alt:
            continue
        records.append(Record(f"{kind[0]}{number}", kind, period, spelling, pos, ref, alt))
    # The contig's two ends, where the anchor rule changes: TT at 1-2.
    records.append(Record("edge_start_del", "deletion", 1, "right", 2, seq[1:3], seq[2]))
    records.append(Record("edge_start_ins", "insertion", 1, "left", 2, seq[1], seq[1] * 2))
    n = len(seq)
    records.append(Record("edge_end_ins", "insertion", 1, "left", n, seq[-1], seq[-1] + "A"))
    records.append(Record("edge_end_del", "deletion", 1, "left", n - 1, seq[-2:], seq[-2]))
    return sorted(records, key=lambda r: r.pos)


def write_vcf(path: Path, contig: str, length: int, records: list[Record],
              samples: tuple[str, ...] = ()) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={contig},length={length}>",
        f'##INFO=<ID={TAG},Number=1,Type=String,Description="Record identifier">',
    ]
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    genotypes = ""
    if samples:
        lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        columns += "\tFORMAT\t" + "\t".join(samples)
        genotypes = "\tGT" + "\t1" * len(samples)
    lines.append(columns)
    for r in records:
        lines.append(
            f"{contig}\t{r.pos}\t{r.ident}\t{r.ref}\t{r.alt}\t.\tPASS\t{TAG}={r.ident}{genotypes}"
        )
    path.write_text("\n".join(lines) + "\n")


def compressed(vcf: Path) -> Path:
    """The same file as `.vcf.gz`, which the other reader takes."""
    out = vcf.with_name(vcf.stem + "_gz.vcf.gz")
    subprocess.run([BCFTOOLS, "view", "-Oz", "-o", str(out), str(vcf)],
                   check=True, capture_output=True)
    return out


def bcftools_norm(vcf: Path) -> dict[str, tuple[int, str, str]]:
    proc = subprocess.run(
        [BCFTOOLS, "norm", "-f", str(FASTA), "-m", "-", str(vcf)],
        check=True, capture_output=True, text=True,
    )
    out = {}
    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        out[fields[2]] = (int(fields[1]), fields[3], fields[4])
    return out


def get_mnv(work: Path, vcf: Path, flags: list[str],
            sample: str | None = None) -> dict[str, set[tuple[int, str, str]]]:
    """What get_MNV wrote for each record, read back out of its VCF: a record
    keeps its INFO, and with it the identifier, under --keep-original-info.
    With `sample`, the run is a `--sample all` cohort and that is the sample
    whose output is read."""
    if sample is None or sample == COHORT[0]:
        proc = subprocess.run(
            [GET_MNV, "--vcf", str(vcf), "--fasta", str(FASTA), "--genes", str(GENES),
             "--convert", "--keep-original-info", "--output-dir", str(work), *flags],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            raise SystemExit(f"get_mnv failed on {vcf.name}:\n{proc.stderr[-3000:]}")
    stem = vcf.name.split(".")[0]
    name = f"{stem}.MNV.vcf" if sample is None else f"{stem}.sample_{sample}.MNV.vcf"
    out: dict[str, set[tuple[int, str, str]]] = {}
    for line in (work / name).read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        ident = next(
            (kv[len(TAG) + 1:] for kv in fields[7].split(";") if kv.startswith(f"{TAG}=")),
            None,
        )
        if ident is not None:
            out.setdefault(ident, set()).add((int(fields[1]), fields[3], fields[4]))
    return out


def compare(records, oracle, got) -> list[str]:
    problems = []
    for r in records:
        want = oracle.get(r.ident)
        seen = got.get(r.ident, set())
        if want is None:
            problems.append(f"{r.ident}: bcftools wrote no record for {r.pos} {r.ref}>{r.alt}")
        elif seen != {want}:
            shown = ", ".join(f"{p} {a}>{b}" for p, a, b in sorted(seen)) or "nothing"
            problems.append(
                f"{r.ident} ({r.kind}, period {r.period}, {r.spelling}): "
                f"{r.pos} {r.ref}>{r.alt} -> bcftools {want[0]} {want[1]}>{want[2]}, "
                f"get_MNV {shown}"
            )
    return problems


def main() -> int:
    if shutil.which(BCFTOOLS) is None:
        print("bcftools not found on PATH; skipping the left-alignment comparison.")
        return 0

    contig, seq = read_genome(FASTA)
    rng = random.Random(SEED)
    records = make_records(seq, rng)

    work = HERE / "work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    plain = work / "indels.vcf"
    write_vcf(plain, contig, len(seq), records)
    cohort = work / "cohort.vcf"
    write_vcf(cohort, contig, len(seq), records, COHORT)

    oracle = bcftools_norm(plain)
    problems = []
    for vcf in (plain, compressed(plain)):
        for problem in compare(records, oracle, get_mnv(work, vcf, ["--left-align"])):
            problems.append(f"[{vcf.name}] {problem}")
    for vcf in (cohort, compressed(cohort)):
        for sample in COHORT:
            got = get_mnv(work, vcf, ["--left-align", "--sample", "all"], sample)
            for problem in compare(records, oracle, got):
                problems.append(f"[{vcf.name} {sample}] {problem}")
    unaligned = len(compare(records, oracle, get_mnv(work, plain, [])))

    moved = [r for r in records if oracle.get(r.ident, (r.pos,))[0] != r.pos]
    shifts = [r.pos - oracle[r.ident][0] for r in moved]
    kinds = Counter(r.kind for r in records)
    periods = Counter(
        str(r.period) if r.period < 7 else "7-100" for r in records if r.kind != "complex"
    )
    print(
        f"{len(records)} records against bcftools norm -f: "
        f"{kinds['deletion']} deletions, {kinds['insertion']} insertions, "
        f"{kinds['complex']} complex"
    )
    print("  through both VCF readers, alone and as a two-sample cohort with --sample all")
    print("  repeat unit of the indels: " + ", ".join(
        f"{period} bp x{periods[period]}" for period in ["1", "2", "3", "4", "5", "6", "7-100"]))
    print(f"  {len(moved)} moved by bcftools, by up to {max(shifts, default=0)} bp")
    print(f"  without --left-align, get_MNV would disagree on {unaligned}")
    print()

    if unaligned == 0:
        print("FAIL: nothing here needed left-aligning, so the comparison proves nothing.")
        return 1
    if not problems:
        print("OK: get_MNV --left-align writes every record as bcftools norm does.")
        return 0
    print(f"{len(problems)} disagreement(s):")
    for problem in problems[:40]:
        print(f"  {problem}")
    if len(problems) > 40:
        print(f"  ... and {len(problems) - 40} more")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
