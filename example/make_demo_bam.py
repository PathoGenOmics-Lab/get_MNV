#!/usr/bin/env python3
"""Regenerate the tiny demo BAM bundled with the example dataset.

The bundled ``example/`` ships a VCF + reference + gene table, but no
alignment, so the GUI read viewer (the IGV-style pileup) has nothing to draw
out of the box. This script synthesises a *tiny* coordinate-sorted, indexed BAM
covering a single real MNV call so the read visualisation works with the
distributed example.

Target locus (from ``G35894.var.snp.vcf`` / the MNV TSV):
    gene Rv2036 (+ strand), codon 93, genomic 2282375-2282377
    reference codon GTT -> GCC  (Val93Ala)
    MNV = two adjacent SNVs: 2282376 T>C and 2282377 T>C

Every synthetic read carries both ALT bases (the call is homozygous, FREQ 100%),
matches the reference everywhere else, and is a clean 100M alignment with high
base/mapping qualities. Reads are staggered (IGV-style) but all span the MNV.

Reproducible: fixed RNG seed, reference bases pulled from the bundled FASTA.

Usage:
    python3 make_demo_bam.py           # writes G35894.demo.bam (+ .bai)
Set SAMTOOLS to point at a samtools binary if it is not on PATH.
"""
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")
FASTA = os.path.join(HERE, "MTB_ancestor.fas")
OUT_BAM = os.path.join(HERE, "G35894.demo.bam")

CHROM = "MTB_anc"
CONTIG_LEN = 4411532

# MNV: genomic position -> ALT base (1-based genomic coordinates).
MNV = {2282376: "C", 2282377: "C"}

# Window of reference we carve reads from (1-based, inclusive).
WIN_START = 2282280
WIN_END = 2282480

READ_LEN = 100
N_READS = 24
SEED = 35894


def faidx_region(start: int, end: int) -> str:
    """Return the reference bases for CHROM:start-end (1-based inclusive)."""
    out = subprocess.run(
        [SAMTOOLS, "faidx", FASTA, f"{CHROM}:{start}-{end}"],
        check=True, capture_output=True, text=True,
    ).stdout
    return "".join(line.strip() for line in out.splitlines() if not line.startswith(">"))


def main() -> int:
    # Ensure the FASTA is indexed (the viewer needs the .fai too).
    subprocess.run([SAMTOOLS, "faidx", FASTA], check=True)

    window = list(faidx_region(WIN_START, WIN_END))
    assert len(window) == WIN_END - WIN_START + 1, "unexpected window length"
    # Sanity-check the reference codon before mutating.
    codon = "".join(window[2282375 - WIN_START: 2282378 - WIN_START])
    assert codon == "GTT", f"reference codon changed: {codon!r} (expected GTT)"

    # Apply the MNV ALT bases in window coordinates.
    for gpos, alt in MNV.items():
        window[gpos - WIN_START] = alt
    mut_codon = "".join(window[2282375 - WIN_START: 2282378 - WIN_START])
    assert mut_codon == "GCC", f"mutated codon wrong: {mut_codon!r} (expected GCC)"
    sample = "".join(window)

    rng = random.Random(SEED)
    # Read starts spread across the window but always spanning both MNV bases.
    lo, hi = WIN_START, max(MNV)  # start <= last MNV pos guarantees coverage
    starts = sorted(rng.randint(lo, hi) for _ in range(N_READS))

    sam_lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        f"@SQ\tSN:{CHROM}\tLN:{CONTIG_LEN}",
        "@RG\tID:demo\tSM:G35894\tPL:ILLUMINA",
        f"@PG\tID:make_demo_bam\tPN:make_demo_bam\tVN:1",
    ]
    for i, start in enumerate(starts, 1):
        seq = sample[start - WIN_START: start - WIN_START + READ_LEN]
        qual = "I" * len(seq)  # Phred 40
        flag = 0 if i % 2 else 16  # alternate strands
        sam_lines.append(
            "\t".join([
                f"demo_read_{i:03d}", str(flag), CHROM, str(start),
                "60", f"{len(seq)}M", "*", "0", "0", seq, qual, "RG:Z:demo",
            ])
        )

    sort = subprocess.run(
        [SAMTOOLS, "sort", "-O", "bam", "-o", OUT_BAM, "-"],
        input="\n".join(sam_lines) + "\n", text=True,
    )
    if sort.returncode != 0:
        return sort.returncode
    subprocess.run([SAMTOOLS, "index", OUT_BAM], check=True)
    print(f"wrote {OUT_BAM} ({N_READS} reads) + .bai; indexed {FASTA}.fai")
    return 0


if __name__ == "__main__":
    sys.exit(main())
