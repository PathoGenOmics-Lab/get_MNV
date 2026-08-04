#!/usr/bin/env python3
"""Generate the small, self-contained genomes used by the differential suite.

The fixtures are committed so the comparison is hermetic and the inputs can be
read by eye. Re-run this only to change or extend them:

    python3 tests/differential/make_fixtures.py

Every fixture ships one GFF3 that BOTH get_MNV and `bcftools csq` read, so a
disagreement can never be blamed on the two tools seeing different annotation.
That is why the CDS rows carry `Parent`/`ID` (which bcftools csq needs) *and*
`gene=` (which get_MNV uses for the gene name); each tool ignores what it does
not recognise.
"""

from __future__ import annotations

import random
from pathlib import Path

FIXTURES = Path(__file__).resolve().parent / "fixtures"
STOPS = {"TAA", "TAG", "TGA"}


def sense_codons(rng: random.Random, count: int) -> list[str]:
    """`count` random codons, none of them a stop."""
    out = []
    while len(out) < count:
        codon = "".join(rng.choice("ACGT") for _ in range(3))
        if codon not in STOPS:
            out.append(codon)
    return out


def write_fasta(path: Path, name: str, sequence: str) -> None:
    lines = [f">{name}"]
    lines += [sequence[i : i + 60] for i in range(0, len(sequence), 60)]
    path.write_text("\n".join(lines) + "\n")


def write_vcf(path: Path, contig: str, length: int, records: list[tuple]) -> None:
    header = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={contig},length={length}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    body = [f"{contig}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t." for pos, ref, alt in records]
    path.write_text("\n".join(header + body) + "\n")


def other_base(base: str) -> str:
    return "ACGT".replace(base.upper(), "")[0]


def gff3(
    contig: str,
    gene: str,
    span: tuple[int, int],
    segments: list[tuple[int, int]],
    exons: bool = True,
) -> str:
    """A GFF3 both tools accept: gene + mRNA + CDS (and optionally exon) rows.

    `exons=False` mirrors how an unspliced CDS is annotated in practice: the
    NCBI GFF3 for SARS-CoV-2 gives ORF1ab a gene and a two-segment CDS with a
    ribosomal-slippage exception, and no exon features at all. Writing exon rows
    for those segments would be inventing exons that do not exist, and it
    changes what `bcftools csq` reports, so the fixture must not do it.
    """
    start, end = span
    rows = [
        "##gff-version 3",
        f"##sequence-region {contig} 1 {end + 200}",
        f"{contig}\tfixture\tgene\t{start}\t{end}\t.\t+\t.\t"
        f"ID=gene:{gene};biotype=protein_coding;Name={gene}",
        f"{contig}\tfixture\tmRNA\t{start}\t{end}\t.\t+\t.\t"
        f"ID=transcript:{gene};Parent=gene:{gene};biotype=protein_coding",
    ]
    for seg_start, seg_end in segments:
        # bcftools csq needs exon rows to place splice sites; get_MNV analyses
        # CDS rows and ignores the exons.
        if exons:
            rows.append(
                f"{contig}\tfixture\texon\t{seg_start}\t{seg_end}\t.\t+\t.\t"
                f"Parent=transcript:{gene}"
            )
        rows.append(
            f"{contig}\tfixture\tCDS\t{seg_start}\t{seg_end}\t.\t+\t0\t"
            f"ID=cds:{gene};Parent=transcript:{gene};gene={gene}"
        )
    return "\n".join(rows) + "\n"


def build_slippage() -> None:
    """A ribosomal-slippage join: `join(101..400, 400..702)`.

    The two CDS segments overlap by the re-read base, exactly how SARS-CoV-2
    ORF1ab is annotated. There is no intron, so no splice site and no exon-exon
    junction, and the spliced CDS reads base 400 twice.
    """
    rng = random.Random(42)
    # 201 codons: ATG + 199 sense codons + TAA. Segment 1 contributes 300 nt and
    # segment 2 the remaining 303, with base 400 counted in both.
    while True:
        cds = "ATG" + "".join(sense_codons(rng, 199)) + "TAA"
        # The genome holds one physical base at 400, so the re-read base must be
        # the same nucleotide in both of its CDS roles.
        if cds[299] == cds[300]:
            break

    genome = [rng.choice("ACGT") for _ in range(900)]
    for i in range(300):
        genome[100 + i] = cds[i]  # 101..400
    for i in range(302):
        genome[400 + i] = cds[301 + i]  # 401..702
    sequence = "".join(genome)

    write_fasta(FIXTURES / "slippage.fa", "v1", sequence)
    (FIXTURES / "slippage.gff").write_text(
        gff3("v1", "orf1ab", (101, 702), [(101, 400), (400, 702)], exons=False)
    )
    # Plain coding bases either side of the slippage boundary. None of these is a
    # splice site, because the segments are not separated by an intron.
    positions = [396, 399, 401, 402, 405, 450]
    write_vcf(
        FIXTURES / "slippage.vcf",
        "v1",
        len(sequence),
        [(p, sequence[p - 1], other_base(sequence[p - 1])) for p in positions],
    )


def build_spliced() -> None:
    """A genuinely spliced two-exon transcript: exon 101..250, intron, exon 401..553.

    The positive control for the slippage fixture: here the splice terms *must*
    fire, and both tools must agree on which positions carry them.
    """
    rng = random.Random(7)
    exon1_len, exon2_len = 150, 153  # 303 nt = 101 codons
    cds = "ATG" + "".join(sense_codons(rng, 99)) + "TAA"
    assert len(cds) == exon1_len + exon2_len

    genome = [rng.choice("ACGT") for _ in range(800)]
    for i in range(exon1_len):
        genome[100 + i] = cds[i]  # 101..250
    for i in range(exon2_len):
        genome[400 + i] = cds[exon1_len + i]  # 401..553
    # Canonical GT..AG intron so the annotation is biologically plausible.
    genome[250], genome[251] = "G", "T"  # intron first two bases (251, 252)
    genome[398], genome[399] = "A", "G"  # intron last two bases (399, 400)
    sequence = "".join(genome)

    write_fasta(FIXTURES / "spliced.fa", "v1", sequence)
    (FIXTURES / "spliced.gff").write_text(
        gff3("v1", "geneS", (101, 553), [(101, 250), (401, 553)])
    )
    positions = [
        248,  # exon 1, third-from-last base: splice region
        250,  # exon 1 last base: splice region
        251,  # intron +1: essential donor
        252,  # intron +2: essential donor
        255,  # intron +5: intronic splice region
        320,  # deep intron: no consequence
        399,  # intron -2: essential acceptor
        400,  # intron -1: essential acceptor
        401,  # exon 2 first base: splice region
        500,  # deep in exon 2: plain coding change
    ]
    write_vcf(
        FIXTURES / "spliced.vcf",
        "v1",
        len(sequence),
        [(p, sequence[p - 1], other_base(sequence[p - 1])) for p in positions],
    )


def main() -> None:
    FIXTURES.mkdir(parents=True, exist_ok=True)
    build_slippage()
    build_spliced()
    for path in sorted(FIXTURES.iterdir()):
        print(f"wrote {path.relative_to(FIXTURES.parent)} ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
