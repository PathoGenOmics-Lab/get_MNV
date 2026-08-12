"""
Mini framework para crear escenarios MNV/indel sintéticos, ejecutar
get_mnv y comparar la salida TSV con expectativas pre-definidas.

Mini-genoma (contig "chr_test", 1300 bp):

  geneA  pos 1-300     hebra +     codones ATG + 98*GCT + TAA
  filler pos 301-400   A's
  geneB  pos 401-700   hebra -     genomic = RC(geneA_cds): TTA + 98*AGC + CAT
  filler pos 701-800   A's
  geneC  pos 801-1200  hebra +     multi-exón
    exon 1   pos 801-900   100 bp  ATG + 32*GCT + G   (codones 1-33 + base 1 codon 34)
    intron   pos 901-1000  100 bp  Ts
    exon 2   pos 1001-1200 200 bp  CT + 65*GCT + TAA  (bases 2-3 codon 34 + codones 35-100)
  filler pos 1201-1300 A's

geneB codon math (hebra -):
  codon N ocupa pos (701-3N) a (703-3N) en genomico
  codon 1   = pos 698-700  (CAT en + strand, ATG en - strand)
  codon 10  = pos 671-673  (AGC en + strand, GCT en - strand)
  codon 100 = pos 401-403  (TTA en + strand, TAA en - strand)
  En mRNA: base 1 codon N = RC(pos 703-3N), base 2 = RC(pos 702-3N), base 3 = RC(pos 701-3N)

geneC codon math (spliced):
  codon 1  = pos 801-803 (en exon 1)
  codon 33 = pos 897-899 (final de exon 1, antes del orphan base en pos 900)
  codon 34 = pos 900 (exon 1) + pos 1001-1002 (exon 2)  *** cruza junction ***
  codon 35 = pos 1003-1005 (primer codon completo de exon 2)
  codon 100 = pos 1198-1200 (final de exon 2, stop)
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

# Repo root = two levels up from this file (tests/scenarios/framework.py)
REPO_ROOT = Path(__file__).resolve().parents[2]

# External tools. Override via env vars (SAMTOOLS, BGZIP, TABIX) when not on PATH.
SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")
BGZIP = os.environ.get("BGZIP", "bgzip")
TABIX = os.environ.get("TABIX", "tabix")


def _default_get_mnv() -> str:
    """Locate the get_mnv binary in standard build locations."""
    for relative in ("target/debug/get_mnv", "target/release/get_mnv", "dist/get_mnv"):
        candidate = REPO_ROOT / relative
        if candidate.exists():
            return str(candidate)
    return "get_mnv"  # fall back to PATH


GET_MNV = os.environ.get("GET_MNV", _default_get_mnv())

CONTIG = "chr_test"
CONTIG_LEN = 1300

CONTIG2 = "chr_test2"
CONTIG2_LEN = 600  # geneD (+ strand) pos 1-300, filler 301-600


def reverse_complement(s: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp[b] for b in reversed(s))


def build_reference_seq() -> str:
    # geneA (+) pos 1-300
    cds_template = "ATG" + "GCT" * 98 + "TAA"
    assert len(cds_template) == 300
    geneA = cds_template
    filler1 = "A" * 100
    # geneB (-) pos 401-700 = RC del CDS estandar
    geneB_minus = reverse_complement(cds_template)
    assert len(geneB_minus) == 300
    filler2 = "A" * 100
    # geneC exon 1 pos 801-900 = ATG + 32*GCT + G (codones 1-33 + base 1 de codon 34)
    geneC_exon1 = "ATG" + "GCT" * 32 + "G"
    assert len(geneC_exon1) == 100
    intron = "T" * 100  # pos 901-1000
    # geneC exon 2 pos 1001-1200 = CT (bases 2-3 codon 34) + 65*GCT + TAA
    geneC_exon2 = "CT" + "GCT" * 65 + "TAA"
    assert len(geneC_exon2) == 200
    filler3 = "A" * 100
    seq = geneA + filler1 + geneB_minus + filler2 + geneC_exon1 + intron + geneC_exon2 + filler3
    assert len(seq) == CONTIG_LEN, f"len = {len(seq)}, expected {CONTIG_LEN}"
    return seq


def build_reference_seq2() -> str:
    # geneD (+) pos 1-300
    cds = "ATG" + "GCT" * 98 + "TAA"
    # geneE (+) pos 301-321: ATG GTA (AAA)x4 TAA. Deleting the G at pos 304
    # frameshifts so codon 2 reads TAA (premature stop at protein position 2);
    # used to test that downstream variants are flagged as past the stop.
    gene_e = "ATG" + "GTA" + "AAA" * 4 + "TAA"  # 21 bp
    filler = gene_e + "A" * (300 - len(gene_e))
    seq = cds + filler
    assert len(seq) == CONTIG2_LEN
    return seq


REFERENCE_SEQ = build_reference_seq()
REFERENCE_SEQ2 = build_reference_seq2()


def codon_pos_plus(codon_index: int, gene_start: int = 1) -> tuple[int, int, int]:
    """Posiciones genomicas 1-based de un codon en un gen + strand."""
    start = gene_start + (codon_index - 1) * 3
    return start, start + 1, start + 2


def codon_pos_minus_strand(codon_index: int, gene_end: int) -> tuple[int, int, int]:
    """Posiciones genomicas (no ordenadas por mRNA) de un codon en gen - strand.

    En mRNA, codon N en gen - strand ocupa (gene_end - 3*(N-1) - 2) a (gene_end - 3*(N-1)).
    Devuelve (low_pos, mid_pos, high_pos) en genomico. high_pos = base 1 del codon (mRNA).
    """
    high = gene_end - 3 * (codon_index - 1)
    return high - 2, high - 1, high


# Default GFF: gene records para geneA, geneB (chr_test) y geneD (chr_test2).
GFF_GENE_ONLY = (
    "##gff-version 3\n"
    f"{CONTIG}\tsynth\tgene\t1\t300\t.\t+\t.\tID=gene-geneA;Name=geneA\n"
    f"{CONTIG}\tsynth\tgene\t401\t700\t.\t-\t.\tID=gene-geneB;Name=geneB\n"
    f"{CONTIG2}\tsynth\tgene\t1\t300\t.\t+\t.\tID=gene-geneD;Name=geneD\n"
    f"{CONTIG2}\tsynth\tgene\t301\t321\t.\t+\t.\tID=gene-geneE;Name=geneE\n"
)

# GFF con CDS multi-exón para geneC (para --gff-features CDS).
# Phase: exon 1 tiene 100 bp = 33 codones + 1 base orphan.
# El siguiente codon (#34) empieza a 2 bases dentro de exon 2 -> phase=2.
GFF_CDS_MULTIEXON = (
    "##gff-version 3\n"
    f"{CONTIG}\tsynth\tgene\t801\t1200\t.\t+\t.\tID=gene-geneC;Name=geneC\n"
    f"{CONTIG}\tsynth\tmRNA\t801\t1200\t.\t+\t.\tID=mrna-geneC;Parent=gene-geneC;Name=geneC\n"
    f"{CONTIG}\tsynth\tCDS\t801\t900\t.\t+\t0\tID=cds-geneC-e1;Parent=mrna-geneC;Name=geneC\n"
    f"{CONTIG}\tsynth\tCDS\t1001\t1200\t.\t+\t2\tID=cds-geneC-e2;Parent=mrna-geneC;Name=geneC\n"
)


# ---------------------------------------------------------------------------
# Modelo de operaciones sobre lecturas
# ---------------------------------------------------------------------------


@dataclass
class Op:
    """Edicion aplicada a una lectura sobre la secuencia de referencia que cubre.

    - kind="snv":  sustituir base en `pos` por `seq` (1 base).
    - kind="ins":  insertar `seq` justo despues de `pos` (estilo VCF anchor).
    - kind="del":  borrar `length` bases empezando en `pos`.
    - kind="skip": saltar `length` bases de referencia sin emitir seq (CIGAR N, p.ej. intron).
                   Empieza en `pos`.
    """

    kind: str
    pos: int
    seq: str = ""
    length: int = 0


@dataclass
class ReadGroup:
    name_prefix: str
    start: int          # 1-based reference start
    length: int         # numero de bases de referencia que consume la lectura
    ops: list[Op] = field(default_factory=list)
    count: int = 20
    strand: str = "+"   # "+" -> FLAG 0, "-" -> FLAG 16
    chrom: str = CONTIG  # contig al que se alinea (chr_test por defecto)
    # Si se define, el grupo se emite como pares: cada lectura lleva el mismo
    # QNAME que su mate, de modo que ambas son una sola molecula.
    mate_start: int | None = None
    mate_length: int = 0
    mate_ops: list[Op] = field(default_factory=list)
    mapq: int = 60          # calidad de mapeo emitida en el SAM
    base_quality: str = "I" # caracter de calidad por base ('I' = Q40)


@dataclass
class VcfRecord:
    pos: int
    ref: str
    alt: str
    chrom: str = CONTIG   # por defecto chr_test; usar CONTIG2 para multi-contig
    af: float | None = None  # si se define, se emite AF= en INFO (rellena original_freq)
    # Genotipo emitido en FORMAT. Con '|' declara fase; el default '1/1' no.
    genotype: str = "1/1"
    phase_set: int | None = None  # si se define, se emite PS en FORMAT


@dataclass
class IvarRecord:
    """Registro iVar TSV. ALT prefijado con + para inserciones, - para deleciones."""
    pos: int
    ref: str
    alt: str
    chrom: str = CONTIG
    pass_flag: str = "TRUE"
    total_dp: int = 30
    alt_dp: int = 20
    alt_freq: float = 1.0


@dataclass
class ExpectedRow:
    """Una fila esperada en la salida TSV de get_mnv.

    Solo se comprueban los campos no None.
    """

    positions: str | None = None
    gene: str | None = None
    reference_bases: str | None = None
    base_changes: str | None = None
    aa_changes: str | None = None
    variant_type: str | None = None
    change_type: str | None = None
    event_class: str | None = None
    event_components: str | None = None
    reference_codon: str | None = None
    snp_codon: str | None = None
    mnv_codon: str | None = None
    snp_reads: str | None = None
    mnv_reads: str | None = None
    snp_forward_reads: str | None = None
    snp_reverse_reads: str | None = None
    mnv_forward_reads: str | None = None
    mnv_reverse_reads: str | None = None
    total_reads: str | None = None
    snp_frequencies: str | None = None
    mnv_frequencies: str | None = None
    event_reads: str | None = None
    event_frequency: str | None = None
    mnv_phasing_support: str | None = None
    mnv_phasing_reads: str | None = None
    frameshift_phasing: str | None = None
    declared_phase: str | None = None
    codon_ld: str | None = None


@dataclass
class Scenario:
    name: str
    description: str
    variants: list[VcfRecord]
    reads: list[ReadGroup]
    expected: list[ExpectedRow]
    expected_row_count: int | None = None
    gff_content: str | None = None         # override (default = GFF_GENE_ONLY)
    gff_features: str | None = None         # --gff-features VALUE
    extra_cli_args: list[str] | None = None # flags extra para get_mnv
    ivar_records: list[IvarRecord] | None = None  # si presente, usa --tsv en vez de --vcf


# ---------------------------------------------------------------------------
# Construccion de archivos
# ---------------------------------------------------------------------------


def write_fasta(path: Path) -> None:
    with path.open("w") as fh:
        fh.write(f">{CONTIG}\n")
        for i in range(0, CONTIG_LEN, 60):
            fh.write(REFERENCE_SEQ[i : i + 60] + "\n")
        fh.write(f">{CONTIG2}\n")
        for i in range(0, CONTIG2_LEN, 60):
            fh.write(REFERENCE_SEQ2[i : i + 60] + "\n")
    subprocess.run(
        [SAMTOOLS, "faidx", str(path)], check=True, capture_output=True
    )


def reference_for(chrom: str) -> str:
    return REFERENCE_SEQ if chrom == CONTIG else REFERENCE_SEQ2


def contig_len_for(chrom: str) -> int:
    return CONTIG_LEN if chrom == CONTIG else CONTIG2_LEN


def write_gff(path: Path, content: str | None = None) -> None:
    path.write_text(content or GFF_GENE_ONLY)


def write_vcf(path: Path, records: Iterable[VcfRecord]) -> None:
    header = (
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={CONTIG},length={CONTIG_LEN}>\n"
        f"##contig=<ID={CONTIG2},length={CONTIG2_LEN}>\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n'
        '##INFO=<ID=AF,Number=A,Type=Float,Description="allele frequency">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    with path.open("w") as fh:
        fh.write(header)
        for v in records:
            info = "DP=30"
            if v.af is not None:
                info += f";AF={v.af:.4f}"
            if v.phase_set is None:
                fmt, sample = "GT", v.genotype
            else:
                fmt, sample = "GT:PS", f"{v.genotype}:{v.phase_set}"
            fh.write(
                f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t100\tPASS\t{info}\t{fmt}\t{sample}\n"
            )


def write_ivar_tsv(path: Path, records: Iterable[IvarRecord]) -> None:
    cols = ["REGION", "POS", "REF", "ALT", "REF_DP", "ALT_DP", "ALT_FREQ", "TOTAL_DP", "PASS"]
    with path.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in records:
            ref_dp = max(0, r.total_dp - r.alt_dp)
            row = [
                r.chrom,
                str(r.pos),
                r.ref,
                r.alt,
                str(ref_dp),
                str(r.alt_dp),
                f"{r.alt_freq:.4f}",
                str(r.total_dp),
                r.pass_flag,
            ]
            fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Construccion de BAM
# ---------------------------------------------------------------------------


def _build_read_seq_cigar(start: int, length: int, ops: list[Op], chrom: str) -> tuple[str, str]:
    """Devuelve (SEQ, CIGAR) para una lectura que consume `length` bases de
    referencia de `chrom` empezando en `start` (1-based), aplicando las operaciones."""
    ref = reference_for(chrom)
    snv = {o.pos: o.seq for o in ops if o.kind == "snv"}
    ins_after = {o.pos: o.seq for o in ops if o.kind == "ins"}
    del_at = {o.pos: o.length for o in ops if o.kind == "del"}
    skip_at = {o.pos: o.length for o in ops if o.kind == "skip"}

    seq_chars: list[str] = []
    cigar: list[tuple[int, str]] = []

    def push(n: int, op_char: str) -> None:
        if n <= 0:
            return
        if cigar and cigar[-1][1] == op_char:
            cigar[-1] = (cigar[-1][0] + n, op_char)
        else:
            cigar.append((n, op_char))

    cur = start
    end = start + length - 1
    while cur <= end:
        if cur in skip_at:
            push(skip_at[cur], "N")
            cur += skip_at[cur]
            continue
        if cur in del_at:
            push(del_at[cur], "D")
            cur += del_at[cur]
            continue
        base = snv.get(cur, ref[cur - 1])
        seq_chars.append(base)
        push(1, "M")
        if cur in ins_after:
            ins_seq = ins_after[cur]
            seq_chars.extend(list(ins_seq))
            push(len(ins_seq), "I")
        cur += 1

    cigar_str = "".join(f"{n}{c}" for n, c in cigar)
    return "".join(seq_chars), cigar_str


def write_bam(path: Path, work_dir: Path, reads: list[ReadGroup]) -> None:
    sam_path = work_dir / "reads.sam"
    with sam_path.open("w") as fh:
        fh.write("@HD\tVN:1.6\tSO:coordinate\n")
        fh.write(f"@SQ\tSN:{CONTIG}\tLN:{CONTIG_LEN}\n")
        fh.write(f"@SQ\tSN:{CONTIG2}\tLN:{CONTIG2_LEN}\n")
        fh.write("@RG\tID:test\tSM:sample\n")
        for grp in reads:
            chrom = getattr(grp, "chrom", CONTIG)
            seq, cigar = _build_read_seq_cigar(grp.start, grp.length, grp.ops, chrom)
            qual = grp.base_quality * len(seq)
            paired = getattr(grp, "mate_start", None) is not None
            if paired:
                mate_seq, mate_cigar = _build_read_seq_cigar(
                    grp.mate_start, grp.mate_length, grp.mate_ops, chrom
                )
                mate_qual = grp.base_quality * len(mate_seq)
                # 1 paired, 2 proper pair, 32 mate reverse, 64 first, 128 last.
                segments = [
                    (1 + 2 + 32 + 64, grp.start, cigar, seq, qual, grp.mate_start),
                    (1 + 2 + 16 + 128, grp.mate_start, mate_cigar, mate_seq, mate_qual, grp.start),
                ]
            else:
                segments = [
                    (0 if grp.strand == "+" else 16, grp.start, cigar, seq, qual, None)
                ]
            for i in range(grp.count):
                name = f"{grp.name_prefix}_{i:03d}"
                for flag, start, cig, sq, ql, mate_pos in segments:
                    fh.write(
                        "\t".join(
                            [
                                name,
                                str(flag),
                                chrom,
                                str(start),
                                str(grp.mapq),
                                cig,
                                "=" if mate_pos is not None else "*",
                                str(mate_pos) if mate_pos is not None else "0",
                                "0",
                                sq,
                                ql,
                                "RG:Z:test",
                            ]
                        )
                        + "\n"
                    )
    subprocess.run(
        [SAMTOOLS, "sort", "-O", "bam", "-o", str(path), str(sam_path)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [SAMTOOLS, "index", str(path)], check=True, capture_output=True
    )


# ---------------------------------------------------------------------------
# Ejecucion + comparacion
# ---------------------------------------------------------------------------


def run_get_mnv(
    work_dir: Path,
    variant_input: Path,
    fasta: Path,
    gff: Path,
    bam: Path | None,
    gff_features: str | None = None,
    extra_args: list[str] | None = None,
    input_flag: str = "--vcf",
) -> Path:
    cmd = [
        GET_MNV,
        input_flag,
        str(variant_input),
        "--fasta",
        str(fasta),
        "--gff",
        str(gff),
    ]
    if gff_features is not None:
        cmd += ["--gff-features", gff_features]
    if bam is not None:
        cmd += ["--bam", str(bam)]
    if extra_args:
        cmd += list(extra_args)
    proc = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True)
    log_path = work_dir / "get_mnv.log"
    log_path.write_text(
        "CMD: " + " ".join(cmd) + "\n" + proc.stdout + "\n--- stderr ---\n" + proc.stderr
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"get_mnv failed (rc={proc.returncode}). See {log_path}\n{proc.stderr}"
        )
    out = work_dir / (variant_input.stem + ".MNV.tsv")
    if not out.exists():
        raise RuntimeError(f"get_mnv did not produce expected output {out}")
    return out


def parse_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        parts = line.split("\t")
        rows.append(dict(zip(header, parts)))
    return header, rows


_FIELD_BY_ATTR = {
    "positions": "Positions",
    "gene": "Gene",
    "reference_bases": "Reference Bases",
    "base_changes": "Base Changes",
    "aa_changes": "AA Changes",
    "variant_type": "Variant Type",
    "change_type": "Change Type",
    "event_class": "Event Class",
    "event_components": "Event Components",
    "reference_codon": "Reference Codon",
    "snp_codon": "SNP Codon",
    "mnv_codon": "MNV Codon",
    "snp_reads": "SNP Reads",
    "mnv_reads": "MNV Reads",
    "snp_forward_reads": "SNP Forward Reads",
    "snp_reverse_reads": "SNP Reverse Reads",
    "mnv_forward_reads": "MNV Forward Reads",
    "mnv_reverse_reads": "MNV Reverse Reads",
    "total_reads": "Total Reads",
    "snp_frequencies": "SNP Frequencies",
    "mnv_frequencies": "MNV Frequencies",
    "event_reads": "Event Reads",
    "event_frequency": "Event Frequency",
    "mnv_phasing_support": "MNV Phasing Support",
    "mnv_phasing_reads": "MNV Phasing Reads",
    "frameshift_phasing": "Frameshift Phasing",
    "declared_phase": "Declared Phase",
    "codon_ld": "Haplotype LD",
}


def _find_row(rows: list[dict[str, str]], exp: ExpectedRow) -> dict[str, str] | None:
    for row in rows:
        ok = True
        for attr, col in _FIELD_BY_ATTR.items():
            want = getattr(exp, attr)
            if want is None:
                continue
            if row.get(col, "") != want:
                ok = False
                break
        if ok:
            return row
    return None


def compare(expected: list[ExpectedRow], rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    for i, exp in enumerate(expected):
        match = _find_row(rows, exp)
        if match is None:
            wanted = {
                col: getattr(exp, attr)
                for attr, col in _FIELD_BY_ATTR.items()
                if getattr(exp, attr) is not None
            }
            errors.append(
                f"  fila esperada #{i + 1} no encontrada: {wanted}"
            )
    return errors


def run_scenario(scenario: Scenario, base_work: Path) -> tuple[bool, list[str], Path]:
    work = base_work / scenario.name
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    fasta = work / "ref.fasta"
    gff = work / "genes.gff3"
    bam = work / "reads.bam"

    write_fasta(fasta)
    write_gff(gff, scenario.gff_content)
    write_bam(bam, work, scenario.reads)

    if scenario.ivar_records is not None:
        variant_input = work / "variants.tsv"
        write_ivar_tsv(variant_input, scenario.ivar_records)
        input_flag = "--tsv"
    else:
        variant_input = work / "variants.vcf"
        write_vcf(variant_input, scenario.variants)
        input_flag = "--vcf"

    out = run_get_mnv(
        work, variant_input, fasta, gff, bam,
        gff_features=scenario.gff_features,
        extra_args=scenario.extra_cli_args,
        input_flag=input_flag,
    )
    _, rows = parse_tsv(out)
    errors = compare(scenario.expected, rows)
    if scenario.expected_row_count is not None and len(rows) != scenario.expected_row_count:
        errors.append(
            f"  numero de filas: esperado {scenario.expected_row_count}, obtenido {len(rows)}"
        )
    return (not errors), errors, out
