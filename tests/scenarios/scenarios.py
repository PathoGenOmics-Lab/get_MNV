"""Escenarios de test con expectativas pre-definidas.

Mini-genoma chr_test (1300 bp):
  geneA  pos 1-300     hebra +     ATG + 98*GCT + TAA
  geneB  pos 401-700   hebra -     ATG + 98*GCT + TAA en mRNA (genomico es RC)
  geneC  pos 801-1200  hebra +     multi-exon (exon1: 801-900, intron: 901-1000, exon2: 1001-1200)

geneB - strand math:
  codon 10 mRNA = GCT, ocupa pos genomicas 671-673 (genomic ref = AGC)
  Base 1 mRNA = RC(pos 673), base 2 = RC(pos 672), base 3 = RC(pos 671)
  Para cambiar base 1 mRNA, modificar pos 673 (RC del nuevo nucleotido)

geneC spliced math:
  exon 1 spliced offsets 1-100  (genomic 801-900)
  exon 2 spliced offsets 101-300 (genomic 1001-1200)
  codon 34 cruza junction: base 1 = pos 900 (exon 1), bases 2-3 = pos 1001-1002 (exon 2)
  codon 50 = pos 1048-1050 (exon 2)
"""

from framework import (
    CONTIG2,
    GFF_CDS_MULTIEXON,
    GFF_GENE_ONLY,
    ExpectedRow,
    IvarRecord,
    Op,
    ReadGroup,
    Scenario,
    VcfRecord,
)


# ---------------------------------------------------------------------------
# 1. SNP no-sinonimo simple
# ---------------------------------------------------------------------------
# Codon 10 (pos 28-30) GCT -> ACT.  Ala10Thr.
scenario_snp_simple = Scenario(
    name="01_snp_simple",
    description="SNP no-sinonimo simple: pos 28 G>A en codon 10 (GCT->ACT, Ala10Thr), 20/20 reads alt",
    variants=[VcfRecord(pos=28, ref="G", alt="A")],
    reads=[
        ReadGroup(
            name_prefix="r_alt",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=28, seq="A")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28",
            gene="geneA",
            reference_bases="G",
            base_changes="A",
            aa_changes="Ala10Thr",
            variant_type="SNP",
            change_type="Non-synonymous",
            reference_codon="GCT",
            snp_codon="ACT",
            snp_reads="20",
            total_reads="20",
            snp_frequencies="1.0000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 2a. SNP/MNV: dos SNV en mismo codon (no adyacentes), todas las lecturas
#    llevan ambas SNV. get_mnv reporta como SNP/MNV porque la VCF tiene
#    dos records SNV independientes.
# ---------------------------------------------------------------------------
# Codon 10 GCT -> TCA  (pos 28 G>T, pos 30 T>A).  Ala10Ser.
scenario_snp_mnv_full = Scenario(
    name="02a_snp_mnv_full_phasing",
    description="SNP/MNV: pos28 G>T + pos30 T>A juntos en TODAS las lecturas (haplotipo MNV completo)",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_mnv",
            start=1,
            length=150,
            ops=[
                Op(kind="snv", pos=28, seq="T"),
                Op(kind="snv", pos=30, seq="A"),
            ],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            gene="geneA",
            variant_type="SNP/MNV",
            aa_changes="Ala10Ser",
            reference_codon="GCT",
            mnv_codon="TCA",
            change_type="Non-synonymous",
            # ningun read lleva una sola SNV -> SNP support = 0
            snp_reads="0, 0",
            mnv_reads="20",
            total_reads="20",
            snp_frequencies="0.0000, 0.0000",
            mnv_frequencies="1.0000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 2b. VCF MNP -> get_mnv lo decompone como SNP/MNV
# ---------------------------------------------------------------------------
# Un solo registro VCF MNP (REF=GC, ALT=TA en pos 28-29) se descompone
# internamente en dos SNVs (pos 28 G>T y pos 29 C>A). Ambas caen en codon 10
# GCT -> TAT (Ala10Tyr). get_mnv reporta como SNP/MNV con event_class=mnv.
scenario_mnp_decomposed = Scenario(
    name="02b_vcf_mnp_decomposed",
    description="VCF MNP REF=GC ALT=TA en pos 28-29: decompuesto a SNP/MNV (Ala10Tyr), event_class=mnv",
    variants=[VcfRecord(pos=28, ref="GC", alt="TA")],
    reads=[
        ReadGroup(
            name_prefix="r_mnp",
            start=1,
            length=150,
            ops=[
                Op(kind="snv", pos=28, seq="T"),
                Op(kind="snv", pos=29, seq="A"),
            ],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 29",
            gene="geneA",
            variant_type="SNP/MNV",
            aa_changes="Ala10Tyr",
            reference_codon="GCT",
            mnv_codon="TAT",
            change_type="Non-synonymous",
            event_class="mnv",
            event_components="SNV:28:G>T | SNV:29:C>A",
            mnv_reads="20",
            mnv_frequencies="1.0000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 3. SNP/MNV mixto: solo un tercio de las lecturas llevan ambas SNV
# ---------------------------------------------------------------------------
# pos28 G>A + pos30 T>G en codon 10 GCT:
#   solo pos28 -> ACT = Thr
#   solo pos30 -> GCG = Ala (sinonimo)
#   ambas      -> ACG = Thr
# 10 reads con ambas, 10 con pos28 solo, 10 con pos30 solo. Total 30.
scenario_snp_mnv_mixed = Scenario(
    name="03_snp_mnv_mixed",
    description="SNP/MNV mixto: 10 reads con pos28+pos30, 10 solo pos28, 10 solo pos30 (30 total)",
    variants=[
        VcfRecord(pos=28, ref="G", alt="A"),
        VcfRecord(pos=30, ref="T", alt="G"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_both",
            start=1,
            length=150,
            ops=[
                Op(kind="snv", pos=28, seq="A"),
                Op(kind="snv", pos=30, seq="G"),
            ],
            count=10,
        ),
        ReadGroup(
            name_prefix="r_snv1",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=28, seq="A")],
            count=10,
        ),
        ReadGroup(
            name_prefix="r_snv2",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=30, seq="G")],
            count=10,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            gene="geneA",
            variant_type="SNP/MNV",
            reference_codon="GCT",
            mnv_codon="ACG",
            snp_reads="10, 10",
            mnv_reads="10",
            total_reads="30",
            snp_frequencies="0.3333, 0.3333",
            mnv_frequencies="0.3333",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 4. Insercion in-frame en CDS (3 bp = 1 codon)
# ---------------------------------------------------------------------------
# Insertar GCT despues de pos 30.  Anhade un nuevo residuo Ala.
scenario_ins_inframe = Scenario(
    name="04_ins_inframe_cds",
    description="Insercion in-frame de GCT (1 Ala) entre pos 30 y 31 del CDS, 20 lecturas",
    variants=[VcfRecord(pos=30, ref="T", alt="TGCT")],
    reads=[
        ReadGroup(
            name_prefix="r_ins",
            start=1,
            length=147,  # 147 ref bases consumidas; SEQ = 150 con la insercion
            ops=[Op(kind="ins", pos=30, seq="GCT")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            gene="geneA",
            reference_bases="T",
            base_changes="TGCT",
            variant_type="INDEL",
            change_type="In-frame Indel",
            event_class="insertion",
            event_components="INS:30:+GCT",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 5. Delecion frameshift (1 bp) en CDS
# ---------------------------------------------------------------------------
# Borrar pos 31 (primera base de codon 11).  Frameshift desde codon 11.
# VCF: POS=30 REF=TG ALT=T (anchor T en pos 30, delete G en pos 31)
scenario_del_frameshift = Scenario(
    name="05_del_frameshift_cds",
    description="Delecion frameshift de 1bp en pos 31, 20 lecturas (Ala11Leufs)",
    variants=[VcfRecord(pos=30, ref="TG", alt="T")],
    reads=[
        ReadGroup(
            name_prefix="r_del",
            start=1,
            length=151,
            ops=[Op(kind="del", pos=31, length=1)],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            gene="geneA",
            reference_bases="TG",
            base_changes="T",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            event_class="deletion",
            event_components="DEL:31:G",
            aa_changes="Ala11Leufs",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 6. Haplotipo combinado: insercion + SNV cercana en las mismas lecturas
# ---------------------------------------------------------------------------
# pos 30 ins GCT  +  pos 36 G>A (cambia codon 12 GCT -> ACT = Ala12Thr)
# Las 20 lecturas llevan AMBAS modificaciones en cis.
scenario_indel_plus_snv = Scenario(
    name="06_indel_plus_snv_haplotype",
    description="Haplotipo combinado: insercion in-frame en pos 30 + SNV pos 36 G>A (Ala12Thr) en cis",
    variants=[
        VcfRecord(pos=30, ref="T", alt="TGCT"),
        VcfRecord(pos=36, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_complex",
            start=1,
            length=147,
            ops=[
                Op(kind="ins", pos=30, seq="GCT"),
                Op(kind="snv", pos=36, seq="A"),
            ],
            count=20,
        ),
    ],
    expected=[
        # ambos eventos deben aparecer en la TSV
        ExpectedRow(
            positions="30",
            variant_type="INDEL",
            event_class="insertion",
            event_components="INS:30:+GCT",
            event_reads="20",
        ),
        ExpectedRow(
            positions="36",
            variant_type="SNP",
            reference_codon="GCT",
            snp_codon="GCA",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 7. Frameshift downstream propagation: del frameshift + SNV downstream
# ---------------------------------------------------------------------------
# Del de 1 bp en pos 31 (frameshift desde codon 11).
# SNV pos 39 T>A (en codon 13 GCT). Las 20 lecturas llevan ambos cambios.
# Distancia indel<->SNV = 8 bp (> 3 bp) -> NO se combinan en complex_indel,
# pero el frameshift SI debe propagarse al SNV downstream (marcador "(fs)").
scenario_fs_plus_downstream_snv = Scenario(
    name="07_fs_del_plus_downstream_snv",
    description="Frameshift del en pos 31 + SNV downstream pos 39 a >3bp: el SNV downstream debe llevar marcador (fs)",
    variants=[
        VcfRecord(pos=30, ref="TG", alt="T"),
        VcfRecord(pos=39, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_fs_snv",
            start=1,
            length=151,
            ops=[
                Op(kind="del", pos=31, length=1),
                Op(kind="snv", pos=39, seq="A"),
            ],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            event_class="deletion",
            event_components="DEL:31:G",
            aa_changes="Ala11Leufs",
            event_reads="20",
            event_frequency="1.0000",
        ),
        ExpectedRow(
            positions="39",
            variant_type="SNP",
            change_type="Synonymous (frameshift)",
            aa_changes="Ala13Ala (fs)",
            reference_codon="GCT",
            snp_codon="GCA",
            snp_reads="20",
            snp_frequencies="1.0000",
        ),
    ],
    expected_row_count=2,
)


# ---------------------------------------------------------------------------
# 8. In-frame insertion DENTRO del codon + MNV en mismo codon (todo en cis)
# ---------------------------------------------------------------------------
# Insertar GCT despues de pos 29 (entre pos 29=C y pos 30=T): cae DENTRO
# del codon 10. Mas dos SNV pos 28 G>T y pos 30 T>A en mismo codon.
# Las 20 lecturas llevan todo en cis: en la muestra hay UNA sola especie
# molecular. Get_mnv lee los haplotipos de las moleculas, asi que emite el
# haplotipo completo y no sus subconjuntos, que serian moleculas inexistentes:
#   1) complex_indel SNV:28 + INS:29 + SNV:30     -> Ala10delinsSerLeu
#   2) SNP/MNV pos28+pos30 marcado Indel overlap  -> AA Unknown
#   3) insertion INS:29 sola (el registro del VCF) -> Ala10_Ala11insLeu
scenario_inframe_ins_with_mnv = Scenario(
    name="08_inframe_ins_inside_codon_with_mnv",
    description="Ins in-frame DENTRO de codon 10 + MNV pos28+pos30 todo en cis: 3 filas (el haplotipo completo, MNV-overlap, ins sola); los subconjuntos no son moleculas reales",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=29, ref="C", alt="CGCT"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_combo",
            start=1,
            length=147,
            ops=[
                Op(kind="snv", pos=28, seq="T"),
                Op(kind="ins", pos=29, seq="GCT"),
                Op(kind="snv", pos=30, seq="A"),
            ],
            count=20,
        ),
    ],
    expected=[
        # El haplotipo completo, que es lo que llevan las 20 moleculas
        ExpectedRow(
            event_class="complex_indel",
            event_components="SNV:28:G>T | INS:29:+GCT | SNV:30:T>A",
            variant_type="INDEL",
            change_type="In-frame Indel",
            aa_changes="Ala10delinsSerLeu",
            event_reads="20",
            event_frequency="1.0000",
        ),
        # MNV solapada: pos28+pos30, Indel overlap, AA Unknown
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            change_type="Indel overlap",
            aa_changes="Unknown",
            mnv_codon="TCA",
            mnv_reads="20",
            mnv_frequencies="1.0000",
        ),
        # insertion sola: el registro del VCF anotado por su cuenta
        ExpectedRow(
            positions="29",
            variant_type="INDEL",
            event_class="insertion",
            event_components="INS:29:+GCT",
            aa_changes="Ala10_Ala11insLeu",
            change_type="In-frame Indel",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=3,
)


# ---------------------------------------------------------------------------
# 9. Frameshift del + SNV solapado en mismo codon (todo en cis)
# ---------------------------------------------------------------------------
# Borrar 2 bp en pos 29-30 (los dos ultimos del codon 10): frameshift.
# Mas SNV pos 28 G>T (la unica base no borrada del codon 10).
# Las 20 lecturas llevan todo. Get_mnv emite 3 filas:
#   1) SNP pos28 marcado Indel overlap / AA Unknown
#   2) INDEL solo (DEL:29-30:CT)  -> Event Reads = 0 (nadie lleva SOLO la del)
#   3) complex_indel SNV:28 + DEL:29-30 -> Ala10Cysfs, freq 1.0
scenario_fs_del_with_mnv = Scenario(
    name="09_fs_del_with_snv_overlap",
    description="Frameshift del 2bp consume codon 10 + SNV pos28 en cis: SNP-overlap, INDEL solo (sin reads), complex_indel con frameshift",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=28, ref="GCT", alt="G"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_fs_mnv",
            start=1,
            length=152,
            ops=[
                Op(kind="snv", pos=28, seq="T"),
                Op(kind="del", pos=29, length=2),
            ],
            count=20,
        ),
    ],
    expected=[
        # SNV solapando deletion -> Indel overlap, AA Unknown
        ExpectedRow(
            positions="28",
            variant_type="SNP",
            change_type="Indel overlap",
            aa_changes="Unknown",
            event_class="snp",
            event_components="SNV:28:G>T",
            snp_reads="20",
            snp_frequencies="1.0000",
        ),
        # Indel solo: no hay lecturas con SOLO la del (todas llevan SNV tambien)
        ExpectedRow(
            event_class="deletion",
            event_components="DEL:29-30:CT",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            aa_changes="Ala10Glyfs",
            event_reads="0",
            event_frequency="0.0000",
        ),
        # complex_indel: SNV + del en cis con frameshift
        ExpectedRow(
            event_class="complex_indel",
            event_components="SNV:28:G>T | DEL:29-30:CT",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            aa_changes="Ala10Cysfs",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=3,
)


# ---------------------------------------------------------------------------
# HEBRA NEGATIVA (geneB)
# ---------------------------------------------------------------------------

# 10. SNV non-sin en geneB (- strand)
# Cambio en pos 673 C>T (genomic). mRNA base 1 codon 10: RC(C)=G -> RC(T)=A
# Codon 10 mRNA: GCT -> ACT = Thr. Ala10Thr.
scenario_minus_snp_simple = Scenario(
    name="10_minus_snp_simple",
    description="SNP non-sin en geneB (hebra -): pos 673 C>T -> codon 10 GCT->ACT (Ala10Thr); los codones se reportan en la orientacion del transcrito",
    variants=[VcfRecord(pos=673, ref="C", alt="T")],
    reads=[
        ReadGroup(
            name_prefix="r_minus_snp",
            start=600,
            length=100,
            ops=[Op(kind="snv", pos=673, seq="T")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="673",
            gene="geneB",
            reference_bases="C",
            base_changes="T",
            aa_changes="Ala10Thr",
            variant_type="SNP",
            change_type="Non-synonymous",
            reference_codon="GCT",  # orientacion de transcrito: traduce a Ala
            snp_codon="ACT",        # traduce a Thr, el aminoacido reportado
            snp_reads="20",
            total_reads="20",
            snp_frequencies="1.0000",
        ),
    ],
    expected_row_count=1,
)


# 11. MNV en geneB codon 10 (- strand)
# 2 SNVs en pos 671 (A>T) y pos 673 (C>A) -> mRNA bases 3 (RC(A)=T->RC(T)=A) y 1 (RC(C)=G->RC(A)=T)
# mRNA codon 10: GCT -> TCA = Ser. Ala10Ser.
scenario_minus_mnv = Scenario(
    name="11_minus_mnv_same_codon",
    description="MNV en geneB codon 10: pos 671 A>T + pos 673 C>A en cis (codon GCT->TCA, Ala10Ser); los codones se reportan en la orientacion del transcrito",
    variants=[
        VcfRecord(pos=671, ref="A", alt="T"),
        VcfRecord(pos=673, ref="C", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_minus_mnv",
            start=600,
            length=100,
            ops=[
                Op(kind="snv", pos=671, seq="T"),
                Op(kind="snv", pos=673, seq="A"),
            ],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="671, 673",
            gene="geneB",
            variant_type="SNP/MNV",
            aa_changes="Ala10Ser",
            reference_codon="GCT",  # orientacion de transcrito: traduce a Ala
            mnv_codon="TCA",        # traduce a Ser, el aminoacido reportado
            change_type="Non-synonymous",
            mnv_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
)


# 12. Deleción frameshift en geneB
# Delete 1bp en pos 671 (anchor en pos 670). El frameshift afecta codon 10
# hacia adelante en mRNA = posiciones genomicas menores.
scenario_minus_fs_del = Scenario(
    name="12_minus_fs_del",
    description="Frameshift del 1bp en pos 671 (geneB hebra -): codon 10 con frameshift",
    variants=[VcfRecord(pos=670, ref="CA", alt="C")],
    reads=[
        ReadGroup(
            name_prefix="r_minus_fs",
            start=600,
            length=100,
            ops=[Op(kind="del", pos=671, length=1)],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="670",
            gene="geneB",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            event_class="deletion",
            event_components="DEL:671:A",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=1,
)


# ---------------------------------------------------------------------------
# MULTI-EXON (geneC, requires --gff-features CDS)
# ---------------------------------------------------------------------------

# 13. SNV no-sin en exon 2 (no junction)
# Codon 50 = pos 1048-1050 (exon 2). Cambio pos 1048 G>A: GCT -> ACT = Thr.
# Ala50Thr.
scenario_multiexon_snp_exon2 = Scenario(
    name="13_multiexon_snp_exon2",
    description="SNP no-sin en exon 2 de geneC: pos 1048 G>A en codon 50 (GCT->ACT, Ala50Thr)",
    variants=[VcfRecord(pos=1048, ref="G", alt="A")],
    reads=[
        ReadGroup(
            name_prefix="r_mx_e2",
            start=1001,
            length=100,
            ops=[Op(kind="snv", pos=1048, seq="A")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="1048",
            gene="geneC",
            base_changes="A",
            aa_changes="Ala50Thr",
            variant_type="SNP",
            change_type="Non-synonymous",
            reference_codon="GCT",
            snp_codon="ACT",
            snp_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
    gff_content=GFF_CDS_MULTIEXON,
    gff_features="CDS",
)


# 14. SNV en codon que cruza exon junction
# Codon 34 = pos 900 (exon 1) + pos 1001-1002 (exon 2). Cambio pos 900 G>A:
# mRNA codon 34 GCT -> ACT = Thr. Ala34Thr.
# Read needs to span both exons via skip CIGAR (RNA-seq style).
scenario_multiexon_junction_snp = Scenario(
    name="14_multiexon_junction_snp",
    description="SNP en codon junction de geneC: pos 900 G>A (base 1 codon 34 GCT->ACT, Ala34Thr); reads N-CIGAR",
    variants=[VcfRecord(pos=900, ref="G", alt="A")],
    reads=[
        ReadGroup(
            name_prefix="r_mx_junc",
            start=851,
            length=200,  # 50 bp exon1 + 100 bp intron + 50 bp exon2
            ops=[
                Op(kind="snv", pos=900, seq="A"),
                Op(kind="skip", pos=901, length=100),
            ],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="900",
            gene="geneC",
            base_changes="A",
            aa_changes="Ala34Thr",
            variant_type="SNP",
            change_type="Non-synonymous",
            reference_codon="GCT",
            snp_codon="ACT",
            snp_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
    gff_content=GFF_CDS_MULTIEXON,
    gff_features="CDS",
)


# 15. MNV cruzando exon junction
# Codon 34 con cambios en pos 900 G>A (exon 1) + pos 1002 T>C (exon 2):
# mRNA codon 34: GCT -> ACC = Thr. Posiciones MUY separadas en genomico
# (902 bp apart) pero adyacentes en spliced CDS.
scenario_multiexon_junction_mnv = Scenario(
    name="15_multiexon_junction_mnv",
    description="MNV cruzando junction: pos 900 G>A + pos 1002 T>C en cis (mRNA codon 34 GCT->ACC, Ala34Thr)",
    variants=[
        VcfRecord(pos=900, ref="G", alt="A"),
        VcfRecord(pos=1002, ref="T", alt="C"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_mx_mnv",
            start=851,
            length=200,
            ops=[
                Op(kind="snv", pos=900, seq="A"),
                Op(kind="skip", pos=901, length=100),
                Op(kind="snv", pos=1002, seq="C"),
            ],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="900, 1002",
            gene="geneC",
            variant_type="SNP/MNV",
            aa_changes="Ala34Thr",
            reference_codon="GCT",
            mnv_codon="ACC",
            change_type="Non-synonymous",
            mnv_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
    gff_content=GFF_CDS_MULTIEXON,
    gff_features="CDS",
)


# ---------------------------------------------------------------------------
# EDGE CASES
# ---------------------------------------------------------------------------

# 16. Variant en VCF pero 0 reads en BAM en esa posicion
# Variant en pos 28 (codon 10 geneA), pero todas las lecturas estan en pos 500-600
# (region intergenica). Verifica como reporta get_mnv un evento sin soporte BAM.
scenario_no_bam_coverage = Scenario(
    name="16_no_bam_coverage",
    description="VCF con SNV en pos 28 (geneA) pero 0 reads en BAM en esa zona",
    variants=[VcfRecord(pos=28, ref="G", alt="A")],
    reads=[
        # 20 reads en pos 500-600 (lejos de la variant)
        ReadGroup(
            name_prefix="r_far",
            start=500,
            length=100,
            ops=[],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28",
            gene="geneA",
            variant_type="SNP",
            aa_changes="Ala10Thr",
            snp_reads="0",
            total_reads="0",
            snp_frequencies="0.0000",
        ),
    ],
    expected_row_count=1,
)


# 17. Filtro --min-snp-frequency descarta variant de baja frecuencia
# 2 reads con la SNV + 18 reads ref -> freq 10%. Con --min-snp-frequency 0.5
# debe ser descartada del output.
scenario_min_snp_frequency_filter = Scenario(
    name="17_min_snp_frequency_filter",
    description="SNV con freq 10% (2/20 reads) filtrada por --min-snp-frequency 0.5 -> 0 filas",
    variants=[VcfRecord(pos=28, ref="G", alt="A")],
    reads=[
        ReadGroup(
            name_prefix="r_alt",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="A")],
            count=2,
        ),
        ReadGroup(
            name_prefix="r_ref",
            start=1,
            length=100,
            ops=[],
            count=18,
        ),
    ],
    expected=[],   # esperamos NINGUNA fila por el filtro
    expected_row_count=0,
    extra_cli_args=["--min-snp-frequency", "0.5"],
)


# ---------------------------------------------------------------------------
# AA EFFECTS EDGE CASES
# ---------------------------------------------------------------------------

# 18. Stop-gained via MNV (geneA codon 50 GCT -> TAA stop)
# GCT no llega a stop con 2 SNVs; necesita 3: pos148 G>T + pos149 C>A + pos150 T>A -> TAA
scenario_stop_gained = Scenario(
    name="18_stop_gained_via_mnv",
    description="Stop-gained via MNV de 3 SNVs: codon 50 GCT->TAA (Ala50Ter)",
    variants=[
        VcfRecord(pos=148, ref="G", alt="T"),
        VcfRecord(pos=149, ref="C", alt="A"),
        VcfRecord(pos=150, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_stop",
            start=100,
            length=100,
            ops=[
                Op(kind="snv", pos=148, seq="T"),
                Op(kind="snv", pos=149, seq="A"),
                Op(kind="snv", pos=150, seq="A"),
            ],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="148, 149, 150",
            gene="geneA",
            variant_type="SNP/MNV",
            reference_codon="GCT",
            mnv_codon="TAA",
            change_type="Stop gained",
            mnv_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
)


# 19. Start codon alterado: ATG -> ACG (Met1Thr)
# NOTA: get_mnv reporta como "Non-synonymous" estandar; no tiene categoria
# especifica "Start lost" para SNVs en el codon inicial. La AA si es correcta.
scenario_start_lost = Scenario(
    name="19_start_codon_altered",
    description="SNV en el codon iniciador ATG: pos 2 T>C -> ACG (Met1Thr) -> change_type='Start lost'",
    variants=[VcfRecord(pos=2, ref="T", alt="C")],
    reads=[
        ReadGroup(
            name_prefix="r_start",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=2, seq="C")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="2",
            gene="geneA",
            variant_type="SNP",
            reference_codon="ATG",
            snp_codon="ACG",
            aa_changes="Met1Thr",
            change_type="Start lost",
            snp_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
)


# 20. Stop-lost: codon 100 TAA -> CAA (pos 298 T>C, Gln en vez de stop)
scenario_stop_lost = Scenario(
    name="20_stop_lost",
    description="Stop-lost: codon 100 TAA -> CAA (pos 298 T>C, Ter100Gln)",
    variants=[VcfRecord(pos=298, ref="T", alt="C")],
    reads=[
        ReadGroup(
            name_prefix="r_stoplost",
            start=200,
            length=100,
            ops=[Op(kind="snv", pos=298, seq="C")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="298",
            gene="geneA",
            variant_type="SNP",
            reference_codon="TAA",
            snp_codon="CAA",
            change_type="Stop lost",
            snp_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
)


# ---------------------------------------------------------------------------
# ALELOS COMPLEJOS
# ---------------------------------------------------------------------------

# 21. Variant en intron (geneC pos 950, dentro del intron 901-1000)
# Sin CDS feature en el intron -> intergenic con --gff-features CDS
scenario_intron_variant = Scenario(
    name="21_intron_variant",
    description=(
        "SNV en el intron de geneC (pos 950 T>A), lejos de los sitios de splice: "
        "se reporta contra geneC como intron_variant, no como intergenic"
    ),
    variants=[VcfRecord(pos=950, ref="T", alt="A")],
    reads=[
        ReadGroup(
            name_prefix="r_intron",
            start=900,
            length=100,
            ops=[Op(kind="snv", pos=950, seq="A")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="950",
            gene="geneC",
            variant_type="SNP",
            change_type="Unknown",
        ),
    ],
    expected_row_count=1,
    gff_content=GFF_CDS_MULTIEXON,
    gff_features="CDS",
)


# 22. Multiallelic SNV con --split-multiallelic
# Con el fix de merge_snp_into_groups, get_mnv emite UNA fila por cada alt
# alternativo en la misma posicion. Las 2 SNV (G>A y G>T en pos 28) producen
# 2 filas independientes con AA y soporte BAM cada una.
scenario_multiallelic = Scenario(
    name="22_multiallelic_split",
    description="Multiallelic SNV pos 28 G>A,T con --split-multiallelic: 2 filas SNV independientes (G>A Ala10Thr + G>T Ala10Ser)",
    # El genotipo lleva los dos alelos: 1/2. Con el 1/1 por defecto la muestra
    # solo llevaria el primero y el segundo no se anotaria, que es lo correcto
    # para ese genotipo pero no lo que este escenario quiere medir.
    variants=[VcfRecord(pos=28, ref="G", alt="A,T", genotype="1/2")],
    reads=[
        ReadGroup(
            name_prefix="r_alt_a",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="A")],
            count=10,
        ),
        ReadGroup(
            name_prefix="r_alt_t",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="T")],
            count=10,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28",
            gene="geneA",
            base_changes="A",
            variant_type="SNP",
            snp_codon="ACT",
            aa_changes="Ala10Thr",
            snp_reads="10",
            total_reads="20",
            snp_frequencies="0.5000",
        ),
        ExpectedRow(
            positions="28",
            gene="geneA",
            base_changes="T",
            variant_type="SNP",
            snp_codon="TCT",
            aa_changes="Ala10Ser",
            snp_reads="10",
            total_reads="20",
            snp_frequencies="0.5000",
        ),
    ],
    expected_row_count=2,
    extra_cli_args=["--split-multiallelic"],
)


# 23. Delins: REF=GCT ALT=GA en pos 28
# - anchor base G en pos 28 conservada
# - C en pos 29 sustituida por A
# - T en pos 30 borrada
# Net: -1 bp + substitucion = frameshift
scenario_delins = Scenario(
    name="23_delins_subst_plus_del",
    description="Delins REF=GCT ALT=GA en pos 28 (subst pos 29 C>A + del pos 30): frameshift",
    variants=[VcfRecord(pos=28, ref="GCT", alt="GA")],
    reads=[
        ReadGroup(
            name_prefix="r_delins",
            start=1,
            length=150,
            ops=[
                Op(kind="snv", pos=29, seq="A"),
                Op(kind="del", pos=30, length=1),
            ],
            count=20,
        ),
    ],
    expected=[
        # get_mnv decompone el delins en componentes; comprobamos que sale al menos un INDEL frameshift
        ExpectedRow(
            positions="28",
            gene="geneA",
            variant_type="INDEL",
            change_type="Frameshift Indel",
        ),
    ],
)


# 24. Inserción in-frame grande (9 bp = 3 codones)
# Insertar GCTGCTGCT después de pos 30. Inserta 3 nuevos Ala in-frame.
scenario_large_ins = Scenario(
    name="24_large_inframe_insertion",
    description="Insercion in-frame de 9 bp (GCTGCTGCT = 3 Ala) entre pos 30 y 31",
    variants=[VcfRecord(pos=30, ref="T", alt="TGCTGCTGCT")],
    reads=[
        ReadGroup(
            name_prefix="r_largein",
            start=1,
            length=141,
            ops=[Op(kind="ins", pos=30, seq="GCTGCTGCT")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            gene="geneA",
            variant_type="INDEL",
            change_type="In-frame Indel",
            event_class="insertion",
            event_components="INS:30:+GCTGCTGCT",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=1,
)


# 25. Deleción in-frame grande (6 bp = 2 codones)
# Borrar pos 31-36 (codones 11+12 GCT GCT). VCF: pos=30 REF=TGCTGCT ALT=T
scenario_large_del = Scenario(
    name="25_large_inframe_deletion",
    description="Delecion in-frame de 6 bp (2 Ala consecutivos) en pos 31-36",
    variants=[VcfRecord(pos=30, ref="TGCTGCT", alt="T")],
    reads=[
        ReadGroup(
            name_prefix="r_largedel",
            start=1,
            length=156,
            ops=[Op(kind="del", pos=31, length=6)],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            gene="geneA",
            variant_type="INDEL",
            change_type="In-frame Indel",
            event_class="deletion",
            event_components="DEL:31-36:GCTGCT",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=1,
)


# ---------------------------------------------------------------------------
# MULTI-CONTIG + iVar TSV INPUT
# ---------------------------------------------------------------------------

# 26. Variant en chr_test2 (multi-contig): SNV en geneD codon 10
scenario_multicontig = Scenario(
    name="26_multicontig",
    description="Multi-contig: SNV en chr_test2 geneD codon 10 (GCT->ACT, Ala10Thr) + SNV en chr_test geneA",
    variants=[
        VcfRecord(pos=28, ref="G", alt="A", chrom="chr_test"),
        VcfRecord(pos=28, ref="G", alt="A", chrom="chr_test2"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_c1",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="A")],
            count=20,
        ),
        ReadGroup(
            name_prefix="r_c2",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="A")],
            count=20,
            chrom="chr_test2",
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28",
            gene="geneA",
            aa_changes="Ala10Thr",
            variant_type="SNP",
        ),
        ExpectedRow(
            positions="28",
            gene="geneD",
            aa_changes="Ala10Thr",
            variant_type="SNP",
        ),
    ],
    expected_row_count=2,
)


# 27. iVar TSV input - SNV simple
scenario_ivar_snv = Scenario(
    name="27_ivar_tsv_snv",
    description="iVar TSV input: SNV pos 28 G>A (Ala10Thr)",
    variants=[],   # se usa ivar_records
    reads=[
        ReadGroup(
            name_prefix="r_ivar",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="A")],
            count=20,
        ),
    ],
    ivar_records=[
        IvarRecord(pos=28, ref="G", alt="A", total_dp=20, alt_dp=20, alt_freq=1.0),
    ],
    expected=[
        ExpectedRow(
            positions="28",
            gene="geneA",
            base_changes="A",
            aa_changes="Ala10Thr",
            variant_type="SNP",
            snp_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
)


# 28. iVar TSV input - insertion (notacion +SEQ)
scenario_ivar_insertion = Scenario(
    name="28_ivar_tsv_insertion",
    description="iVar TSV: insercion in-frame +GCT en pos 30 (notacion +SEQ)",
    variants=[],
    reads=[
        ReadGroup(
            name_prefix="r_ivarin",
            start=1,
            length=147,
            ops=[Op(kind="ins", pos=30, seq="GCT")],
            count=20,
        ),
    ],
    ivar_records=[
        IvarRecord(pos=30, ref="T", alt="+GCT", total_dp=20, alt_dp=20, alt_freq=1.0),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            gene="geneA",
            variant_type="INDEL",
            change_type="In-frame Indel",
            event_class="insertion",
            event_components="INS:30:+GCT",
        ),
    ],
    expected_row_count=1,
)


# 29. iVar TSV input - deletion (notacion -SEQ)
# iVar ancla la delecion en la base ANTERIOR al hueco: POS=30 (ancla, T),
# ALT=-G borra la G en pos 31 (primera base del codon 11). Convierte a la
# variante VCF (30, TG, T) -> DEL:31:G, frameshift desde el codon 11.
scenario_ivar_deletion = Scenario(
    name="29_ivar_tsv_deletion",
    description="iVar TSV: delecion frameshift -G (ancla POS=30, borra pos 31; notacion -SEQ)",
    variants=[],
    reads=[
        ReadGroup(
            name_prefix="r_ivardel",
            start=1,
            length=151,
            ops=[Op(kind="del", pos=31, length=1)],
            count=20,
        ),
    ],
    ivar_records=[
        IvarRecord(pos=30, ref="T", alt="-G", total_dp=20, alt_dp=20, alt_freq=1.0),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            gene="geneA",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            event_class="deletion",
            event_components="DEL:31:G",
        ),
    ],
    expected_row_count=1,
)


# ---------------------------------------------------------------------------
# REFINAMIENTOS DE INDEL (rama indels): stop por indel in-frame + gate de freq
# ---------------------------------------------------------------------------

# 30. Stop GANADO por insercion in-frame.
# geneA codon 10 = pos 28-30 (GCT). Insertamos TAA justo despues de pos 30
# (entre codon 10 y codon 11). La insercion es de 3 bp -> in-frame, y forma
# un codon de stop completo nuevo. indel_stop_effect cuenta los '*': la
# proteina ref tiene 1 stop (terminal), la alt tiene 2 -> "Stop gained".
# (Sin este refinamiento se etiquetaria como el generico "In-frame Indel".)
scenario_stop_gained_inframe_ins = Scenario(
    name="30_stop_gained_inframe_ins",
    description="Insercion in-frame de TAA tras pos 30 crea stop prematuro -> Change Type 'Stop gained' (no 'In-frame Indel')",
    variants=[VcfRecord(pos=30, ref="T", alt="TTAA")],
    reads=[
        ReadGroup(
            name_prefix="r_stopgain_ins",
            start=1,
            length=147,
            ops=[Op(kind="ins", pos=30, seq="TAA")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            gene="geneA",
            variant_type="INDEL",
            event_class="insertion",
            event_components="INS:30:+TAA",
            change_type="Stop gained",
            aa_changes="Ala10_Ala11ins*",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=1,
)


# 31. Stop PERDIDO por delecion in-frame.
# geneA stop = TAA en pos 298-300. Borramos esas 3 bases (in-frame).
# VCF: anchor en pos 297 (T, ultima base del codon 99 GCT@295-297),
# REF=TTAA (pos 297-300), ALT=T -> borra TAA@298-300.
# La proteina ref tiene 1 stop, la alt tiene 0 -> "Stop lost".
scenario_stop_lost_inframe_del = Scenario(
    name="31_stop_lost_inframe_del",
    description="Delecion in-frame de 3bp que elimina el stop TAA (298-300) -> Change Type 'Stop lost'",
    variants=[VcfRecord(pos=297, ref="TTAA", alt="T")],
    reads=[
        ReadGroup(
            name_prefix="r_stoplost_del",
            start=200,
            length=104,  # cubre 200-303; D en 298-300 flanqueada por M
            ops=[Op(kind="del", pos=298, length=3)],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="297",
            gene="geneA",
            variant_type="INDEL",
            event_class="deletion",
            event_components="DEL:298-300:TAA",
            change_type="Stop lost",
            aa_changes="*100del",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=1,
)


# 32. Gate de frameshift --frameshift-min-freq 0.0: COMPORTAMIENTO OPT-IN.
# Misma estructura que el escenario 07: del frameshift de 1bp en pos 31 +
# SNV downstream pos 39. El VCF declara AF=0.20 para la del y AF=0.95 para
# el SNV. Con --frameshift-min-freq 0.0 (opt-in) cualquier indel upstream
# propaga el frameshift -> el SNV downstream lleva marcador "(fs)". El
# default ahora es 0.5, que SUPRIME esta propagacion (ver escenario 33).
scenario_fs_gate_default = Scenario(
    name="32_fs_gate_zero_propagates",
    description="Del fs upstream AF=0.20 + SNV downstream con --frameshift-min-freq 0.0 (opt-in): el frameshift se propaga desde cualquier indel -> SNV con '(fs)'. El default ahora es 0.5 (ver 33)",
    variants=[
        VcfRecord(pos=30, ref="TG", alt="T", af=0.20),
        VcfRecord(pos=39, ref="T", alt="A", af=0.95),
    ],
    reads=[
        # La del esta en 4 de 20 lecturas: un 20% real, coherente con el AF=0.20
        # que declara el VCF. Antes las 20 lecturas la llevaban, asi que el BAM
        # decia 100% mientras el VCF decia 20% y el escenario se contradecia.
        # Las 4 llevan tambien el SNV: con 4 lecturas cis de 20 informativas la
        # regla de trans (cis <= 10%) no se dispara, asi que lo que se prueba
        # aqui sigue siendo la puerta de frecuencia y no la evidencia de trans.
        ReadGroup(
            name_prefix="r_gate_def_with_del",
            start=1,
            length=151,
            ops=[
                Op(kind="del", pos=31, length=1),
                Op(kind="snv", pos=39, seq="A"),
            ],
            count=4,
        ),
        ReadGroup(
            name_prefix="r_gate_def_no_del",
            start=1,
            length=151,
            ops=[
                Op(kind="snv", pos=39, seq="A"),
            ],
            count=16,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            event_class="deletion",
            event_components="DEL:31:G",
        ),
        ExpectedRow(
            positions="39",
            variant_type="SNP",
            change_type="Synonymous (frameshift)",
            aa_changes="Ala13Ala (fs)",
            reference_codon="GCT",
            snp_codon="GCA",
        ),
    ],
    expected_row_count=2,
    extra_cli_args=["--frameshift-min-freq", "0.0"],
)


# 33. Gate de frameshift --frameshift-min-freq 0.5: SUPRIME la propagacion.
# Mismos inputs que el escenario 32. La del upstream tiene AF=0.20 < 0.5,
# asi que NO supera el gate y NO propaga el frameshift: el SNV downstream
# (AF=0.95, casi seguro en otra molecula) se anota como sinonimo normal,
# SIN marcador "(fs)". La del sigue siendo "Frameshift Indel" por si misma.
scenario_fs_gate_suppressed = Scenario(
    name="33_fs_gate_high_freq_suppressed",
    description="Mismos inputs que 32 pero con --frameshift-min-freq 0.5: la del AF=0.20 no supera el gate -> SNV downstream SIN '(fs)'",
    variants=[
        VcfRecord(pos=30, ref="TG", alt="T", af=0.20),
        VcfRecord(pos=39, ref="T", alt="A", af=0.95),
    ],
    reads=[
        # La del esta en 4 de 20 lecturas: un 20% real, coherente con el AF=0.20
        # que declara el VCF. Antes las 20 lecturas la llevaban, asi que el BAM
        # decia 100% mientras el VCF decia 20% y el escenario se contradecia.
        # Las 4 llevan tambien el SNV: con 4 lecturas cis de 20 informativas la
        # regla de trans (cis <= 10%) no se dispara, asi que lo que se prueba
        # aqui sigue siendo la puerta de frecuencia y no la evidencia de trans.
        ReadGroup(
            name_prefix="r_gate_sup_with_del",
            start=1,
            length=151,
            ops=[
                Op(kind="del", pos=31, length=1),
                Op(kind="snv", pos=39, seq="A"),
            ],
            count=4,
        ),
        ReadGroup(
            name_prefix="r_gate_sup_no_del",
            start=1,
            length=151,
            ops=[
                Op(kind="snv", pos=39, seq="A"),
            ],
            count=16,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            event_class="deletion",
            event_components="DEL:31:G",
        ),
        ExpectedRow(
            positions="39",
            variant_type="SNP",
            change_type="Synonymous",
            aa_changes="Ala13Ala",
            reference_codon="GCT",
            snp_codon="GCA",
        ),
    ],
    expected_row_count=2,
    extra_cli_args=["--frameshift-min-freq", "0.5"],
)



# ---------------------------------------------------------------------------
# 58. Gate de frameshift SIN AF declarada en el VCF: manda la frecuencia del BAM.
# ---------------------------------------------------------------------------
# Muchos llamadores no escriben AF (el propio VCF de ejemplo del repo, estilo
# VarScan, no lo hace). El gate solo miraba la frecuencia declarada, asi que con
# un VCF sin AF pasaba cualquier indel a cualquier umbral y --frameshift-min-freq
# quedaba inerte justo donde habia evidencia de lecturas. Mismos datos que el 33
# pero sin AF: la del esta en 4 de 20 lecturas (0.20 < 0.5), y el gate debe
# suprimir la propagacion leyendo eso del BAM.
scenario_fs_gate_no_declared_af = Scenario(
    name="58_fs_gate_reads_decide_without_af",
    description="Del fs upstream en 4/20 lecturas y VCF SIN AF, con --frameshift-min-freq 0.5: el gate usa la frecuencia del BAM (0.20) -> SNV downstream SIN '(fs)'",
    variants=[
        VcfRecord(pos=30, ref="TG", alt="T"),
        VcfRecord(pos=39, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_gate_noaf_with_del",
            start=1,
            length=151,
            ops=[
                Op(kind="del", pos=31, length=1),
                Op(kind="snv", pos=39, seq="A"),
            ],
            count=4,
        ),
        ReadGroup(
            name_prefix="r_gate_noaf_no_del",
            start=1,
            length=151,
            ops=[
                Op(kind="snv", pos=39, seq="A"),
            ],
            count=16,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            event_class="deletion",
            event_components="DEL:31:G",
        ),
        ExpectedRow(
            positions="39",
            variant_type="SNP",
            change_type="Synonymous",
            aa_changes="Ala13Ala",
            reference_codon="GCT",
            snp_codon="GCA",
        ),
    ],
    expected_row_count=2,
    extra_cli_args=["--frameshift-min-freq", "0.5"],
)


# ---------------------------------------------------------------------------
# 59. Un indel FUERA de todo feature conserva su soporte de lecturas.
# ---------------------------------------------------------------------------
# El contador intergenico solo trataba sustituciones de una posicion, asi que un
# indel fuera de cualquier gen no llegaba a ningun contador y la fila salia
# afirmando "Event Reads = 0" a "Event Depth = 0" para un alelo que llevaban
# todas las lecturas. Misma insercion que el escenario 04 y las mismas lecturas;
# lo unico que cambia es que el gen esta lejos, asi que el locus queda
# intergenico. El soporte tiene que ser identico: 20 de 20.
GFF_GENE_FAR_AWAY = (
    "##gff-version 3\n"
    "chr_test\tsynth\tgene\t400\t700\t.\t+\t.\tID=gene-geneZ;Name=geneZ\n"
)

scenario_intergenic_indel_support = Scenario(
    name="59_intergenic_indel_keeps_read_support",
    description="Insercion fuera de todo gen con 20 lecturas que la llevan: el soporte de evento se cuenta igual que dentro de un gen (20/20), no se reporta 0 sobre profundidad 0",
    variants=[
        VcfRecord(pos=30, ref="T", alt="TGCT"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_interg_ins",
            start=1,
            length=151,
            ops=[Op(kind="ins", pos=30, seq="GCT")],
            count=20,
        ),
    ],
    gff_content=GFF_GENE_FAR_AWAY,
    expected=[
        ExpectedRow(
            positions="30",
            gene="intergenic",
            variant_type="INDEL",
            event_class="insertion",
            event_components="INS:30:+GCT",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=1,
)


# ---------------------------------------------------------------------------
# 60. Lecturas que no deben contarse: duplicado, secundaria, suplementaria y QC fail.
# ---------------------------------------------------------------------------
# 10 lecturas limpias llevan el ALT, y otras 40 lo llevan tambien pero marcadas
# con uno de los cuatro flags que un contador debe ignorar. El soporte tiene que
# ser 10, no 50. QC fail (0x200) era el unico de los cuatro que no se miraba, asi
# que get_MNV contaba lecturas que su propio oraculo, bcftools mpileup, descarta.
scenario_uncountable_flags = Scenario(
    name="60_flagged_reads_are_not_counted",
    description="Lecturas duplicadas, secundarias, suplementarias y QC-fail que llevan el ALT no cuentan como soporte: 10 limpias -> 10, no 50",
    variants=[VcfRecord(pos=28, ref="G", alt="T")],
    reads=[
        ReadGroup(name_prefix="r_clean", start=1, length=100,
                  ops=[Op(kind="snv", pos=28, seq="T")], count=10),
        ReadGroup(name_prefix="r_dup", start=1, length=100,
                  ops=[Op(kind="snv", pos=28, seq="T")], count=10, extra_flags=0x400),
        ReadGroup(name_prefix="r_secondary", start=1, length=100,
                  ops=[Op(kind="snv", pos=28, seq="T")], count=10, extra_flags=0x100),
        ReadGroup(name_prefix="r_supplementary", start=1, length=100,
                  ops=[Op(kind="snv", pos=28, seq="T")], count=10, extra_flags=0x800),
        ReadGroup(name_prefix="r_qcfail", start=1, length=100,
                  ops=[Op(kind="snv", pos=28, seq="T")], count=10, extra_flags=0x200),
    ],
    expected=[
        ExpectedRow(
            positions="28",
            variant_type="SNP",
            snp_reads="10",
            total_reads="10",
        ),
    ],
    expected_row_count=1,
)



# ---------------------------------------------------------------------------
# 61. Insercion escrita contra la base SIGUIENTE, con soporte de lecturas.
# ---------------------------------------------------------------------------
# VCF permite `31 G>AAAG`, que inserta AAA entre la 30 y la 31, igual que
# `30 T>TAAA`. El ancla cae entonces antes del POS del registro, y todo lo que
# observa un alelo en get_MNV lee hacia delante desde POS: el contador exacto y
# el descubrimiento de haplotipos no podian ver esa union, asi que la fila salia
# con 0 lecturas sobre 20 de profundidad. Se re-ancla el registro al entrar. Se
# inserta AAA, y no una copia del repeat GCT, para que la forma derecha no sea
# ambigua.
scenario_right_anchored_insertion = Scenario(
    name="61_right_anchored_insertion_is_counted",
    description="Insercion escrita contra la base siguiente (31 G>AAAG): se re-ancla a 30 T>TAAA y sus 20 lecturas se cuentan",
    variants=[VcfRecord(pos=31, ref="G", alt="AAAG")],
    reads=[
        ReadGroup(
            name_prefix="r_right_ins",
            start=1,
            length=151,
            ops=[Op(kind="ins", pos=30, seq="AAA")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            reference_bases="T",
            base_changes="TAAA",
            variant_type="INDEL",
            event_class="insertion",
            event_components="INS:30:+AAA",
            event_reads="20",
        ),
    ],
    expected_row_count=1,
)


# ---------------------------------------------------------------------------
# 62. Sustitucion multibase FUERA de todo gen: tambien se cuenta.
# ---------------------------------------------------------------------------
# El contador intergenico solo trataba filas SNP e INDEL, asi que una fila MNV
# fuera de un CDS no llegaba a ningun contador: salia sin recuentos, y un umbral
# leia esa ausencia como cero y borraba la fila. Aqui las 20 lecturas llevan las
# dos bases, asi que la fila tiene que reportar 20 y sobrevivir a --mnv 5.
GFF_GENE_ELSEWHERE = (
    "##gff-version 3\n"
    "chr_test\tsynth\tgene\t400\t700\t.\t+\t.\tID=gene-geneZ;Name=geneZ\n"
)

scenario_intergenic_mnv_counted = Scenario(
    name="62_intergenic_mnv_is_counted",
    description="Sustitucion multibase fuera de todo gen con 20 lecturas que la llevan entera: se cuenta (MNV Reads 20) y sobrevive a --mnv 5",
    variants=[VcfRecord(pos=28, ref="GC", alt="TA")],
    reads=[
        ReadGroup(
            name_prefix="r_interg_mnv",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="T"), Op(kind="snv", pos=29, seq="A")],
            count=20,
        ),
    ],
    gff_content=GFF_GENE_ELSEWHERE,
    expected=[
        ExpectedRow(
            positions="28, 29",
            gene="intergenic",
            variant_type="MNV",
            mnv_reads="20",
            total_reads="20",
        ),
    ],
    expected_row_count=1,
    extra_cli_args=["--mnv", "5"],
)


# ---------------------------------------------------------------------------
# 63. Union por deslizamiento ribosomal: la base compartida pertenece a 2 codones.
# ---------------------------------------------------------------------------
# En un join(a..b, b..c) (la forma de ORF1ab de SARS-CoV-2, que el modelo admite)
# dos filas CDS comparten una base y el ribosoma la lee dos veces: ocupa dos
# posiciones del transcrito. Se tomaba solo la primera, asi que el segundo codon
# se quedaba con la base de referencia y su consecuencia no se reportaba. Aqui la
# segunda es una perdida de codon de parada, que es justo lo que se perdia.
GFF_SLIPPAGE_JOIN = (
    "##gff-version 3\n"
    "chr_test\tsynth\tCDS\t1\t30\t.\t+\t0\tID=cds-a;Parent=t1;Name=orfSlip\n"
    "chr_test\tsynth\tCDS\t30\t60\t.\t+\t0\tID=cds-b;Parent=t1;Name=orfSlip\n"
)

scenario_slippage_shared_base = Scenario(
    name="63_slippage_shared_base_hits_both_codons",
    description="CDS join(1..30, 30..60): una sustitucion en la base compartida se anota en los dos codones que ocupa, no solo en el primero",
    variants=[VcfRecord(pos=30, ref="T", alt="A")],
    reads=[],
    gff_content=GFF_SLIPPAGE_JOIN,
    gff_features="CDS",
    expected=[
        ExpectedRow(positions="30", reference_codon="GCT", snp_codon="GCA"),
    ],
    expected_row_count=2,
)

# ---------------------------------------------------------------------------
# 34. Eucariotas: GFF con features CDS pero SIN --gff-features -> auto-CDS.
# ---------------------------------------------------------------------------
# Mismo SNV que el escenario 13 (pos 1048 G>A, geneC codon 50 -> Ala50Thr) pero
# sin pasar --gff-features. Por defecto se usaria gene,pseudogene y geneC se
# anotaria sobre el span gen 801-1200 (con intron) dando un resultado erroneo;
# ahora get_mnv detecta que el GFF tiene CDS y auto-selecciona el modelo CDS
# (phase/splice-aware), produciendo la anotacion multi-exon correcta.
scenario_eukaryote_autocds = Scenario(
    name="34_eukaryote_autocds_default",
    description="GFF con CDS pero sin --gff-features: auto-CDS -> anotacion multi-exon correcta (pos 1048 Ala50Thr en geneC)",
    variants=[VcfRecord(pos=1048, ref="G", alt="A")],
    reads=[
        ReadGroup(
            name_prefix="r_autocds",
            start=1001,
            length=100,
            ops=[Op(kind="snv", pos=1048, seq="A")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="1048",
            gene="geneC",
            aa_changes="Ala50Thr",
            variant_type="SNP",
            change_type="Non-synonymous",
            reference_codon="GCT",
            snp_codon="ACT",
        ),
    ],
    expected_row_count=1,
    gff_content=GFF_CDS_MULTIEXON,
    # gff_features intencionadamente SIN definir -> ejercita el auto-CDS.
)


# ---------------------------------------------------------------------------
# 35. Frameshift mas alla del stop prematuro (geneE en chr_test2).
# ---------------------------------------------------------------------------
# geneE = ATG GTA (AAA)x4 TAA en pos 301-321. Borrar la G en pos 304 desplaza
# el marco: el codon 2 pasa a leer TAA -> stop prematuro en la posicion 2 de la
# proteina. Una SNV downstream en pos 313 (codon 5) queda DESPUES del stop, asi
# que no se traduce: get_mnv debe marcarla "downstream of premature stop" en vez
# de "(fs)". Distancia del>SNV = 9 bp (>3) -> filas separadas, no complex_indel.
scenario_frameshift_past_stop = Scenario(
    name="35_frameshift_past_premature_stop",
    description="geneE: del frameshift pos 304 crea stop prematuro (codon 2); la SNV downstream pos 313 (codon 5) se marca 'downstream of premature stop', no '(fs)'",
    variants=[
        VcfRecord(pos=303, ref="GG", alt="G", chrom="chr_test2"),
        VcfRecord(pos=313, ref="A", alt="C", chrom="chr_test2"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_fs_ptc",
            start=301,
            length=60,
            ops=[
                Op(kind="del", pos=304, length=1),
                Op(kind="snv", pos=313, seq="C"),
            ],
            count=20,
            chrom="chr_test2",
        ),
    ],
    expected=[
        ExpectedRow(
            positions="303",
            gene="geneE",
            variant_type="INDEL",
            change_type="Frameshift Indel",
            event_class="deletion",
        ),
        ExpectedRow(
            positions="313",
            gene="geneE",
            variant_type="SNP",
            change_type="Downstream of premature stop",
            aa_changes="downstream of premature stop",
        ),
    ],
    expected_row_count=2,
)


# ---------------------------------------------------------------------------
# 36. Dos haplotipos locales reales conviviendo (control positivo)
# ---------------------------------------------------------------------------
# Mismos tres variantes que el escenario 08, pero repartidos en DOS especies
# moleculares distintas y ninguna molecula lleva las tres:
#   12 lecturas: SNV:28 + INS:29        (referencia en pos 30)
#    8 lecturas: INS:29 + SNV:30        (referencia en pos 28)
# Los dos haplotipos existen y deben salir los dos, con sus recuentos propios
# (12/20 y 8/20). El triple NO debe salir: ninguna molecula lo lleva. Es el
# caso que la enumeracion de subconjuntos no podia distinguir del escenario 08,
# porque alli emitia los mismos dos subconjuntos a partir de una sola especie.
scenario_two_real_haplotypes = Scenario(
    name="36_two_coexisting_local_haplotypes",
    description="Dos haplotipos locales reales (SNV:28+INS:29 en 12 lecturas, INS:29+SNV:30 en 8): salen ambos con sus frecuencias y el triple no sale",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=29, ref="C", alt="CGCT"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_hap_left",
            start=1,
            length=147,
            ops=[
                Op(kind="snv", pos=28, seq="T"),
                Op(kind="ins", pos=29, seq="GCT"),
            ],
            count=12,
        ),
        ReadGroup(
            name_prefix="r_hap_right",
            start=1,
            length=147,
            ops=[
                Op(kind="ins", pos=29, seq="GCT"),
                Op(kind="snv", pos=30, seq="A"),
            ],
            count=8,
        ),
    ],
    expected=[
        ExpectedRow(
            event_class="complex_indel",
            event_components="SNV:28:G>T | INS:29:+GCT",
            variant_type="INDEL",
            event_reads="12",
            event_frequency="0.6000",
        ),
        ExpectedRow(
            event_class="complex_indel",
            event_components="INS:29:+GCT | SNV:30:T>A",
            variant_type="INDEL",
            event_reads="8",
            event_frequency="0.4000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 37. Trans probado por las lecturas: se publica la evidencia
# ---------------------------------------------------------------------------
# Delecion frameshift de 1 bp en pos 4 (codon 2) y SNV en pos 28 (codon 10), en
# MOLECULAS DISTINTAS: 10 lecturas llevan solo la delecion y 10 solo la SNV.
# El frameshift no debe propagarse al codon 10, porque su molecula no lleva la
# delecion. Sin AF en el VCF el gate de frecuencia deja pasar la delecion, asi
# que la unica razon para no propagarlo son las lecturas, y la columna
# Frameshift Phasing lo dice: trans, indel en pos 3, 0 de 10 lecturas en cis.
scenario_trans_evidence_reported = Scenario(
    name="37_frameshift_trans_evidence",
    description="Del frameshift pos 4 y SNV pos 28 en moleculas distintas: no se propaga el frameshift y la columna publica trans:3:0/10",
    variants=[
        VcfRecord(pos=3, ref="GG", alt="G"),
        VcfRecord(pos=28, ref="G", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_del_only",
            start=1,
            length=100,
            ops=[Op(kind="del", pos=4, length=1)],
            count=10,
        ),
        ReadGroup(
            name_prefix="r_snv_only",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="A")],
            count=10,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28",
            gene="geneA",
            variant_type="SNP",
            aa_changes="Ala10Thr",
            change_type="Non-synonymous",
            frameshift_phasing="trans:3:0/10",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 38. Un haplotipo de una sola lectura no se emite por defecto
# ---------------------------------------------------------------------------
# 19 lecturas llevan solo la insercion y 1 lleva insercion + SNV. Esa unica
# lectura es indistinguible de un error de secuenciacion, y desde que los
# haplotipos se leen de las moleculas cada error acuñaria su propia fila. Con
# el default (--phased-indel-min-reads 2) el complex_indel no sale; bajandolo
# a 1 si sale, con su unica lectura.
scenario_singleton_haplotype = Scenario(
    name="38_singleton_haplotype_suppressed",
    description="Haplotipo ins+SNV visto en 1 sola lectura: suprimido con el default de 2 lecturas",
    variants=[
        VcfRecord(pos=30, ref="T", alt="TGCT"),
        VcfRecord(pos=31, ref="G", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_ins_only",
            start=1,
            length=147,
            ops=[Op(kind="ins", pos=30, seq="GCT")],
            count=19,
        ),
        ReadGroup(
            name_prefix="r_ins_snv",
            start=1,
            length=147,
            ops=[Op(kind="ins", pos=30, seq="GCT"), Op(kind="snv", pos=31, seq="A")],
            count=1,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            variant_type="INDEL",
            event_class="insertion",
            event_components="INS:30:+GCT",
        ),
    ],
    expected_row_count=2,
)


# El mismo dato con el umbral bajado a 1: el haplotipo de una lectura si sale.
scenario_singleton_haplotype_opt_in = Scenario(
    name="39_singleton_haplotype_opt_in",
    description="Mismos datos que 38 con --phased-indel-min-reads 1: el haplotipo de 1 lectura si se emite",
    variants=scenario_singleton_haplotype.variants,
    reads=scenario_singleton_haplotype.reads,
    expected=[
        ExpectedRow(
            event_class="complex_indel",
            event_components="INS:30:+GCT | SNV:31:G>A",
            variant_type="INDEL",
            event_reads="1",
        ),
    ],
    expected_row_count=3,
    extra_cli_args=["--phased-indel-min-reads", "1"],
)


# ---------------------------------------------------------------------------
# 40. Fase declarada en el VCF (estilo WhatsHap/HiPhase/Clair3)
# ---------------------------------------------------------------------------
# Dos SNV del codon 10 con GT faseado y el mismo PS: el llamador ya dijo que
# van en la misma molecula. get_mnv no vuelve a fasear, pero lo publica en la
# columna Declared Phase, y las lecturas (que llevan ambas) no lo contradicen.
scenario_declared_phase_cis = Scenario(
    name="40_declared_phase_cis",
    description="VCF faseado con GT 1|0 y PS=42 en ambas SNV del codon: Declared Phase = cis:42",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T", genotype="1|0", phase_set=42),
        VcfRecord(pos=30, ref="T", alt="A", genotype="1|0", phase_set=42),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_cis",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=28, seq="T"), Op(kind="snv", pos=30, seq="A")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            declared_phase="cis:42",
            mnv_phasing_support="1.0000",
        ),
    ],
)


# 41. Fase declarada en haplotipos opuestos: el llamador dice trans
scenario_declared_phase_trans = Scenario(
    name="41_declared_phase_trans",
    description="VCF faseado con 1|0 y 0|1 en el mismo PS: Declared Phase = trans:42",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T", genotype="1|0", phase_set=42),
        VcfRecord(pos=30, ref="T", alt="A", genotype="0|1", phase_set=42),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_left",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=28, seq="T")],
            count=10,
        ),
        ReadGroup(
            name_prefix="r_right",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=30, seq="A")],
            count=10,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            declared_phase="trans:42",
            mnv_phasing_support="0.0000",
        ),
    ],
)


# 42. El VCF declara cis y las lecturas lo desmienten
# Ninguna de las 20 lecturas que cruzan el codon lleva las dos variantes, asi
# que la afirmacion del llamador no se sostiene. No es cuestion de grado: es
# cero. La columna lo marca en vez de repetir la afirmacion sin mas.
scenario_declared_phase_contradicted = Scenario(
    name="42_declared_phase_contradicted",
    description="El VCF declara cis pero ninguna lectura lleva ambas: Declared Phase lo marca como contradicho",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T", genotype="1|0", phase_set=42),
        VcfRecord(pos=30, ref="T", alt="A", genotype="1|0", phase_set=42),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_left",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=28, seq="T")],
            count=10,
        ),
        ReadGroup(
            name_prefix="r_right",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=30, seq="A")],
            count=10,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            declared_phase="cis:42|contradicted-by-reads",
            mnv_phasing_support="0.0000",
        ),
    ],
)


# 43. GT sin fasear: el llamador no afirma nada y get_mnv tampoco
scenario_unphased_declares_nothing = Scenario(
    name="43_unphased_genotype_declares_nothing",
    description="GT 1/1 (sin fasear) con PS presente: Declared Phase = '-', porque '/' no es una afirmacion de fase",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T", genotype="1/1", phase_set=42),
        VcfRecord(pos=30, ref="T", alt="A", genotype="1/1", phase_set=42),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_cis",
            start=1,
            length=150,
            ops=[Op(kind="snv", pos=28, seq="T"), Op(kind="snv", pos=30, seq="A")],
            count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            declared_phase="-",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 44. Un codon que cruza el intron, resuelto por el par de mates
# ---------------------------------------------------------------------------
# El codon 34 de geneC son las pos 900 (exon 1) y 1002 (exon 2), 102 pb aparte.
# Ninguna lectura de 40 pb alcanza las dos, pero un par si: mate 1 cubre la
# 900 y mate 2 la 1002. Son la misma molecula, asi que la pregunta de fase si
# tiene respuesta, y sale 1.0000 sobre 20 moleculas. Contando los mates por
# separado la respuesta seria '-' (nadie pudo mirar).
scenario_mate_pair_spans_junction = Scenario(
    name="44_mate_pair_resolves_junction_codon",
    description="Codon 34 partido por el intron: ninguna lectura sola lo cruza pero el par de mates si, y la fase se resuelve",
    variants=[
        VcfRecord(pos=900, ref="G", alt="A"),
        VcfRecord(pos=1002, ref="T", alt="C"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_pair",
            start=861,
            length=40,
            ops=[Op(kind="snv", pos=900, seq="A")],
            count=20,
            mate_start=1001,
            mate_length=40,
            mate_ops=[Op(kind="snv", pos=1002, seq="C")],
        ),
    ],
    expected=[
        ExpectedRow(
            positions="900, 1002",
            gene="geneC",
            variant_type="SNP/MNV",
            aa_changes="Ala34Thr",
            mnv_phasing_support="1.0000",
            mnv_phasing_reads="20",
        ),
    ],
    gff_content=GFF_CDS_MULTIEXON,
    gff_features="CDS",
)


# 45. Los mismos datos contando los mates por separado: nadie pudo mirar
scenario_mates_separately_cannot_answer = Scenario(
    name="45_mates_separately_cannot_answer",
    description="Mismos datos que 44 con --count-mates-separately: ninguna lectura cruza el codon, asi que la fase es desconocida",
    variants=scenario_mate_pair_spans_junction.variants,
    reads=scenario_mate_pair_spans_junction.reads,
    expected=[
        ExpectedRow(
            positions="900, 1002",
            gene="geneC",
            variant_type="SNP/MNV",
            mnv_phasing_support="-",
            mnv_phasing_reads="-",
        ),
    ],
    gff_content=GFF_CDS_MULTIEXON,
    gff_features="CDS",
    extra_cli_args=["--count-mates-separately"],
)


# ---------------------------------------------------------------------------
# 46. Mates solapados: una molecula, no dos observaciones
# ---------------------------------------------------------------------------
# Los dos mates cubren el codon entero, asi que cada molecula se ve dos veces.
# La profundidad son 20 moleculas, no 40 observaciones: contar las dos mitades
# del mismo fragmento infla la cobertura y hace que la frecuencia parezca
# sostenida por el doble de material del que hay.
scenario_overlapping_mates_counted_once = Scenario(
    name="46_overlapping_mates_counted_once",
    description="Mates que solapan el codon: 20 moleculas de profundidad, no 40 observaciones",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_overlap",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="T"), Op(kind="snv", pos=30, seq="A")],
            count=20,
            mate_start=1,
            mate_length=100,
            mate_ops=[Op(kind="snv", pos=28, seq="T"), Op(kind="snv", pos=30, seq="A")],
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            mnv_reads="20",
            total_reads="20",
            mnv_frequencies="1.0000",
            # Cada molecula se leyo por las dos hebras, asi que cuenta en las
            # dos. Acreditarla solo a forward dejaria el brazo reverse a cero en
            # todo dato paired-end, que se lee como sesgo de hebra total y hace
            # saltar --min-mnv-strand sobre datos perfectamente equilibrados.
            mnv_forward_reads="20",
            mnv_reverse_reads="20",
        ),
    ],
)


# 47. Mates que se contradicen en el solape
# Un mate ve la ALT en pos 28 y el otro ve la referencia. Uno de los dos se
# equivoca y no hay forma de saber cual, asi que la molecula no ha observado
# esa posicion. El codon queda sin lecturas que puedan responder por el.
scenario_mates_disagree = Scenario(
    name="47_mates_disagree_in_overlap",
    description="Los mates discrepan en pos 28: la molecula no observa esa posicion y la fase queda sin responder",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_conflict",
            start=1,
            length=100,
            ops=[Op(kind="snv", pos=28, seq="T"), Op(kind="snv", pos=30, seq="A")],
            count=20,
            mate_start=1,
            mate_length=100,
            # Sin la edicion en pos 28: este mate ve la referencia alli.
            mate_ops=[Op(kind="snv", pos=30, seq="A")],
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            mnv_reads="0",
            total_reads="0",
            mnv_phasing_support="-",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 48. Mates que se contradicen sobre un indel
# ---------------------------------------------------------------------------
# Los dos mates cubren la pos 30. Uno lleva la insercion y el otro lee la
# referencia alli. Uno de los dos se equivoca y no hay forma de saber cual, asi
# que la molecula no puede respaldar la llamada. El solape de mates es
# justamente donde aparecen los artefactos de realineamiento de indels, y sin
# esta regla un solo mate con un artefacto bastaba para acreditar la molecula
# entera: el resultado era identico al de 20 moleculas de acuerdo.
scenario_indel_mates_disagree = Scenario(
    name="48_indel_mates_disagree",
    description="Mates que discrepan sobre la insercion en pos 30: la molecula no respalda la llamada",
    variants=[VcfRecord(pos=30, ref="T", alt="TGCT")],
    reads=[
        ReadGroup(
            name_prefix="r_conflict",
            start=1,
            length=100,
            ops=[Op(kind="ins", pos=30, seq="GCT")],
            count=20,
            mate_start=1,
            mate_length=100,
            mate_ops=[],
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            variant_type="INDEL",
            event_class="insertion",
            event_reads="0",
            event_frequency="0.0000",
        ),
    ],
)


# 49. El mismo dato con los dos mates de acuerdo: si respalda
scenario_indel_mates_agree = Scenario(
    name="49_indel_mates_agree",
    description="Mismos fragmentos con la insercion en los dos mates: 20 moleculas la respaldan",
    variants=[VcfRecord(pos=30, ref="T", alt="TGCT")],
    reads=[
        ReadGroup(
            name_prefix="r_agree",
            start=1,
            length=100,
            ops=[Op(kind="ins", pos=30, seq="GCT")],
            count=20,
            mate_start=1,
            mate_length=100,
            mate_ops=[Op(kind="ins", pos=30, seq="GCT")],
        ),
    ],
    expected=[
        ExpectedRow(
            positions="30",
            variant_type="INDEL",
            event_class="insertion",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 50. Una lectura que termina en la base ancla de una insercion no la niega
# ---------------------------------------------------------------------------
# Ventana: delecion en 28, SNV en 30, insercion anclada en 31.
#   20 moleculas llevan las tres.
#   10 lecturas terminan justo en la base 31: ensenan la delecion y la SNV, y no
#      pueden ver si la insercion viene detras.
# Esas 10 son igual de compatibles con llevar la insercion. Tratarlas como si
# la negaran acuñaba una fila delecion+SNV que ninguna molecula demostro, en
# contra de la regla que la propia herramienta documenta: un subconjunto solo se
# emite si hay lecturas que lo ensenan.
scenario_read_ending_on_insertion_anchor = Scenario(
    name="50_read_ending_on_insertion_anchor",
    description="Lecturas que acaban en la base ancla de una insercion no cuentan como evidencia de que falta",
    variants=[
        VcfRecord(pos=28, ref="GC", alt="G"),
        VcfRecord(pos=30, ref="T", alt="A"),
        VcfRecord(pos=31, ref="G", alt="GGCT"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_full",
            start=1,
            length=100,
            ops=[
                Op(kind="del", pos=29, length=1),
                Op(kind="snv", pos=30, seq="A"),
                Op(kind="ins", pos=31, seq="GCT"),
            ],
            count=20,
        ),
        ReadGroup(
            name_prefix="r_stops",
            start=1,
            length=31,
            ops=[Op(kind="del", pos=29, length=1), Op(kind="snv", pos=30, seq="A")],
            count=10,
        ),
    ],
    expected=[
        ExpectedRow(
            event_class="complex_indel",
            event_components="DEL:29:C | SNV:30:T>A | INS:31:+GCT",
            event_reads="20",
        ),
    ],
    # Las tres filas del VCF mas el unico haplotipo que las lecturas ensenan.
    # Sin la regla salia una quinta: DEL:29:C | SNV:30:T>A.
    expected_row_count=4,
)


# ---------------------------------------------------------------------------
# 51. Dos sustituciones en el mismo codon que nunca comparten molecula
# ---------------------------------------------------------------------------
# 20 moleculas llevan la de pos 28 y otras 20 la de pos 30, ninguna las dos.
# El ratio de co-ocurrencia dice 0.0000, que es cierto pero no distingue "no
# hay evidencia" de "se excluyen activamente". D' = -1 con un p diminuto dice lo
# segundo: en una poblacion haploide eso son dos linajes en competencia, la
# firma de una infeccion mixta, no un MNV.
scenario_competing_haplotypes = Scenario(
    name="51_competing_haplotypes",
    description="Dos sustituciones del mismo codon que se excluyen: D' = -1, la firma de dos linajes en competencia",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(name_prefix="r_first", start=1, length=100,
                  ops=[Op(kind="snv", pos=28, seq="T")], count=20),
        ReadGroup(name_prefix="r_second", start=1, length=100,
                  ops=[Op(kind="snv", pos=30, seq="A")], count=20),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            mnv_phasing_support="0.0000",
            codon_ld="-1.0000",
        ),
    ],
)


# 52. Dos sustituciones comunes que solo coinciden
# Cada una esta en el 90% de las moleculas y se encuentran juntas en el 81%,
# que es exactamente lo que predice el azar. El ratio dice 0.9000 y se lee como
# ligamiento fuerte; D' = 0 dice que no hay ningun exceso sobre las frecuencias.
scenario_coincidence_not_linkage = Scenario(
    name="52_coincidence_is_not_linkage",
    description="Dos sustituciones comunes e independientes: el ratio dice 0.9000 y D' dice 0.0000",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(name_prefix="r_both", start=1, length=100,
                  ops=[Op(kind="snv", pos=28, seq="T"), Op(kind="snv", pos=30, seq="A")],
                  count=81),
        ReadGroup(name_prefix="r_first", start=1, length=100,
                  ops=[Op(kind="snv", pos=28, seq="T")], count=9),
        ReadGroup(name_prefix="r_second", start=1, length=100,
                  ops=[Op(kind="snv", pos=30, seq="A")], count=9),
        ReadGroup(name_prefix="r_neither", start=1, length=100, ops=[], count=1),
    ],
    expected=[
        ExpectedRow(
            positions="28, 30",
            variant_type="SNP/MNV",
            mnv_phasing_support="0.9000",
            codon_ld="0.0000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 53. La hebra de un mate no se presta a la posicion que leyo el otro
# ---------------------------------------------------------------------------
# Dos variantes en geneA. Cada fragmento es mate1 forward sobre la pos 28 (que
# lee referencia) y mate2 reverse sobre la pos 200 (que lleva la ALT). La 200
# solo la han leido mates reverse, asi que su brazo forward tiene que estar
# vacio. Con banderas de hebra a nivel de fragmento salia 20/20 y cualquier
# sitio parecia equilibrado, con lo que el sesgo de hebra no podia rechazar nada.
scenario_strand_is_per_position = Scenario(
    name="53_strand_is_per_position",
    description="Mates lejanos: la hebra que leyo una variante no cuenta para la que leyo la otra",
    variants=[
        VcfRecord(pos=28, ref="G", alt="A"),
        VcfRecord(pos=200, ref="C", alt="T"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_pair", start=1, length=60, ops=[], count=20,
            mate_start=160, mate_length=60,
            mate_ops=[Op(kind="snv", pos=200, seq="T")],
        ),
    ],
    expected=[
        ExpectedRow(
            positions="200",
            variant_type="SNP",
            snp_reads="20",
            snp_forward_reads="0",
            snp_reverse_reads="20",
        ),
    ],
)


# 54. Una delecion mas larga no es soporte exacto de la corta
# El VCF declara borrar solo la pos 29. Las lecturas borran 29 y 30. Dentro del
# tramo REF las dos son indistinguibles, asi que hacia falta mirar una base mas
# alla para ver donde termina la delecion de la lectura.
scenario_longer_deletion_is_not_support = Scenario(
    name="54_longer_deletion_is_not_support",
    description="Lecturas con una delecion de 2 bp no respaldan la delecion de 1 bp declarada",
    variants=[VcfRecord(pos=28, ref="GC", alt="G")],
    reads=[
        ReadGroup(
            name_prefix="r_longer_del", start=1, length=100,
            ops=[Op(kind="del", pos=29, length=2)], count=20,
        ),
    ],
    expected=[
        ExpectedRow(
            positions="28",
            variant_type="INDEL",
            event_class="deletion",
            event_reads="0",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 56. Un haplotipo dentro de otro se cuenta por sus propias moleculas
# ---------------------------------------------------------------------------
# 12 moleculas llevan las tres variantes y 8 solo las dos primeras. Las dos
# combinaciones son reales y salen las dos, pero la del par son 8 moleculas, no
# 20: las 12 del triple no son evidencia de una molecula con solo dos.
# El conteo exacto empareja el alelo compuesto sobre SU tramo, asi que una
# molecula que lleva eso y ademas otra cosa fuera del tramo tambien casa, y la
# fila del par salia con las 20 lecturas del triple al 100% de frecuencia.
scenario_nested_haplotype_counts = Scenario(
    name="56_nested_haplotype_counts_its_own_molecules",
    description="Par contenido en un triple: 8 moleculas para el par y 12 para el triple, no 20 para ambos",
    variants=[
        VcfRecord(pos=28, ref="G", alt="T"),
        VcfRecord(pos=29, ref="C", alt="CGCT"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_all_three", start=1, length=147,
            ops=[
                Op(kind="snv", pos=28, seq="T"),
                Op(kind="ins", pos=29, seq="GCT"),
                Op(kind="snv", pos=30, seq="A"),
            ],
            count=12,
        ),
        ReadGroup(
            name_prefix="r_first_two", start=1, length=147,
            ops=[Op(kind="snv", pos=28, seq="T"), Op(kind="ins", pos=29, seq="GCT")],
            count=8,
        ),
    ],
    expected=[
        ExpectedRow(
            event_class="complex_indel",
            event_components="SNV:28:G>T | INS:29:+GCT | SNV:30:T>A",
            event_reads="12",
            event_frequency="0.6000",
        ),
        ExpectedRow(
            event_class="complex_indel",
            event_components="SNV:28:G>T | INS:29:+GCT",
            event_reads="8",
            event_frequency="0.4000",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 57. Un haplotipo de indel real que aun asi es una casualidad
# ---------------------------------------------------------------------------
# 5 moleculas llevan insercion + SNV, 15 solo la insercion, 10 solo la SNV y 20
# ninguna. La combinacion existe y sale con sus 5 moleculas, pero cada variante
# esta al 40% y al 30%, asi que encontrarlas juntas en el 10% es lo que predice
# el azar: D' = -0.1667 con p = 0.75. El recuento dice cuantas moleculas son esa
# combinacion; D' dice si tiene algo que ver una con otra.
scenario_indel_haplotype_by_chance = Scenario(
    name="57_indel_haplotype_by_chance",
    description="Haplotipo indel+SNV en 5 moleculas cuyas variantes son independientes: D' cerca de 0 y p no significativo",
    variants=[
        VcfRecord(pos=29, ref="C", alt="CGCT"),
        VcfRecord(pos=30, ref="T", alt="A"),
    ],
    reads=[
        ReadGroup(name_prefix="r_both", start=1, length=147,
                  ops=[Op(kind="ins", pos=29, seq="GCT"), Op(kind="snv", pos=30, seq="A")],
                  count=5),
        ReadGroup(name_prefix="r_ins", start=1, length=147,
                  ops=[Op(kind="ins", pos=29, seq="GCT")], count=15),
        ReadGroup(name_prefix="r_snv", start=1, length=147,
                  ops=[Op(kind="snv", pos=30, seq="A")], count=10),
        ReadGroup(name_prefix="r_ref", start=1, length=147, ops=[], count=20),
    ],
    expected=[
        ExpectedRow(
            event_class="complex_indel",
            event_components="INS:29:+GCT | SNV:30:T>A",
            event_reads="5",
            codon_ld="-0.1667",
        ),
    ],
)


# ---------------------------------------------------------------------------
# 64: un MNV intergenico que las lecturas llevan entero sobrevive a un umbral
# de SNP. Sus lecturas de SNP son 0 por construccion, porque ninguna lectura
# lleva una base sola, asi que juzgarlo por ese lado lo borraba del VCF mientras
# el TSV lo conservaba con sus lecturas de haplotipo.
scenario_intergenic_mnv_threshold = Scenario(
    name="64_intergenic_mnv_survives_snp_threshold",
    description="MNV intergenico 750-751 con 12/20 lecturas ligadas y --snp 5: TSV y VCF lo conservan los dos",
    variants=[VcfRecord(pos=750, ref="AA", alt="TG")],
    reads=[
        ReadGroup(
            name_prefix="r_hap",
            start=700,
            length=100,
            ops=[Op(kind="snv", pos=750, seq="T"), Op(kind="snv", pos=751, seq="G")],
            count=12,
        ),
        ReadGroup(name_prefix="r_ref", start=700, length=100, ops=[], count=8),
    ],
    extra_cli_args=["--snp", "5"],
    expected=[
        ExpectedRow(
            gene="intergenic",
            positions="750, 751",
            variant_type="MNV",
            mnv_reads="12",
        ),
    ],
    expected_row_count=1,
    expected_vcf_positions=[750, 751],
)


# ---------------------------------------------------------------------------
# 65: un indel intergenico que 20 de 20 lecturas llevan sobrevive a un umbral
# de MNV. El VCF lo juzgaba con soporte 0 sobre profundidad 0 aunque la misma
# linea escribia ER=20 y EFREQ=1.0000, asi que el TSV lo conservaba y el VCF lo
# tiraba en cuanto habia un umbral activo.
scenario_intergenic_indel_threshold = Scenario(
    name="65_intergenic_indel_survives_mnv_threshold",
    description="Insercion fuera de todo gen con 20/20 lecturas y --mnv 5: TSV y VCF la conservan los dos",
    variants=[VcfRecord(pos=30, ref="T", alt="TGCT")],
    reads=[
        ReadGroup(
            name_prefix="r_interg_ins",
            start=1,
            length=151,
            ops=[Op(kind="ins", pos=30, seq="GCT")],
            count=20,
        ),
    ],
    gff_content=GFF_GENE_FAR_AWAY,
    extra_cli_args=["--mnv", "5"],
    expected=[
        ExpectedRow(
            positions="30",
            gene="intergenic",
            variant_type="INDEL",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=1,
    expected_vcf_positions=[30],
)


ALL_SCENARIOS = [
    scenario_snp_simple,
    scenario_snp_mnv_full,
    scenario_mnp_decomposed,
    scenario_snp_mnv_mixed,
    scenario_ins_inframe,
    scenario_del_frameshift,
    scenario_indel_plus_snv,
    scenario_fs_plus_downstream_snv,
    scenario_inframe_ins_with_mnv,
    scenario_fs_del_with_mnv,
    scenario_minus_snp_simple,
    scenario_minus_mnv,
    scenario_minus_fs_del,
    scenario_multiexon_snp_exon2,
    scenario_multiexon_junction_snp,
    scenario_multiexon_junction_mnv,
    scenario_no_bam_coverage,
    scenario_min_snp_frequency_filter,
    # nuevos AA edge cases:
    scenario_stop_gained,
    scenario_start_lost,
    scenario_stop_lost,
    # alelos complejos:
    scenario_intron_variant,
    scenario_multiallelic,
    scenario_delins,
    scenario_large_ins,
    scenario_large_del,
    # multi-contig + iVar:
    scenario_multicontig,
    scenario_ivar_snv,
    scenario_ivar_insertion,
    scenario_ivar_deletion,
    # refinamientos de indel (rama indels):
    scenario_stop_gained_inframe_ins,
    scenario_stop_lost_inframe_del,
    scenario_fs_gate_default,
    scenario_fs_gate_suppressed,
    # eucariotas (auto-CDS):
    scenario_eukaryote_autocds,
    # frameshift mas alla del stop prematuro:
    scenario_frameshift_past_stop,
    # haplotipos locales leidos de las moleculas:
    scenario_two_real_haplotypes,
    # evidencia de fase publicada y umbral de haplotipo:
    scenario_trans_evidence_reported,
    scenario_singleton_haplotype,
    scenario_singleton_haplotype_opt_in,
    # fase declarada por el llamador de entrada:
    scenario_declared_phase_cis,
    scenario_declared_phase_trans,
    scenario_declared_phase_contradicted,
    scenario_unphased_declares_nothing,
    # fase a nivel de fragmento:
    scenario_mate_pair_spans_junction,
    scenario_mates_separately_cannot_answer,
    scenario_overlapping_mates_counted_once,
    scenario_mates_disagree,
    scenario_indel_mates_disagree,
    scenario_indel_mates_agree,
    scenario_read_ending_on_insertion_anchor,
    scenario_competing_haplotypes,
    scenario_coincidence_not_linkage,
    scenario_strand_is_per_position,
    scenario_longer_deletion_is_not_support,
    scenario_nested_haplotype_counts,
    scenario_indel_haplotype_by_chance,
    scenario_fs_gate_no_declared_af,
    scenario_intergenic_indel_support,
    scenario_uncountable_flags,
    scenario_right_anchored_insertion,
    scenario_intergenic_mnv_counted,
    scenario_slippage_shared_base,
    scenario_intergenic_mnv_threshold,
    scenario_intergenic_indel_threshold,
]


# Por defecto los escenarios simples no exigen numero de filas exacto.
# Para los escenarios donde get_mnv emite filas adicionales por exploracion
# de haplotipos locales, comprobamos numero exacto via expected_row_count.
scenario_snp_simple.expected_row_count = 1
scenario_snp_mnv_full.expected_row_count = 1
scenario_mnp_decomposed.expected_row_count = 1
scenario_snp_mnv_mixed.expected_row_count = 1
scenario_ins_inframe.expected_row_count = 1
scenario_del_frameshift.expected_row_count = 1
scenario_indel_plus_snv.expected_row_count = 2
