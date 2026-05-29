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
# Las 20 lecturas llevan todo en cis. Distancias 1-2 bp -> ventana local.
# Get_mnv explora subconjuntos del haplotipo local y emite 5 filas:
#   1) complex_indel SNV:28 + INS:29              -> Ala10delinsSerLeu
#   2) complex_indel SNV:28 + INS:29 + SNV:30     -> Ala10delinsSerLeu (todo)
#   3) SNP/MNV pos28+pos30 marcado Indel overlap  -> AA Unknown
#   4) insertion INS:29 sola                       -> 10_11insLeu
#   5) complex_indel INS:29 + SNV:30              -> 10_11insLeu
scenario_inframe_ins_with_mnv = Scenario(
    name="08_inframe_ins_inside_codon_with_mnv",
    description="Ins in-frame DENTRO de codon 10 + MNV pos28+pos30 todo en cis: 5 filas (2 complex_indel, MNV-overlap, ins solo, complex_indel ins+SNV)",
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
        # complex_indel SNV:28 + INS:29
        ExpectedRow(
            event_class="complex_indel",
            event_components="SNV:28:G>T | INS:29:+GCT",
            variant_type="INDEL",
            change_type="In-frame Indel",
            aa_changes="Ala10delinsSerLeu",
            event_reads="20",
            event_frequency="1.0000",
        ),
        # complex_indel completo SNV:28 + INS:29 + SNV:30
        ExpectedRow(
            event_class="complex_indel",
            event_components="SNV:28:G>T | INS:29:+GCT | SNV:30:T>A",
            variant_type="INDEL",
            change_type="In-frame Indel",
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
        # insertion sola
        ExpectedRow(
            positions="29",
            variant_type="INDEL",
            event_class="insertion",
            event_components="INS:29:+GCT",
            aa_changes="10_11insLeu",
            change_type="In-frame Indel",
            event_reads="20",
            event_frequency="1.0000",
        ),
        # complex_indel INS:29 + SNV:30 (sin pos28)
        ExpectedRow(
            event_class="complex_indel",
            event_components="INS:29:+GCT | SNV:30:T>A",
            variant_type="INDEL",
            change_type="In-frame Indel",
            event_reads="20",
            event_frequency="1.0000",
        ),
    ],
    expected_row_count=5,
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
    description="SNP non-sin en geneB (hebra -): pos 673 C>T -> mRNA codon 10 GCT->ACT (Ala10Thr). NOTA: codones mostrados en +strand genomico",
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
            reference_codon="AGC",   # genomic (RC=GCT mRNA Ala)
            snp_codon="AGT",          # genomic (RC=ACT mRNA Thr)
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
    description="MNV en geneB codon 10: pos 671 A>T + pos 673 C>A en cis (mRNA GCT->TCA, Ala10Ser); codones genomicos",
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
            reference_codon="AGC",   # genomic (RC=GCT mRNA)
            mnv_codon="TGA",          # genomic (RC=TCA mRNA Ser)
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
    description="SNV en intron de geneC (pos 950 T>A): debe reportarse como intergenic",
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
            gene="intergenic",
            variant_type="SNP",
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
    variants=[VcfRecord(pos=28, ref="G", alt="A,T")],
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
# iVar: POS=31, REF=G, ALT=-G  ->  borrar pos 31 (G). Anchor en pos 30.
scenario_ivar_deletion = Scenario(
    name="29_ivar_tsv_deletion",
    description="iVar TSV: delecion frameshift -G en pos 31 (notacion -SEQ)",
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
        IvarRecord(pos=31, ref="G", alt="-G", total_dp=20, alt_dp=20, alt_freq=1.0),
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
            aa_changes="10_11ins*",
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


# 32. Gate de frameshift --frameshift-min-freq: COMPORTAMIENTO POR DEFECTO.
# Misma estructura que el escenario 07: del frameshift de 1bp en pos 31 +
# SNV downstream pos 39. El VCF declara AF=0.20 para la del y AF=0.95 para
# el SNV. Por defecto (--frameshift-min-freq 0.0) cualquier indel upstream
# propaga el frameshift -> el SNV downstream lleva marcador "(fs)".
# Comparar con el escenario 33 (mismos inputs, flag 0.5).
scenario_fs_gate_default = Scenario(
    name="32_fs_gate_default_propagates",
    description="Del fs upstream AF=0.20 + SNV downstream: por defecto el frameshift se propaga -> SNV con '(fs)'",
    variants=[
        VcfRecord(pos=30, ref="TG", alt="T", af=0.20),
        VcfRecord(pos=39, ref="T", alt="A", af=0.95),
    ],
    reads=[
        ReadGroup(
            name_prefix="r_gate_def",
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
        ReadGroup(
            name_prefix="r_gate_sup",
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
