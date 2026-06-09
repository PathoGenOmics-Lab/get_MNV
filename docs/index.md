# get_MNV

**Codon-level Multi-Nucleotide Variant detection and indel annotation from VCF or iVar TSV.**
Pure Rust · no C dependencies · cross-platform (macOS, Linux, Windows).

get_MNV finds cases where two or more SNVs fall in the **same codon** and should
be interpreted together. These combined changes can produce a different amino
acid effect than the individual SNVs alone.

<div style="text-align: center;" markdown>
![MNV amino acid reclassification](assets/get_mnv_aa.png){ width="620" }
</div>

## What it takes

| Input | Format |
|---|---|
| Variant calls | VCF or iVar `variants.tsv` |
| Reference sequence | FASTA |
| Gene annotation | GFF/GFF3/GTF or a simple TSV |
| Aligned reads *(optional)* | BAM, to count SNP/MNV/indel event support |

It writes annotated variants as TSV, VCF, or both.

## Main features

- Groups SNVs by codon and reports **SNP**, **MNV**, or **SNP/MNV** calls.
- Recalculates amino acid changes from the **full codon haplotype**.
- Decomposes `REF/ALT` alleles into SNV, MNV, insertion, deletion, delins, and
  complex indel components.
- Reads VCF and iVar TSV, including iVar `+SEQ` / `-SEQ` indel notation.
- Uses BAM reads (when provided) for SNP/MNV support, exact indel event support,
  and strand bias.
- Supports 9 NCBI genetic code tables.
- Includes a desktop GUI for drag-and-drop analysis.

## Installation

=== "Desktop GUI"

    Download the latest release for your platform from the
    [Releases page](https://github.com/PathoGenOmics-Lab/get_MNV/releases/latest).

    !!! note "macOS users"
        The app is not signed with an Apple Developer certificate. On first
        launch, right-click the app → **Open** → click **Open** in the dialog.

=== "Command line"

    ```bash
    conda install -c bioconda get_mnv
    ```

    Or build from source:

    ```bash
    git clone https://github.com/PathoGenOmics-Lab/get_MNV.git
    cd get_MNV
    cargo install --path .
    ```

## Quick start

```bash
get_mnv \
  --vcf variants.vcf \
  --fasta reference.fasta \
  --gff genes.gff3
```

A ready-to-run *M. tuberculosis* dataset (reference, genes, VCF, and a tiny demo
BAM for the read viewer) lives in the
[`example/`](https://github.com/PathoGenOmics-Lab/get_MNV/tree/indels/example)
folder of the repository.

!!! tip "Next steps"
    - [Usage](usage.md) — common commands and what each argument means.
    - [Input Formats](input-formats.md) — VCF, iVar TSV, FASTA, and annotations.
    - [Output Formats](output-formats.md) — TSV/VCF/BCF/JSON columns and fields.
    - [Indel & MNV Semantics](indel-mnv-semantics.md) — boundary rules and tuning.
    - [Troubleshooting](troubleshooting.md) — fixes for common input mismatches.

## Citation

If you use get_MNV, please cite it via its
[Zenodo DOI](https://doi.org/10.5281/zenodo.13907423).

## License

get_MNV is released under the **AGPL-3.0** license.
