# Usage

## Basic command

```bash
get_mnv --vcf <VCF_FILE> --fasta <FASTA_FILE> (--genes <GENES_FILE> | --gff <GFF_FILE>)
```

## Required arguments

| Argument | Description |
|----------|-------------|
| `-v, --vcf <VCF_FILE>` | VCF file containing the SNVs |
| `-f, --fasta <FASTA_FILE>` | FASTA file with the reference sequence |
| `-g, --genes <GENES_FILE>` | Gene annotation file in TSV format (if `--gff` is not used) |
| `--gff <GFF_FILE>` | Gene annotation file in GFF/GFF3 format (if `--genes` is not used) |

## Optional arguments

### Input options

| Argument | Description |
|----------|-------------|
| `-b, --bam <BAM_FILE>` | BAM file with aligned reads |
| `--sample <SAMPLE>` | Sample name for FORMAT metrics (default: first sample). Use `all` to process every sample |
| `--chrom <CHROM>` | Comma-separated contig selection (default: all contigs in VCF) |
| `--gff-features <FEATURES>` | Comma-separated GFF feature types (default: `gene,pseudogene`) |

### Quality filters

| Argument | Default | Description |
|----------|---------|-------------|
| `-q, --quality <QUALITY>` | 20 | Minimum Phred quality score |
| `--mapq <MAPQ>` | 0 | Minimum mapping quality (MAPQ) |
| `--snp <N>` | 0 | Minimum supporting reads for SNP |
| `--mnv <N>` | 0 | Minimum supporting reads for MNV |
| `--min-snp-strand <N>` | 0 | Minimum per-strand reads for SNP |
| `--min-mnv-strand <N>` | 0 | Minimum per-strand reads for MNV |
| `--min-strand-bias-p <P>` | 0 | Fisher exact strand-bias p-value threshold |
| `--strict` | off | Fail if original metrics (ODP/OFREQ) are missing |

### Output options

| Argument | Description |
|----------|-------------|
| `--convert` | Output VCF instead of TSV |
| `--both` | Output both TSV and VCF |
| `--vcf-gz` | Write BGZF-compressed `.MNV.vcf.gz` |
| `--index-vcf-gz` | Build Tabix index for `.vcf.gz` |
| `--bcf` | Also write BCF output |
| `--emit-filtered` | Keep filtered variants with FILTER tags instead of skipping |
| `--strand-bias-info` | Add strand-bias p-values to VCF INFO |
| `--keep-original-info` | Preserve original VCF INFO fields in output |
| `--exclude-intergenic` | Exclude variants outside annotated genes |
| `--summary-json <FILE>` | Write structured run summary (JSON) |
| `--error-json <FILE>` | Write structured error details (JSON) |
| `--run-manifest <FILE>` | Write reproducibility manifest with checksums |

### Performance

| Argument | Description |
|----------|-------------|
| `--threads <N>` | Number of worker threads (default: Rayon auto) |
| `--dry-run` | Validate inputs without writing output files |
| `--normalize-alleles` | Normalize REF/ALT alleles before processing |
| `--split-multiallelic` | Split multiallelic VCF records into individual ALT alleles |

## Examples

### Basic MNV detection

```bash
get_mnv --vcf variants.vcf --fasta reference.fasta --gff genes.gff3
```

### With BAM reads and quality filters

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --quality 30 \
  --mapq 20
```

### VCF output with strand-bias filtering

```bash
get_mnv \
  --vcf variants.vcf \
  --bam reads.bam \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --both \
  --vcf-gz \
  --index-vcf-gz \
  --emit-filtered \
  --strand-bias-info \
  --min-strand-bias-p 0.05
```

### Multi-sample processing

```bash
get_mnv \
  --vcf multisample.vcf \
  --fasta reference.fasta \
  --gff genes.gff3 \
  --sample all \
  --threads 8
```
