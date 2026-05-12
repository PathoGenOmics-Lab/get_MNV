// TypeScript types matching Rust structs from get_mnv

export interface AnalysisConfig {
  vcfFile: string;
  inputFormat: "auto" | "vcf" | "tsv";
  bamFile?: string;
  fastaFile: string;
  genesFile: string;
  gffFeatures: string[];
  sample?: string;
  chrom?: string;
  minQuality: number;
  minMapq: number;
  threads?: number;
  minSnpReads: number;
  minSnpFrequency: number;
  minMnvReads: number;
  minMnvFrequency: number;
  minSnpStrandReads: number;
  minMnvStrandReads: number;
  minStrandBiasP: number;
  normalizeAlleles: boolean;
  splitMultiallelic: boolean;
  strict: boolean;
  emitFiltered: boolean;
  strandBiasInfo: boolean;
  keepOriginalInfo: boolean;
  excludeIntergenic: boolean;
  translationTable: number;
  outputTsv: boolean;
  outputVcf: boolean;
  vcfGz: boolean;
  outputDir?: string;
  outputPrefix?: string;
}

export interface TsvData {
  headers: string[];
  rows: string[][];
}

export interface BamVariantSite {
  position: number;
  referenceBase: string;
  altBase: string;
}

export interface BamReadView {
  name: string;
  strand: string;
  support: "mnv" | "partial" | "reference" | "other";
  start: number;
  end: number;
  mapq: number;
  bases: string[];
}

export interface BamSupportCounts {
  total: number;
  mnv: number;
  partial: number;
  reference: number;
  other: number;
}

export interface BamViewResponse {
  chrom: string;
  displayStart: number;
  displayEnd: number;
  reference: string;
  sites: BamVariantSite[];
  reads: BamReadView[];
  counts: BamSupportCounts;
  totalReads: number;
  truncated: boolean;
  /** Per-position depth from ALL reads (not just displayed subset). */
  coverage: number[];
}

export interface ContigSummary {
  contig: string;
  snp_records_in_vcf: number;
  mapped_genes: number;
  produced_variants: number;
  snp_variants: number;
  mnv_variants: number;
  snp_mnv_variants: number;
  indel_variants: number;
  intergenic_variants: number;
  region_cache_hits: number;
  region_cache_misses: number;
}

export interface GlobalSummary {
  contig_count: number;
  snp_records_in_vcf: number;
  mapped_genes: number;
  produced_variants: number;
  snp_variants: number;
  mnv_variants: number;
  snp_mnv_variants: number;
  indel_variants: number;
  intergenic_variants: number;
  region_cache_hits: number;
  region_cache_misses: number;
}

export interface RunTimings {
  parse_inputs_ms: number;
  process_ms: number;
  emit_ms: number;
  total_ms: number;
}

export interface RunInputs {
  vcf: string;
  fasta: string;
  annotation: string;
  bam?: string;
  checksums: InputChecksums;
}

export interface InputChecksums {
  vcf_sha256: string;
  fasta_sha256: string;
  annotation_sha256: string;
  bam_sha256?: string;
}

export interface RunSummary {
  schema_version: string;
  sample?: string;
  dry_run: boolean;
  bam_provided: boolean;
  translation_table: number;
  inputs: RunInputs;
  output_tsv?: string;
  output_vcf?: string;
  output_bcf?: string;
  contigs: ContigSummary[];
  timings: RunTimings;
  global: GlobalSummary;
}

export interface SampleEntry {
  id: string;
  name: string;
  vcfPath: string;
  bamPath?: string;
  result?: RunSummary;
  tsvData?: TsvData;
  status: "pending" | "running" | "done" | "error";
  error?: string;
  runConfig?: AnalysisConfig;
}

export type VariantType = "Snp" | "Mnv" | "SnpMnv" | "Indel";

export type ChangeType =
  | "Unknown"
  | "Synonymous"
  | "NonSynonymous"
  | "StopGained"
  | "StopLost"
  | "IndelOverlap"
  | "FrameshiftSynonymous"
  | "FrameshiftNonSynonymous"
  | "FrameshiftStopGained"
  | "FrameshiftStopLost"
  | "FrameshiftUnknown"
  | "FrameshiftIndel"
  | "InFrameIndel";

export const DEFAULT_CONFIG: AnalysisConfig = {
  vcfFile: "",
  inputFormat: "auto",
  fastaFile: "",
  genesFile: "",
  gffFeatures: ["gene", "pseudogene"],
  minQuality: 20,
  minMapq: 20,
  minSnpReads: 2,
  minSnpFrequency: 0,
  minMnvReads: 2,
  minMnvFrequency: 0,
  minSnpStrandReads: 0,
  minMnvStrandReads: 0,
  minStrandBiasP: 0,
  normalizeAlleles: true,
  splitMultiallelic: true,
  strict: false,
  emitFiltered: false,
  strandBiasInfo: false,
  keepOriginalInfo: false,
  excludeIntergenic: false,
  translationTable: 11,
  outputTsv: true,
  outputVcf: false,
  vcfGz: false,
};

export interface ConfigPreset {
  name: string;
  description: string;
  config: Partial<AnalysisConfig>;
}

export const BUILT_IN_PRESETS: ConfigPreset[] = [
  {
    name: "Default",
    description: "Balanced sensitivity and specificity",
    config: {},
  },
  {
    name: "Strict",
    description: "High confidence MNVs only",
    config: { minQuality: 30, minMapq: 30, minMnvReads: 4, minSnpReads: 4, strict: true },
  },
  {
    name: "Exploratory",
    description: "Maximally sensitive, more false positives",
    config: { minQuality: 10, minMapq: 10, minSnpReads: 1, minMnvReads: 1, emitFiltered: true },
  },
  {
    name: "Publication",
    description: "Strict + strand bias + VCF output",
    config: { minQuality: 30, minMapq: 30, minMnvReads: 4, minSnpReads: 4, strict: true, strandBiasInfo: true, outputVcf: true },
  },
];
