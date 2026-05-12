import { useState, useRef, useCallback } from "react";
import { createPortal } from "react-dom";
import { open } from "@tauri-apps/plugin-dialog";
import { type AnalysisConfig, DEFAULT_CONFIG, BUILT_IN_PRESETS } from "../types";

interface ParameterFormProps {
  config: AnalysisConfig;
  onChange: (config: AnalysisConfig) => void;
  /** Is the loaded annotation file a GFF/GFF3/GTF? */
  isGff: boolean;
  /** Available feature types detected from the GFF file */
  availableFeatures: string[];
}

/* ── Tooltip icon with hover popover (rendered via portal so it's never clipped) ── */
function Tip({ text }: { text: string }) {
  const [show, setShow] = useState(false);
  const [pos, setPos] = useState({ x: 0, y: 0 });
  const iconRef = useRef<HTMLSpanElement>(null);

  const handleEnter = useCallback(() => {
    if (iconRef.current) {
      const rect = iconRef.current.getBoundingClientRect();
      // getBoundingClientRect is viewport-relative, which is correct for position:fixed tooltips
      setPos({ x: rect.left + rect.width / 2, y: rect.top });
    }
    setShow(true);
  }, []);

  return (
    <>
      <span
        className="tooltip-icon"
        ref={iconRef}
        onMouseEnter={handleEnter}
        onMouseLeave={() => setShow(false)}
      >
        ?
      </span>
      {show &&
        createPortal(
          <span
            className="tooltip-popover tooltip-popover--fixed"
            style={{ left: pos.x, top: pos.y }}
          >
            {text}
          </span>,
          document.body
        )}
    </>
  );
}

/* ── Slider + number input combo ── */
interface SliderFieldProps {
  label: string;
  tip: string;
  value: number;
  min: number;
  max: number;
  step?: number;
  onChange: (v: number) => void;
}

function SliderField({ label, tip, value, min, max, step = 1, onChange }: SliderFieldProps) {
  const pct = max > min ? ((value - min) / (max - min)) * 100 : 0;
  const isDecimal = step < 1;
  // Format display value: always use dot as decimal separator
  const displayValue = isDecimal ? value.toFixed(String(step).split(".")[1]?.length ?? 2) : String(value);

  return (
    <div className="param-slider-row">
      <div className="param-slider-label">
        <span>{label}</span>
        <Tip text={tip} />
      </div>
      <div className="param-slider-control">
        <input
          type="range"
          className="param-range"
          min={min}
          max={max}
          step={step}
          value={value}
          onChange={(e) => onChange(Number(e.target.value))}
          style={{ "--pct": `${pct}%` } as React.CSSProperties}
        />
        <input
          type="text"
          inputMode="decimal"
          className="param-number"
          value={displayValue}
          onChange={(e) => {
            const raw = e.target.value.replace(",", ".");
            const v = Number(raw);
            if (!isNaN(v) && raw !== "" && raw !== ".") onChange(Math.max(min, Math.min(max, v)));
          }}
        />
      </div>
    </div>
  );
}

/* ── Toggle switch ── */
interface ToggleFieldProps {
  label: string;
  tip: string;
  checked: boolean;
  onChange: (v: boolean) => void;
}

function ToggleField({ label, tip, checked, onChange }: ToggleFieldProps) {
  return (
    <label className="param-toggle-row">
      <div className="param-toggle-label">
        <span>{label}</span>
        <Tip text={tip} />
      </div>
      <button
        type="button"
        role="switch"
        aria-checked={checked}
        className={`toggle-switch${checked ? " toggle-switch--on" : ""}`}
        onClick={() => onChange(!checked)}
      >
        <span className="toggle-knob" />
      </button>
    </label>
  );
}

/* ── Collapsed/expanded group ── */
function ParamGroup({
  title,
  accent,
  children,
  defaultOpen = true,
}: {
  title: string;
  accent?: string;
  children: React.ReactNode;
  defaultOpen?: boolean;
}) {
  const [open, setOpen] = useState(defaultOpen);

  return (
    <div className={`param-group${open ? " param-group--open" : ""}`}>
      <button
        type="button"
        className="param-group-header"
        onClick={() => setOpen(!open)}
      >
        {accent && <span className="param-group-accent" style={{ background: accent }} />}
        <h4>{title}</h4>
        <svg
          className="param-group-chevron"
          width="14"
          height="14"
          viewBox="0 0 14 14"
          fill="none"
          stroke="currentColor"
          strokeWidth="2"
          strokeLinecap="round"
        >
          <path d="M4 5.5l3 3 3-3" />
        </svg>
      </button>
      {open && <div className="param-group-body">{children}</div>}
    </div>
  );
}

export default function ParameterForm({ config, onChange, isGff, availableFeatures }: ParameterFormProps) {
  const update = <K extends keyof AnalysisConfig>(
    key: K,
    value: AnalysisConfig[K]
  ) => {
    onChange({ ...config, [key]: value });
  };

  const updateMany = (updates: Partial<AnalysisConfig>) => {
    onChange({ ...config, ...updates });
  };

  const toggleFeature = (feature: string) => {
    const current = config.gffFeatures;
    const next = current.includes(feature)
      ? current.filter((f) => f !== feature)
      : [...current, feature];
    update("gffFeatures", next);
  };

  // Detect which preset (if any) matches the current config
  const activePreset = BUILT_IN_PRESETS.find((p) => {
    const merged = { ...DEFAULT_CONFIG, ...p.config };
    // Compare all non-file, non-output-path fields
    return (
      config.minQuality === merged.minQuality &&
      config.minMapq === merged.minMapq &&
      config.minSnpReads === merged.minSnpReads &&
      config.minSnpFrequency === merged.minSnpFrequency &&
      config.minMnvReads === merged.minMnvReads &&
      config.minMnvFrequency === merged.minMnvFrequency &&
      config.minSnpStrandReads === merged.minSnpStrandReads &&
      config.minMnvStrandReads === merged.minMnvStrandReads &&
      config.minStrandBiasP === merged.minStrandBiasP &&
      config.normalizeAlleles === merged.normalizeAlleles &&
      config.splitMultiallelic === merged.splitMultiallelic &&
      config.strict === merged.strict &&
      config.emitFiltered === merged.emitFiltered &&
      config.strandBiasInfo === merged.strandBiasInfo &&
      config.keepOriginalInfo === merged.keepOriginalInfo &&
      config.translationTable === merged.translationTable &&
      config.outputTsv === merged.outputTsv &&
      config.outputVcf === merged.outputVcf &&
      config.vcfGz === merged.vcfGz
    );
  });

  const applyPreset = (name: string) => {
    const preset = BUILT_IN_PRESETS.find((p) => p.name === name);
    if (!preset) return;
    // Merge preset config on top of defaults, preserving user-selected files
    onChange({
      ...DEFAULT_CONFIG,
      ...preset.config,
      // Preserve user-selected files and environment settings
      vcfFile: config.vcfFile,
      inputFormat: config.inputFormat,
      fastaFile: config.fastaFile,
      genesFile: config.genesFile,
      bamFile: config.bamFile,
      gffFeatures: config.gffFeatures,
      sample: config.sample,
      chrom: config.chrom,
      threads: config.threads,
      keepOriginalInfo: config.keepOriginalInfo,
      translationTable: config.translationTable,
      excludeIntergenic: config.excludeIntergenic,
      vcfGz: config.vcfGz,
      outputDir: config.outputDir,
      outputPrefix: config.outputPrefix,
    });
  };

  // Infer strand filter method from config values
  type StrandMethod = "none" | "reads" | "bias";
  const hasStrandReads = config.minSnpStrandReads > 0 || config.minMnvStrandReads > 0;
  const hasStrandBias = config.minStrandBiasP > 0;
  const strandMethod: StrandMethod = hasStrandBias ? "bias" : hasStrandReads ? "reads" : "none";

  const setStrandMethod = (method: StrandMethod) => {
    switch (method) {
      case "none":
        updateMany({ minSnpStrandReads: 0, minMnvStrandReads: 0, minStrandBiasP: 0 });
        break;
      case "reads":
        // Switch to reads: zero out bias, set sensible defaults if both were 0
        updateMany({
          minStrandBiasP: 0,
          minSnpStrandReads: config.minSnpStrandReads || 1,
          minMnvStrandReads: config.minMnvStrandReads || 1,
        });
        break;
      case "bias":
        // Switch to bias: zero out strand reads, set sensible default if bias was 0
        updateMany({
          minSnpStrandReads: 0,
          minMnvStrandReads: 0,
          minStrandBiasP: config.minStrandBiasP || 0.05,
        });
        break;
    }
  };

  return (
    <div className="parameter-form">
      {/* ── Preset selector ── */}
      <div className="preset-bar">
        <span className="preset-label">Preset</span>
        <div className="preset-chips">
          {BUILT_IN_PRESETS.map((p) => (
            <button
              key={p.name}
              type="button"
              className={`preset-chip${activePreset?.name === p.name ? " preset-chip--active" : ""}`}
              onClick={() => applyPreset(p.name)}
              title={p.description}
            >
              {p.name}
            </button>
          ))}
          {!activePreset && <span className="preset-chip preset-chip--custom">Custom</span>}
        </div>
      </div>

      {/* ── Two-column grid for main groups ── */}
      <div className="param-grid">
        {/* Left column */}
        <ParamGroup title="Quality Filters" accent="#149389">
          <SliderField
            label="Min Phred quality"
            tip="Minimum base quality score (Phred scale). Variants below this are filtered. 20 ≈ 99% accuracy."
            value={config.minQuality}
            min={0}
            max={60}
            onChange={(v) => update("minQuality", v)}
          />
          <SliderField
            label="Min MAPQ"
            tip="Minimum mapping quality for aligned reads. Higher values require more confident alignments. 0 = no filter."
            value={config.minMapq}
            min={0}
            max={60}
            onChange={(v) => update("minMapq", v)}
          />
        </ParamGroup>

        {/* Right column */}
        <ParamGroup title="Read Support" accent="#E7216A">
          <SliderField
            label="Min SNP reads"
            tip="Minimum number of reads supporting a SNP call. 0 = accept any read depth."
            value={config.minSnpReads}
            min={0}
            max={50}
            onChange={(v) => update("minSnpReads", v)}
          />
          <SliderField
            label="Min MNV reads"
            tip="Minimum reads where multiple SNPs co-occur on the same read, confirming they form an MNV."
            value={config.minMnvReads}
            min={0}
            max={50}
            onChange={(v) => update("minMnvReads", v)}
          />
          <SliderField
            label="Min SNP frequency"
            tip="Minimum BAM-derived allele frequency for SNP records. 0 = disabled; requires a BAM file when enabled."
            value={config.minSnpFrequency}
            min={0}
            max={1}
            step={0.01}
            onChange={(v) => update("minSnpFrequency", v)}
          />
          <SliderField
            label="Min MNV frequency"
            tip="Minimum BAM-derived haplotype frequency for MNV records. 0 = disabled; requires a BAM file when enabled."
            value={config.minMnvFrequency}
            min={0}
            max={1}
            step={0.01}
            onChange={(v) => update("minMnvFrequency", v)}
          />

          <div className="param-divider" />

          {/* ── Strand filtering method selector ── */}
          <div className="param-slider-label" style={{ marginBottom: "0.3rem" }}>
            <span>Strand filtering</span>
            <Tip text="Choose how to filter variants by strand support. 'Min reads' requires a minimum number of reads on each strand. 'Strand bias' uses Fisher's exact test to detect strand imbalance." />
          </div>
          <div className="strand-method-selector">
            <button
              type="button"
              className={`strand-method-btn${strandMethod === "none" ? " strand-method-btn--active" : ""}`}
              onClick={() => setStrandMethod("none")}
            >
              None
            </button>
            <button
              type="button"
              className={`strand-method-btn${strandMethod === "reads" ? " strand-method-btn--active" : ""}`}
              onClick={() => setStrandMethod("reads")}
            >
              Min reads
            </button>
            <button
              type="button"
              className={`strand-method-btn${strandMethod === "bias" ? " strand-method-btn--active" : ""}`}
              onClick={() => setStrandMethod("bias")}
            >
              Strand bias
            </button>
          </div>

          {/* Show sliders only for the selected method */}
          {strandMethod === "reads" && (
            <>
              <SliderField
                label="Min SNP strand reads"
                tip="Minimum reads per strand for SNP calls. Ensures the variant is observed on both forward and reverse strands."
                value={config.minSnpStrandReads}
                min={0}
                max={30}
                onChange={(v) => update("minSnpStrandReads", v)}
              />
              <SliderField
                label="Min MNV strand reads"
                tip="Minimum reads per strand supporting co-occurring SNPs (MNV phasing). 0 = no strand requirement."
                value={config.minMnvStrandReads}
                min={0}
                max={30}
                onChange={(v) => update("minMnvStrandReads", v)}
              />
            </>
          )}
          {strandMethod === "bias" && (
            <SliderField
              label="Strand bias p-value"
              tip="Fisher's exact test p-value threshold for strand bias. Variants below this are flagged. 0 = disabled."
              value={config.minStrandBiasP}
              min={0}
              max={1}
              step={0.01}
              onChange={(v) => update("minStrandBiasP", v)}
            />
          )}
        </ParamGroup>
      </div>

      {/* ── Full-width GFF row — only show when annotation file is GFF format ── */}
      {isGff && (
        <ParamGroup title="GFF Features" accent="#6F398D">
          {availableFeatures.length > 0 && (
            <div className="gff-chips-row">
              <div className="param-slider-label">
                <span>Detected features</span>
                <Tip text="Feature types found in your GFF file. Click to toggle which ones to use as gene regions." />
              </div>
              <div className="gff-chips">
                {availableFeatures.map((f) => (
                  <button
                    key={f}
                    type="button"
                    className={`gff-chip${config.gffFeatures.includes(f) ? " gff-chip--on" : ""}`}
                    onClick={() => toggleFeature(f)}
                  >
                    {f}
                  </button>
                ))}
              </div>
            </div>
          )}
          <div className="param-text-row">
            <div className="param-slider-label">
              <span>Custom types</span>
              <Tip text="Comma-separated GFF3 feature types used as gene regions for codon analysis. Edit freely or use the chips above." />
            </div>
            <input
              type="text"
              className="param-text-input"
              value={config.gffFeatures.join(", ")}
              onChange={(e) =>
                update(
                  "gffFeatures",
                  e.target.value
                    .split(",")
                    .map((s) => s.trim())
                    .filter((s) => s.length > 0)
                )
              }
              placeholder="gene, pseudogene"
            />
          </div>
        </ParamGroup>
      )}

      {/* ── Two-column grid for toggles ── */}
      <div className="param-grid">
        <ParamGroup title="Analysis Options" accent="#149389">
          <ToggleField
            label="Normalize alleles"
            tip="Left-align and trim VCF alleles before analysis. iVar TSV SNVs are already position-based."
            checked={config.normalizeAlleles}
            onChange={(v) => update("normalizeAlleles", v)}
          />
          <ToggleField
            label="Split multiallelic"
            tip="Decompose multi-allelic VCF records into separate biallelic entries. iVar TSV inputs already use one row per allele."
            checked={config.splitMultiallelic}
            onChange={(v) => update("splitMultiallelic", v)}
          />
          <ToggleField
            label="Strict mode"
            tip="Fail on ambiguous or problematic records instead of silently skipping them. Recommended for production."
            checked={config.strict}
            onChange={(v) => update("strict", v)}
          />
          <ToggleField
            label="Emit filtered"
            tip="Include filtered VCF records with FILTER tags such as LowSupport or LowFrequency. Useful for auditing which variants were removed."
            checked={config.emitFiltered}
            onChange={(v) => update("emitFiltered", v)}
          />
          <ToggleField
            label="Strand bias INFO"
            tip="Add strand bias statistics (SB, FS, SOR) to VCF INFO fields for downstream filtering."
            checked={config.strandBiasInfo}
            onChange={(v) => update("strandBiasInfo", v)}
          />
          <ToggleField
            label="Keep original INFO"
            tip="Preserve original INFO fields from VCF input in the output VCF (e.g. SnpEff ANN, CSQ). Requires VCF output."
            checked={config.keepOriginalInfo}
            onChange={(v) => update("keepOriginalInfo", v)}
          />
          <ToggleField
            label="Exclude intergenic"
            tip="Exclude intergenic SNPs (variants outside annotated genes) from the output. By default, intergenic variants are included."
            checked={config.excludeIntergenic}
            onChange={(v) => update("excludeIntergenic", v)}
          />

          <div className="param-divider" />

          <div className="param-text-row">
            <div className="param-slider-label">
              <span>Genetic code</span>
              <Tip text="NCBI translation table for codon-to-amino-acid translation. Table 11 (Bacterial/Archaeal) is standard for prokaryotes. Use table 1 for nuclear eukaryotic genes or mitochondrial tables (2–5) for organellar genomes." />
            </div>
            <select
              className="param-select"
              value={config.translationTable}
              onChange={(e) => update("translationTable", Number(e.target.value))}
            >
              <option value={11}>11 — Bacterial, Archaeal & Plant Plastid</option>
              <option value={1}>1 — Standard (Eukaryotic Nuclear)</option>
              <option value={2}>2 — Vertebrate Mitochondrial</option>
              <option value={3}>3 — Yeast Mitochondrial</option>
              <option value={4}>4 — Mold, Protozoan & Coelenterate Mito</option>
              <option value={5}>5 — Invertebrate Mitochondrial</option>
              <option value={6}>6 — Ciliate, Dasycladacean & Hexamita Nuclear</option>
              <option value={12}>12 — Alternative Yeast Nuclear</option>
              <option value={25}>25 — SR1 & Gracilibacteria</option>
            </select>
          </div>
        </ParamGroup>

        <ParamGroup title="Output & Performance" accent="#E7216A">
          <ToggleField
            label="TSV output"
            tip="Generate a tab-separated file with annotated variants, amino acid changes, and MNV classifications."
            checked={config.outputTsv}
            onChange={(v) => {
              // Prevent disabling both outputs — at least one must be active
              if (!v && !config.outputVcf) return;
              update("outputTsv", v);
            }}
          />
          <ToggleField
            label="VCF output"
            tip="Generate a VCF file with MNV annotations in INFO fields. Compatible with downstream VCF tools."
            checked={config.outputVcf}
            onChange={(v) => {
              // Prevent disabling both outputs — at least one must be active
              if (!v && !config.outputTsv) return;
              if (!v && config.vcfGz) {
                // Disable both VCF and VCF-GZ in one update
                onChange({ ...config, outputVcf: false, vcfGz: false });
              } else {
                update("outputVcf", v);
              }
            }}
          />
          {config.outputVcf && (
            <ToggleField
              label="Compress VCF (BGZF)"
              tip="Write BGZF-compressed .vcf.gz output with Tabix index. Recommended for large datasets and IGV compatibility."
              checked={config.vcfGz}
              onChange={(v) => update("vcfGz", v)}
            />
          )}

          <div className="param-divider" />

          {/* ── Output location ── */}
          <div className="param-text-row">
            <div className="param-slider-label">
              <span>Output directory</span>
              <Tip text="Directory where output files will be written. Leave empty to use the same directory as the variant input file." />
            </div>
            <div className="output-dir-row">
              <input
                type="text"
                className="param-text-input output-dir-input"
                value={config.outputDir ?? ""}
                onChange={(e) =>
                  update("outputDir", e.target.value || undefined)
                }
                placeholder={
                  config.vcfFile
                    ? config.vcfFile.split(/[\\/]/).slice(0, -1).join("/") || "/"
                    : "Same as variant file"
                }
              />
              <button
                type="button"
                className="output-dir-browse"
                onClick={async () => {
                  const dir = await open({ directory: true, multiple: false });
                  if (dir && typeof dir === "string") {
                    update("outputDir", dir);
                  }
                }}
                title="Browse for directory"
              >
                <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
                  <path d="M22 19a2 2 0 01-2 2H4a2 2 0 01-2-2V5a2 2 0 012-2h5l2 3h9a2 2 0 012 2z" />
                </svg>
              </button>
            </div>
          </div>

          <div className="param-text-row">
            <div className="param-slider-label">
              <span>Output prefix</span>
              <Tip text="Filename prefix for output files (e.g. 'my_run' produces 'my_run.MNV.tsv'). Leave empty to derive from the variant input filename." />
            </div>
            <input
              type="text"
              className="param-text-input"
              value={config.outputPrefix ?? ""}
              onChange={(e) =>
                update("outputPrefix", e.target.value || undefined)
              }
              placeholder={
                config.vcfFile
                  ? config.vcfFile.split(/[\\/]/).pop()?.replace(/\.(vcf\.gz|vcf|tsv|txt)$/i, "") ?? "auto"
                  : "Derived from variants"
              }
            />
          </div>

          <div className="param-divider" />

          <div className="param-slider-row">
            <div className="param-slider-label">
              <span>Threads</span>
              <Tip text="Number of parallel threads for contig processing. Leave at 0 to auto-detect all available CPU cores." />
            </div>
            <div className="param-slider-control">
              <input
                type="range"
                className="param-range"
                min={0}
                max={32}
                step={1}
                value={config.threads ?? 0}
                onChange={(e) => {
                  const v = Number(e.target.value);
                  update("threads", v === 0 ? undefined : v);
                }}
                style={{ "--pct": `${((config.threads ?? 0) / 32) * 100}%` } as React.CSSProperties}
              />
              <input
                type="text"
                inputMode="numeric"
                className="param-number"
                value={config.threads ?? 0}
                onChange={(e) => {
                  const v = parseInt(e.target.value, 10);
                  if (!isNaN(v)) update("threads", v === 0 ? undefined : Math.max(0, Math.min(32, v)));
                }}
                placeholder="auto"
              />
            </div>
            {(config.threads ?? 0) === 0 && (
              <span className="param-auto-badge">auto</span>
            )}
          </div>
        </ParamGroup>
      </div>
    </div>
  );
}
