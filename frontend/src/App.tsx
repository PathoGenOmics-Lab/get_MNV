import { useState, useEffect, useCallback, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { getCurrentWebview } from "@tauri-apps/api/webview";
import { listen } from "@tauri-apps/api/event";
import { open } from "@tauri-apps/plugin-shell";
import "./App.css";
import FileSelector from "./components/FileSelector";
import ParameterForm from "./components/ParameterForm";
import Summary from "./components/Summary";
import ResultsTable from "./components/ResultsTable";
import VariantTable from "./components/VariantTable";
import BamViewer from "./components/BamViewer";
import type { AnalysisConfig, RunSummary, TsvData, SampleEntry } from "./types";
import { DEFAULT_CONFIG } from "./types";

type Tab = "analysis" | "results";

/** Map file extension to config field */
const EXT_MAP: Record<string, keyof AnalysisConfig> = {
  vcf: "vcfFile",
  "vcf.gz": "vcfFile",
  fasta: "fastaFile",
  fas: "fastaFile",
  fa: "fastaFile",
  fna: "fastaFile",
  gff: "genesFile",
  gff3: "genesFile",
  gtf: "genesFile",
  bed: "genesFile",
  tsv: "genesFile",
  csv: "genesFile",
  txt: "genesFile",
  bam: "bamFile",
};

/** Keyword patterns in filenames → config field (fallback when extension fails) */
const NAME_HINTS: [RegExp, keyof AnalysisConfig][] = [
  [/vcf|variant|snp/i, "vcfFile"],
  [/\.fa$|fasta|refseq|reference|genome/i, "fastaFile"],
  [/gene|annot|gff|gtf|feature|cds/i, "genesFile"],
  [/\.bam$|align|mapping/i, "bamFile"],
];

function classifyFile(path: string): keyof AnalysisConfig | null {
  const lower = path.toLowerCase();

  // 1. Compound extensions first
  if (lower.endsWith(".vcf.gz")) return "vcfFile";
  if (lower.endsWith(".gff.gz") || lower.endsWith(".gtf.gz")) return "genesFile";

  // 2. Simple extension lookup
  const ext = lower.split(".").pop() ?? "";
  const byExt = EXT_MAP[ext];
  if (byExt) return byExt;

  // 3. Keyword fallback on filename (handles files without extensions)
  const filename = lower.split(/[\\/]/).pop() ?? "";
  for (const [pattern, field] of NAME_HINTS) {
    if (pattern.test(filename)) return field;
  }

  return null;
}

/** Check if the annotation file is GFF/GFF3/GTF format (vs plain genes TSV) */
function isGffFile(path: string): boolean {
  const lower = path.toLowerCase();
  return /\.(gff3?|gtf)(\.gz)?$/.test(lower);
}

/** Extract filename stem — removes path and known bioinformatics extensions */
function filenameStem(path: string): string {
  return (path.split(/[\\/]/).pop() ?? "").replace(/\.(vcf\.gz|vcf|bam)$/i, "");
}

/** Extract the core sample ID (first dot-segment) from a filename */
function sampleIdFromPath(path: string): string {
  return filenameStem(path).split(".")[0];
}

/**
 * Match a BAM file to a VCF using multi-level fallback:
 * 1. Exact stem match (e.g., "sample1.vcf" ↔ "sample1.bam")
 * 2. Prefix match (e.g., "G35894.var.snp.vcf" ↔ "G35894.bam" — one stem starts with the other)
 * 3. Sample ID match (e.g., "MIP00022.MTB_anc.ann.vcf" ↔ "MIP00022.MTB_anc.final.bam" — first segment)
 */
function matchBamToVcf(vcfPath: string, bamPaths: string[]): string | undefined {
  if (bamPaths.length === 0) return undefined;
  const vcfStem = filenameStem(vcfPath);

  // 1. Exact stem
  const exact = bamPaths.find((b) => filenameStem(b) === vcfStem);
  if (exact) return exact;

  // 2. One stem is prefix of the other (handles "G35894.var.snp" vs "G35894")
  const prefix = bamPaths.find((b) => {
    const bamStem = filenameStem(b);
    return vcfStem.startsWith(bamStem + ".") || bamStem.startsWith(vcfStem + ".");
  });
  if (prefix) return prefix;

  // 3. First segment match (sample ID: "MIP00022" vs "MIP00022")
  const vcfId = sampleIdFromPath(vcfPath);
  if (vcfId) {
    const idMatch = bamPaths.find((b) => sampleIdFromPath(b) === vcfId);
    if (idMatch) return idMatch;
  }

  return undefined;
}

/** Codon background — each letter illuminates randomly with MNV colors */
const CODON_SEQ = "ATG CTC ATC CGT TCT ATG GCA CTG TCA GAT CTC ATG ATC CGT TCT GCA CTG";

// Three animation classes with different durations for pseudo-random feel
const MNV_CLASSES = ["codon-mnv-m", "codon-mnv-n", "codon-mnv-v"] as const;

// Simple hash-like function to create pseudo-random but deterministic delays
function pseudoRandom(i: number): number {
  return ((i * 7 + 3) * 13) % 37;
}

function CodonBackground() {
  const letters = CODON_SEQ.split("").map((ch, i) => {
    if (ch === " ") return <span key={i} className="codon-space">{"\u00A0"}</span>;
    // Cycle through MNV colors pseudo-randomly based on position
    const colorClass = MNV_CLASSES[pseudoRandom(i) % 3];
    // Pseudo-random delay between 0–12s so they don't sync
    const delay = (pseudoRandom(i) / 37) * 12;
    return (
      <span
        key={i}
        className={`codon-nt ${colorClass}`}
        style={{ animationDelay: `${delay.toFixed(1)}s` }}
      >
        {ch}
      </span>
    );
  });

  return <div className="codon-bg">{letters}{letters}</div>;
}

/** Animated crab mascot — uses native emoji for best look */
function CrabMascot() {
  return (
    <div className="crab-traverse-track">
      <span className="crab-emoji">🦀</span>
    </div>
  );
}

function DnaIcon() {
  return (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
      <path d="M2 15c6.667-6 13.333 0 20-6" />
      <path d="M9 22c1.798-1.998 2.518-3.995 2.807-5.993" />
      <path d="M15 2c-1.798 1.998-2.518 3.995-2.807 5.993" />
      <path d="M17 6H3" />
      <path d="M21 18H7" />
    </svg>
  );
}

function ChartIcon() {
  return (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
      <path d="M3 3v18h18" />
      <path d="M7 16l4-8 4 4 6-6" />
    </svg>
  );
}

function SunIcon() {
  return (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
      <circle cx="12" cy="12" r="5" />
      <path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42" />
    </svg>
  );
}

function MoonIcon() {
  return (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
      <path d="M21 12.79A9 9 0 1111.21 3 7 7 0 0021 12.79z" />
    </svg>
  );
}

function GitHubIcon() {
  return (
    <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
      <path d="M12 0C5.37 0 0 5.37 0 12c0 5.31 3.435 9.795 8.205 11.385.6.105.825-.255.825-.57 0-.285-.015-1.23-.015-2.235-3.015.555-3.795-.735-4.035-1.41-.135-.345-.72-1.41-1.23-1.695-.42-.225-1.02-.78-.015-.795.945-.015 1.62.87 1.845 1.23 1.08 1.815 2.805 1.305 3.495.99.105-.78.42-1.305.765-1.605-2.67-.3-5.46-1.335-5.46-5.925 0-1.305.465-2.385 1.23-3.225-.12-.3-.54-1.53.12-3.18 0 0 1.005-.315 3.3 1.23.96-.27 1.98-.405 3-.405s2.04.135 3 .405c2.295-1.56 3.3-1.23 3.3-1.23.66 1.65.24 2.88.12 3.18.765.84 1.23 1.905 1.23 3.225 0 4.605-2.805 5.625-5.475 5.925.435.375.81 1.095.81 2.22 0 1.605-.015 2.895-.015 3.3 0 .315.225.69.825.57A12.02 12.02 0 0024 12c0-6.63-5.37-12-12-12z" />
    </svg>
  );
}

function getInitialTheme(): "light" | "dark" {
  try {
    const stored = localStorage.getItem("get_mnv_theme");
    if (stored === "dark" || stored === "light") return stored;
  } catch { /* ignore */ }
  return "light";
}

function getInitialCrab(): boolean {
  try {
    const stored = localStorage.getItem("get_mnv_crab");
    if (stored === "false") return false;
  } catch { /* ignore */ }
  return true;
}

interface ProgressEvent {
  phase: string;
  contig: string | null;
  current: number;
  total: number;
}

function App() {
  const [tab, setTab] = useState<Tab>("analysis");
  const [config, setConfig] = useState<AnalysisConfig>(DEFAULT_CONFIG);
  const [running, setRunning] = useState(false);
  const [justFinished, setJustFinished] = useState(false);
  const [samples, setSamples] = useState<SampleEntry[]>([]);
  const [activeSampleId, setActiveSampleId] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [showValidation, setShowValidation] = useState(false);
  const [dragActive, setDragActive] = useState(false);
  const [justAssigned, setJustAssigned] = useState<Set<string>>(new Set());
  const [gffAvailableFeatures, setGffAvailableFeatures] = useState<string[]>([]);
  const [theme, setTheme] = useState<"light" | "dark">(getInitialTheme);
  const [showCrab, setShowCrab] = useState(getInitialCrab);
  const [appVersion, setAppVersion] = useState("v1.1.1");
  const [batchProgress, setBatchProgress] = useState<{ current: number; total: number } | null>(null);
  const [runProgress, setRunProgress] = useState<ProgressEvent | null>(null);

  const configRef = useRef(config);
  useEffect(() => { configRef.current = config; }, [config]);

  // Derived state
  const activeSample = samples.find((s) => s.id === activeSampleId) ?? null;
  const activeResult = activeSample?.result ?? null;
  const activeTsvData = activeSample?.tsvData ?? null;
  const anyDone = samples.some((s) => s.status === "done");

  // Sync theme to DOM and localStorage
  useEffect(() => {
    document.documentElement.setAttribute("data-theme", theme);
    try { localStorage.setItem("get_mnv_theme", theme); } catch { /* ignore */ }
  }, [theme]);

  const toggleTheme = useCallback(() => {
    setTheme((prev) => (prev === "light" ? "dark" : "light"));
  }, []);

  const toggleCrab = useCallback(() => {
    setShowCrab((prev) => {
      const next = !prev;
      try { localStorage.setItem("get_mnv_crab", String(next)); } catch { /* ignore */ }
      return next;
    });
  }, []);

  // Fetch version from Tauri backend
  useEffect(() => {
    invoke<string>("get_core_version")
      .then((v) => setAppVersion(`v${v}`))
      .catch(() => { /* keep fallback */ });
  }, []);

  // Measure header + tabs and set CSS variables for fixed layout
  const headerRef = useRef<HTMLElement>(null);
  const tabsRef = useRef<HTMLElement>(null);
  useEffect(() => {
    const measure = () => {
      const hh = headerRef.current?.offsetHeight ?? 52;
      const th = tabsRef.current?.offsetHeight ?? 38;
      document.documentElement.style.setProperty("--header-height", `${hh}px`);
      document.documentElement.style.setProperty("--tabs-height", `${th}px`);
    };
    measure();
    window.addEventListener("resize", measure);
    return () => window.removeEventListener("resize", measure);
  }, []);

  const genesIsGff = isGffFile(config.genesFile);

  const filesReady =
    samples.length > 0 &&
    config.fastaFile.length > 0 &&
    config.genesFile.length > 0;

  // Clear validation errors once all required fields are filled
  useEffect(() => {
    if (filesReady) setShowValidation(false);
  }, [filesReady]);

  const fileCount = [
    samples.length > 0 ? "vcf" : "",
    config.fastaFile,
    config.genesFile,
  ].filter((f) => f.length > 0).length;

  const flashAssigned = useCallback((fields: string[]) => {
    setJustAssigned(new Set(fields));
    setTimeout(() => setJustAssigned(new Set()), 800);
  }, []);

  // VCF FileSelector change handler — single VCF
  const handleVcfChange = useCallback((path: string) => {
    if (path) {
      setSamples([{
        id: crypto.randomUUID(),
        name: sampleIdFromPath(path) || filenameStem(path),
        vcfPath: path,
        status: "pending",
      }]);
    } else {
      setSamples([]);
    }
    setConfig((prev) => ({ ...prev, vcfFile: path }));
  }, []);

  // Remove a sample from the list
  const removeSample = useCallback((id: string) => {
    setSamples((prev) => {
      const next = prev.filter((s) => s.id !== id);
      // Sync config.vcfFile with first remaining sample
      if (next.length > 0) {
        setConfig((c) => ({ ...c, vcfFile: next[0].vcfPath }));
      } else {
        setConfig((c) => ({ ...c, vcfFile: "" }));
      }
      return next;
    });
    setActiveSampleId((prev) => prev === id ? null : prev);
  }, []);

  // Auto-detect GFF feature types when a GFF file is loaded
  useEffect(() => {
    let cancelled = false;

    if (!config.genesFile || !isGffFile(config.genesFile)) {
      Promise.resolve().then(() => { if (!cancelled) setGffAvailableFeatures([]); });
    } else {
      invoke<string[]>("get_gff_features", { path: config.genesFile })
        .then((features) => { if (!cancelled) setGffAvailableFeatures(features); })
        .catch(() => { if (!cancelled) setGffAvailableFeatures([]); });
    }

    return () => { cancelled = true; };
  }, [config.genesFile]);

  // Global Tauri drag-drop listener — supports multiple VCFs
  useEffect(() => {
    let unlisten: (() => void) | undefined;

    getCurrentWebview()
      .onDragDropEvent((event) => {
        const payload = event.payload;
        if (payload.type === "enter") {
          setDragActive(true);
        } else if (payload.type === "leave") {
          setDragActive(false);
        } else if (payload.type === "drop") {
          setDragActive(false);

          const vcfPaths: string[] = [];
          const bamPaths: string[] = [];
          const updates: Partial<AnalysisConfig> = {};
          const assigned: string[] = [];

          for (const path of payload.paths) {
            const field = classifyFile(path);
            if (field === "vcfFile") {
              vcfPaths.push(path);
            } else if (field === "bamFile") {
              bamPaths.push(path);
            } else if (field) {
              (updates as Record<string, unknown>)[field] = path;
              assigned.push(field);
            }
          }

          // Build sample entries from VCFs (deduplicated)
          if (vcfPaths.length > 0) {
            // Deduplicate within the drop itself
            const uniqueVcfPaths = [...new Set(vcfPaths)];

            setSamples((prev) => {
              // Filter out VCFs that already exist as samples
              const existingPaths = new Set(prev.map((s) => s.vcfPath));
              const freshVcfPaths = uniqueVcfPaths.filter((p) => !existingPaths.has(p));
              if (freshVcfPaths.length === 0) return prev; // nothing new

              const usedBams = new Set<string>();
              const newSamples: SampleEntry[] = freshVcfPaths.map((vcf) => {
                // Match BAM using multi-level fallback (exact stem → prefix → sample ID)
                const availableBams = bamPaths.filter((b) => !usedBams.has(b));
                const matchedBam = matchBamToVcf(vcf, availableBams);
                if (matchedBam) usedBams.add(matchedBam);
                return {
                  id: crypto.randomUUID(),
                  name: sampleIdFromPath(vcf) || filenameStem(vcf),
                  vcfPath: vcf,
                  bamPath: matchedBam,
                  status: "pending" as const,
                };
              });

              // If there's only one BAM and no stem match happened, assign it to all new samples
              if (bamPaths.length === 1 && !newSamples.some((s) => s.bamPath)) {
                for (const s of newSamples) s.bamPath = bamPaths[0];
              }

              // Append to existing done/running samples
              const keepPrev = prev.filter((s) => s.status === "done" || s.status === "running");
              return [...keepPrev, ...newSamples];
            });
            updates.vcfFile = vcfPaths[0];
            assigned.push("vcfFile");

            // If a single BAM was dropped alongside VCFs, also set config.bamFile
            if (bamPaths.length === 1) {
              updates.bamFile = bamPaths[0];
              assigned.push("bamFile");
            }
          } else if (bamPaths.length > 0) {
            // Only BAMs dropped (no VCFs) — assign to config + existing samples
            updates.bamFile = bamPaths[0];
            assigned.push("bamFile");
            if (bamPaths.length === 1) {
              setSamples((prev) => prev.map((s) => ({ ...s, bamPath: bamPaths[0] })));
            }
          }

          if (assigned.length > 0) {
            setConfig((prev) => ({ ...prev, ...updates }));
            flashAssigned(assigned);
            // Auto-generate .fai when FASTA is assigned via drag-and-drop
            if (updates.fastaFile) {
              invoke("ensure_fasta_index", { fastaPath: updates.fastaFile }).catch(() => {});
            }
          }
        }
      })
      .then((fn) => {
        unlisten = fn;
      });

    return () => {
      unlisten?.();
    };
  }, [flashAssigned]);

  // ── Batch analysis runner ──
  const handleRunAll = useCallback(async () => {
    const toRun = samples.filter((s) => s.status === "pending" || s.status === "error");
    if (toRun.length === 0) return;

    setRunning(true);
    setError(null);
    setBatchProgress({ current: 0, total: toRun.length });

    const unlisten = await listen<ProgressEvent>("analysis-progress", (event) => {
      setRunProgress(event.payload);
    });

    try {
      for (let i = 0; i < toRun.length; i++) {
        const sample = toRun[i];
        setBatchProgress({ current: i + 1, total: toRun.length });
        setRunProgress(null);

        setSamples((prev) =>
          prev.map((s) => (s.id === sample.id ? { ...s, status: "running" } : s))
        );

        const sampleConfig: AnalysisConfig = {
          ...configRef.current,
          vcfFile: sample.vcfPath,
          bamFile: sample.bamPath ?? configRef.current.bamFile,
        };

        try {
          const result = await invoke<RunSummary>("run_analysis", { config: sampleConfig });

          let tsvData: TsvData | undefined;
          if (result.output_tsv) {
            try {
              tsvData = await invoke<TsvData>("read_tsv_file", { path: result.output_tsv });
            } catch { /* ignore */ }
          }

          setSamples((prev) =>
            prev.map((s) =>
              s.id === sample.id
                ? { ...s, status: "done" as const, result, tsvData, runConfig: sampleConfig }
                : s
            )
          );

          // Auto-activate first completed sample
          setActiveSampleId((prev) => prev ?? sample.id);
        } catch (err) {
          setSamples((prev) =>
            prev.map((s) =>
              s.id === sample.id ? { ...s, status: "error", error: String(err) } : s
            )
          );
        }
      }
    } finally {
      unlisten();
      setRunning(false);
      setBatchProgress(null);
      setRunProgress(null);

      // Brief success flash on the run button
      setJustFinished(true);
      setTimeout(() => setJustFinished(false), 1800);

      // Switch to results if any completed
      setSamples((prev) => {
        const firstDone = prev.find((s) => s.status === "done");
        if (firstDone) {
          setActiveSampleId((id) => id ?? firstDone.id);
          setTab("results");
        }
        return prev;
      });
    }
  }, [samples]);

  // Re-run a single sample (e.g., after error or completed)
  const handleRerunSample = useCallback((sampleId: string) => {
    setSamples((prev) =>
      prev.map((s) => (s.id === sampleId ? { ...s, status: "pending", error: undefined, result: undefined, tsvData: undefined } : s))
    );
  }, []);

  // Re-run all completed/error samples — resets them to pending, then auto-triggers batch run
  const rerunRequestedRef = useRef(false);

  const handleRerunAll = useCallback(() => {
    setSamples((prev) =>
      prev.map((s) =>
        s.status === "done" || s.status === "error"
          ? { ...s, status: "pending" as const, error: undefined, result: undefined, tsvData: undefined }
          : s
      )
    );
    rerunRequestedRef.current = true;
  }, []);

  // Auto-trigger batch run after re-run resets samples to pending
  useEffect(() => {
    if (rerunRequestedRef.current && !running && samples.some((s) => s.status === "pending")) {
      rerunRequestedRef.current = false;
      handleRunAll();
    }
  }, [samples, running, handleRunAll]);

  // Clear all samples
  const clearAllSamples = useCallback(() => {
    setSamples([]);
    setActiveSampleId(null);
    setConfig((prev) => ({ ...prev, vcfFile: "", bamFile: undefined }));
  }, []);

  // Build progress status text
  const progressStatusText = (() => {
    if (!running) return "";
    const parts: string[] = [];
    if (batchProgress && batchProgress.total > 1) {
      const currentSample = samples.find((s) => s.status === "running");
      parts.push(`Sample ${batchProgress.current}/${batchProgress.total}${currentSample ? ` — ${currentSample.name}` : ""}`);
    }
    if (runProgress) {
      if (runProgress.phase === "parsing") {
        parts.push("Parsing inputs\u2026");
      } else if (runProgress.phase === "processing" && runProgress.contig) {
        parts.push(`Processing ${runProgress.contig}\u2026 (${runProgress.current}/${runProgress.total})`);
      } else if (runProgress.phase === "complete") {
        parts.push("Finalizing output\u2026");
      }
    } else {
      parts.push("Starting\u2026");
    }
    return parts.join(" · ");
  })();

  const hasDeterminateProgress = runProgress && runProgress.total > 0;
  const progressPct = hasDeterminateProgress
    ? Math.round((runProgress.current / runProgress.total) * 100)
    : 0;

  // Count pending + error samples (runnable)
  const runnableCount = samples.filter((s) => s.status === "pending" || s.status === "error").length;

  // Ctrl+Enter / Cmd+Enter shortcut to run analysis
  // Use a ref to always capture the latest handleRunAll without re-registering the listener
  const handleRunAllRef = useRef(handleRunAll);
  handleRunAllRef.current = handleRunAll;
  const runStateRef = useRef({ running, filesReady, runnableCount });
  runStateRef.current = { running, filesReady, runnableCount };

  useEffect(() => {
    const handler = (e: KeyboardEvent) => {
      if ((e.ctrlKey || e.metaKey) && e.key === "Enter") {
        e.preventDefault();
        const { running: r, filesReady: f, runnableCount: rc } = runStateRef.current;
        if (!r && f && rc > 0) {
          handleRunAllRef.current();
        }
      }
    };
    document.addEventListener("keydown", handler);
    return () => document.removeEventListener("keydown", handler);
  }, []);

  return (
    <div className="app">
      <header className="app-header" ref={headerRef}>
        <CodonBackground />
        {showCrab && <CrabMascot />}
        <div className="app-header-left">
          <div className="app-header-title">
            <span className="header-get">get</span>
            <span className="header-underscore">_</span>
            <span className="header-m">M</span>
            <span className="header-n">N</span>
            <span className="header-v-letter">V</span>
          </div>
          <span className="app-version">{appVersion}</span>
        </div>
        <div className="app-header-right">
          <p className="app-subtitle">Multi-Nucleotide Variant Analysis</p>
          <button
            className="github-link"
            onClick={() => open("https://github.com/PathoGenOmics-Lab/get_MNV")}
            title="View on GitHub"
            aria-label="Open GitHub repository"
          >
            <GitHubIcon />
          </button>
          <button
            className="crab-toggle"
            onClick={toggleCrab}
            title={showCrab ? "Hide crab mascot" : "Show crab mascot"}
            aria-label="Toggle crab mascot"
          >
            🦀
          </button>
          <button
            className="theme-toggle"
            onClick={toggleTheme}
            title={theme === "light" ? "Switch to dark mode" : "Switch to light mode"}
            aria-label="Toggle dark mode"
          >
            {theme === "light" ? <MoonIcon /> : <SunIcon />}
          </button>
        </div>
      </header>

      <nav className="tabs" ref={tabsRef}>
        <button
          className={tab === "analysis" ? "tab active" : "tab"}
          onClick={() => setTab("analysis")}
        >
          <span className="tab-icon"><DnaIcon /></span>
          Analysis
        </button>
        <button
          className={tab === "results" ? "tab active" : "tab"}
          onClick={() => setTab("results")}
          disabled={!anyDone}
        >
          <span className="tab-icon"><ChartIcon /></span>
          Results
          {anyDone && activeTsvData && (
            <span className="tab-badge">{activeTsvData.rows.length.toLocaleString("en-US")}</span>
          )}
        </button>
      </nav>

      <main className="content">
        {tab === "analysis" && (
          <div className="analysis-layout">
            {/* ══════ LEFT: Files + Run ══════ */}
            <div className="analysis-main">
              <section className="step-section">
                <div className="step-header">
                  <div className={`step-badge step-badge-1${filesReady ? " step-badge--done" : ""}`}>
                    {filesReady ? "\u2713" : "1"}
                  </div>
                  <h3 className="step-title">Input Files</h3>
                  <span className="step-subtitle">
                    {fileCount}/3 required{config.bamFile ? " + BAM" : ""}
                    {samples.length > 1 ? ` · ${samples.length} samples` : ""}
                  </span>
                  {(fileCount > 0 || samples.length > 0) && (
                    <button
                      className="clear-all-btn"
                      onClick={() => {
                        clearAllSamples();
                        setConfig((prev) => ({
                          ...prev,
                          vcfFile: "",
                          fastaFile: "",
                          genesFile: "",
                          bamFile: undefined,
                        }));
                      }}
                      title="Clear all input files"
                    >
                      <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                        <path d="M3 6h18" />
                        <path d="M8 6V4a2 2 0 012-2h4a2 2 0 012 2v2" />
                        <path d="M19 6l-1 14a2 2 0 01-2 2H8a2 2 0 01-2-2L5 6" />
                      </svg>
                      Clear all
                    </button>
                  )}
                </div>
                <div className="files-grid">
                  <FileSelector
                    label={samples.length > 1 ? `VCF variants (${samples.length})` : "VCF variants"}
                    value={samples.length === 1 ? samples[0].vcfPath : samples.length > 1 ? `${samples.length} samples loaded` : ""}
                    onChange={handleVcfChange}
                    filters={[{ name: "VCF", extensions: ["vcf", "vcf.gz"] }]}
                    required
                    isDragActive={dragActive}
                    justAssigned={justAssigned.has("vcfFile")}
                    fileType="VCF"
                    validationError={showValidation && samples.length === 0}
                  />
                  <FileSelector
                    label="FASTA reference"
                    value={config.fastaFile}
                    onChange={(path) => {
                      setConfig({ ...config, fastaFile: path });
                      if (path) {
                        invoke("ensure_fasta_index", { fastaPath: path }).catch(() => {});
                      }
                    }}
                    filters={[{ name: "FASTA", extensions: ["fasta", "fas", "fa", "fna"] }]}
                    required
                    isDragActive={dragActive}
                    justAssigned={justAssigned.has("fastaFile")}
                    fileType="FASTA"
                    validationError={showValidation && !config.fastaFile}
                  />
                  <FileSelector
                    label="Gene annotation"
                    value={config.genesFile}
                    onChange={(path) => setConfig({ ...config, genesFile: path })}
                    filters={[{ name: "Annotation", extensions: ["gff", "gff3", "tsv", "txt"] }]}
                    required
                    isDragActive={dragActive}
                    justAssigned={justAssigned.has("genesFile")}
                    fileType="Annotation"
                    validationError={showValidation && !config.genesFile}
                  />
                  <FileSelector
                    label="BAM alignment"
                    value={config.bamFile ?? ""}
                    onChange={(path) =>
                      setConfig({
                        ...config,
                        bamFile: path.length > 0 ? path : undefined,
                      })
                    }
                    filters={[{ name: "BAM", extensions: ["bam"] }]}
                    isDragActive={dragActive}
                    justAssigned={justAssigned.has("bamFile")}
                    fileType="BAM"
                  />
                </div>

                {/* ── Sample list (when multiple VCFs) ── */}
                {samples.length > 1 && (
                  <div className="sample-list">
                    <div className="sample-list-header">
                      <span className="sample-list-title">Samples ({samples.length})</span>
                    </div>
                    {samples.map((s) => (
                      <div key={s.id} className={`sample-list-item sample-list-item--${s.status}`}>
                        <span className="sample-status-dot" />
                        <span className="sample-name" title={s.vcfPath}>{s.name}</span>
                        {s.bamPath && <span className="sample-bam-tag">BAM</span>}
                        {(s.status === "error" || s.status === "done") && (
                          <button
                            className="sample-retry-btn"
                            onClick={() => handleRerunSample(s.id)}
                            title={s.status === "error" ? "Retry this sample" : "Re-run this sample"}
                          >
                            ↻
                          </button>
                        )}
                        {s.status !== "running" && (
                          <button
                            className="sample-remove-btn"
                            onClick={() => removeSample(s.id)}
                            title="Remove sample"
                          >
                            ×
                          </button>
                        )}
                      </div>
                    ))}
                  </div>
                )}
              </section>

              {/* Run */}
              <section className="step-section">
                <div className="step-header">
                  <div className="step-badge step-badge-2">2</div>
                  <h3 className="step-title">Run Analysis</h3>
                </div>
                <div className="analysis-runner">
                  <button
                    className={`run-button${justFinished ? " run-button--success" : ""}`}
                    disabled={running}
                    onClick={() => {
                      if (!filesReady || samples.length === 0) {
                        setShowValidation(true);
                        return;
                      }
                      setShowValidation(false);
                      (runnableCount > 0 ? handleRunAll : handleRerunAll)();
                    }}
                  >
                    {justFinished ? (
                      <>
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
                          <path d="M20 6L9 17l-5-5" />
                        </svg>
                        Done!
                      </>
                    ) : running ? (
                      <>
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" strokeWidth="2.5" strokeLinecap="round" style={{ animation: "spin 0.8s linear infinite" }}>
                          <defs>
                            <linearGradient id="mnv-spin" x1="0" y1="0" x2="1" y2="1">
                              <stop offset="0%" stopColor="#149389" />
                              <stop offset="50%" stopColor="#E7216A" />
                              <stop offset="100%" stopColor="#6F398D" />
                            </linearGradient>
                          </defs>
                          <path d="M12 2a10 10 0 0 1 10 10" stroke="url(#mnv-spin)" />
                          <style>{`@keyframes spin { to { transform: rotate(360deg); } }`}</style>
                        </svg>
                        Translating codons...
                      </>
                    ) : runnableCount > 0 ? (
                      <>
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="currentColor">
                          <path d="M8 5.14v14l11-7-11-7z" />
                        </svg>
                        {samples.length > 1
                          ? `Run Analysis (${runnableCount} sample${runnableCount !== 1 ? "s" : ""})`
                          : "Run Analysis"}
                      </>
                    ) : (
                      <>
                        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
                          <path d="M1 4v6h6" />
                          <path d="M3.51 15a9 9 0 1 0 2.13-9.36L1 10" />
                        </svg>
                        {samples.length > 1
                          ? `Re-run Analysis (${samples.length} sample${samples.length !== 1 ? "s" : ""})`
                          : "Re-run Analysis"}
                      </>
                    )}
                  </button>
                  {running && (
                    <>
                      <div className="progress-bar">
                        {hasDeterminateProgress ? (
                          <div
                            className="progress-bar-fill"
                            style={{ width: `${progressPct}%` }}
                          />
                        ) : (
                          <div className="progress-bar-indeterminate" />
                        )}
                      </div>
                      <p className="progress-status">{progressStatusText}</p>
                      {batchProgress && batchProgress.total > 1 && (
                        <div className="batch-progress-bar">
                          <div
                            className="batch-progress-fill"
                            style={{ width: `${Math.round((batchProgress.current / batchProgress.total) * 100)}%` }}
                          />
                        </div>
                      )}
                    </>
                  )}
                  {!filesReady && !running && (
                    <p className={`hint${showValidation ? " hint--error" : ""}`}>
                      {showValidation
                        ? "⚠ Required files are missing — select them above"
                        : "Select VCF, FASTA, and gene annotation files to enable analysis"}
                    </p>
                  )}
                  {filesReady && !running && (
                    <p className="hint hint--shortcut">
                      <kbd>{navigator.platform?.includes("Mac") ? "⌘" : "Ctrl"}</kbd>+<kbd>Enter</kbd> to run
                    </p>
                  )}
                </div>
              </section>

              {error && (
                <div className="error-banner">
                  <strong>Error:</strong> {error}
                </div>
              )}
            </div>

            {/* ══════ RIGHT: Parameters sidebar ══════ */}
            <aside className="analysis-sidebar">
              <div className="sidebar-header">
                <div className="step-badge step-badge-2" style={{ width: 22, height: 22, fontSize: "0.6rem" }}>
                  <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round">
                    <path d="M12 15.5A3.5 3.5 0 1012 8.5a3.5 3.5 0 000 7z" />
                    <path d="M19.4 15a1.65 1.65 0 00.33 1.82l.06.06a2 2 0 01-2.83 2.83l-.06-.06a1.65 1.65 0 00-1.82-.33 1.65 1.65 0 00-1 1.51V21a2 2 0 01-4 0v-.09A1.65 1.65 0 009 19.4a1.65 1.65 0 00-1.82.33l-.06.06a2 2 0 01-2.83-2.83l.06-.06A1.65 1.65 0 004.68 15a1.65 1.65 0 00-1.51-1H3a2 2 0 010-4h.09A1.65 1.65 0 004.6 9a1.65 1.65 0 00-.33-1.82l-.06-.06a2 2 0 012.83-2.83l.06.06A1.65 1.65 0 009 4.68a1.65 1.65 0 001-1.51V3a2 2 0 014 0v.09a1.65 1.65 0 001 1.51 1.65 1.65 0 001.82-.33l.06-.06a2 2 0 012.83 2.83l-.06.06A1.65 1.65 0 0019.4 9a1.65 1.65 0 001.51 1H21a2 2 0 010 4h-.09a1.65 1.65 0 00-1.51 1z" />
                  </svg>
                </div>
                <h3 className="step-title">Parameters</h3>
                <button
                  className="reset-defaults-btn"
                  onClick={() =>
                    setConfig((prev) => ({
                      ...DEFAULT_CONFIG,
                      vcfFile: prev.vcfFile,
                      fastaFile: prev.fastaFile,
                      genesFile: prev.genesFile,
                      bamFile: prev.bamFile,
                    }))
                  }
                  title="Reset all parameters to defaults"
                >
                  <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                    <path d="M3 12a9 9 0 019-9 9.75 9.75 0 016.74 2.74L21 8" />
                    <path d="M21 3v5h-5" />
                    <path d="M21 12a9 9 0 01-9 9 9.75 9.75 0 01-6.74-2.74L3 16" />
                    <path d="M3 21v-5h5" />
                  </svg>
                  Defaults
                </button>
              </div>
              <ParameterForm
                config={config}
                onChange={setConfig}
                isGff={genesIsGff}
                availableFeatures={gffAvailableFeatures}
              />
            </aside>
          </div>
        )}

        {tab === "results" && anyDone && (
          <div className="results-panel">
            {/* ── Sample selector (when multiple samples) ── */}
            {samples.filter((s) => s.status === "done").length > 1 && (
              <div className="sample-selector">
                <span className="sample-selector__label">Sample:</span>
                {samples
                  .filter((s) => s.status === "done")
                  .map((s) => (
                    <button
                      key={s.id}
                      className={`sample-tab${s.id === activeSampleId ? " sample-tab--active" : ""}`}
                      onClick={() => setActiveSampleId(s.id)}
                    >
                      {s.name}
                    </button>
                  ))}
              </div>
            )}

            {activeResult && (
              <>
                <Summary result={activeResult} />
                <ResultsTable
                  outputTsv={activeResult.output_tsv}
                  outputVcf={activeResult.output_vcf}
                  outputBcf={activeResult.output_bcf}
                />
              </>
            )}

            {activeSample && activeSample.status === "done" && activeTsvData && (
              <>
                {activeResult?.inputs.bam && activeSample.runConfig && (
                  <BamViewer
                    bamPath={activeResult.inputs.bam}
                    fastaPath={activeResult.inputs.fasta}
                    data={activeTsvData}
                    minMapq={activeSample.runConfig.minMapq}
                    minBaseQuality={activeSample.runConfig.minQuality}
                  />
                )}
                <section className="step-section">
                  <div className="step-header">
                    <div className="step-badge step-badge-2">
                      <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
                        <rect x="3" y="3" width="18" height="18" rx="2" />
                        <path d="M3 9h18M3 15h18M9 3v18M15 3v18" />
                      </svg>
                    </div>
                    <h3 className="step-title">Variant Data</h3>
                    <span className="step-subtitle">
                      {activeTsvData.rows.length.toLocaleString("en-US")} rows
                    </span>
                  </div>
                  <VariantTable data={activeTsvData} />
                </section>
              </>
            )}
          </div>
        )}
      </main>

      {/* Global drop overlay */}
      {dragActive && (
        <div className="drop-overlay">
          <div className="drop-overlay-content">
            <svg width="48" height="48" viewBox="0 0 48 48" fill="none">
              <path d="M24 32V8m0 0l-8 8m8-8l8 8" stroke="#EF781A" strokeWidth="3" strokeLinecap="round" strokeLinejoin="round" />
              <path d="M8 34v4a4 4 0 004 4h24a4 4 0 004-4v-4" stroke="#EF781A" strokeWidth="3" strokeLinecap="round" />
            </svg>
            <p>Drop files to auto-assign by type</p>
            <span className="drop-overlay-sub">VCF · FASTA · GFF/TSV · BAM</span>
          </div>
        </div>
      )}
    </div>
  );
}

export default App;
