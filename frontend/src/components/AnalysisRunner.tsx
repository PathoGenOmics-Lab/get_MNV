import { useState, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import type { AnalysisConfig, RunSummary } from "../types";

interface ProgressEvent {
  phase: string;
  contig: string | null;
  current: number;
  total: number;
}

interface AnalysisRunnerProps {
  config: AnalysisConfig;
  running: boolean;
  onStart: () => void;
  onComplete: (result: RunSummary) => void;
  onError: (error: string) => void;
}

function PlayIcon() {
  return (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="currentColor">
      <path d="M8 5.14v14l11-7-11-7z" />
    </svg>
  );
}

function SpinnerIcon() {
  return (
    <svg width="18" height="18" viewBox="0 0 24 24" fill="none" strokeWidth="2.5" strokeLinecap="round" style={{ animation: "spin 0.8s linear infinite" }}>
      {/* MNV tricolor spinner arc */}
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
  );
}

export default function AnalysisRunner({
  config,
  running,
  onStart,
  onComplete,
  onError,
}: AnalysisRunnerProps) {
  const [progress, setProgress] = useState<ProgressEvent | null>(null);

  const canRun =
    config.vcfFile.length > 0 &&
    config.fastaFile.length > 0 &&
    config.genesFile.length > 0;

  // Listen for progress events from the Rust backend
  useEffect(() => {
    if (!running) return;
    const unlistenPromise = listen<ProgressEvent>("analysis-progress", (event) => {
      setProgress(event.payload);
    });
    return () => {
      unlistenPromise.then((fn) => fn());
    };
  }, [running]);

  const handleRun = async () => {
    setProgress(null);
    onStart();
    try {
      const result = await invoke<RunSummary>("run_analysis", { config });
      onComplete(result);
    } catch (err) {
      onError(String(err));
    }
  };

  // Determine progress bar state
  const hasDeterminateProgress = progress && progress.total > 0;
  const progressPct = hasDeterminateProgress
    ? Math.round((progress.current / progress.total) * 100)
    : 0;

  // Build status text
  let statusText = "Parsing inputs\u2026";
  if (progress) {
    if (progress.phase === "parsing") {
      statusText = "Parsing inputs\u2026";
    } else if (progress.phase === "processing" && progress.contig) {
      statusText = `Processing ${progress.contig}\u2026 (${progress.current}/${progress.total})`;
    } else if (progress.phase === "complete") {
      statusText = "Finalizing output\u2026";
    }
  }

  return (
    <div className="analysis-runner">
      <button
        className="run-button"
        disabled={!canRun || running}
        onClick={handleRun}
      >
        {running ? (
          <>
            <SpinnerIcon />
            Translating codons...
          </>
        ) : (
          <>
            <PlayIcon />
            Run Analysis
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
          <p className="progress-status">{statusText}</p>
        </>
      )}
      {!canRun && !running && (
        <p className="hint">
          Select VCF, FASTA, and gene annotation files to enable analysis
        </p>
      )}
    </div>
  );
}
