import { open } from "@tauri-apps/plugin-shell";

interface ResultsTableProps {
  outputTsv?: string;
  outputVcf?: string;
}

/** Extract just the filename from a full path */
function basename(path: string): string {
  return path.split(/[\\/]/).pop() ?? path;
}

/** Extract the directory from a full path */
function dirname(path: string): string {
  const parts = path.split(/[\\/]/);
  parts.pop();
  return parts.join("/") || "/";
}

/** Open the containing folder in the system file manager */
async function revealInFinder(path: string) {
  try {
    const dir = dirname(path);
    await open(dir);
  } catch {
    // Silently ignore if shell plugin isn't available
  }
}

function TsvIcon() {
  return (
    <svg width="28" height="28" viewBox="0 0 32 32" fill="none">
      <rect x="3" y="2" width="26" height="28" rx="3" fill="#e6f7f5" stroke="#149389" strokeWidth="1.5" />
      <line x1="10" y1="12" x2="22" y2="12" stroke="#149389" strokeWidth="1.5" strokeLinecap="round" />
      <line x1="10" y1="17" x2="22" y2="17" stroke="#149389" strokeWidth="1.5" strokeLinecap="round" opacity="0.6" />
      <line x1="10" y1="22" x2="18" y2="22" stroke="#149389" strokeWidth="1.5" strokeLinecap="round" opacity="0.35" />
    </svg>
  );
}

function VcfIcon() {
  return (
    <svg width="28" height="28" viewBox="0 0 32 32" fill="none">
      <rect x="3" y="2" width="26" height="28" rx="3" fill="#fef2f2" stroke="#E52421" strokeWidth="1.5" />
      <path d="M10 13l3 7 3-7" stroke="#E52421" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" />
      <circle cx="22" cy="13" r="2" fill="#E52421" opacity="0.4" />
      <line x1="10" y1="23" x2="22" y2="23" stroke="#E52421" strokeWidth="1.5" strokeLinecap="round" opacity="0.35" />
    </svg>
  );
}

function RevealIcon() {
  return (
    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
      <path d="M18 13v6a2 2 0 01-2 2H5a2 2 0 01-2-2V8a2 2 0 012-2h6" />
      <polyline points="15 3 21 3 21 9" />
      <line x1="10" y1="14" x2="21" y2="3" />
    </svg>
  );
}

export default function ResultsTable({
  outputTsv,
  outputVcf,
}: ResultsTableProps) {
  if (!outputTsv && !outputVcf) {
    return (
      <div className="results-table">
        <p className="hint">Run an analysis to see results here.</p>
      </div>
    );
  }

  return (
    <div className="results-table">
      <div className="output-header">
        <h4>Output Files</h4>
      </div>
      <div className="output-files">
        {outputTsv && (
          <div className="output-file-card output-file-card--tsv">
            <div className="output-file-icon">
              <TsvIcon />
            </div>
            <div className="output-file-info">
              <div className="output-file-top">
                <span className="output-file-badge output-file-badge--tsv">TSV</span>
                <span className="output-file-name">{basename(outputTsv)}</span>
              </div>
              <span className="output-file-dir">{dirname(outputTsv)}</span>
            </div>
            <button
              className="output-file-reveal"
              onClick={() => revealInFinder(outputTsv)}
              title="Open containing folder"
            >
              <RevealIcon />
              <span>Reveal</span>
            </button>
          </div>
        )}
        {outputVcf && (
          <div className="output-file-card output-file-card--vcf">
            <div className="output-file-icon">
              <VcfIcon />
            </div>
            <div className="output-file-info">
              <div className="output-file-top">
                <span className="output-file-badge output-file-badge--vcf">VCF</span>
                <span className="output-file-name">{basename(outputVcf)}</span>
              </div>
              <span className="output-file-dir">{dirname(outputVcf)}</span>
            </div>
            <button
              className="output-file-reveal"
              onClick={() => revealInFinder(outputVcf)}
              title="Open containing folder"
            >
              <RevealIcon />
              <span>Reveal</span>
            </button>
          </div>
        )}
      </div>
    </div>
  );
}
