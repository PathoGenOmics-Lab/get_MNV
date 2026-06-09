// Standalone demo entry: renders the real BamViewer component with a real
// get_bam_view fixture (scenario 02a: chr_test MNV 28 G>T / 30 T>A, GCT->TCA
// Ala10Ser, 20 supporting reads). Used to verify the read visualization renders
// well without launching the full Tauri desktop app.
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import BamViewer from "../src/components/BamViewer";
import type { TsvData } from "../src/types";
import sampleTsv from "./sample_tsv.json";
import "../src/index.css";
import "../src/App.css";

const root = document.getElementById("root");
if (!root) throw new Error("missing #root");

// Demo aid: once the (mocked) reads have rendered, scroll the reference row to
// the top so a fresh-page screenshot frames the read pileup.
if (new URLSearchParams(location.search).get("reads") === "1") {
  const tryScroll = () => {
    const ref = document.querySelector(".bam-grid-row--reference");
    if (ref) {
      ref.scrollIntoView({ block: "start" });
      window.scrollBy(0, -8);
    } else {
      setTimeout(tryScroll, 100);
    }
  };
  setTimeout(tryScroll, 300);
}

createRoot(root).render(
  <StrictMode>
    <div className="app" style={{ padding: 16 }}>
      <h2 style={{ margin: "4px 0 12px" }}>get_MNV — Genomic Track Viewer (demo)</h2>
      <BamViewer
        bamPath="scenario02a/reads.bam"
        fastaPath="scenario02a/ref.fasta"
        data={sampleTsv as TsvData}
        minMapq={0}
        minBaseQuality={20}
      />
    </div>
  </StrictMode>,
);
