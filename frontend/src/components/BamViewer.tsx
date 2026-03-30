import { useCallback, useEffect, useMemo, useRef, useState, type CSSProperties } from "react";
import { invoke } from "@tauri-apps/api/core";
import type { BamReadView, BamVariantSite, BamViewResponse, TsvData } from "../types";

interface BamViewerProps {
  bamPath: string;
  fastaPath: string;
  data: TsvData;
  minMapq: number;
  minBaseQuality: number;
}

interface ViewerLocus {
  id: string;
  chrom: string;
  gene: string;
  positions: number[];
  refBases: string[];
  altBases: string[];
  variantType: string;
  changeType: string;
  aaChanges: string;
  mnvReads: string;
  mnvFrequency: string;
  refCodon: string;
  mnvCodon: string;
  snpCodons: string;
  snpAaChanges: string;
}

interface WindowRange {
  start: number;
  end: number;
}

interface LaidOutLocus extends ViewerLocus {
  start: number;
  end: number;
  lane: number;
}

/* ── Constants ────────────────────────────────────────── */
const BUFFER_PADDING = 80; // bases fetched on each side of the locus
const VIEWER_MAX_READS = 80;
const DEFAULT_CELL_SIZE = 14;
const MIN_CELL_SIZE = 6;
const MAX_CELL_SIZE = 24;

/* ── Pure helpers ─────────────────────────────────────── */

function headerIndex(headers: string[], label: string): number {
  return headers.findIndex((h) => h.trim().toLowerCase() === label.trim().toLowerCase());
}

function splitField(value: string): string[] {
  return value.split(",").map((s) => s.trim()).filter((s) => s.length > 0);
}

/** Parse an AA-change string like "Gly92Asp" → { refAa: "Gly", pos: "92", altAa: "Asp" } */
function parseAaChange(aa: string): { refAa: string; pos: string; altAa: string } | null {
  const m = aa.match(/^([A-Z][a-z]{0,2}|\*)(\d+)([A-Z][a-z]{0,2}|\*)$/);
  if (!m) return null;
  return { refAa: m[1], pos: m[2], altAa: m[3] };
}

function parsePositions(value: string): number[] {
  return splitField(value)
    .map((s) => Number.parseInt(s, 10))
    .filter((n) => Number.isFinite(n) && n > 0);
}

function buildLoci(data: TsvData): ViewerLocus[] {
  const ci = headerIndex(data.headers, "Chromosome");
  const gi = headerIndex(data.headers, "Gene");
  const pi = headerIndex(data.headers, "Positions");
  const ri = headerIndex(data.headers, "Reference Bases");
  const ai = headerIndex(data.headers, "Base Changes");
  const vi = headerIndex(data.headers, "Variant Type");
  const cti = headerIndex(data.headers, "Change Type");
  const aai = headerIndex(data.headers, "AA Changes");
  const mri = headerIndex(data.headers, "MNV Reads");
  const mfi = headerIndex(data.headers, "MNV Frequencies");
  const rci = headerIndex(data.headers, "Reference Codon");
  const mci = headerIndex(data.headers, "MNV Codon");
  const sci = headerIndex(data.headers, "SNP Codon");
  const sai = headerIndex(data.headers, "SNP AA Changes");

  if ([ci, gi, pi, ri, ai, vi].some((i) => i < 0)) return [];

  return data.rows
    .map((row) => {
      const vt = row[vi] ?? "";
      if (!vt.toLowerCase().includes("mnv")) return null;

      // Skip loci with 0 MNV reads — nothing meaningful to visualise
      const mnvReadCount = mri >= 0 ? Number.parseInt(row[mri] ?? "0", 10) : 0;
      if (!mnvReadCount || mnvReadCount <= 0) return null;

      const positions = parsePositions(row[pi] ?? "");
      const refBases = splitField(row[ri] ?? "");
      const altBases = splitField(row[ai] ?? "");
      if (positions.length < 2 || positions.length !== refBases.length || positions.length !== altBases.length) return null;

      const chrom = row[ci] ?? "";
      return {
        id: `${chrom}:${positions.join("-")}:${altBases.join("")}`,
        chrom,
        gene: row[gi] ?? "Unknown",
        positions,
        refBases,
        altBases,
        variantType: vt,
        changeType: cti >= 0 ? row[cti] ?? "" : "",
        aaChanges: aai >= 0 ? row[aai] ?? "" : "",
        mnvReads: mri >= 0 ? row[mri] ?? "" : "",
        mnvFrequency: mfi >= 0 ? row[mfi] ?? "" : "",
        refCodon: rci >= 0 ? row[rci] ?? "" : "",
        mnvCodon: mci >= 0 ? row[mci] ?? "" : "",
        snpCodons: sci >= 0 ? row[sci] ?? "" : "",
        snpAaChanges: sai >= 0 ? row[sai] ?? "" : "",
      };
    })
    .filter((l): l is ViewerLocus => l !== null);
}

function locusBounds(locus: ViewerLocus): WindowRange {
  return { start: Math.min(...locus.positions), end: Math.max(...locus.positions) };
}

function countLabel(n: number): string {
  return n.toLocaleString("en-US");
}

function supportLabel(s: BamReadView["support"]): string {
  switch (s) {
    case "mnv": return "MNV";
    case "partial": return "Partial";
    case "reference": return "Ref";
    case "other": return "Other";
    default: return s;
  }
}

function locusTitle(l: ViewerLocus): string {
  return `${l.gene} · ${l.chrom}:${l.positions.join(", ")}`;
}

function focusSiteMap(sites: BamVariantSite[]): Map<number, BamVariantSite> {
  return new Map(sites.map((s) => [s.position, s]));
}

function nucleotideClass(base: string): string {
  switch (base.toUpperCase()) {
    case "A": return "bam-cell--nt-a";
    case "T": return "bam-cell--nt-t";
    case "G": return "bam-cell--nt-g";
    case "C": return "bam-cell--nt-c";
    default: return "";
  }
}

function tickStepForCellSize(cellSize: number): number {
  const rawStep = Math.ceil(55 / cellSize);
  return [1, 2, 5, 10, 15, 20, 25, 50].find((s) => s >= rawStep) ?? rawStep;
}

/** Reverse complement a short DNA string */
function reverseComplement(seq: string): string {
  const comp: Record<string, string> = { A: "T", T: "A", G: "C", C: "G", N: "N" };
  return seq.split("").reverse().map((b) => comp[b.toUpperCase()] ?? "N").join("");
}

/**
 * Find the 0-based index in the reference array where the codon starts.
 * Tries 3 candidate offsets from min(positions) and verifies the extracted
 * triplet matches refCodon — first in forward orientation, then reverse
 * complement (for genes on the minus strand). Returns -1 if not found.
 */
function findCodonStartIndex(
  refCodon: string,
  positions: number[],
  reference: string,
  displayStart: number,
): number {
  if (!refCodon || refCodon.length !== 3 || positions.length === 0) return -1;
  const minPos = Math.min(...positions);
  const maxPos = Math.max(...positions);
  const rcCodon = reverseComplement(refCodon);

  for (const offset of [0, -1, -2]) {
    const candidateStart = minPos + offset;
    const candidateEnd = candidateStart + 2;
    if (candidateEnd < maxPos) continue;
    const idx = candidateStart - displayStart;
    if (idx < 0 || idx + 3 > reference.length) continue;
    const extracted = reference.substring(idx, idx + 3).toUpperCase();
    const matches = extracted === refCodon.toUpperCase() || extracted === rcCodon.toUpperCase();
    if (!matches) continue;
    const allWithin = positions.every(
      (p) => p >= candidateStart && p <= candidateEnd,
    );
    if (allWithin) return idx;
  }
  return -1;
}

/** Map change type to a CSS class for color-coded annotation */
function changeTypeColorClass(changeType: string): string {
  const ct = changeType.toLowerCase();
  if (ct.includes("stopgained") || ct.includes("stop_gained")) return "bam-codon--stop-gained";
  if (ct.includes("stoplost") || ct.includes("stop_lost")) return "bam-codon--stop-lost";
  if (ct.includes("nonsynonymous") || ct.includes("non_synonymous")) return "bam-codon--nonsynonymous";
  if (ct.includes("synonymous")) return "bam-codon--synonymous";
  if (ct.includes("frameshift")) return "bam-codon--frameshift";
  if (ct.includes("indel")) return "bam-codon--indel";
  return "bam-codon--unknown";
}

function layoutVisibleLoci(loci: ViewerLocus[], chrom: string, range: WindowRange): LaidOutLocus[] {
  const overlapping = loci
    .filter((l) => l.chrom === chrom)
    .map((l) => ({ ...l, ...locusBounds(l) }))
    .filter((l) => l.start <= range.end && l.end >= range.start)
    .sort((a, b) => a.start - b.start || a.end - b.end);

  const laneEnds: number[] = [];
  return overlapping.map((l) => {
    let lane = laneEnds.findIndex((end) => end < l.start);
    if (lane < 0) { lane = laneEnds.length; laneEnds.push(l.end); }
    else { laneEnds[lane] = l.end; }
    return { ...l, lane };
  });
}

/* ── BamCell ──────────────────────────────────────────── */

function BamCell({ value, position, referenceBase, site }: {
  value: string;
  position: number;
  referenceBase: string;
  site?: BamVariantSite;
}) {
  const cls = ["bam-cell"];
  let text = value || "";

  const uc = value?.toUpperCase() ?? "";
  if (!value) {
    cls.push("bam-cell--empty");
  } else if (site) {
    cls.push("bam-cell--focus");
    if (uc === site.altBase.toUpperCase()) cls.push("bam-cell--focus-alt");
    else if (uc === site.referenceBase.toUpperCase()) cls.push("bam-cell--focus-ref");
    else if (value === "-") cls.push("bam-cell--focus-gap");
    else cls.push("bam-cell--focus-other");
  } else if (value === "-") {
    cls.push("bam-cell--gap");
  } else if (uc === referenceBase.toUpperCase()) {
    cls.push("bam-cell--match");
    text = ""; // IGV-style: match = colored bar, no text
  } else {
    cls.push("bam-cell--mismatch");
    const nc = nucleotideClass(value);
    if (nc) cls.push(nc);
  }

  return (
    <span className={cls.join(" ")} title={`${position}: ${value || "–"}`}>
      {text}
    </span>
  );
}

/* ── Main component ───────────────────────────────────── */

export default function BamViewer({ bamPath, fastaPath, data, minMapq, minBaseQuality }: BamViewerProps) {
  const [search, setSearch] = useState("");
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [view, setView] = useState<BamViewResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [expanded, setExpanded] = useState(false);
  const [cellSize, setCellSize] = useState(DEFAULT_CELL_SIZE);
  const [visibleSupport, setVisibleSupport] = useState<Set<string>>(
    () => new Set(["mnv", "partial", "reference", "other"]),
  );
  const gridWrapRef = useRef<HTMLDivElement>(null);
  const prevCellSizeRef = useRef(DEFAULT_CELL_SIZE);
  /* ── Loci ── */
  const loci = useMemo(() => buildLoci(data), [data]);

  const filteredLoci = useMemo(() => {
    const q = search.trim().toLowerCase();
    if (!q) return loci;
    return loci.filter((l) =>
      [l.chrom, l.gene, l.variantType, l.aaChanges, l.positions.join(","),
       l.refCodon, l.mnvCodon, l.snpAaChanges, l.changeType].join(" ").toLowerCase().includes(q),
    );
  }, [loci, search]);

  // Auto-select first locus
  useEffect(() => {
    if (filteredLoci.length === 0) { setSelectedId(null); setView(null); return; }
    if (!selectedId || !filteredLoci.some((l) => l.id === selectedId)) {
      setSelectedId(filteredLoci[0].id);
    }
  }, [filteredLoci, selectedId]);

  const selectedLocus = useMemo(
    () => filteredLoci.find((l) => l.id === selectedId) ?? null,
    [filteredLoci, selectedId],
  );

  const selectedIndex = selectedLocus
    ? filteredLoci.findIndex((l) => l.id === selectedLocus.id)
    : -1;

  /* ── Keyboard navigation ── */
  const shiftLocus = useCallback((delta: number) => {
    if (selectedIndex < 0 || filteredLoci.length === 0) return;
    const next = selectedIndex + delta;
    if (next >= 0 && next < filteredLoci.length) setSelectedId(filteredLoci[next].id);
  }, [selectedIndex, filteredLoci]);

  useEffect(() => {
    if (loci.length === 0) return;
    const h = (e: KeyboardEvent) => {
      // Don't capture when focus is in an input/textarea (e.g. search box)
      if (e.target instanceof HTMLInputElement || e.target instanceof HTMLTextAreaElement) return;
      switch (e.key) {
        case "Escape":
          if (expanded) setExpanded(false);
          break;
        case "ArrowUp":
        case "ArrowLeft":
          e.preventDefault();
          shiftLocus(-1);
          break;
        case "ArrowDown":
        case "ArrowRight":
          e.preventDefault();
          shiftLocus(1);
          break;
      }
    };
    document.addEventListener("keydown", h);
    return () => document.removeEventListener("keydown", h);
  }, [loci.length, expanded, shiftLocus]);

  /* ── Scroll sidebar to active locus (keyboard nav) ── */
  useEffect(() => {
    if (!selectedId) return;
    requestAnimationFrame(() => {
      const active = document.querySelector(".bam-locus-item--active");
      active?.scrollIntoView({ block: "nearest", behavior: "smooth" });
    });
  }, [selectedId]);

  /* ── Fetch: ONE call per locus change, wide window ── */
  // Depends on `selectedId` + `loci` (stable references) instead of `selectedLocus`
  // (object ref that changes on every search keystroke). This avoids redundant
  // backend calls when the user types in the search box.
  useEffect(() => {
    const locus = loci.find((l) => l.id === selectedId) ?? null;
    if (!locus) { setView(null); setError(null); setLoading(false); return; }
    let cancelled = false;
    // Clear previous view immediately to avoid showing stale reads from
    // a different locus while the new data loads.
    setView(null);
    setLoading(true);
    setError(null);

    const bounds = locusBounds(locus);
    const windowStart = Math.max(1, bounds.start - BUFFER_PADDING);
    const windowEnd = bounds.end + BUFFER_PADDING;

    invoke<BamViewResponse>("get_bam_view", {
      request: {
        bamPath,
        fastaPath,
        chrom: locus.chrom,
        positions: locus.positions,
        refBases: locus.refBases,
        altBases: locus.altBases,
        minMapq,
        minBaseQuality,
        maxReads: VIEWER_MAX_READS,
        windowStart,
        windowEnd,
      },
    })
      .then((r) => {
        if (cancelled) return;
        setView(r);
        setLoading(false);
      })
      .catch((e) => { if (!cancelled) { setError(String(e)); setLoading(false); } });

    return () => { cancelled = true; };
    // `view` intentionally omitted — only used for conditional loading bar,
    // including it would cause an infinite fetch loop.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [bamPath, fastaPath, minBaseQuality, minMapq, selectedId, loci]);

  /* ── Auto-scroll to center MNV when data loads ── */
  // Uses `selectedId` (string) instead of `selectedLocus` (object) to avoid
  // spurious scroll resets when typing in search. Guards against stale view
  // data by checking that the locus positions fall within the display range.
  useEffect(() => {
    if (!view || !selectedId || !gridWrapRef.current) return;
    const locus = loci.find((l) => l.id === selectedId);
    if (!locus) return;
    const bounds = locusBounds(locus);
    // Only scroll if the view data actually covers this locus
    if (bounds.start < view.displayStart || bounds.end > view.displayEnd) return;
    const center = Math.round((bounds.start + bounds.end) / 2);
    const offsetPx = (center - view.displayStart) * cellSize;
    const vw = gridWrapRef.current.clientWidth;
    requestAnimationFrame(() => {
      if (gridWrapRef.current) {
        gridWrapRef.current.scrollLeft = Math.max(0, offsetPx - vw / 2 + 60);
      }
    });
    // `loci` and `cellSize` intentionally omitted — loci is only for lookup
    // (stable per ID), cellSize zoom is handled by the zoom compensation effect.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [view, selectedId]);

  /* ── Compensate scroll position on zoom ── */
  useEffect(() => {
    if (!gridWrapRef.current || prevCellSizeRef.current === cellSize) return;
    const ratio = cellSize / prevCellSizeRef.current;
    gridWrapRef.current.scrollLeft = Math.round(gridWrapRef.current.scrollLeft * ratio);
    prevCellSizeRef.current = cellSize;
  }, [cellSize]);

  /* ── Drag-to-scroll (grab & pan) ── */
  const dragState = useRef<{ active: boolean; startX: number; startY: number; scrollLeft: number; scrollTop: number }>({
    active: false, startX: 0, startY: 0, scrollLeft: 0, scrollTop: 0,
  });

  useEffect(() => {
    const el = gridWrapRef.current;
    if (!el) return;

    const onMouseDown = (e: MouseEvent) => {
      // Only start drag on primary button, skip if clicking interactive elements
      if (e.button !== 0) return;
      const target = e.target as HTMLElement;
      if (target.closest("button, a, input, select, [role='button']")) return;
      dragState.current = {
        active: true,
        startX: e.clientX,
        startY: e.clientY,
        scrollLeft: el.scrollLeft,
        scrollTop: el.scrollTop,
      };
      el.classList.add("bam-grid-wrap--dragging");
    };

    const onMouseMove = (e: MouseEvent) => {
      if (!dragState.current.active) return;
      e.preventDefault();
      const dx = e.clientX - dragState.current.startX;
      const dy = e.clientY - dragState.current.startY;
      el.scrollLeft = dragState.current.scrollLeft - dx;
      el.scrollTop = dragState.current.scrollTop - dy;
    };

    const onMouseUp = () => {
      if (!dragState.current.active) return;
      dragState.current.active = false;
      el.classList.remove("bam-grid-wrap--dragging");
    };

    el.addEventListener("mousedown", onMouseDown);
    window.addEventListener("mousemove", onMouseMove);
    window.addEventListener("mouseup", onMouseUp);
    return () => {
      el.removeEventListener("mousedown", onMouseDown);
      window.removeEventListener("mousemove", onMouseMove);
      window.removeEventListener("mouseup", onMouseUp);
    };
    // Re-attach when view changes (grid mounts/unmounts conditionally)
  }, [view]);

  /* ── Derived data ── */
  const siteMap = useMemo(
    () => view ? focusSiteMap(view.sites) : new Map<number, BamVariantSite>(),
    [view],
  );
  const referenceBases = useMemo(
    () => view ? view.reference.split("") : [],
    [view],
  );
  const windowPositions = useMemo(
    () => view ? Array.from({ length: view.reference.length }, (_, i) => view.displayStart + i) : [],
    [view],
  );
  const trackWidth = windowPositions.length * cellSize;
  const tickStep = tickStepForCellSize(cellSize);
  const trackStyle = { "--bam-cell-size": `${cellSize}px` } as CSSProperties;

  // Coverage — use real per-position depth from backend (all reads, not just displayed)
  const coverageData = useMemo(() => {
    if (!view) return [];
    return view.coverage ?? [];
  }, [view]);
  const maxCoverage = useMemo(() => Math.max(1, ...coverageData), [coverageData]);

  // Use real counts from ALL reads (backend computes before truncation)
  const displayCounts = useMemo(() => {
    if (!view) return { mnv: 0, partial: 0, reference: 0, other: 0 };
    return {
      mnv: view.counts.mnv,
      partial: view.counts.partial,
      reference: view.counts.reference,
      other: view.counts.other,
    };
  }, [view]);

  // Filtered reads by support type
  const filteredReads = useMemo(
    () => view ? view.reads.filter((r) => visibleSupport.has(r.support)) : [],
    [view, visibleSupport],
  );

  // Codon annotation for the selected locus
  const codonAnnotation = useMemo(() => {
    if (!selectedLocus || !view || !selectedLocus.refCodon) return null;
    const startIdx = findCodonStartIndex(
      selectedLocus.refCodon,
      selectedLocus.positions,
      view.reference,
      view.displayStart,
    );
    if (startIdx < 0) return null;
    return {
      startIdx,
      startPos: view.displayStart + startIdx,
      refCodon: selectedLocus.refCodon,
      mnvCodon: selectedLocus.mnvCodon,
      aaChange: selectedLocus.aaChanges,
      changeType: selectedLocus.changeType,
      snpCodons: selectedLocus.snpCodons,
      snpAaChanges: selectedLocus.snpAaChanges,
      variantType: selectedLocus.variantType,
    };
  }, [selectedLocus, view]);

  const toggleSupport = useCallback((type: string) => {
    setVisibleSupport((prev) => {
      const next = new Set(prev);
      if (next.has(type)) next.delete(type);
      else next.add(type);
      return next;
    });
  }, []);

  // Visible loci for variants track
  const visibleLoci = useMemo(
    () => selectedLocus && view
      ? layoutVisibleLoci(loci, selectedLocus.chrom, { start: view.displayStart, end: view.displayEnd })
      : [],
    [loci, selectedLocus, view],
  );
  const visibleLocusLanes = visibleLoci.reduce((m, l) => Math.max(m, l.lane + 1), 0);

  /* ── Render ── */
  if (loci.length === 0) return null;

  return (
    <section className={`step-section bam-viewer${expanded ? " bam-viewer--expanded" : ""}`}>
      <div className="step-header">
        <div className="step-badge step-badge-3">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
            <path d="M4 6h16M4 12h16M4 18h16" />
            <path d="M9 4v16M15 4v16" opacity="0.45" />
          </svg>
        </div>
        <h3 className="step-title">Genomic Track Viewer</h3>
        <span className="step-subtitle">{countLabel(loci.length)} MNV loci</span>
      </div>

        <div className="bam-viewer-shell">
          {/* ── Locus sidebar ── */}
          <aside className="bam-locus-panel">
            <div className="bam-panel-head">
              <h4>Loci</h4>
              <span>{countLabel(filteredLoci.length)}</span>
            </div>
            <div className="vt-search-wrapper bam-locus-search">
              <svg className="vt-search-icon" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
                <circle cx="11" cy="11" r="8" /><path d="M21 21l-4.35-4.35" />
              </svg>
              <input type="text" className="vt-search" placeholder="Search gene, contig..." value={search} onChange={(e) => setSearch(e.target.value)} />
            </div>
            <div className="bam-locus-list">
              {filteredLoci.map((l) => (
                <button
                  key={l.id}
                  type="button"
                  className={`bam-locus-item${selectedId === l.id ? " bam-locus-item--active" : ""}`}
                  onClick={() => setSelectedId(l.id)}
                >
                  <div className="bam-locus-item-top">
                    <span className="bam-locus-gene">{l.gene}</span>
                    <span className="bam-locus-type">{l.variantType}</span>
                  </div>
                  <div className="bam-locus-item-mid">{l.chrom}:{l.positions.join(", ")}</div>
                  <div className="bam-locus-item-bottom">
                    <span>{l.refBases.join(",")} → {l.altBases.join(",")}</span>
                    {l.mnvReads && <span>×{l.mnvReads}</span>}
                  </div>
                </button>
              ))}
            </div>
          </aside>

          {/* ── Stage ── */}
          <div className="bam-stage">
            <div className="bam-stage-toolbar">
              <div className="bam-stage-title">
                <h4>{selectedLocus ? locusTitle(selectedLocus) : "Select a locus"}</h4>
                {selectedLocus && (
                  <span className="bam-stage-subtitle">
                    {selectedLocus.aaChanges || selectedLocus.changeType || ""}
                  </span>
                )}
              </div>
              <div className="bam-stage-actions">
                <button type="button" className="bam-nav-btn" onClick={() => shiftLocus(-1)} disabled={selectedIndex <= 0}>
                  <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round"><path d="M15 18l-6-6 6-6" /></svg>
                  Prev
                </button>
                <button type="button" className="bam-nav-btn" onClick={() => shiftLocus(1)} disabled={selectedIndex < 0 || selectedIndex >= filteredLoci.length - 1}>
                  Next
                  <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round"><path d="M9 18l6-6-6-6" /></svg>
                </button>
                <label className="bam-zoom-slider-label" title="Zoom level">
                  <svg width="11" height="11" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
                    <circle cx="11" cy="11" r="8" /><path d="M21 21l-4.35-4.35" /><path d="M8 11h6" />
                  </svg>
                  <input
                    type="range"
                    className="bam-zoom-slider"
                    min={MIN_CELL_SIZE}
                    max={MAX_CELL_SIZE}
                    value={cellSize}
                    onChange={(e) => setCellSize(Number(e.target.value))}
                  />
                  <svg width="11" height="11" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
                    <circle cx="11" cy="11" r="8" /><path d="M21 21l-4.35-4.35" /><path d="M8 11h6M11 8v6" />
                  </svg>
                </label>
                <button type="button" className="bam-nav-btn bam-expand-btn" onClick={() => setExpanded(!expanded)} title={expanded ? "Exit fullscreen" : "Expand to fullscreen"}>
                  {expanded ? (
                    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round">
                      <path d="M8 3v3a2 2 0 01-2 2H3m18 0h-3a2 2 0 01-2-2V3m0 18v-3a2 2 0 012-2h3M3 16h3a2 2 0 012 2v3" />
                    </svg>
                  ) : (
                    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round">
                      <path d="M15 3h6v6M9 21H3v-6M21 3l-7 7M3 21l7-7" />
                    </svg>
                  )}
                </button>
              </div>
            </div>

            {/* ── Focus strip + counts ── */}
            {selectedLocus && (
              <div className="bam-focus-strip">
                {selectedLocus.positions.map((pos, i) => (
                  <span key={pos} className="bam-focus-chip">
                    {pos} {selectedLocus.refBases[i]}&gt;{selectedLocus.altBases[i]}
                  </span>
                ))}
                {selectedLocus.refCodon && (
                  <>
                    <span className={`bam-focus-chip bam-focus-chip--codon ${changeTypeColorClass(selectedLocus.changeType)}`}>
                      <span className="bam-track-tag bam-track-tag--mnv">MNV</span>
                      {selectedLocus.refCodon} → {selectedLocus.mnvCodon || "?"}{" "}
                      ({selectedLocus.aaChanges || selectedLocus.changeType})
                    </span>
                    {selectedLocus.snpAaChanges && (() => {
                      const snpCodonList = selectedLocus.snpCodons ? selectedLocus.snpCodons.split(",").map((s) => s.trim()) : [];
                      const snpAaList = selectedLocus.snpAaChanges.split(",").map((s) => s.trim());
                      return snpAaList.map((aa, idx) => {
                        const parsed = parseAaChange(aa);
                        return (
                          <span key={idx} className="bam-focus-chip bam-focus-chip--snp-aa">
                            <span className="bam-track-tag bam-track-tag--snp">SNP</span>
                            {snpCodonList[idx] && (
                              <>
                                <span className="bam-snp-badge bam-snp-badge--ref bam-snp-badge--inline">REF</span>
                                <span className="bam-focus-mono">{selectedLocus.refCodon}</span>
                                {" → "}
                                <span className="bam-snp-badge bam-snp-badge--alt bam-snp-badge--inline">ALT</span>
                                <span className="bam-focus-mono">{snpCodonList[idx]}</span>
                                {" "}
                              </>
                            )}
                            {parsed && (
                              <span className="bam-focus-aa-change">
                                ({parsed.refAa}→{parsed.altAa})
                              </span>
                            )}
                          </span>
                        );
                      });
                    })()}
                  </>
                )}
                {view && (
                  <span className="bam-focus-summary">
                    {filteredReads.length !== view.reads.length
                      ? `${countLabel(filteredReads.length)} shown · `
                      : ""}
                    {countLabel(view.reads.length)} of {countLabel(view.totalReads)} reads
                    {view.truncated ? " (truncated)" : ""}
                  </span>
                )}
              </div>
            )}

            {/* ── Loading — only on initial load, not when switching loci ── */}
            {loading && !view && (
              <div className="vt-loading bam-loading">
                <div className="progress-bar" style={{ maxWidth: 240, margin: "0 auto" }}>
                  <div className="progress-bar-indeterminate" />
                </div>
                <p className="hint">Loading reads…</p>
              </div>
            )}

            {error && !view && (
              <div className="error-banner"><strong>Error:</strong> {error}</div>
            )}

            {/* ── Stale-data error hint (view exists but last fetch failed) ── */}
            {error && view && (
              <div className="bam-stale-error" title={error}>
                ⚠ Failed to load new data — showing previous view
              </div>
            )}

            {/* ── Tracks — native horizontal scroll ── */}
            {view && (
              <>
                <div className="bam-count-grid">
                  <button type="button" className={`bam-count-card${visibleSupport.has("mnv") ? "" : " bam-count-card--off"}`} onClick={() => toggleSupport("mnv")} title="Toggle MNV reads">
                    <strong>{countLabel(displayCounts.mnv)}</strong><span>MNV</span>
                  </button>
                  <button type="button" className={`bam-count-card${visibleSupport.has("partial") ? "" : " bam-count-card--off"}`} onClick={() => toggleSupport("partial")} title="Toggle Partial reads">
                    <strong>{countLabel(displayCounts.partial)}</strong><span>Partial</span>
                  </button>
                  <button type="button" className={`bam-count-card${visibleSupport.has("reference") ? "" : " bam-count-card--off"}`} onClick={() => toggleSupport("reference")} title="Toggle Reference reads">
                    <strong>{countLabel(displayCounts.reference)}</strong><span>Ref</span>
                  </button>
                  <button type="button" className={`bam-count-card${visibleSupport.has("other") ? "" : " bam-count-card--off"}`} onClick={() => toggleSupport("other")} title="Toggle Other reads">
                    <strong>{countLabel(displayCounts.other)}</strong><span>Other</span>
                  </button>
                </div>

                <div className="bam-grid-wrap" ref={gridWrapRef}>
                  <div className="bam-grid">
                    {/* Ruler */}
                    <div className="bam-track-shell">
                      <div className="bam-grid-label bam-grid-label--track">
                        <span className="bam-read-name">Ruler</span>
                      </div>
                      <div className="bam-scale-track">
                        <div className="bam-scale-track-inner" style={{ width: `${trackWidth}px` }}>
                          {windowPositions.map((pos, i) => {
                            const isFocus = siteMap.has(pos);
                            const isEdge = i === 0 || i === windowPositions.length - 1;
                            const isTick = i % tickStep === 0;
                            if (!isFocus && !isEdge && !isTick) return null;
                            // Skip regular ticks too close to a focus position to avoid label overlap.
                            // Label width ≈ digits × 5px; need at least that many cells clearance.
                            if (!isFocus && (isEdge || isTick)) {
                              const minGap = Math.ceil(35 / cellSize);
                              for (const [focusPos] of siteMap) {
                                const focusIdx = focusPos - view.displayStart;
                                if (Math.abs(i - focusIdx) > 0 && Math.abs(i - focusIdx) < minGap) return null;
                              }
                            }
                            return (
                              <span key={pos} className={`bam-scale-tick${isFocus ? " bam-scale-tick--focus" : ""}`} style={{ left: `${i * cellSize}px` }}>
                                <span className="bam-scale-mark" />
                                <span className="bam-scale-label">{pos}</span>
                              </span>
                            );
                          })}
                        </div>
                      </div>
                    </div>

                    {/* Variants */}
                    {visibleLoci.length > 0 && (
                      <div className="bam-track-shell">
                        <div className="bam-grid-label bam-grid-label--track">
                          <span className="bam-read-name">Variants</span>
                        </div>
                        <div className="bam-feature-stage">
                          <div className="bam-feature-canvas" style={{ width: `${trackWidth}px`, height: `${Math.max(1, visibleLocusLanes) * 28}px` }}>
                            {visibleLoci.map((l) => {
                              const left = (l.start - view.displayStart) * cellSize;
                              const w = Math.max(cellSize, (l.end - l.start + 1) * cellSize);
                              return (
                                <button
                                  key={l.id}
                                  type="button"
                                  className={`bam-locus-feature${l.id === selectedId ? " bam-locus-feature--active" : ""}`}
                                  style={{ left: `${left}px`, width: `${w}px`, top: `${l.lane * 28}px` }}
                                  title={`${l.gene} · ${l.chrom}:${l.positions.join(",")} ${l.aaChanges || l.changeType || ""}`}
                                  onClick={() => setSelectedId(l.id)}
                                >
                                  <span className="bam-locus-feature-gene">{l.gene}</span>
                                  <span className="bam-locus-feature-coords">
                                    {l.positions.length <= 3
                                      ? l.positions.join(", ")
                                      : `${l.positions[0]}..${l.positions[l.positions.length - 1]}`}
                                  </span>
                                </button>
                              );
                            })}
                          </div>
                        </div>
                      </div>
                    )}

                    {/* Coverage */}
                    <div className="bam-track-shell bam-track-shell--coverage">
                      <div className="bam-grid-label bam-grid-label--track">
                        <span className="bam-read-name">Coverage</span>
                        <span className="bam-read-meta">max {maxCoverage}×</span>
                      </div>
                      <div className="bam-coverage-track" style={{ width: `${trackWidth + 12}px` }}>
                        {coverageData.map((depth, i) => {
                          const pos = view.displayStart + i;
                          const site = siteMap.get(pos);
                          return (
                            <span
                              key={pos}
                              className={`bam-coverage-bar${site ? " bam-coverage-bar--variant" : ""}`}
                              style={{ width: `${cellSize}px`, height: `${(depth / maxCoverage) * 100}%` }}
                              title={`${pos}: ${depth}×`}
                            />
                          );
                        })}
                      </div>
                    </div>

                    {/* Reference */}
                    <div className="bam-grid-row bam-grid-row--reference">
                      <div className="bam-grid-label">
                        <span className="bam-read-name">Reference</span>
                      </div>
                      <div className="bam-track" style={trackStyle}>
                        {referenceBases.map((base, i) => {
                          const pos = view.displayStart + i;
                          const site = siteMap.get(pos);
                          return (
                            <span key={pos} className={`bam-cell bam-cell--ref-base ${nucleotideClass(base)}${site ? " bam-cell--focus" : ""}`} title={`${pos}: ${base}`}>
                              {base}
                            </span>
                          );
                        })}
                      </div>
                    </div>

                    {/* Codon / AA annotation — MNV interpretation */}
                    {codonAnnotation && (
                      <div className={`bam-grid-row bam-grid-row--codon ${changeTypeColorClass(codonAnnotation.changeType)}`}>
                        <div className="bam-grid-label bam-grid-label--codon">
                          <span className="bam-read-name">
                            <span className="bam-track-tag bam-track-tag--mnv">MNV</span>
                            Codon
                          </span>
                          <span className="bam-codon-change-label">
                            {codonAnnotation.refCodon} → {codonAnnotation.mnvCodon || "?"}
                          </span>
                          {codonAnnotation.aaChange && (
                            <span className="bam-aa-change-label">{codonAnnotation.aaChange}</span>
                          )}
                        </div>
                        <div className="bam-track" style={trackStyle}>
                          {windowPositions.map((pos, i) => {
                            const inCodon = i >= codonAnnotation.startIdx && i < codonAnnotation.startIdx + 3;
                            if (!inCodon) {
                              return <span key={pos} className="bam-cell bam-cell--codon-empty" />;
                            }
                            const codonOffset = i - codonAnnotation.startIdx;
                            const refBase = codonAnnotation.refCodon[codonOffset] ?? "";
                            const mnvBase = codonAnnotation.mnvCodon?.[codonOffset] ?? "";
                            const changed = mnvBase !== "" && refBase.toUpperCase() !== mnvBase.toUpperCase();
                            const isVariant = selectedLocus!.positions.includes(pos);
                            return (
                              <span
                                key={pos}
                                className={[
                                  "bam-cell",
                                  "bam-cell--codon",
                                  changed ? "bam-cell--codon-changed" : "bam-cell--codon-ref",
                                  isVariant ? "bam-cell--codon-variant" : "",
                                  codonOffset === 0 ? "bam-cell--codon-start" : "",
                                  codonOffset === 2 ? "bam-cell--codon-end" : "",
                                ].filter(Boolean).join(" ")}
                                title={`Codon pos ${codonOffset + 1}: ${refBase}${changed ? ` → ${mnvBase}` : ""}`}
                              >
                                {changed ? mnvBase : refBase}
                              </span>
                            );
                          })}
                        </div>
                      </div>
                    )}

                    {/* SNP individual effects (only for SNP/MNV variants) */}
                    {codonAnnotation && codonAnnotation.variantType.toLowerCase().includes("snp") && codonAnnotation.snpCodons && (
                      <div className="bam-grid-row bam-grid-row--snp-codons">
                        <div className="bam-grid-label bam-grid-label--codon">
                          <span className="bam-read-name">
                            <span className="bam-track-tag bam-track-tag--snp">SNP</span>
                            Codons
                          </span>
                          {codonAnnotation.snpCodons && codonAnnotation.snpAaChanges && (() => {
                            const snpCodonList = codonAnnotation.snpCodons.split(",").map((s) => s.trim());
                            const snpAaList = codonAnnotation.snpAaChanges.split(",").map((s) => s.trim());
                            const firstParsed = parseAaChange(snpAaList[0] ?? "");
                            return (
                              <div className="bam-snp-detail-list">
                                <span className="bam-snp-detail-row bam-snp-detail-row--ref">
                                  <span className="bam-snp-badge bam-snp-badge--ref">REF</span>
                                  <span className="bam-snp-codon-mono">{codonAnnotation.refCodon}</span>
                                  {firstParsed && <span className="bam-snp-aa-part">({firstParsed.refAa})</span>}
                                </span>
                                {snpCodonList.map((codon, idx) => {
                                  const parsed = parseAaChange(snpAaList[idx] ?? "");
                                  return (
                                    <span key={idx} className="bam-snp-detail-row bam-snp-detail-row--alt">
                                      <span className="bam-snp-badge bam-snp-badge--alt">ALT</span>
                                      <span className="bam-snp-codon-mono">{codon}</span>
                                      {parsed ? (
                                        <span className="bam-snp-aa-part">({parsed.altAa})</span>
                                      ) : snpAaList[idx] ? (
                                        <span className="bam-snp-aa-part">({snpAaList[idx]})</span>
                                      ) : null}
                                    </span>
                                  );
                                })}
                              </div>
                            );
                          })()}
                        </div>
                        <div className="bam-track" style={trackStyle}>
                          {windowPositions.map((pos, i) => {
                            const inCodon = i >= codonAnnotation.startIdx && i < codonAnnotation.startIdx + 3;
                            if (!inCodon) {
                              return <span key={pos} className="bam-cell bam-cell--codon-empty" />;
                            }
                            const codonOffset = i - codonAnnotation.startIdx;
                            const isVariantPos = selectedLocus!.positions.includes(pos);
                            if (isVariantPos) {
                              const snpIdx = selectedLocus!.positions.indexOf(pos);
                              const snpCodonList = codonAnnotation.snpCodons.split(",").map((s) => s.trim());
                              const snpCodon = snpCodonList[snpIdx] ?? "";
                              const snpBase = snpCodon[codonOffset] ?? "";
                              const snpAaList = codonAnnotation.snpAaChanges.split(",").map((s) => s.trim());
                              const parsed = parseAaChange(snpAaList[snpIdx] ?? "");
                              const refBase = codonAnnotation.refCodon[codonOffset] ?? "";
                              return (
                                <span
                                  key={pos}
                                  className="bam-cell bam-cell--codon bam-cell--snp-highlight"
                                  title={`SNP${snpIdx + 1}: REF ${codonAnnotation.refCodon}${parsed ? ` (${parsed.refAa})` : ""} → ALT ${snpCodon}${parsed ? ` (${parsed.altAa})` : snpAaList[snpIdx] ? ` (${snpAaList[snpIdx]})` : ""} | Base: ${refBase}→${snpBase}`}
                                >
                                  {snpBase || selectedLocus!.altBases[snpIdx] || "?"}
                                </span>
                              );
                            }
                            const refBase = codonAnnotation.refCodon[codonOffset] ?? "";
                            return (
                              <span key={pos} className="bam-cell bam-cell--codon bam-cell--codon-ref">
                                {refBase}
                              </span>
                            );
                          })}
                        </div>
                      </div>
                    )}

                    {/* Read pileup — compact IGV-style */}
                    {filteredReads.map((read) => (
                      <div
                        key={`${read.name}-${read.start}-${read.end}`}
                        className={`bam-grid-row bam-grid-row--read bam-grid-row--strand-${read.strand === "+" ? "fwd" : "rev"}`}
                        title={`${read.name} | ${read.strand} | MAPQ ${read.mapq} | ${supportLabel(read.support)}`}
                      >
                        <div className="bam-grid-label">
                          <span className={`bam-support-pill bam-support-pill--${read.support}`}>
                            {supportLabel(read.support)}
                          </span>
                          <span className="bam-strand-indicator">{read.strand}</span>
                        </div>
                        <div className="bam-track" style={trackStyle}>
                          {read.bases.map((value, i) => (
                            <BamCell
                              key={`${read.name}-${view.displayStart + i}`}
                              value={value}
                              position={view.displayStart + i}
                              referenceBase={referenceBases[i] ?? ""}
                              site={siteMap.get(view.displayStart + i)}
                            />
                          ))}
                        </div>
                      </div>
                    ))}
                  </div>
                </div>
              </>
            )}
          </div>
        </div>
    </section>
  );
}
