import { useCallback, useEffect, useMemo, useRef, useState, type CSSProperties } from "react";
import { invoke } from "@tauri-apps/api/core";
import type { BamReadView, BamVariantSite, BamViewColumn, BamViewResponse, TsvData } from "../types";

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
  eventClass: string;
  eventComponents: string;
  eventReads: string;
  eventFrequency: string;
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
const SUPPORT_TYPES = ["mnv", "partial", "reference", "other"] as const;
const INDEL_EVENT_CLASSES = new Set(["insertion", "deletion", "delins", "complex_indel", "symbolic"]);

/* ── Pure helpers ─────────────────────────────────────── */

function headerIndex(headers: string[], label: string): number {
  return headers.findIndex((h) => h.trim().toLowerCase() === label.trim().toLowerCase());
}

function splitField(value: string): string[] {
  return value.split(",").map((s) => s.trim()).filter((s) => s.length > 0);
}

function nonPlaceholder(value: string | undefined): string {
  const trimmed = (value ?? "").trim();
  return trimmed === "-" ? "" : trimmed;
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

function parseFirstCount(value: string): number {
  const raw = splitField(value)[0] ?? value.trim();
  const parsed = Number.parseInt(raw, 10);
  return Number.isFinite(parsed) ? parsed : 0;
}

function isIndelEventClass(eventClass: string): boolean {
  return INDEL_EVENT_CLASSES.has(eventClass.trim().toLowerCase());
}

function isRenderableLocus(variantType: string, eventClass: string): boolean {
  const vt = variantType.toLowerCase();
  return vt.includes("mnv") || vt.includes("indel") || isIndelEventClass(eventClass);
}

function locusEvidenceReads(locus: Pick<ViewerLocus, "eventReads" | "mnvReads">): string {
  return locus.eventReads && locus.eventReads !== "-" ? locus.eventReads : locus.mnvReads;
}

function locusEvidenceFrequency(locus: Pick<ViewerLocus, "eventFrequency" | "mnvFrequency">): string {
  return locus.eventFrequency && locus.eventFrequency !== "-" ? locus.eventFrequency : locus.mnvFrequency;
}

function locusTypeLabel(locus: Pick<ViewerLocus, "eventClass" | "variantType">): string {
  if (isIndelEventClass(locus.eventClass)) return locus.eventClass.replace("_", " ").toUpperCase();
  return locus.variantType || "Variant";
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
  const eci = headerIndex(data.headers, "Event Class");
  const ecpi = headerIndex(data.headers, "Event Components");
  const eri = headerIndex(data.headers, "Event Reads");
  const efi = headerIndex(data.headers, "Event Frequency");
  const rci = headerIndex(data.headers, "Reference Codon");
  const mci = headerIndex(data.headers, "MNV Codon");
  const sci = headerIndex(data.headers, "SNP Codon");
  const sai = headerIndex(data.headers, "SNP AA Changes");

  if ([ci, gi, pi, ri, ai, vi].some((i) => i < 0)) return [];

  return data.rows
    .map((row) => {
      const vt = row[vi] ?? "";
      const eventClass = eci >= 0 ? nonPlaceholder(row[eci]) : "";
      if (!isRenderableLocus(vt, eventClass)) return null;

      // Skip loci with zero BAM support when read-derived columns are present.
      const eventReadCount = eri >= 0 ? parseFirstCount(row[eri] ?? "0") : 0;
      const mnvReadCount = mri >= 0 ? parseFirstCount(row[mri] ?? "0") : 0;
      if ((eri >= 0 || mri >= 0) && Math.max(eventReadCount, mnvReadCount) <= 0) return null;

      const positions = parsePositions(row[pi] ?? "");
      const refBases = splitField(row[ri] ?? "");
      const altBases = splitField(row[ai] ?? "");
      if (positions.length === 0 || positions.length !== refBases.length || positions.length !== altBases.length) return null;

      const chrom = row[ci] ?? "";
      return {
        id: `${chrom}:${positions.join("-")}:${refBases.join("")}>${altBases.join("")}:${eventClass || vt}`,
        chrom,
        gene: (row[gi] && row[gi].trim()) || `${row[ci] ?? "?"}:${parsePositions(row[pi] ?? "").join("-")}`,
        positions,
        refBases,
        altBases,
        variantType: vt,
        changeType: cti >= 0 ? nonPlaceholder(row[cti]) : "",
        aaChanges: aai >= 0 ? nonPlaceholder(row[aai]) : "",
        mnvReads: mri >= 0 ? nonPlaceholder(row[mri]) : "",
        mnvFrequency: mfi >= 0 ? nonPlaceholder(row[mfi]) : "",
        eventClass,
        eventComponents: ecpi >= 0 ? nonPlaceholder(row[ecpi]) : "",
        eventReads: eri >= 0 ? nonPlaceholder(row[eri]) : "",
        eventFrequency: efi >= 0 ? nonPlaceholder(row[efi]) : "",
        refCodon: rci >= 0 ? nonPlaceholder(row[rci]) : "",
        mnvCodon: mci >= 0 ? nonPlaceholder(row[mci]) : "",
        snpCodons: sci >= 0 ? nonPlaceholder(row[sci]) : "",
        snpAaChanges: sai >= 0 ? nonPlaceholder(row[sai]) : "",
      };
    })
    .filter((l): l is ViewerLocus => l !== null);
}

function locusBounds(locus: ViewerLocus): WindowRange {
  const starts = locus.positions;
  const ends = locus.positions.map((p, i) => p + Math.max(1, locus.refBases[i]?.length ?? 1) - 1);
  return { start: Math.min(...starts), end: Math.max(...ends) };
}

function countLabel(n: number): string {
  return n.toLocaleString("en-US");
}

function supportLabel(s: BamReadView["support"]): string {
  switch (s) {
    case "mnv": return "ALT";
    case "partial": return "SNP/partial";
    case "reference": return "Ref";
    case "other": return "Other";
    default: return s;
  }
}

function supportButtonTitle(visibleSupport: Set<string>, type: string, label: string): string {
  return visibleSupport.size === 1 && visibleSupport.has(type)
    ? "Show all reads"
    : `Show only ${label} reads`;
}

function locusTitle(l: ViewerLocus): string {
  return `${l.gene} · ${l.chrom}:${l.positions.join(", ")} · ${locusTypeLabel(l)}`;
}

function expectedInsertionAfterSite(site: BamVariantSite): { anchor: number; sequence: string } | null {
  const refChars = Array.from(site.referenceBase);
  const altChars = Array.from(site.altBase);
  if (altChars.length <= refChars.length) return null;

  const minLen = Math.min(refChars.length, altChars.length);
  let prefix = 0;
  while (prefix < minLen && refChars[prefix].toUpperCase() === altChars[prefix].toUpperCase()) {
    prefix += 1;
  }

  const insertedLen = altChars.length - refChars.length;
  const anchor = prefix === 0 ? site.position : site.position + prefix - 1;
  return { anchor, sequence: altChars.slice(prefix, prefix + insertedLen).join("") };
}

function siteForColumn(column: BamViewColumn, sites: BamVariantSite[]): BamVariantSite | undefined {
  if (column.kind === "ins") {
    return sites.find((site) => {
      const insertion = expectedInsertionAfterSite(site);
      return Boolean(
        insertion
        && insertion.anchor === column.position
        && (column.insertionIndex ?? 0) >= 1
        && (column.insertionIndex ?? 0) <= insertion.sequence.length,
      );
    });
  }

  return sites.find((site) => {
    const spanEnd = site.position + Math.max(1, site.referenceBase.length) - 1;
    return column.position >= site.position && column.position <= spanEnd;
  });
}

function expectedBaseForColumn(column: BamViewColumn, site?: BamVariantSite): string | null {
  if (!site) return null;
  if (column.kind === "ins") {
    const insertion = expectedInsertionAfterSite(site);
    if (!insertion || insertion.anchor !== column.position) return null;
    return insertion.sequence[(column.insertionIndex ?? 1) - 1] ?? null;
  }
  const offset = column.position - site.position;
  if (site.altBase.length === site.referenceBase.length) return site.altBase[offset] ?? null;
  return null;
}

function refColumnIndex(columns: BamViewColumn[], position: number): number {
  return columns.findIndex((column) => column.kind === "ref" && column.position === position);
}

function columnRangeForLocus(columns: BamViewColumn[], locus: ViewerLocus): { left: number; width: number } {
  const bounds = locusBounds(locus);
  let left = columns.findIndex((column) => column.kind === "ref" && column.position >= bounds.start);
  if (left < 0) left = 0;

  let right = left;
  columns.forEach((column, idx) => {
    if (column.position >= bounds.start && column.position <= bounds.end) {
      right = Math.max(right, idx);
    }
  });

  return {
    left,
    width: Math.max(1, right - left + 1),
  };
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
 *
 * Strategy: the SNP positions define which codon positions are mutated.
 * A codon is 3 consecutive bases. We try all 3 possible reading frames
 * (offset 0, -1, -2 from min position) that would contain ALL SNP positions.
 * We then verify by matching against refCodon or its reverse complement
 * (for minus-strand genes). If no match, we fall back to the first frame
 * that geometrically contains all positions — the codon info from the TSV
 * is authoritative, so we show it even if the genomic sequence looks different
 * (e.g., complex indel context, IUPAC ambiguity).
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

  let fallbackIdx = -1;

  for (const offset of [0, -1, -2]) {
    const candidateStart = minPos + offset;
    const candidateEnd = candidateStart + 2;
    if (candidateEnd < maxPos) continue; // doesn't span all positions
    const idx = candidateStart - displayStart;
    if (idx < 0 || idx + 3 > reference.length) continue;

    const allWithin = positions.every(
      (p) => p >= candidateStart && p <= candidateEnd,
    );
    if (!allWithin) continue;

    // Remember first valid frame as fallback
    if (fallbackIdx < 0) fallbackIdx = idx;

    // Try to verify against reference
    const extracted = reference.substring(idx, idx + 3).toUpperCase();
    if (extracted === refCodon.toUpperCase() || extracted === rcCodon.toUpperCase()) {
      return idx; // verified match
    }
  }

  // No exact match found — use geometric fallback (TSV codon is authoritative)
  return fallbackIdx;
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

function BamCell({ value, column, site, expectedBase, isReadStart, isReadEnd, strand }: {
  value: string;
  column: BamViewColumn;
  site?: BamVariantSite;
  expectedBase?: string | null;
  isReadStart?: boolean;
  isReadEnd?: boolean;
  strand?: string;
}) {
  const cls = ["bam-cell"];
  let text = value || "";
  const isInsertionColumn = column.kind === "ins";

  const uc = value?.toUpperCase() ?? "";
  if (!value || value === " ") {
    cls.push("bam-cell--empty");
    if (isInsertionColumn) cls.push("bam-cell--insertion");
  } else if (isInsertionColumn) {
    cls.push("bam-cell--insertion-base");
    if (site) {
      cls.push("bam-cell--focus");
      if (expectedBase && uc === expectedBase.toUpperCase()) cls.push("bam-cell--focus-alt");
      else cls.push("bam-cell--focus-other");
    }
    const nc = nucleotideClass(value);
    if (nc) cls.push(nc);
  } else if (site) {
    cls.push("bam-cell--focus");
    if (expectedBase && uc === expectedBase.toUpperCase()) cls.push("bam-cell--focus-alt");
    else if (uc === site.altBase.toUpperCase()) cls.push("bam-cell--focus-alt");
    else if (uc === site.referenceBase.toUpperCase()) cls.push("bam-cell--focus-ref");
    else if (value === "-") cls.push("bam-cell--focus-gap");
    else cls.push("bam-cell--focus-other");
  } else if (value === "-") {
    cls.push("bam-cell--gap");
  } else if (uc === column.referenceBase.toUpperCase()) {
    cls.push("bam-cell--match");
    text = ""; // IGV-style: match = colored bar, no text
  } else {
    cls.push("bam-cell--mismatch");
    const nc = nucleotideClass(value);
    if (nc) cls.push(nc);
  }

  // Read boundary markers
  if (isReadStart) cls.push("bam-cell--read-start");
  if (isReadEnd) cls.push("bam-cell--read-end");

  // Arrow direction: forward reads point right (▶ at end), reverse point left (◀ at start)
  const isFwd = strand === "+";
  const showArrow = (isFwd && isReadEnd) || (!isFwd && isReadStart);
  const titlePosition = isInsertionColumn
    ? `${column.position}+${column.insertionIndex ?? ""}`
    : `${column.position}`;

  return (
    <span className={cls.join(" ")} title={`${titlePosition}: ${value || "–"}`}>
      {text}
      {showArrow && <span className={`bam-cell-arrow bam-cell-arrow--${isFwd ? "fwd" : "rev"}`}>{isFwd ? "▸" : "◂"}</span>}
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
    () => new Set(SUPPORT_TYPES),
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

  // Auto-select first locus.
  //
  // The rule is right about this one and it is not fixed here. Correcting the
  // selection in an effect costs a render: the list paints with a selection
  // that is not in it, then paints again. The answer is to derive the effective
  // selection while rendering and keep `selectedId` for what the user clicked,
  // which means changing the ten places that read `selectedId`, among them the
  // dependency list of the fetch below, whose exact contents are load-bearing
  // and are explained in its own comment.
  //
  // That is a change to how this viewer chooses what to show, and nothing in
  // the repository tests that today. It wants its own pull request, with tests
  // for the selection first, not a passing remark in a dependency bump.
  /* eslint-disable react-hooks/set-state-in-effect -- deferred, see above */
  useEffect(() => {
    if (filteredLoci.length === 0) { setSelectedId(null); setView(null); return; }
    if (!selectedId || !filteredLoci.some((l) => l.id === selectedId)) {
      setSelectedId(filteredLoci[0].id);
    }
  }, [filteredLoci, selectedId]);
  /* eslint-enable react-hooks/set-state-in-effect */

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
  //
  // The rule reports one position per effect, so it looks like it objects to the
  // early return. It does not: removing that line only moves the report to
  // `setLoading(true)` below. What it objects to is a fetch written as an
  // effect, and every fetch written as an effect announces that it has started
  // before it can await anything. Satisfying it here means moving this
  // component onto Suspense or a fetching library, which is an architecture
  // decision and not a lint fix, so this stays as it is deliberately.
  /* eslint-disable react-hooks/set-state-in-effect -- a fetch effect has to announce that it started */
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
  }, [bamPath, fastaPath, minBaseQuality, minMapq, selectedId, loci]);
  /* eslint-enable react-hooks/set-state-in-effect */

  /* ── Auto-scroll to center variant when data loads ── */
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
    const centerIdx = refColumnIndex(view.columns, center);
    const offsetPx = Math.max(0, centerIdx) * cellSize;
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
  const displayColumns = useMemo(
    () => view ? view.columns : [],
    [view],
  );
  const focusColumnIndexes = useMemo(
    () => displayColumns
      .map((column, idx) => (siteForColumn(column, view?.sites ?? []) ? idx : -1))
      .filter((idx) => idx >= 0),
    [displayColumns, view],
  );
  const trackWidth = displayColumns.length * cellSize;
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
    let startIdx = findCodonStartIndex(
      selectedLocus.refCodon,
      selectedLocus.positions,
      view.reference,
      view.displayStart,
    );
    // Fallback: center the codon on the first SNP position even if we can't
    // verify it against the reference (e.g., minus-strand, IUPAC, edge cases).
    if (startIdx < 0) {
      const minPos = Math.min(...selectedLocus.positions);
      startIdx = Math.max(0, minPos - view.displayStart);
      // Clamp so codon doesn't overflow the reference
      if (startIdx + 3 > view.reference.length) {
        startIdx = Math.max(0, view.reference.length - 3);
      }
    }
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
      if (prev.size === 1 && prev.has(type)) {
        return new Set(SUPPORT_TYPES);
      }
      return new Set([type]);
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
        <span className="step-subtitle">{countLabel(loci.length)} variant loci</span>
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
                    <span className="bam-locus-type">{locusTypeLabel(l)}</span>
                  </div>
                  <div className="bam-locus-item-mid">{l.chrom}:{l.positions.join(", ")}</div>
                  <div className="bam-locus-item-bottom">
                    <span>{l.refBases.join(",")} → {l.altBases.join(",")}</span>
                    {locusEvidenceReads(l) && <span>×{locusEvidenceReads(l)}</span>}
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
                  <span key={`${pos}-${selectedLocus.refBases[i]}-${selectedLocus.altBases[i]}`} className="bam-focus-chip">
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
                    {locusEvidenceFrequency(selectedLocus) && locusEvidenceFrequency(selectedLocus) !== "-" ? `${locusEvidenceFrequency(selectedLocus)} freq · ` : ""}
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
                  <button type="button" className={`bam-count-card${visibleSupport.has("mnv") ? "" : " bam-count-card--off"}`} onClick={() => toggleSupport("mnv")} title={supportButtonTitle(visibleSupport, "mnv", "ALT/event")}>
                    <strong>{countLabel(displayCounts.mnv)}</strong><span>ALT</span>
                  </button>
                  <button type="button" className={`bam-count-card${visibleSupport.has("partial") ? "" : " bam-count-card--off"}`} onClick={() => toggleSupport("partial")} title={supportButtonTitle(visibleSupport, "partial", "SNP/partial")}>
                    <strong>{countLabel(displayCounts.partial)}</strong><span>SNP/partial</span>
                  </button>
                  <button type="button" className={`bam-count-card${visibleSupport.has("reference") ? "" : " bam-count-card--off"}`} onClick={() => toggleSupport("reference")} title={supportButtonTitle(visibleSupport, "reference", "Reference")}>
                    <strong>{countLabel(displayCounts.reference)}</strong><span>Ref</span>
                  </button>
                  <button type="button" className={`bam-count-card${visibleSupport.has("other") ? "" : " bam-count-card--off"}`} onClick={() => toggleSupport("other")} title={supportButtonTitle(visibleSupport, "other", "Other")}>
                    <strong>{countLabel(displayCounts.other)}</strong><span>Other</span>
                  </button>
                </div>

                <div className="bam-legend">
                  <span className="bam-legend-item"><span className="bam-legend-swatch bam-legend-swatch--match" />Match (same as ref)</span>
                  <span className="bam-legend-item"><span className="bam-legend-swatch bam-legend-swatch--a" />A</span>
                  <span className="bam-legend-item"><span className="bam-legend-swatch bam-legend-swatch--t" />T</span>
                  <span className="bam-legend-item"><span className="bam-legend-swatch bam-legend-swatch--g" />G</span>
                  <span className="bam-legend-item"><span className="bam-legend-swatch bam-legend-swatch--c" />C</span>
                  <span className="bam-legend-item"><span className="bam-legend-swatch bam-legend-swatch--gap" />Deletion</span>
                  <span className="bam-legend-item"><span className="bam-legend-swatch bam-legend-swatch--insertion" />Insertion</span>
                  <span className="bam-legend-item"><span className="bam-legend-swatch bam-legend-swatch--variant" />Variant site</span>
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
                          {displayColumns.map((column, i) => {
                            const site = siteForColumn(column, view.sites);
                            const isFocus = Boolean(site);
                            const isInsertion = column.kind === "ins";
                            const isEdge = i === 0 || i === displayColumns.length - 1;
                            const isTick = i % tickStep === 0;
                            if (isInsertion && !isFocus) return null;
                            if (!isFocus && !isEdge && !isTick) return null;
                            // Skip regular ticks too close to a focus position to avoid label overlap.
                            // Label width ≈ digits × 5px; need at least that many cells clearance.
                            if (!isFocus && (isEdge || isTick)) {
                              const minGap = Math.ceil(35 / cellSize);
                              for (const focusIdx of focusColumnIndexes) {
                                if (Math.abs(i - focusIdx) > 0 && Math.abs(i - focusIdx) < minGap) return null;
                              }
                            }
                            return (
                              <span key={column.key} className={`bam-scale-tick${isFocus ? " bam-scale-tick--focus" : ""}${isInsertion ? " bam-scale-tick--insertion" : ""}`} style={{ left: `${i * cellSize}px` }}>
                                <span className="bam-scale-mark" />
                                <span className="bam-scale-label">{isInsertion ? column.label : column.position}</span>
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
                              const columnRange = columnRangeForLocus(displayColumns, l);
                              const left = columnRange.left * cellSize;
                              const w = Math.max(cellSize, columnRange.width * cellSize);
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
                          const column = displayColumns[i];
                          if (!column) return null;
                          const site = siteForColumn(column, view.sites);
                          const label = column.kind === "ins" ? `${column.position}+${column.insertionIndex ?? ""}` : `${column.position}`;
                          return (
                            <span
                              key={column.key}
                              className={`bam-coverage-bar${site ? " bam-coverage-bar--variant" : ""}${column.kind === "ins" ? " bam-coverage-bar--insertion" : ""}`}
                              style={{ width: `${cellSize}px`, height: `${(depth / maxCoverage) * 100}%` }}
                              title={`${label}: ${depth}×`}
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
                        {displayColumns.map((column) => {
                          const base = column.referenceBase;
                          const site = siteForColumn(column, view.sites);
                          const isInsertion = column.kind === "ins";
                          const titlePosition = isInsertion ? `${column.position}+${column.insertionIndex ?? ""}` : `${column.position}`;
                          return (
                            <span key={column.key} className={`bam-cell bam-cell--ref-base ${isInsertion ? " bam-cell--ref-insertion" : nucleotideClass(base)}${site ? " bam-cell--focus" : ""}`} title={`${titlePosition}: ${isInsertion ? "insertion slot" : base}`}>
                              {isInsertion ? "+" : base}
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
                          {displayColumns.map((column) => {
                            const inCodon = column.kind === "ref" && column.position >= codonAnnotation.startPos && column.position < codonAnnotation.startPos + 3;
                            if (!inCodon) {
                              return <span key={column.key} className={`bam-cell bam-cell--codon-empty${column.kind === "ins" ? " bam-cell--insertion" : ""}`} />;
                            }
                            const codonOffset = column.position - codonAnnotation.startPos;
                            const refBase = codonAnnotation.refCodon[codonOffset] ?? "";
                            const mnvBase = codonAnnotation.mnvCodon?.[codonOffset] ?? "";
                            const changed = mnvBase !== "" && refBase.toUpperCase() !== mnvBase.toUpperCase();
                            const isVariant = selectedLocus!.positions.includes(column.position);
                            return (
                              <span
                                key={column.key}
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
                          {displayColumns.map((column) => {
                            const inCodon = column.kind === "ref" && column.position >= codonAnnotation.startPos && column.position < codonAnnotation.startPos + 3;
                            if (!inCodon) {
                              return <span key={column.key} className={`bam-cell bam-cell--codon-empty${column.kind === "ins" ? " bam-cell--insertion" : ""}`} />;
                            }
                            const codonOffset = column.position - codonAnnotation.startPos;
                            const isVariantPos = selectedLocus!.positions.includes(column.position);
                            if (isVariantPos) {
                              const snpIdx = selectedLocus!.positions.indexOf(column.position);
                              const snpCodonList = codonAnnotation.snpCodons.split(",").map((s) => s.trim());
                              const snpCodon = snpCodonList[snpIdx] ?? "";
                              const snpBase = snpCodon[codonOffset] ?? "";
                              const snpAaList = codonAnnotation.snpAaChanges.split(",").map((s) => s.trim());
                              const parsed = parseAaChange(snpAaList[snpIdx] ?? "");
                              const refBase = codonAnnotation.refCodon[codonOffset] ?? "";
                              return (
                                <span
                                  key={column.key}
                                  className="bam-cell bam-cell--codon bam-cell--snp-highlight"
                                  title={`SNP${snpIdx + 1}: REF ${codonAnnotation.refCodon}${parsed ? ` (${parsed.refAa})` : ""} → ALT ${snpCodon}${parsed ? ` (${parsed.altAa})` : snpAaList[snpIdx] ? ` (${snpAaList[snpIdx]})` : ""} | Base: ${refBase}→${snpBase}`}
                                >
                                  {snpBase || selectedLocus!.altBases[snpIdx] || "?"}
                                </span>
                              );
                            }
                            const refBase = codonAnnotation.refCodon[codonOffset] ?? "";
                            return (
                              <span key={column.key} className="bam-cell bam-cell--codon bam-cell--codon-ref">
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
                          {(() => {
                            // Find first and last covered base indices for edge markers
                            let firstIdx = -1, lastIdx = -1;
                            for (let j = 0; j < read.bases.length; j++) {
                              if (read.bases[j] && read.bases[j] !== " ") {
                                if (firstIdx < 0) firstIdx = j;
                                lastIdx = j;
                              }
                            }
                            return displayColumns.map((column, i) => (
                              <BamCell
                                key={`${read.name}-${column.key}`}
                                value={read.bases[i] ?? " "}
                                column={column}
                                site={siteForColumn(column, view.sites)}
                                expectedBase={expectedBaseForColumn(column, siteForColumn(column, view.sites))}
                                isReadStart={i === firstIdx}
                                isReadEnd={i === lastIdx}
                                strand={read.strand}
                              />
                            ));
                          })()}
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
