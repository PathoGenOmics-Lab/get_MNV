import { useState, useMemo, useCallback, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
import { save } from "@tauri-apps/plugin-dialog";
import type { TsvData } from "../types";

interface VariantTableProps {
  data: TsvData;
}

const PAGE_SIZE = 50;
const MAX_DROPDOWN_OPTIONS = 30; // columns with ≤ this many unique values get a dropdown

/** Map variant type strings to CSS class for color-coded badges */
function variantTypeClass(val: string): string {
  const v = val.toLowerCase();
  if (v === "mnv") return "vt-badge vt-badge--mnv";
  if (v === "snpmnv" || v === "snp/mnv") return "vt-badge vt-badge--snpmnv";
  if (v === "snp") return "vt-badge vt-badge--snp";
  if (v === "indel") return "vt-badge vt-badge--indel";
  return "";
}

/** Map change type strings to CSS class */
function changeTypeClass(val: string): string {
  const v = val.toLowerCase();
  if (v.includes("synonymous") && !v.includes("non")) return "ct-badge ct-badge--syn";
  if (v.includes("nonsynonymous") || v.includes("non-synonymous") || v.includes("non_synonymous"))
    return "ct-badge ct-badge--nonsyn";
  if (v.includes("stopgained") || v.includes("stop_gained")) return "ct-badge ct-badge--stop";
  if (v.includes("stoplost") || v.includes("stop_lost")) return "ct-badge ct-badge--stoplost";
  if (v.includes("frameshift")) return "ct-badge ct-badge--frameshift";
  if (v.includes("indel")) return "ct-badge ct-badge--indel";
  return "";
}

/** Is this column a "special" column that gets a badge? */
function isVariantTypeCol(header: string): boolean {
  return header.toLowerCase() === "variant type";
}

/** Check if a comma-separated read-count field has any non-zero value */
function hasReads(row: string[], colIdx: number): boolean {
  if (colIdx < 0) return false;
  const v = row[colIdx] ?? "";
  return v.split(",").some((s) => { const n = parseInt(s.trim(), 10); return Number.isFinite(n) && n > 0; });
}

/** Reclassify SNP/MNV variant type based on actual read support */
function reclassifyVariantType(raw: string, row: string[], snpReadsIdx: number, mnvReadsIdx: number): string {
  const v = raw.toLowerCase();
  if (v === "snp/mnv" || v === "snpmnv") {
    const hasMnv = mnvReadsIdx >= 0 ? hasReads(row, mnvReadsIdx) : true;
    const hasSnp = snpReadsIdx >= 0 ? hasReads(row, snpReadsIdx) : true;
    if (hasMnv && !hasSnp) return "MNV";
    if (hasSnp && !hasMnv) return "SNP";
  }
  return raw;
}
function isChangeTypeCol(header: string): boolean {
  return header.toLowerCase() === "change type";
}
function isCodonCol(header: string): boolean {
  const h = header.toLowerCase();
  return h.includes("codon") || h === "base changes";
}
/** Classify codon columns by their origin: only MNV or SNP get a tag */
function codonColType(header: string): "mnv" | "snp" | false {
  const h = header.toLowerCase();
  if (h === "mnv codon") return "mnv";
  if (h === "snp codon") return "snp";
  return false;
}
function isAACol(header: string): "mnv" | "snp" | false {
  const h = header.toLowerCase();
  if (h === "aa changes") return "mnv";
  if (h === "snp aa changes") return "snp";
  return false;
}
function isNumericCol(header: string): boolean {
  const h = header.toLowerCase();
  return (
    h.includes("reads") ||
    h.includes("frequencies") ||
    h === "total reads" ||
    h === "positions"
  );
}

/** Determine if a column is better served by a dropdown (categorical) vs text input */
function isCategoricalCol(header: string): boolean {
  const h = header.toLowerCase();
  return (
    h === "variant type" ||
    h === "change type" ||
    h === "chromosome"
  );
}

export default function VariantTable({ data }: VariantTableProps) {
  const [search, setSearch] = useState("");
  const [sortCol, setSortCol] = useState<number | null>(null);
  const [sortAsc, setSortAsc] = useState(true);
  const [page, setPage] = useState(0);
  // Per-column filters: colIdx → filter value
  const [colFilters, setColFilters] = useState<Record<number, string>>({});
  const [expanded, setExpanded] = useState(false);

  // ESC key closes fullscreen
  useEffect(() => {
    if (!expanded) return;
    const handler = (e: KeyboardEvent) => {
      if (e.key === "Escape") setExpanded(false);
    };
    document.addEventListener("keydown", handler);
    return () => document.removeEventListener("keydown", handler);
  }, [expanded]);

  const { headers, rows } = data;

  // Column indices for read-count gating (only show MNV/SNP tags when reads > 0)
  const mnvReadsIdx = useMemo(() => headers.findIndex((h) => h.trim().toLowerCase() === "mnv reads"), [headers]);
  const snpReadsIdx = useMemo(() => headers.findIndex((h) => h.trim().toLowerCase() === "snp reads"), [headers]);
  const vtColIdx = useMemo(() => headers.findIndex((h) => h.trim().toLowerCase() === "variant type"), [headers]);

  /** Get the effective (reclassified) cell value for a given row and column */
  const effectiveCellValue = useCallback((row: string[], colIdx: number): string => {
    const raw = row[colIdx] ?? "";
    if (colIdx === vtColIdx) return reclassifyVariantType(raw, row, snpReadsIdx, mnvReadsIdx);
    return raw;
  }, [vtColIdx, snpReadsIdx, mnvReadsIdx]);

  // Compute unique values per column for dropdowns (using effective/reclassified values)
  const columnUniques = useMemo(() => {
    const result: Record<number, string[]> = {};
    headers.forEach((h, i) => {
      // Only compute for columns that could be dropdowns
      if (isCategoricalCol(h) || !isNumericCol(h)) {
        const unique = new Set<string>();
        for (const row of rows) {
          const v = effectiveCellValue(row, i);
          if (v && v !== "-") unique.add(v);
        }
        if (unique.size <= MAX_DROPDOWN_OPTIONS) {
          result[i] = Array.from(unique).sort();
        }
      }
    });
    return result;
  }, [headers, rows, effectiveCellValue]);

  const activeFilterCount = useMemo(
    () => Object.values(colFilters).filter((v) => v.length > 0).length,
    [colFilters],
  );

  const setColFilter = useCallback((colIdx: number, value: string) => {
    setColFilters((prev) => {
      const next = { ...prev };
      if (value === "") {
        delete next[colIdx];
      } else {
        next[colIdx] = value;
      }
      return next;
    });
    setPage(0);
  }, []);

  const clearAllFilters = useCallback(() => {
    setColFilters({});
    setSearch("");
    setPage(0);
  }, []);

  // Filter rows — global search + per-column filters (uses reclassified values)
  const filtered = useMemo(() => {
    let result = rows;

    // Per-column filters
    const activeFilters = Object.entries(colFilters).filter(([, v]) => v.length > 0);
    if (activeFilters.length > 0) {
      result = result.filter((row) =>
        activeFilters.every(([colStr, filterVal]) => {
          const col = parseInt(colStr, 10);
          const cellVal = effectiveCellValue(row, col);
          // If this is a dropdown column, match exactly
          if (columnUniques[col]) {
            return cellVal === filterVal;
          }
          // Otherwise text substring match (case-insensitive)
          return cellVal.toLowerCase().includes(filterVal.toLowerCase());
        }),
      );
    }

    // Global text search
    if (search.trim()) {
      const q = search.toLowerCase();
      result = result.filter((r) =>
        r.some((_, i) => effectiveCellValue(r, i).toLowerCase().includes(q)),
      );
    }
    return result;
  }, [rows, search, colFilters, columnUniques, effectiveCellValue]);

  // Sort rows
  const sorted = useMemo(() => {
    if (sortCol === null) return filtered;
    const col = sortCol;
    const isNum = isNumericCol(headers[col]);
    return [...filtered].sort((a, b) => {
      const va = effectiveCellValue(a, col);
      const vb = effectiveCellValue(b, col);
      if (isNum) {
        const na = parseFloat(va) || 0;
        const nb = parseFloat(vb) || 0;
        return sortAsc ? na - nb : nb - na;
      }
      return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
    });
  }, [filtered, sortCol, sortAsc, headers, effectiveCellValue]);

  // Pagination
  const totalPages = Math.max(1, Math.ceil(sorted.length / PAGE_SIZE));
  const safePage = Math.min(page, totalPages - 1);
  const pageRows = sorted.slice(safePage * PAGE_SIZE, (safePage + 1) * PAGE_SIZE);

  const handleSort = (col: number) => {
    if (sortCol === col) {
      setSortAsc(!sortAsc);
    } else {
      setSortCol(col);
      setSortAsc(true);
    }
  };

  const handleExport = async () => {
    const filePath = await save({
      filters: [
        { name: "TSV", extensions: ["tsv"] },
        { name: "CSV", extensions: ["csv"] },
      ],
      defaultPath: "variants_filtered.tsv",
    });
    if (!filePath) return;
    const sep = filePath.toLowerCase().endsWith(".csv") ? "," : "\t";
    const lines = [
      headers.join(sep),
      ...sorted.map((row) =>
        headers.map((_, ci) => effectiveCellValue(row, ci)).join(sep)
      ),
    ];
    await invoke("write_text_file", { path: filePath, content: lines.join("\n") });
  };

  /** Render a cell with appropriate formatting */
  const renderCell = (value: string, colIdx: number, row: string[]) => {
    if (isVariantTypeCol(headers[colIdx])) {
      const displayValue = reclassifyVariantType(value, row, snpReadsIdx, mnvReadsIdx);
      const cls = variantTypeClass(displayValue);
      return cls ? <span className={cls}>{displayValue}</span> : displayValue;
    }
    if (isChangeTypeCol(headers[colIdx])) {
      const cls = changeTypeClass(value);
      return cls ? <span className={cls}>{value}</span> : value;
    }
    const aaType = isAACol(headers[colIdx]);
    if (aaType && value && value !== "-") {
      // Only show MNV/SNP tag if there are actual supporting reads
      const showTag =
        (aaType === "mnv" && hasReads(row, mnvReadsIdx)) ||
        (aaType === "snp" && hasReads(row, snpReadsIdx));
      const parts = value.split(/,\s*/);
      return (
        <span className="vt-aa-cell">
          {showTag && (
            <span className={`vt-aa-tag vt-aa-tag--${aaType}`}>{aaType === "mnv" ? "MNV" : "SNP"}</span>
          )}
          {parts.map((p, i) => (
            <span key={i} className={showTag ? `vt-aa-value vt-aa-value--${aaType}` : "vt-aa-value"}>
              {p}
            </span>
          ))}
        </span>
      );
    }
    if (isCodonCol(headers[colIdx])) {
      const cType = codonColType(headers[colIdx]);
      if (cType && value && value !== "-") {
        // Only show MNV/SNP tag if there are actual supporting reads
        const showTag =
          (cType === "mnv" && hasReads(row, mnvReadsIdx)) ||
          (cType === "snp" && hasReads(row, snpReadsIdx));
        const parts = value.split(/,\s*/);
        return (
          <span className="vt-aa-cell">
            {showTag && (
              <span className={`vt-codon-tag vt-codon-tag--${cType}`}>
                {cType === "mnv" ? "MNV" : "SNP"}
              </span>
            )}
            {parts.map((p, i) => (
              <span key={i} className={showTag ? `vt-codon vt-codon--${cType}` : "vt-codon"}>{p}</span>
            ))}
          </span>
        );
      }
      return <span className="vt-codon">{value}</span>;
    }
    return value;
  };

  /** Render per-column filter input */
  const renderColumnFilter = (colIdx: number) => {
    const uniqueVals = columnUniques[colIdx];
    const currentVal = colFilters[colIdx] ?? "";

    // Dropdown for categorical / low-cardinality columns
    if (uniqueVals) {
      return (
        <select
          className="vt-col-filter-select"
          value={currentVal}
          onChange={(e) => setColFilter(colIdx, e.target.value)}
          title={`Filter ${headers[colIdx]}`}
        >
          <option value="">All</option>
          {uniqueVals.map((v) => (
            <option key={v} value={v}>
              {v}
            </option>
          ))}
        </select>
      );
    }

    // Text input for other columns
    return (
      <input
        type="text"
        className="vt-col-filter-input"
        placeholder="Filter..."
        value={currentVal}
        onChange={(e) => setColFilter(colIdx, e.target.value)}
        title={`Filter ${headers[colIdx]}`}
      />
    );
  };

  return (
    <div className={`variant-table-container${expanded ? " variant-table-container--expanded" : ""}`}>
      {/* ── Toolbar ── */}
      <div className="vt-toolbar">
        <div className="vt-search-wrapper">
          <svg className="vt-search-icon" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
            <circle cx="11" cy="11" r="8" />
            <path d="M21 21l-4.35-4.35" />
          </svg>
          <input
            type="text"
            className="vt-search"
            placeholder="Search all columns..."
            value={search}
            onChange={(e) => { setSearch(e.target.value); setPage(0); }}
          />
          {search && (
            <button className="vt-search-clear" onClick={() => { setSearch(""); setPage(0); }} title="Clear search">
              <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round"><path d="M18 6L6 18M6 6l12 12" /></svg>
            </button>
          )}
        </div>
        {activeFilterCount > 0 && (
          <button className="vt-clear-filters-btn" onClick={clearAllFilters} title="Clear all filters">
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round"><path d="M18 6L6 18M6 6l12 12" /></svg>
            Clear {activeFilterCount} filter{activeFilterCount !== 1 ? "s" : ""}
          </button>
        )}
        <span className="vt-count">
          {sorted.length.toLocaleString("en-US")} variant{sorted.length !== 1 ? "s" : ""}
          {sorted.length !== rows.length && (
            <span className="vt-count-total"> of {rows.length.toLocaleString("en-US")}</span>
          )}
        </span>
        <button
          className="vt-export-btn"
          onClick={handleExport}
          title="Export filtered data as TSV or CSV"
        >
          <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round">
            <path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4" />
            <polyline points="7 10 12 15 17 10" />
            <line x1="12" y1="15" x2="12" y2="3" />
          </svg>
          Export
        </button>
        <button
          className="vt-expand-btn"
          onClick={() => setExpanded(!expanded)}
          title={expanded ? "Exit fullscreen" : "Expand to fullscreen"}
        >
          {expanded ? (
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
              <path d="M8 3v3a2 2 0 01-2 2H3m18 0h-3a2 2 0 01-2-2V3m0 18v-3a2 2 0 012-2h3M3 16h3a2 2 0 012 2v3" />
            </svg>
          ) : (
            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round">
              <path d="M15 3h6v6M9 21H3v-6M21 3l-7 7M3 21l7-7" />
            </svg>
          )}
        </button>
      </div>

      {/* ── Table ── */}
      <div className="vt-scroll">
        <table className="vt-table">
          <thead>
            <tr>
              {headers.map((h, i) => (
                <th
                  key={i}
                  className={`vt-th${sortCol === i ? " vt-th--sorted" : ""}${isNumericCol(h) ? " vt-th--numeric" : ""}`}
                  onClick={() => handleSort(i)}
                >
                  <span className="vt-th-label">{h}</span>
                  {sortCol === i && (
                    <span className="vt-sort-arrow">{sortAsc ? "\u25B2" : "\u25BC"}</span>
                  )}
                </th>
              ))}
            </tr>
            {/* ── Per-column filter row ── */}
            <tr className="vt-filter-row">
              {headers.map((_, i) => (
                <th key={i} className="vt-filter-cell">
                  {renderColumnFilter(i)}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {pageRows.length === 0 ? (
              <tr>
                <td colSpan={headers.length} className="vt-empty">
                  {search || activeFilterCount > 0 ? "No variants match your filters" : "No variants found"}
                </td>
              </tr>
            ) : (
              pageRows.map((row, ri) => (
                <tr key={ri} className={ri % 2 === 1 ? "vt-row--alt" : ""}>
                  {headers.map((_, ci) => (
                    <td
                      key={ci}
                      className={`${isNumericCol(headers[ci]) ? "vt-td--numeric" : ""}${isCodonCol(headers[ci]) ? " vt-td--codon" : ""}`}
                    >
                      {renderCell(row[ci] ?? "", ci, row)}
                    </td>
                  ))}
                </tr>
              ))
            )}
          </tbody>
        </table>
      </div>

      {/* ── Pagination ── */}
      {totalPages > 1 && (
        <div className="vt-pagination">
          <button
            className="vt-page-btn"
            disabled={safePage === 0}
            onClick={() => setPage(0)}
            title="First page"
          >
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round"><path d="M11 17l-5-5 5-5M18 17l-5-5 5-5" /></svg>
          </button>
          <button
            className="vt-page-btn"
            disabled={safePage === 0}
            onClick={() => setPage(safePage - 1)}
            title="Previous page"
          >
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round"><path d="M15 18l-6-6 6-6" /></svg>
          </button>
          <span className="vt-page-info">
            Page {safePage + 1} of {totalPages}
          </span>
          <button
            className="vt-page-btn"
            disabled={safePage >= totalPages - 1}
            onClick={() => setPage(safePage + 1)}
            title="Next page"
          >
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round"><path d="M9 18l6-6-6-6" /></svg>
          </button>
          <button
            className="vt-page-btn"
            disabled={safePage >= totalPages - 1}
            onClick={() => setPage(totalPages - 1)}
            title="Last page"
          >
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round"><path d="M13 17l5-5-5-5M6 17l5-5-5-5" /></svg>
          </button>
        </div>
      )}
    </div>
  );
}
