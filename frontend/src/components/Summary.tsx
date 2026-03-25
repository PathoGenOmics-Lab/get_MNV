import type { RunSummary } from "../types";

interface SummaryProps {
  result: RunSummary;
}

/** Format numbers with English locale (dots for thousands, dots for decimals) */
const fmtN = (n: number) => n.toLocaleString("en-US");

function formatMs(ms: number): string {
  if (ms < 1000) return `${ms.toFixed(0)} ms`;
  return `${(ms / 1000).toLocaleString("en-US", { minimumFractionDigits: 2, maximumFractionDigits: 2 })} s`;
}

/** Compute percentage safely */
function pct(part: number, total: number): string {
  if (total === 0) return "—";
  return `${((part / total) * 100).toLocaleString("en-US", { maximumFractionDigits: 1 })}%`;
}

export default function Summary({ result }: SummaryProps) {
  const { global: g, timings, contigs, sample, bam_provided } = result;

  return (
    <div className="summary">
      <div className="step-header" style={{ marginBottom: "0.75rem" }}>
        <div className="step-badge step-badge--done">{"\u2713"}</div>
        <h3 className="step-title">Analysis Complete</h3>
        <span className="step-subtitle">{formatMs(timings.total_ms)}</span>
      </div>

      {/* ── Quick stats banner ── */}
      <div className="summary-quick-stats">
        <div className="quick-stat">
          <span className="quick-stat-value quick-stat-value--m">{fmtN(g.produced_variants)}</span>
          <span className="quick-stat-label">Variants</span>
        </div>
        <div className="quick-stat-divider" />
        <div className="quick-stat">
          <span className="quick-stat-value quick-stat-value--n">{fmtN(g.snp_records_in_vcf)}</span>
          <span className="quick-stat-label">VCF Records</span>
        </div>
        <div className="quick-stat-divider" />
        <div className="quick-stat">
          <span className="quick-stat-value quick-stat-value--v">{fmtN(g.mapped_genes)}</span>
          <span className="quick-stat-label">Genes Mapped</span>
        </div>
        <div className="quick-stat-divider" />
        <div className="quick-stat">
          <span className="quick-stat-value">{fmtN(g.contig_count)}</span>
          <span className="quick-stat-label">Contigs</span>
        </div>
      </div>

      <div className="summary-grid">
        {/* Teal accent — M for variants */}
        <div className="summary-card summary-card--variants">
          <h4>
            <span className="summary-card-icon summary-card-icon--m">
              <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round"><circle cx="12" cy="12" r="3"/><path d="M12 2v4m0 12v4m-7.07-3.93l2.83-2.83m8.48-8.48l2.83-2.83M2 12h4m12 0h4M4.93 4.93l2.83 2.83m8.48 8.48l2.83 2.83"/></svg>
            </span>
            Variant Breakdown
          </h4>
          <table>
            <tbody>
              <tr className="stat-highlight">
                <td>Total produced</td>
                <td>{fmtN(g.produced_variants)}</td>
              </tr>
              <tr>
                <td>SNP</td>
                <td>
                  <span className="stat-with-pct">
                    {fmtN(g.snp_variants)}
                    <span className="stat-pct">{pct(g.snp_variants, g.produced_variants)}</span>
                  </span>
                </td>
              </tr>
              <tr>
                <td>MNV</td>
                <td>
                  <span className="stat-with-pct">
                    {fmtN(g.mnv_variants)}
                    <span className="stat-pct">{pct(g.mnv_variants, g.produced_variants)}</span>
                  </span>
                </td>
              </tr>
              <tr>
                <td>SNP/MNV</td>
                <td>
                  <span className="stat-with-pct">
                    {fmtN(g.snp_mnv_variants)}
                    <span className="stat-pct">{pct(g.snp_mnv_variants, g.produced_variants)}</span>
                  </span>
                </td>
              </tr>
              <tr>
                <td>Indel</td>
                <td>
                  <span className="stat-with-pct">
                    {fmtN(g.indel_variants)}
                    <span className="stat-pct">{pct(g.indel_variants, g.produced_variants)}</span>
                  </span>
                </td>
              </tr>
            </tbody>
          </table>
        </div>

        {/* Magenta accent — N for inputs */}
        <div className="summary-card summary-card--inputs">
          <h4>
            <span className="summary-card-icon summary-card-icon--n">
              <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><path d="M14 2v6h6"/></svg>
            </span>
            Input Summary
          </h4>
          <table>
            <tbody>
              <tr>
                <td>VCF records</td>
                <td>{fmtN(g.snp_records_in_vcf)}</td>
              </tr>
              <tr>
                <td>Contigs</td>
                <td>{fmtN(g.contig_count)}</td>
              </tr>
              <tr>
                <td>Mapped genes</td>
                <td>{fmtN(g.mapped_genes)}</td>
              </tr>
              {sample && (
                <tr>
                  <td>Sample</td>
                  <td className="stat-mono">{sample}</td>
                </tr>
              )}
              <tr>
                <td>BAM provided</td>
                <td>{bam_provided ? "Yes" : "No"}</td>
              </tr>
              <tr>
                <td>Cache hit rate</td>
                <td>
                  {g.region_cache_hits + g.region_cache_misses > 0
                    ? pct(g.region_cache_hits, g.region_cache_hits + g.region_cache_misses)
                    : "—"}
                </td>
              </tr>
            </tbody>
          </table>
        </div>

        {/* Purple accent — V for performance */}
        <div className="summary-card summary-card--perf">
          <h4>
            <span className="summary-card-icon summary-card-icon--v">
              <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round"><circle cx="12" cy="12" r="10"/><polyline points="12 6 12 12 16 14"/></svg>
            </span>
            Performance
          </h4>
          <table>
            <tbody>
              <tr>
                <td>Parse inputs</td>
                <td>{formatMs(timings.parse_inputs_ms)}</td>
              </tr>
              <tr>
                <td>Processing</td>
                <td>{formatMs(timings.process_ms)}</td>
              </tr>
              <tr>
                <td>Output</td>
                <td>{formatMs(timings.emit_ms)}</td>
              </tr>
              <tr className="stat-highlight">
                <td><strong>Total</strong></td>
                <td><strong>{formatMs(timings.total_ms)}</strong></td>
              </tr>
            </tbody>
          </table>
          {/* Timing bar visualization */}
          <div className="timing-bar">
            <div
              className="timing-bar-seg timing-bar-seg--parse"
              style={{ flex: timings.parse_inputs_ms || 1 }}
              title={`Parse: ${formatMs(timings.parse_inputs_ms)}`}
            />
            <div
              className="timing-bar-seg timing-bar-seg--process"
              style={{ flex: timings.process_ms || 1 }}
              title={`Process: ${formatMs(timings.process_ms)}`}
            />
            <div
              className="timing-bar-seg timing-bar-seg--emit"
              style={{ flex: timings.emit_ms || 1 }}
              title={`Output: ${formatMs(timings.emit_ms)}`}
            />
          </div>
          <div className="timing-bar-legend">
            <span className="timing-legend-item"><span className="timing-dot timing-dot--parse" /> Parse</span>
            <span className="timing-legend-item"><span className="timing-dot timing-dot--process" /> Process</span>
            <span className="timing-legend-item"><span className="timing-dot timing-dot--emit" /> Output</span>
          </div>
        </div>
      </div>

      {contigs.length > 0 && (
        <div className="contig-table-wrapper">
          <div className="contig-table-header">
            <h4>Per-Contig Breakdown</h4>
            <span className="contig-table-count">{fmtN(contigs.length)} contigs</span>
          </div>
          <div className="contig-table-scroll">
            <table className="contig-table">
              <thead>
                <tr>
                  <th>Contig</th>
                  <th>VCF Records</th>
                  <th>Genes</th>
                  <th>Variants</th>
                  <th>SNP</th>
                  <th>MNV</th>
                  <th>SNP/MNV</th>
                  <th>Indel</th>
                </tr>
              </thead>
              <tbody>
                {contigs.map((c, i) => (
                  <tr key={c.contig} className={i % 2 === 1 ? "contig-row--alt" : ""}>
                    <td className="contig-name">{c.contig}</td>
                    <td>{fmtN(c.snp_records_in_vcf)}</td>
                    <td>{fmtN(c.mapped_genes)}</td>
                    <td className="contig-total">{fmtN(c.produced_variants)}</td>
                    <td>{fmtN(c.snp_variants)}</td>
                    <td>{fmtN(c.mnv_variants)}</td>
                    <td>{fmtN(c.snp_mnv_variants)}</td>
                    <td>{fmtN(c.indel_variants)}</td>
                  </tr>
                ))}
              </tbody>
              {contigs.length > 1 && (
                <tfoot>
                  <tr className="contig-totals-row">
                    <td><strong>Total</strong></td>
                    <td><strong>{fmtN(g.snp_records_in_vcf)}</strong></td>
                    <td><strong>{fmtN(g.mapped_genes)}</strong></td>
                    <td><strong>{fmtN(g.produced_variants)}</strong></td>
                    <td><strong>{fmtN(g.snp_variants)}</strong></td>
                    <td><strong>{fmtN(g.mnv_variants)}</strong></td>
                    <td><strong>{fmtN(g.snp_mnv_variants)}</strong></td>
                    <td><strong>{fmtN(g.indel_variants)}</strong></td>
                  </tr>
                </tfoot>
              )}
            </table>
          </div>
        </div>
      )}
    </div>
  );
}
