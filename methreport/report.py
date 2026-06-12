"""Generate a self-contained interactive HTML report using Plotly + Jinja2."""

from __future__ import annotations

import json
import logging
from datetime import datetime
from importlib import resources
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from methreport import __version__
from methreport.analysis import RegionResult, SampleAnalysis
from methreport.controls import FLAG_HIGH, FLAG_LOW
from methreport.reader import RegionMethylation

log = logging.getLogger(__name__)

# Plotly color palette
COLORS = {
    "unphased": "#4C72B0",
    "hp1": "#DD4444",
    "hp2": "#2196F3",
    "control_band": "rgba(150,150,150,0.25)",
    "control_line": "rgba(100,100,100,0.6)",
    "flag_low": "#E53935",
    "flag_high": "#E53935",
    "normal": "#43A047",
    "na": "#9E9E9E",
    "bg": "#F8F9FA",
    "card_bg": "#FFFFFF",
}

FLAG_COLORS = {
    "NORMAL": "#43A047",
    "LOW": "#E53935",
    "HIGH": "#E53935",
    "NA": "#9E9E9E",
}

FLAG_BADGES = {
    "NORMAL": '<span class="badge badge-normal">NORMAL</span>',
    "LOW": '<span class="badge badge-abnormal">LOW ⚠</span>',
    "HIGH": '<span class="badge badge-abnormal">HIGH ⚠</span>',
    "NA": '<span class="badge badge-na">NO DATA</span>',
}


def _make_region_figure(result: RegionResult) -> go.Figure:
    """Create a Plotly figure for one DMR — unphased + phased subplots."""
    has_hp = result.is_phased
    n_cols = 2 if has_hp else 1
    subplot_titles = [f"{result.region.label} — Unphased"]
    if has_hp:
        subplot_titles.append(f"{result.region.label} — Phased (HP1/HP2)")

    fig = make_subplots(
        rows=1,
        cols=n_cols,
        subplot_titles=subplot_titles,
        horizontal_spacing=0.08,
    )

    # --- Unphased panel ---
    _add_methylation_trace(
        fig, result.unphased, result.reference,
        row=1, col=1,
        color=COLORS["unphased"],
        show_controls=True,
    )

    # --- Phased panel ---
    if has_hp:
        _add_methylation_trace(
            fig, result.hp1, result.reference,
            row=1, col=2,
            color=COLORS["hp1"],
            name="HP1",
            show_controls=True,
        )
        _add_methylation_trace(
            fig, result.hp2, result.reference,
            row=1, col=2,
            color=COLORS["hp2"],
            name="HP2",
            show_controls=False,
        )

    fig.update_yaxes(range=[-5, 105], title_text="Methylation (%)", row=1, col=1)
    if has_hp:
        fig.update_yaxes(range=[-5, 105], row=1, col=2)

    fig.update_layout(
        height=350,
        margin=dict(l=50, r=30, t=50, b=40),
        paper_bgcolor="white",
        plot_bgcolor=COLORS["bg"],
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        font=dict(family="Inter, Arial, sans-serif", size=12),
        hovermode="x unified",
    )
    return fig


def _add_methylation_trace(
    fig: go.Figure,
    rm: RegionMethylation,
    reference,
    row: int,
    col: int,
    color: str,
    name: str | None = None,
    show_controls: bool = True,
) -> None:
    label = name or "Sample"

    # Control reference band (±1 SD)
    if show_controls and not reference.empty:
        ref_pos = reference["position"].tolist()
        upper = reference["upper_1sd"].tolist()
        lower = reference["lower_1sd"].tolist()

        fig.add_trace(
            go.Scatter(
                x=ref_pos + ref_pos[::-1],
                y=upper + lower[::-1],
                fill="toself",
                fillcolor=COLORS["control_band"],
                line=dict(color="rgba(0,0,0,0)"),
                name="Control ±1 SD",
                showlegend=(row == 1 and col == 1),
                hoverinfo="skip",
            ),
            row=row,
            col=col,
        )
        # Control mean line
        fig.add_trace(
            go.Scatter(
                x=ref_pos,
                y=reference["mean"].tolist(),
                mode="lines",
                line=dict(color=COLORS["control_line"], width=1, dash="dash"),
                name="Control mean",
                showlegend=(row == 1 and col == 1),
                hoverinfo="skip",
            ),
            row=row,
            col=col,
        )

    if rm.n_cpg == 0:
        return

    df = rm.to_dataframe().sort_values("position")
    positions = df["position"].tolist()
    methyl_pct = df["methylation_pct"].tolist()
    coverages = df["n_total"].tolist()

    # Scatter + smoothed line
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=methyl_pct,
            mode="markers",
            marker=dict(color=color, size=5, opacity=0.6),
            name=f"{label} CpGs",
            customdata=list(zip(coverages, methyl_pct)),
            hovertemplate=(
                "Pos: %{x:,}<br>"
                "Methylation: %{customdata[1]:.1f}%<br>"
                "Coverage: %{customdata[0]}x<extra></extra>"
            ),
            showlegend=True,
        ),
        row=row,
        col=col,
    )

    # LOWESS-style smoothing via rolling mean (requires ≥5 points)
    if len(methyl_pct) >= 5:
        import pandas as pd
        s = pd.Series(methyl_pct, index=positions)
        window = max(3, len(methyl_pct) // 5)
        smoothed = s.rolling(window=window, center=True, min_periods=1).mean()
        fig.add_trace(
            go.Scatter(
                x=smoothed.index.tolist(),
                y=smoothed.tolist(),
                mode="lines",
                line=dict(color=color, width=2.5),
                name=f"{label} (smoothed)",
                hoverinfo="skip",
                showlegend=True,
            ),
            row=row,
            col=col,
        )


def _build_summary_table_html(analysis: SampleAnalysis) -> str:
    rows_html = []
    for r in analysis.results:
        meth = r.mean_methylation
        meth_str = f"{meth:.1f}%" if not np.isnan(meth) else "—"
        hp1_str = f"{r.mean_methylation_hp1:.1f}%" if not np.isnan(r.mean_methylation_hp1) else "—"
        hp2_str = f"{r.mean_methylation_hp2:.1f}%" if not np.isnan(r.mean_methylation_hp2) else "—"
        badge = FLAG_BADGES.get(r.flag, FLAG_BADGES["NA"])
        flag_cls = "row-abnormal" if r.flag in ("LOW", "HIGH") else ""
        rows_html.append(f"""
        <tr class="{flag_cls}">
            <td><strong>{r.region.label}</strong></td>
            <td>{r.region.disorder}</td>
            <td class="mono">{r.region.chrom}:{r.region.start:,}–{r.region.end:,}</td>
            <td>{r.unphased.n_cpg}</td>
            <td>{r.unphased.mean_coverage:.0f}×</td>
            <td>{meth_str}</td>
            <td>{hp1_str}</td>
            <td>{hp2_str}</td>
            <td>{badge}</td>
        </tr>""")
    return "\n".join(rows_html)


def generate_report(
    analysis: SampleAnalysis,
    out_dir: Path,
    run_metadata: dict | None = None,
) -> Path:
    """Render the full interactive HTML report."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{analysis.sample_id}_methreport.html"

    # Build per-region figures and serialize to JSON for embedding
    region_plots = []
    for result in analysis.results:
        fig = _make_region_figure(result)
        region_plots.append({
            "name": result.region.name,
            "label": result.region.label,
            "disorder": result.region.disorder,
            "flag": result.flag,
            "plot_json": fig.to_json(),
            "n_cpg": result.unphased.n_cpg,
            "mean_cov": f"{result.unphased.mean_coverage:.1f}",
            "mean_meth": f"{result.mean_methylation:.1f}" if not np.isnan(result.mean_methylation) else "N/A",
            "is_phased": result.is_phased,
        })

    summary_rows_html = _build_summary_table_html(analysis)
    n_flagged = len(analysis.flagged_regions)
    n_total = len(analysis.results)
    report_date = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Group regions by disorder for navigation
    disorder_groups: dict[str, list[dict]] = {}
    for p in region_plots:
        disorder_groups.setdefault(p["disorder"], []).append(p)

    # Serialize plot data
    plots_json = json.dumps(region_plots)

    html = _render_html(
        sample_id=analysis.sample_id,
        genome=analysis.genome,
        report_date=report_date,
        version=__version__,
        n_flagged=n_flagged,
        n_total=n_total,
        summary_rows_html=summary_rows_html,
        disorder_groups=disorder_groups,
        plots_json=plots_json,
        run_metadata=run_metadata or {},
    )

    out_path.write_text(html, encoding="utf-8")
    log.info("HTML report → %s", out_path)
    return out_path


def _render_html(
    sample_id: str,
    genome: str,
    report_date: str,
    version: str,
    n_flagged: int,
    n_total: int,
    summary_rows_html: str,
    disorder_groups: dict,
    plots_json: str,
    run_metadata: dict,
) -> str:
    status_color = "#E53935" if n_flagged > 0 else "#43A047"
    status_text = f"{n_flagged} region(s) FLAGGED" if n_flagged > 0 else "All regions NORMAL"

    nav_items = "\n".join(
        f'<li><a href="#disorder-{d.replace("/", "-").replace(" ", "-")}">{d}</a></li>'
        for d in disorder_groups
    )

    disorder_sections = ""
    for disorder, plots in disorder_groups.items():
        anchor = disorder.replace("/", "-").replace(" ", "-")
        cards = ""
        for p in plots:
            flag_cls = "card-abnormal" if p["flag"] in ("LOW", "HIGH") else ""
            badge = FLAG_BADGES.get(p["flag"], FLAG_BADGES["NA"])
            cards += f"""
            <div class="region-card {flag_cls}" id="region-{p['name']}">
                <div class="card-header">
                    <span class="card-title">{p['label']}</span>
                    {badge}
                    <span class="card-meta">{p['n_cpg']} CpGs · {p['mean_cov']}× coverage · {p['mean_meth']}% methylation</span>
                </div>
                <div class="plot-container" id="plot-{p['name']}"></div>
            </div>"""

        disorder_sections += f"""
        <section class="disorder-section" id="disorder-{anchor}">
            <h2 class="disorder-title">{disorder}</h2>
            {cards}
        </section>"""

    meta_rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in run_metadata.items()
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>MethReport — {sample_id}</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>
<style>
  :root {{
    --primary: #1A237E;
    --primary-light: #3949AB;
    --accent: #E53935;
    --normal: #43A047;
    --bg: #F0F2F5;
    --card: #FFFFFF;
    --border: #DDE1E7;
    --text: #1C1C1E;
    --text-muted: #6B7280;
    --mono: 'JetBrains Mono', 'Fira Code', monospace;
  }}

  * {{ box-sizing: border-box; margin: 0; padding: 0; }}

  body {{
    font-family: 'Inter', 'Segoe UI', Arial, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.5;
  }}

  /* Top header */
  .page-header {{
    background: linear-gradient(135deg, var(--primary) 0%, var(--primary-light) 100%);
    color: white;
    padding: 24px 40px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    box-shadow: 0 2px 8px rgba(0,0,0,0.2);
  }}
  .header-title {{ font-size: 1.6rem; font-weight: 700; letter-spacing: -0.5px; }}
  .header-subtitle {{ font-size: 0.85rem; opacity: 0.8; margin-top: 4px; }}
  .header-meta {{ text-align: right; font-size: 0.82rem; opacity: 0.85; }}

  /* Layout */
  .layout {{ display: flex; min-height: calc(100vh - 80px); }}

  .sidebar {{
    width: 220px;
    flex-shrink: 0;
    background: var(--card);
    border-right: 1px solid var(--border);
    padding: 20px 0;
    position: sticky;
    top: 0;
    height: calc(100vh - 80px);
    overflow-y: auto;
  }}
  .sidebar h3 {{
    font-size: 0.7rem;
    font-weight: 700;
    letter-spacing: 1.2px;
    text-transform: uppercase;
    color: var(--text-muted);
    padding: 0 16px 8px;
  }}
  .sidebar ul {{ list-style: none; }}
  .sidebar ul li a {{
    display: block;
    padding: 8px 16px;
    font-size: 0.85rem;
    color: var(--text);
    text-decoration: none;
    border-left: 3px solid transparent;
    transition: all 0.15s;
  }}
  .sidebar ul li a:hover {{
    background: var(--bg);
    border-left-color: var(--primary-light);
    color: var(--primary);
  }}

  .main-content {{ flex: 1; padding: 28px 32px; max-width: calc(100% - 220px); }}

  /* Status banner */
  .status-banner {{
    background: var(--card);
    border: 1px solid var(--border);
    border-left: 5px solid {status_color};
    border-radius: 8px;
    padding: 16px 20px;
    margin-bottom: 24px;
    display: flex;
    align-items: center;
    gap: 12px;
  }}
  .status-icon {{ font-size: 1.5rem; }}
  .status-text {{ font-size: 1rem; font-weight: 600; color: {status_color}; }}
  .status-sub {{ font-size: 0.82rem; color: var(--text-muted); margin-top: 2px; }}

  /* Summary table */
  .section-title {{
    font-size: 1.1rem;
    font-weight: 700;
    margin-bottom: 12px;
    color: var(--primary);
  }}
  .table-wrapper {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 8px;
    overflow: hidden;
    margin-bottom: 32px;
    overflow-x: auto;
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.88rem;
  }}
  thead tr {{ background: var(--primary); color: white; }}
  th {{ padding: 10px 14px; text-align: left; font-weight: 600; white-space: nowrap; }}
  td {{ padding: 9px 14px; border-top: 1px solid var(--border); }}
  tr.row-abnormal {{ background: #FFF5F5; }}
  tr:hover td {{ background: #F5F7FF; }}
  .mono {{ font-family: var(--mono); font-size: 0.8rem; }}

  /* Badges */
  .badge {{
    display: inline-block;
    padding: 2px 8px;
    border-radius: 12px;
    font-size: 0.75rem;
    font-weight: 700;
    white-space: nowrap;
  }}
  .badge-normal {{ background: #E8F5E9; color: #2E7D32; }}
  .badge-abnormal {{ background: #FFEBEE; color: #C62828; }}
  .badge-na {{ background: #F5F5F5; color: #757575; }}

  /* Region cards */
  .disorder-section {{ margin-bottom: 36px; }}
  .disorder-title {{
    font-size: 1rem;
    font-weight: 700;
    color: var(--primary);
    border-bottom: 2px solid var(--primary-light);
    padding-bottom: 6px;
    margin-bottom: 16px;
  }}
  .region-card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 8px;
    margin-bottom: 16px;
    overflow: hidden;
  }}
  .region-card.card-abnormal {{ border-left: 4px solid var(--accent); }}
  .card-header {{
    padding: 12px 16px;
    background: #F8F9FA;
    border-bottom: 1px solid var(--border);
    display: flex;
    align-items: center;
    gap: 10px;
    flex-wrap: wrap;
  }}
  .card-title {{ font-weight: 700; font-size: 0.95rem; }}
  .card-meta {{ font-size: 0.8rem; color: var(--text-muted); margin-left: auto; }}
  .plot-container {{ padding: 8px; }}

  /* Run metadata */
  .meta-table {{ font-size: 0.82rem; }}
  .meta-table td {{ padding: 5px 14px; }}
  .meta-table td:first-child {{ color: var(--text-muted); width: 180px; }}

  /* Responsive */
  @media (max-width: 768px) {{
    .sidebar {{ display: none; }}
    .main-content {{ max-width: 100%; padding: 16px; }}
    .page-header {{ flex-direction: column; gap: 8px; }}
  }}
</style>
</head>
<body>

<header class="page-header">
  <div>
    <div class="header-title">MethReport</div>
    <div class="header-subtitle">DNA Methylation Imprinting Disorder Analysis</div>
  </div>
  <div class="header-meta">
    <div><strong>Sample:</strong> {sample_id}</div>
    <div><strong>Genome:</strong> {genome.upper()}</div>
    <div><strong>Date:</strong> {report_date}</div>
    <div style="margin-top:4px;opacity:0.6">v{version}</div>
  </div>
</header>

<div class="layout">
  <nav class="sidebar">
    <h3>Navigation</h3>
    <ul>
      <li><a href="#summary">Summary Table</a></li>
      <li><a href="#plots">Regional Plots</a></li>
      {nav_items}
      <li><a href="#metadata">Run Info</a></li>
    </ul>
  </nav>

  <main class="main-content">

    <!-- Status Banner -->
    <div class="status-banner">
      <div class="status-icon">{'⚠️' if n_flagged > 0 else '✅'}</div>
      <div>
        <div class="status-text">{status_text}</div>
        <div class="status-sub">{n_total} DMRs analysed · Abnormal threshold: &lt;{FLAG_LOW}% or &gt;{FLAG_HIGH}%</div>
      </div>
    </div>

    <!-- Summary Table -->
    <div id="summary">
      <div class="section-title">Summary Table</div>
      <div class="table-wrapper">
        <table>
          <thead>
            <tr>
              <th>Region</th>
              <th>Disorder</th>
              <th>Coordinates</th>
              <th>CpGs</th>
              <th>Coverage</th>
              <th>Methylation</th>
              <th>HP1 Meth.</th>
              <th>HP2 Meth.</th>
              <th>Status</th>
            </tr>
          </thead>
          <tbody>
            {summary_rows_html}
          </tbody>
        </table>
      </div>
    </div>

    <!-- Regional Plots -->
    <div id="plots">
      <div class="section-title">Regional Methylation Plots</div>
      {disorder_sections}
    </div>

    <!-- Run Metadata -->
    <div id="metadata">
      <div class="section-title">Run Information</div>
      <div class="table-wrapper">
        <table class="meta-table">
          <tbody>
            <tr><td>Sample ID</td><td>{sample_id}</td></tr>
            <tr><td>Reference Genome</td><td>{genome.upper()}</td></tr>
            <tr><td>Report Generated</td><td>{report_date}</td></tr>
            <tr><td>MethReport Version</td><td>{version}</td></tr>
            {meta_rows}
          </tbody>
        </table>
      </div>
    </div>

  </main>
</div>

<script>
const plotData = {plots_json};

plotData.forEach(function(region) {{
  const el = document.getElementById('plot-' + region.name);
  if (!el) return;
  const fig = JSON.parse(region.plot_json);
  Plotly.newPlot(el, fig.data, fig.layout, {{
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    toImageButtonOptions: {{
      format: 'svg',
      filename: region.name,
      width: 1200,
      height: 400,
    }}
  }});
}});
</script>

</body>
</html>"""
