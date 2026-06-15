"""Generate a self-contained interactive HTML report using Plotly."""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from methreport import __version__
from methreport.analysis import RegionResult, SampleAnalysis
from methreport.controls import FLAG_HIGH, FLAG_LOW
from methreport.disorders import DISORDER_INFO, DisorderInfo, get_disorder_info
from methreport.reader import RegionMethylation

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
COLORS = {
    "unphased": "#2563EB",
    "hp1": "#DC2626",
    "hp2": "#0891B2",
    "control_band": "rgba(148,163,184,0.20)",
    "control_line": "rgba(100,116,139,0.55)",
    "plot_bg": "#F8FAFC",
}

FLAG_BADGE_HTML = {
    "NORMAL": '<span class="badge badge-normal">Normal</span>',
    "LOW":    '<span class="badge badge-low">Low ↓</span>',
    "HIGH":   '<span class="badge badge-high">High ↑</span>',
    "NA":     '<span class="badge badge-na">No data</span>',
}

# ---------------------------------------------------------------------------
# Plotly figure construction
# ---------------------------------------------------------------------------

def _make_region_figure(result: RegionResult) -> go.Figure:
    has_hp = result.is_phased
    cols = 2 if has_hp else 1
    titles = ["Unphased", "Phased (HP1 / HP2)"] if has_hp else ["Unphased"]

    fig = make_subplots(
        rows=1, cols=cols,
        subplot_titles=titles,
        horizontal_spacing=0.08,
    )

    _add_trace(fig, result.unphased, result.reference, 1, 1,
               COLORS["unphased"], label="Sample", show_controls=True)

    if has_hp:
        _add_trace(fig, result.hp1, result.reference, 1, 2,
                   COLORS["hp1"], label="HP1", show_controls=True)
        _add_trace(fig, result.hp2, result.reference, 1, 2,
                   COLORS["hp2"], label="HP2", show_controls=False)

    fig.update_yaxes(range=[-5, 105], ticksuffix="%", gridcolor="#E2E8F0",
                     zeroline=False, row=1, col=1)
    if has_hp:
        fig.update_yaxes(range=[-5, 105], ticksuffix="%", gridcolor="#E2E8F0",
                         zeroline=False, row=1, col=2)
    fig.update_xaxes(gridcolor="#E2E8F0", zeroline=False)

    # shaded normal band (40-60%) as layout shape
    for col in range(1, cols + 1):
        fig.add_hrect(y0=40, y1=60, line_width=0,
                      fillcolor="rgba(34,197,94,0.06)",
                      row=1, col=col)

    fig.update_layout(
        height=380,
        margin=dict(l=48, r=24, t=44, b=36),
        paper_bgcolor="white",
        plot_bgcolor=COLORS["plot_bg"],
        legend=dict(orientation="h", yanchor="bottom", y=1.04,
                    xanchor="right", x=1, font_size=11),
        font=dict(family="Inter, system-ui, sans-serif", size=11.5),
        hovermode="x unified",
    )
    return fig


def _add_trace(
    fig: go.Figure,
    rm: RegionMethylation,
    reference,
    row: int,
    col: int,
    color: str,
    label: str,
    show_controls: bool,
) -> None:
    # Control band
    if show_controls and not reference.empty:
        pos = reference["position"].tolist()
        fig.add_trace(
            go.Scatter(
                x=pos + pos[::-1],
                y=reference["upper_1sd"].tolist() + reference["lower_1sd"].tolist()[::-1],
                fill="toself",
                fillcolor=COLORS["control_band"],
                line=dict(color="rgba(0,0,0,0)"),
                name="Control ±1 SD",
                showlegend=(row == 1 and col == 1),
                hoverinfo="skip",
            ),
            row=row, col=col,
        )
        fig.add_trace(
            go.Scatter(
                x=pos,
                y=reference["mean"].tolist(),
                mode="lines",
                line=dict(color=COLORS["control_line"], width=1.2, dash="dot"),
                name="Control mean",
                showlegend=(row == 1 and col == 1),
                hoverinfo="skip",
            ),
            row=row, col=col,
        )

    if rm.n_cpg == 0:
        return

    import pandas as pd
    df = rm.to_dataframe().sort_values("position")
    positions = df["position"].tolist()
    pct = df["methylation_pct"].tolist()
    cov = df["n_total"].tolist()

    fig.add_trace(
        go.Scatter(
            x=positions, y=pct,
            mode="markers",
            marker=dict(color=color, size=5, opacity=0.55,
                        line=dict(width=0.5, color="white")),
            name=f"{label} CpGs",
            customdata=list(zip(cov, pct)),
            hovertemplate=(
                "Pos: %{x:,}<br>"
                "Methylation: %{customdata[1]:.1f}%<br>"
                "Coverage: %{customdata[0]}×<extra></extra>"
            ),
            showlegend=True,
        ),
        row=row, col=col,
    )

    if len(pct) >= 5:
        s = pd.Series(pct, index=positions)
        w = max(3, len(pct) // 5)
        smoothed = s.rolling(window=w, center=True, min_periods=1).mean()
        fig.add_trace(
            go.Scatter(
                x=smoothed.index.tolist(), y=smoothed.tolist(),
                mode="lines",
                line=dict(color=color, width=2.5),
                name=f"{label} (smoothed)",
                hoverinfo="skip",
                showlegend=True,
            ),
            row=row, col=col,
        )


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

def _meth_cell_style(pct: float | None) -> str:
    """Return inline style for a methylation percentage table cell."""
    if pct is None or np.isnan(pct):
        return ""
    if pct < FLAG_LOW:
        return 'style="color:#DC2626;font-weight:700;"'
    if pct > FLAG_HIGH:
        return 'style="color:#DC2626;font-weight:700;"'
    if 40 <= pct <= 60:
        return 'style="color:#059669;font-weight:600;"'
    return ""


def _build_stats_tiles(analysis: SampleAnalysis) -> str:
    n_flagged = len(analysis.flagged_regions)
    n_total = len(analysis.results)
    n_cpg = sum(r.unphased.n_cpg for r in analysis.results)
    is_phased = any(r.is_phased for r in analysis.results)
    phased_txt = "Yes" if is_phased else "No"
    phased_sub = "HP1/HP2 available" if is_phased else "Unphased only"
    flag_cls = "tile-danger" if n_flagged > 0 else "tile-success"
    flag_icon = "⚠" if n_flagged > 0 else "✓"

    return f"""
    <div class="stats-row">
      <div class="stat-tile tile-primary">
        <div class="tile-value">{n_total}</div>
        <div class="tile-label">Regions analysed</div>
        <div class="tile-sub">14 imprinting DMRs</div>
      </div>
      <div class="stat-tile tile-neutral">
        <div class="tile-value">{n_cpg:,}</div>
        <div class="tile-label">CpG sites</div>
        <div class="tile-sub">Above coverage threshold</div>
      </div>
      <div class="stat-tile {flag_cls}">
        <div class="tile-value">{flag_icon} {n_flagged}</div>
        <div class="tile-label">Flagged regions</div>
        <div class="tile-sub">&lt;{FLAG_LOW}% or &gt;{FLAG_HIGH}% methylation</div>
      </div>
      <div class="stat-tile tile-neutral">
        <div class="tile-value">{phased_txt}</div>
        <div class="tile-label">Phased analysis</div>
        <div class="tile-sub">{phased_sub}</div>
      </div>
    </div>"""


def _build_interpretation_section(analysis: SampleAnalysis) -> str:
    flagged = analysis.flagged_regions
    if not flagged:
        return """
    <section id="interpretation" class="section">
      <h2 class="section-heading">Clinical Interpretation</h2>
      <div class="normal-panel">
        <div class="normal-icon">✓</div>
        <div>
          <div class="normal-title">All regions show normal methylation</div>
          <div class="normal-sub">
            Methylation levels at all 14 DMRs are within the expected reference range
            (30–75%). No imprinting disorder is indicated by methylation analysis alone.
            Clinical correlation is always recommended.
          </div>
        </div>
      </div>
    </section>"""

    # Group flagged results by disorder
    by_disorder: dict[str, list[RegionResult]] = {}
    for r in flagged:
        by_disorder.setdefault(r.region.disorder, []).append(r)

    cards_html = ""
    for disorder_key, results in by_disorder.items():
        cards_html += _build_finding_card(disorder_key, results)

    n = len(flagged)
    noun = "region" if n == 1 else "regions"
    return f"""
    <section id="interpretation" class="section">
      <h2 class="section-heading">Clinical Interpretation</h2>
      <div class="abnormal-banner">
        <span class="abnormal-icon">⚠</span>
        <div>
          <div class="abnormal-title">Abnormal methylation detected at {n} {noun}</div>
          <div class="abnormal-sub">
            The findings below may indicate an imprinting disorder. This report is
            intended for use by qualified clinical professionals and does not constitute
            a standalone clinical diagnosis. Further confirmatory testing is recommended.
          </div>
        </div>
      </div>
      {cards_html}
    </section>"""


def _build_finding_card(disorder_key: str, results: list[RegionResult]) -> str:
    info: DisorderInfo | None = get_disorder_info(disorder_key)

    # OMIM badges
    omim_html = ""
    if info:
        omim_html = " ".join(e.badge for e in info.omim)

    # Affected regions row
    region_chips = ""
    for r in results:
        pct = r.mean_methylation
        pct_str = f"{pct:.1f}%" if not np.isnan(pct) else "N/A"
        arrow = "↓" if r.flag == "LOW" else "↑"
        chip_cls = "chip-low" if r.flag == "LOW" else "chip-high"
        region_chips += (
            f'<span class="region-chip {chip_cls}">'
            f'{r.region.label} &nbsp;{arrow} {pct_str}'
            f'</span>'
        )

    if info is None:
        return f"""
      <div class="finding-card">
        <div class="finding-header">
          <div class="finding-title">{disorder_key}</div>
        </div>
        <div class="finding-body">
          <div class="affected-label">Affected regions:</div>
          <div class="region-chips">{region_chips}</div>
        </div>
      </div>"""

    # Key features
    features_html = "".join(f"<li>{f}</li>" for f in info.key_features)

    # Direction-specific interpretation
    flags = {r.flag for r in results}
    interp_html = ""
    if "LOW" in flags:
        interp_html += f"""
          <div class="interp-block interp-low">
            <div class="interp-label">↓ Low methylation interpretation</div>
            <p>{info.interpretation_low}</p>
          </div>"""
    if "HIGH" in flags:
        interp_html += f"""
          <div class="interp-block interp-high">
            <div class="interp-label">↑ High methylation interpretation</div>
            <p>{info.interpretation_high}</p>
          </div>"""

    return f"""
      <div class="finding-card">
        <div class="finding-header">
          <div>
            <div class="finding-title">{info.full_name}</div>
            <div class="finding-subtitle">{info.subtitle} &nbsp;·&nbsp; {info.gene_locus}</div>
          </div>
          <div class="omim-group">{omim_html}</div>
        </div>
        <div class="finding-body">
          <div class="affected-label">Affected regions:</div>
          <div class="region-chips">{region_chips}</div>

          <p class="finding-description">{info.description}</p>

          <div class="finding-columns">
            <div class="finding-col">
              <div class="col-heading">Key clinical features</div>
              <ul class="features-list">{features_html}</ul>
            </div>
            <div class="finding-col">
              <div class="col-heading">Methylation interpretation</div>
              {interp_html}
            </div>
          </div>

          <div class="finding-footer">
            <span class="inherit-badge">Inheritance: {info.inheritance}</span>
          </div>
        </div>
      </div>"""


def _build_summary_table(analysis: SampleAnalysis) -> str:
    rows = []
    for r in analysis.results:
        meth = r.mean_methylation
        hp1 = r.mean_methylation_hp1
        hp2 = r.mean_methylation_hp2
        z   = r.z_score
        meth_str = f"{meth:.1f}%" if not np.isnan(meth) else "—"
        hp1_str  = f"{hp1:.1f}%"  if not np.isnan(hp1)  else "—"
        hp2_str  = f"{hp2:.1f}%"  if not np.isnan(hp2)  else "—"
        z_str    = f"{z:+.2f}"    if not np.isnan(z)     else "—"
        z_style  = ""
        if not np.isnan(z):
            if abs(z) >= 2.0:
                z_style = 'style="color:#DC2626;font-weight:700;"'
            elif abs(z) >= 1.5:
                z_style = 'style="color:#D97706;font-weight:600;"'
        badge    = FLAG_BADGE_HTML.get(r.flag, FLAG_BADGE_HTML["NA"])
        row_cls  = "row-flagged" if r.flag in ("LOW", "HIGH") else ""
        meth_style = _meth_cell_style(meth if not np.isnan(meth) else None)
        # informative CpG tooltip
        info_title = f'title="{r.n_informative} of {r.unphased.n_cpg} CpGs overlap control-defined informative positions"'
        rows.append(f"""
        <tr class="{row_cls}">
          <td class="td-region"><strong>{r.region.label}</strong></td>
          <td>{r.region.disorder}</td>
          <td class="td-mono">{r.region.chrom}:{r.region.start:,}–{r.region.end:,}</td>
          <td class="td-num" {info_title}>{r.unphased.n_cpg} <span class="cpg-info">({r.n_informative})</span></td>
          <td class="td-num">{r.unphased.mean_coverage:.0f}×</td>
          <td class="td-num" {meth_style}>{meth_str}</td>
          <td class="td-num">{hp1_str}</td>
          <td class="td-num">{hp2_str}</td>
          <td class="td-num" {z_style}>{z_str}</td>
          <td>{badge}</td>
        </tr>""")
    return "\n".join(rows)


def _build_disorder_sections(
    disorder_groups: dict[str, list[dict]],
    analysis: SampleAnalysis,
) -> str:
    # Build a quick lookup: disorder → flagged region names
    flagged_names = {r.region.name for r in analysis.flagged_regions}
    html = ""
    for disorder, plots in disorder_groups.items():
        anchor = disorder.replace("/", "-").replace(" ", "-")
        any_flagged = any(p["flag"] in ("LOW", "HIGH") for p in plots)
        section_cls = "disorder-section flagged-section" if any_flagged else "disorder-section"

        info = get_disorder_info(disorder)
        desc_html = ""
        if info:
            omim_badges = " ".join(e.badge for e in info.omim)
            desc_html = f"""
          <div class="disorder-desc">
            <div class="disorder-desc-top">
              <span class="disorder-full-name">{info.full_name}</span>
              <span class="disorder-omim">{omim_badges}</span>
            </div>
            <p class="disorder-desc-text">{info.description}</p>
          </div>"""

        cards_html = ""
        for p in plots:
            flag_cls = "region-card-flagged" if p["flag"] in ("LOW", "HIGH") else ""
            badge = FLAG_BADGE_HTML.get(p["flag"], FLAG_BADGE_HTML["NA"])
            phased_tag = '<span class="phased-tag">Phased</span>' if p["is_phased"] else ""
            cards_html += f"""
          <div class="region-card {flag_cls}" id="region-{p['name']}">
            <div class="region-card-header">
              <div class="region-card-left">
                <span class="region-card-title">{p['label']}</span>
                {badge}
                {phased_tag}
              </div>
              <div class="region-card-meta">
                {p['n_cpg']} CpGs &nbsp;·&nbsp; {p['mean_cov']}× coverage &nbsp;·&nbsp; {p['mean_meth']}% methylation
              </div>
            </div>
            <div class="region-plot" id="plot-{p['name']}"></div>
          </div>"""

        html += f"""
      <section class="{section_cls}" id="disorder-{anchor}">
        <div class="disorder-header-row">
          <h2 class="disorder-name">{disorder}</h2>
          {'<span class="disorder-flag-tag">⚠ Flagged</span>' if any_flagged else ''}
        </div>
        {desc_html}
        {cards_html}
      </section>"""
    return html


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def generate_report(
    analysis: SampleAnalysis,
    out_dir: Path,
    run_metadata: dict | None = None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{analysis.sample_id}_methreport.html"

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

    disorder_groups: dict[str, list[dict]] = {}
    for p in region_plots:
        disorder_groups.setdefault(p["disorder"], []).append(p)

    html = _render_html(
        analysis=analysis,
        disorder_groups=disorder_groups,
        plots_json=json.dumps(region_plots),
        report_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        run_metadata=run_metadata or {},
    )
    out_path.write_text(html, encoding="utf-8")
    log.info("HTML report → %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

def _render_html(
    analysis: SampleAnalysis,
    disorder_groups: dict,
    plots_json: str,
    report_date: str,
    run_metadata: dict,
) -> str:
    sample_id = analysis.sample_id
    genome = analysis.genome
    n_flagged = len(analysis.flagged_regions)
    n_total = len(analysis.results)

    stats_tiles    = _build_stats_tiles(analysis)
    interpretation = _build_interpretation_section(analysis)
    summary_rows   = _build_summary_table(analysis)
    disorder_secs  = _build_disorder_sections(disorder_groups, analysis)

    nav_disorders = "\n".join(
        f'<li><a class="nav-link" href="#disorder-{d.replace("/","-").replace(" ","-")}">'
        f'<span class="nav-dot {'nav-dot-flag' if any(p["flag"] in ("LOW","HIGH") for p in plots) else ''}"></span>'
        f'{d}</a></li>'
        for d, plots in disorder_groups.items()
    )

    meta_rows = "".join(
        f"<tr><td class='meta-key'>{k}</td><td>{v}</td></tr>"
        for k, v in run_metadata.items()
    )

    header_status = (
        f'<div class="header-flag-pill">⚠ {n_flagged} region{"s" if n_flagged!=1 else ""} flagged</div>'
        if n_flagged > 0 else
        '<div class="header-ok-pill">✓ All normal</div>'
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1.0"/>
<title>MethReport — {sample_id}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"/>
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin/>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet"/>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>
<style>
/* ── Reset & base ───────────────────────────────────────────────────────── */
*, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
:root {{
  --navy:       #0F172A;
  --navy-mid:   #1E3A5F;
  --blue:       #2563EB;
  --blue-light: #3B82F6;
  --danger:     #DC2626;
  --danger-bg:  #FEF2F2;
  --success:    #059669;
  --success-bg: #F0FDF4;
  --amber:      #D97706;
  --amber-bg:   #FFFBEB;
  --bg:         #F1F5F9;
  --surface:    #FFFFFF;
  --border:     #E2E8F0;
  --border-mid: #CBD5E1;
  --text:       #0F172A;
  --text-2:     #475569;
  --text-3:     #94A3B8;
  --mono:       'JetBrains Mono', 'Fira Code', monospace;
  --radius:     10px;
  --shadow-sm:  0 1px 3px rgba(15,23,42,0.08), 0 1px 2px rgba(15,23,42,0.05);
  --shadow-md:  0 4px 12px rgba(15,23,42,0.10), 0 2px 4px rgba(15,23,42,0.06);
  --shadow-lg:  0 10px 24px rgba(15,23,42,0.12), 0 4px 8px rgba(15,23,42,0.08);
}}
html {{ scroll-behavior: smooth; }}
body {{
  font-family: 'Inter', system-ui, -apple-system, sans-serif;
  background: var(--bg);
  color: var(--text);
  line-height: 1.6;
  font-size: 14px;
}}

/* ── Page header ────────────────────────────────────────────────────────── */
.page-header {{
  background: linear-gradient(135deg, var(--navy) 0%, var(--navy-mid) 55%, #1D4ED8 100%);
  color: #fff;
  padding: 0 40px;
  height: 72px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  position: sticky;
  top: 0;
  z-index: 100;
  box-shadow: 0 2px 12px rgba(15,23,42,0.35);
}}
.header-brand {{ display: flex; align-items: center; gap: 12px; }}
.header-logo {{
  width: 36px; height: 36px;
  background: rgba(255,255,255,0.15);
  border-radius: 8px;
  display: flex; align-items: center; justify-content: center;
  font-size: 18px;
}}
.header-title-wrap {{ line-height: 1.2; }}
.header-title {{ font-size: 1.15rem; font-weight: 700; letter-spacing: -0.3px; }}
.header-subtitle {{ font-size: 0.72rem; opacity: 0.65; letter-spacing: 0.3px; text-transform: uppercase; }}
.header-right {{ display: flex; align-items: center; gap: 16px; }}
.header-info {{ text-align: right; font-size: 0.78rem; opacity: 0.8; line-height: 1.5; }}
.header-info strong {{ opacity: 1; font-weight: 600; }}
.header-flag-pill {{
  background: rgba(220,38,38,0.85);
  border: 1px solid rgba(255,100,100,0.4);
  color: #fff;
  padding: 4px 12px;
  border-radius: 20px;
  font-size: 0.78rem;
  font-weight: 600;
  white-space: nowrap;
}}
.header-ok-pill {{
  background: rgba(5,150,105,0.75);
  border: 1px solid rgba(52,211,153,0.4);
  color: #fff;
  padding: 4px 12px;
  border-radius: 20px;
  font-size: 0.78rem;
  font-weight: 600;
  white-space: nowrap;
}}

/* ── Layout ─────────────────────────────────────────────────────────────── */
.layout {{ display: flex; min-height: calc(100vh - 72px); }}

.sidebar {{
  width: 232px;
  flex-shrink: 0;
  background: var(--surface);
  border-right: 1px solid var(--border);
  padding: 24px 0 40px;
  position: sticky;
  top: 72px;
  height: calc(100vh - 72px);
  overflow-y: auto;
  scrollbar-width: thin;
}}
.sidebar-section {{ margin-bottom: 8px; }}
.sidebar-label {{
  font-size: 0.68rem;
  font-weight: 700;
  letter-spacing: 1px;
  text-transform: uppercase;
  color: var(--text-3);
  padding: 12px 20px 6px;
}}
.nav-link {{
  display: flex;
  align-items: center;
  gap: 8px;
  padding: 7px 20px;
  font-size: 0.83rem;
  color: var(--text-2);
  text-decoration: none;
  border-left: 3px solid transparent;
  transition: all 0.15s ease;
}}
.nav-link:hover {{
  background: var(--bg);
  border-left-color: var(--blue-light);
  color: var(--blue);
}}
.nav-dot {{
  width: 7px; height: 7px;
  border-radius: 50%;
  background: var(--border-mid);
  flex-shrink: 0;
}}
.nav-dot-flag {{ background: var(--danger); }}
.sidebar-divider {{
  height: 1px;
  background: var(--border);
  margin: 8px 20px;
}}

.main-content {{
  flex: 1;
  padding: 28px 36px 56px;
  min-width: 0;
  max-width: calc(100% - 232px);
}}

/* ── Sections ───────────────────────────────────────────────────────────── */
.section {{ margin-bottom: 36px; }}
.section-heading {{
  font-size: 0.95rem;
  font-weight: 700;
  color: var(--navy);
  text-transform: uppercase;
  letter-spacing: 0.6px;
  margin-bottom: 14px;
  padding-bottom: 8px;
  border-bottom: 2px solid var(--border);
  display: flex;
  align-items: center;
  gap: 8px;
}}
.section-heading::before {{
  content: '';
  display: inline-block;
  width: 3px;
  height: 16px;
  background: var(--blue);
  border-radius: 2px;
}}

/* ── Stats tiles ─────────────────────────────────────────────────────────── */
.stats-row {{
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: 14px;
  margin-bottom: 28px;
}}
.stat-tile {{
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 16px 20px;
  box-shadow: var(--shadow-sm);
  border-top: 3px solid transparent;
}}
.tile-primary  {{ border-top-color: var(--blue); }}
.tile-success  {{ border-top-color: var(--success); }}
.tile-danger   {{ border-top-color: var(--danger); }}
.tile-neutral  {{ border-top-color: var(--border-mid); }}
.tile-value  {{ font-size: 1.65rem; font-weight: 700; line-height: 1.1; margin-bottom: 3px; }}
.tile-primary  .tile-value {{ color: var(--blue); }}
.tile-success  .tile-value {{ color: var(--success); }}
.tile-danger   .tile-value {{ color: var(--danger); }}
.tile-neutral  .tile-value {{ color: var(--navy); }}
.tile-label {{ font-size: 0.82rem; font-weight: 600; color: var(--text); }}
.tile-sub   {{ font-size: 0.72rem; color: var(--text-3); margin-top: 2px; }}

/* ── Interpretation section ──────────────────────────────────────────────── */
.normal-panel {{
  background: var(--success-bg);
  border: 1px solid #86EFAC;
  border-radius: var(--radius);
  padding: 18px 22px;
  display: flex;
  align-items: flex-start;
  gap: 16px;
}}
.normal-icon {{
  font-size: 1.6rem;
  color: var(--success);
  flex-shrink: 0;
  line-height: 1;
  margin-top: 2px;
}}
.normal-title {{ font-weight: 700; font-size: 1rem; color: #065F46; margin-bottom: 4px; }}
.normal-sub   {{ font-size: 0.84rem; color: #047857; }}

.abnormal-banner {{
  background: var(--danger-bg);
  border: 1px solid #FECACA;
  border-radius: var(--radius);
  padding: 16px 20px;
  display: flex;
  align-items: flex-start;
  gap: 14px;
  margin-bottom: 18px;
}}
.abnormal-icon  {{ font-size: 1.4rem; color: var(--danger); flex-shrink: 0; margin-top: 2px; }}
.abnormal-title {{ font-weight: 700; font-size: 0.95rem; color: #991B1B; margin-bottom: 3px; }}
.abnormal-sub   {{ font-size: 0.81rem; color: #B91C1C; line-height: 1.5; }}

/* ── Finding cards ──────────────────────────────────────────────────────── */
.finding-card {{
  background: var(--surface);
  border: 1px solid var(--border);
  border-left: 4px solid var(--danger);
  border-radius: var(--radius);
  margin-bottom: 16px;
  box-shadow: var(--shadow-sm);
  overflow: hidden;
}}
.finding-header {{
  background: #FFFAFA;
  border-bottom: 1px solid var(--border);
  padding: 14px 20px;
  display: flex;
  align-items: flex-start;
  justify-content: space-between;
  gap: 12px;
  flex-wrap: wrap;
}}
.finding-title    {{ font-size: 1.05rem; font-weight: 700; color: var(--navy); line-height: 1.3; }}
.finding-subtitle {{ font-size: 0.78rem; color: var(--text-2); margin-top: 3px; font-family: var(--mono); }}
.omim-group {{ display: flex; flex-wrap: wrap; gap: 6px; align-items: center; }}
.omim-badge {{
  display: inline-block;
  padding: 3px 10px;
  border-radius: 14px;
  font-size: 0.72rem;
  font-weight: 600;
  background: #EFF6FF;
  color: #1D4ED8;
  border: 1px solid #BFDBFE;
  text-decoration: none;
  transition: background 0.15s;
}}
.omim-badge:hover {{ background: #DBEAFE; }}

.finding-body {{ padding: 18px 20px; }}
.affected-label {{ font-size: 0.72rem; font-weight: 700; text-transform: uppercase;
                   letter-spacing: 0.5px; color: var(--text-3); margin-bottom: 8px; }}
.region-chips {{ display: flex; flex-wrap: wrap; gap: 7px; margin-bottom: 16px; }}
.region-chip {{
  padding: 4px 12px;
  border-radius: 16px;
  font-size: 0.78rem;
  font-weight: 600;
  font-family: var(--mono);
}}
.chip-low  {{ background: #FEF2F2; color: var(--danger); border: 1px solid #FECACA; }}
.chip-high {{ background: #FFF7ED; color: var(--amber);  border: 1px solid #FED7AA; }}

.finding-description {{
  font-size: 0.85rem;
  color: var(--text-2);
  line-height: 1.65;
  margin-bottom: 16px;
  padding: 12px 14px;
  background: #F8FAFC;
  border-left: 3px solid var(--border-mid);
  border-radius: 0 6px 6px 0;
}}

.finding-columns {{
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 20px;
  margin-bottom: 14px;
}}
@media (max-width: 860px) {{ .finding-columns {{ grid-template-columns: 1fr; }} }}
.finding-col {{ }}
.col-heading {{
  font-size: 0.72rem;
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: 0.5px;
  color: var(--text-3);
  margin-bottom: 8px;
}}
.features-list {{
  list-style: none;
  display: flex;
  flex-direction: column;
  gap: 5px;
}}
.features-list li {{
  font-size: 0.83rem;
  color: var(--text-2);
  padding-left: 14px;
  position: relative;
}}
.features-list li::before {{
  content: '·';
  position: absolute;
  left: 4px;
  color: var(--blue);
  font-weight: 700;
}}
.interp-block {{
  font-size: 0.82rem;
  line-height: 1.6;
  color: var(--text);
  padding: 10px 14px;
  border-radius: 6px;
  margin-bottom: 8px;
}}
.interp-low  {{ background: #FEF2F2; border-left: 3px solid var(--danger); }}
.interp-high {{ background: #FFF7ED; border-left: 3px solid var(--amber);  }}
.interp-label {{ font-weight: 700; margin-bottom: 4px; font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.4px; }}
.interp-low  .interp-label {{ color: var(--danger); }}
.interp-high .interp-label {{ color: var(--amber);  }}

.finding-footer {{
  padding-top: 10px;
  border-top: 1px solid var(--border);
  margin-top: 2px;
}}
.inherit-badge {{
  display: inline-block;
  font-size: 0.72rem;
  color: var(--text-3);
  background: var(--bg);
  border: 1px solid var(--border);
  padding: 3px 10px;
  border-radius: 12px;
}}

/* ── Summary table ───────────────────────────────────────────────────────── */
.table-card {{
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  overflow: hidden;
  box-shadow: var(--shadow-sm);
  margin-bottom: 28px;
  overflow-x: auto;
}}
table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 0.84rem;
}}
thead tr {{
  background: var(--navy);
  color: #fff;
}}
th {{
  padding: 11px 14px;
  text-align: left;
  font-weight: 600;
  font-size: 0.78rem;
  letter-spacing: 0.3px;
  white-space: nowrap;
}}
td {{ padding: 9px 14px; border-top: 1px solid var(--border); }}
tr.row-flagged {{ background: #FFF5F5; }}
tbody tr:hover td {{ background: #F8FAFF; }}
.td-region {{ font-weight: 600; }}
.td-num    {{ text-align: right; font-variant-numeric: tabular-nums; }}
.td-mono   {{ font-family: var(--mono); font-size: 0.76rem; color: var(--text-2); }}

/* ── Badges ──────────────────────────────────────────────────────────────── */
.badge {{
  display: inline-flex;
  align-items: center;
  padding: 2px 9px;
  border-radius: 12px;
  font-size: 0.72rem;
  font-weight: 700;
  white-space: nowrap;
  letter-spacing: 0.2px;
}}
.badge-normal {{ background: #F0FDF4; color: #065F46; border: 1px solid #86EFAC; }}
.badge-low    {{ background: #FEF2F2; color: #991B1B; border: 1px solid #FECACA; }}
.badge-high   {{ background: #FFF7ED; color: #92400E; border: 1px solid #FED7AA; }}
.badge-na     {{ background: #F8FAFC; color: #64748B; border: 1px solid #CBD5E1; }}

/* ── Disorder sections ───────────────────────────────────────────────────── */
.disorder-section {{ margin-bottom: 40px; }}
.disorder-header-row {{
  display: flex;
  align-items: center;
  gap: 10px;
  margin-bottom: 10px;
}}
.disorder-name {{
  font-size: 1rem;
  font-weight: 700;
  color: var(--navy);
}}
.disorder-flag-tag {{
  font-size: 0.72rem;
  font-weight: 700;
  background: var(--danger-bg);
  color: var(--danger);
  border: 1px solid #FECACA;
  padding: 2px 9px;
  border-radius: 10px;
}}
.disorder-desc {{
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 14px 18px;
  margin-bottom: 14px;
  box-shadow: var(--shadow-sm);
}}
.disorder-desc-top {{
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: 10px;
  margin-bottom: 6px;
}}
.disorder-full-name {{ font-size: 0.85rem; font-weight: 700; color: var(--navy); }}
.disorder-omim {{ display: flex; flex-wrap: wrap; gap: 5px; }}
.disorder-desc-text {{ font-size: 0.81rem; color: var(--text-2); line-height: 1.65; }}

/* ── Region cards ────────────────────────────────────────────────────────── */
.region-card {{
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  margin-bottom: 14px;
  overflow: hidden;
  box-shadow: var(--shadow-sm);
  transition: box-shadow 0.2s;
}}
.region-card:hover {{ box-shadow: var(--shadow-md); }}
.region-card-flagged {{ border-left: 3px solid var(--danger); }}
.region-card-header {{
  padding: 11px 16px;
  background: #FAFBFC;
  border-bottom: 1px solid var(--border);
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 10px;
  flex-wrap: wrap;
}}
.region-card-left {{ display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }}
.region-card-title {{ font-weight: 700; font-size: 0.88rem; font-family: var(--mono); }}
.region-card-meta  {{ font-size: 0.77rem; color: var(--text-3); white-space: nowrap; }}
.phased-tag {{
  font-size: 0.68rem;
  font-weight: 600;
  background: #EFF6FF;
  color: #1D4ED8;
  border: 1px solid #BFDBFE;
  padding: 1px 7px;
  border-radius: 8px;
}}
.region-plot {{ padding: 6px 8px; }}

/* ── Metadata table ──────────────────────────────────────────────────────── */
.meta-table {{ font-size: 0.82rem; }}
.meta-table td {{ padding: 8px 16px; }}
.meta-key {{ color: var(--text-2); width: 200px; font-weight: 500; }}
.cpg-info {{ color: var(--text-3); font-size: 0.78rem; }}

/* ── Footer ──────────────────────────────────────────────────────────────── */
.report-footer {{
  text-align: center;
  font-size: 0.72rem;
  color: var(--text-3);
  padding: 24px 0 8px;
  border-top: 1px solid var(--border);
  margin-top: 40px;
}}
.report-footer a {{ color: var(--blue-light); text-decoration: none; }}

/* ── Responsive ──────────────────────────────────────────────────────────── */
@media (max-width: 900px) {{
  .sidebar {{ display: none; }}
  .main-content {{ max-width: 100%; padding: 20px 18px; }}
  .stats-row {{ grid-template-columns: repeat(2, 1fr); }}
  .page-header {{ padding: 0 18px; }}
}}
@media (max-width: 480px) {{
  .stats-row {{ grid-template-columns: 1fr; }}
  .finding-columns {{ grid-template-columns: 1fr; }}
}}
</style>
</head>
<body>

<!-- ─── Header ──────────────────────────────────────────────────────────── -->
<header class="page-header">
  <div class="header-brand">
    <div class="header-logo">🧬</div>
    <div class="header-title-wrap">
      <div class="header-title">MethReport</div>
      <div class="header-subtitle">Methylation Imprinting Analysis</div>
    </div>
  </div>
  <div class="header-right">
    <div class="header-info">
      <div><strong>Sample:</strong> {sample_id}</div>
      <div><strong>Genome:</strong> {genome.upper()} &nbsp;·&nbsp; <strong>Date:</strong> {report_date} &nbsp;·&nbsp; v{__version__}</div>
    </div>
    {header_status}
  </div>
</header>

<!-- ─── Layout ──────────────────────────────────────────────────────────── -->
<div class="layout">

  <nav class="sidebar">
    <div class="sidebar-section">
      <div class="sidebar-label">Overview</div>
      <a class="nav-link" href="#interpretation">
        <span class="nav-dot {'nav-dot-flag' if n_flagged > 0 else ''}"></span>Interpretation
      </a>
      <a class="nav-link" href="#summary">
        <span class="nav-dot"></span>Summary Table
      </a>
    </div>
    <div class="sidebar-divider"></div>
    <div class="sidebar-section">
      <div class="sidebar-label">Disorders</div>
      <ul style="list-style:none">
        {nav_disorders}
      </ul>
    </div>
    <div class="sidebar-divider"></div>
    <div class="sidebar-section">
      <div class="sidebar-label">Report</div>
      <a class="nav-link" href="#runinfo">
        <span class="nav-dot"></span>Run Information
      </a>
    </div>
  </nav>

  <main class="main-content">

    <!-- Stats tiles -->
    {stats_tiles}

    <!-- Clinical interpretation -->
    {interpretation}

    <!-- Summary table -->
    <section id="summary" class="section">
      <h2 class="section-heading">Summary Table</h2>
      <div class="table-card">
        <table>
          <thead>
            <tr>
              <th>Region</th>
              <th>Disorder</th>
              <th>Coordinates</th>
              <th style="text-align:right" title="Total CpGs (informative CpGs used for z-score)">CpGs (info.)</th>
              <th style="text-align:right">Coverage</th>
              <th style="text-align:right">Methylation</th>
              <th style="text-align:right">HP1</th>
              <th style="text-align:right">HP2</th>
              <th style="text-align:right" title="Mean z-score vs per-position control distribution. |z|≥2 triggers flag.">z-score</th>
              <th>Status</th>
            </tr>
          </thead>
          <tbody>
            {summary_rows}
          </tbody>
        </table>
      </div>
    </section>

    <!-- Regional methylation plots -->
    <section id="plots" class="section">
      <h2 class="section-heading">Regional Methylation Plots</h2>
      {disorder_secs}
    </section>

    <!-- Run information -->
    <section id="runinfo" class="section">
      <h2 class="section-heading">Run Information</h2>
      <div class="table-card">
        <table class="meta-table">
          <tbody>
            <tr><td class="meta-key">Sample ID</td><td>{sample_id}</td></tr>
            <tr><td class="meta-key">Reference Genome</td><td>{genome.upper()}</td></tr>
            <tr><td class="meta-key">Report Generated</td><td>{report_date}</td></tr>
            <tr><td class="meta-key">MethReport Version</td><td>{__version__}</td></tr>
            {meta_rows}
          </tbody>
        </table>
      </div>
    </section>

    <footer class="report-footer">
      Generated by <a href="https://github.com/j-jamshidi/MethReport" target="_blank">MethReport</a> v{__version__} ·
      Region coordinates derived from <a href="https://github.com/carolinehey/NanoImprint" target="_blank">NanoImprint</a> ·
      OMIM links are provided for reference only and do not constitute clinical advice.
    </footer>

  </main>
</div>

<!-- ─── Plotly rendering ─────────────────────────────────────────────────── -->
<script>
(function() {{
  const plots = {plots_json};
  plots.forEach(function(r) {{
    const el = document.getElementById('plot-' + r.name);
    if (!el) return;
    const fig = JSON.parse(r.plot_json);
    Plotly.newPlot(el, fig.data, fig.layout, {{
      responsive: true,
      displayModeBar: true,
      modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'],
      toImageButtonOptions: {{
        format: 'svg',
        filename: r.name,
        width: 1200,
        height: 420,
      }},
      displaylogo: false,
    }});
  }});
}})();
</script>

</body>
</html>"""
