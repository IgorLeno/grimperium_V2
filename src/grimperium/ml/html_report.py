"""Self-contained HTML error analysis report generator.

Receives a pre-computed AnalysisResult — does not read CSV, does not
recalculate any statistics. All data comes from analyze().
"""

from __future__ import annotations

import html as html_lib
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from grimperium.ml.error_analysis import AnalysisResult

_SEVERITY_COLORS: dict[str, str] = {
    "normal": "#2ecc71",
    "elevated": "#f39c12",
    "high": "#e67e22",
    "critical": "#c0392b",
    "extreme": "#8e44ad",
}

_CSS = """
body{font-family:Arial,sans-serif;background:#fff;color:#222;margin:2em;max-width:1400px}
h1{color:#2c3e50}
h2{color:#34495e;border-bottom:2px solid #bdc3c7;padding-bottom:.3em;margin-top:2em}
table{border-collapse:collapse;width:100%;margin-bottom:1.5em;font-size:.9em}
th{background:#34495e;color:#fff;padding:.5em .8em;text-align:left}
td{border:1px solid #ddd;padding:.35em .7em}
tr:nth-child(even){background:#f9f9f9}
.crit td{color:#c0392b;font-weight:bold}
.header-box{background:#ecf0f1;padding:1em 1.5em;border-radius:6px;margin-bottom:1.5em;display:flex;flex-wrap:wrap;gap:2em}
.metric{min-width:8em}
.metric-label{font-size:.8em;color:#7f8c8d;text-transform:uppercase;letter-spacing:.05em}
.metric-value{font-size:1.4em;font-weight:bold;color:#2c3e50}
.warn-box{background:#fef9e7;border-left:4px solid #f39c12;padding:.8em 1.2em;border-radius:0 4px 4px 0}
.warn-box ul{margin:.3em 0;padding-left:1.2em}
.no-data{color:#7f8c8d;font-style:italic}
.chart-links{display:flex;gap:1.5em;flex-wrap:wrap}
.chart-links a{background:#2980b9;color:#fff;padding:.4em .9em;border-radius:4px;text-decoration:none;font-size:.9em}
.chart-links a:hover{background:#1a6fa0}
"""


def generate_html_report(result: AnalysisResult, output_path: Path) -> Path:
    """Write a self-contained HTML report for the given AnalysisResult.

    Args:
        result: Pre-computed AnalysisResult from analyze(). No CSV is read here.
        output_path: Destination file path for the HTML report.

    Returns:
        Resolved path where the HTML was written.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    s = result.summary

    body = "\n".join(
        [
            _build_header(s, now),
            _build_severity_table(s),
            _build_top_errors_table(result.top_n_df),
            _build_outliers_table(result.outliers_df),
            _build_warnings_section(result.warnings),
            _build_chart_links(output_path),
        ]
    )

    content = (
        "<!DOCTYPE html>\n"
        '<html lang="en">\n'
        "<head>\n"
        '<meta charset="UTF-8">\n'
        '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
        "<title>Grimperium Error Analysis Report</title>\n"
        f"<style>{_CSS}</style>\n"
        "</head>\n"
        f"<body>\n{body}\n</body>\n"
        "</html>"
    )

    output_path.write_text(content, encoding="utf-8")
    return output_path


# ── Section builders ──────────────────────────────────────────────────────────


def _metric_block(label: str, value: str) -> str:
    return (
        f'<div class="metric">'
        f'<div class="metric-label">{html_lib.escape(label)}</div>'
        f'<div class="metric-value">{html_lib.escape(value)}</div>'
        f"</div>"
    )


def _build_header(s: dict[str, Any], now: str) -> str:
    mae = f"{float(s.get('mae', 0.0)):.4f}"
    rmse = f"{float(s.get('rmse', 0.0)):.4f}"
    r2 = f"{float(s.get('r2', 0.0)):.4f}"
    n = int(s.get("n_molecules", 0))
    n_out = int(s.get("n_outliers", 0))

    blocks = "\n".join(
        [
            _metric_block("Generated", now),
            _metric_block("Molecules", f"{n:,}"),
            _metric_block("MAE (kcal/mol)", mae),
            _metric_block("RMSE (kcal/mol)", rmse),
            _metric_block("R²", r2),
            _metric_block("Outliers", str(n_out)),
        ]
    )
    return (
        "<h1>Grimperium Error Analysis Report</h1>\n"
        f'<div class="header-box">\n{blocks}\n</div>'
    )


def _build_severity_table(s: dict[str, Any]) -> str:
    rows = ""
    for sev in ("normal", "elevated", "high", "critical", "extreme"):
        count = int(s.get(f"n_{sev}", 0))
        color = _SEVERITY_COLORS.get(sev, "#333")
        rows += (
            f"<tr>"
            f'<td style="color:{color};font-weight:bold">{html_lib.escape(sev)}</td>'
            f"<td>{count:,}</td>"
            f"</tr>\n"
        )
    return (
        "<h2>Severity Distribution</h2>\n"
        "<table><tr><th>Severity</th><th>Count</th></tr>\n"
        f"{rows}"
        "</table>"
    )


def _build_top_errors_table(top_n_df: pd.DataFrame) -> str:
    if top_n_df.empty:
        return '<h2>Top 20 Largest Errors</h2>\n<p class="no-data">No data.</p>'

    display_cols = [
        c
        for c in ["rank", "mol_id", "abs_error", "signed_error", "severity"]
        if c in top_n_df.columns
    ]
    display = top_n_df.head(20)
    header = "".join(f"<th>{html_lib.escape(c)}</th>" for c in display_cols)

    rows = ""
    for _, row in display.iterrows():
        sev = str(row.get("severity", ""))
        css = "crit" if sev in ("critical", "extreme") else ""
        cells = _row_cells(row, display_cols)
        rows += f'<tr class="{css}">{cells}</tr>\n'

    return (
        "<h2>Top 20 Largest Errors</h2>\n"
        f"<table><tr>{header}</tr>\n{rows}</table>"
    )


def _build_outliers_table(outliers_df: pd.DataFrame) -> str:
    if outliers_df.empty:
        return (
            "<h2>All Outliers</h2>\n"
            '<p class="no-data">No outliers detected.</p>'
        )

    priority = [
        "mol_id",
        "smiles",
        "H298_cbs",
        "H298_predicted",
        "abs_error",
        "signed_error",
        "severity",
        "outlier_score",
        "iqr_outlier",
        "extreme_iqr_outlier",
        "threshold_outlier",
        "abs_error_zscore",
    ]
    ordered = [c for c in priority if c in outliers_df.columns]
    rest = [c for c in outliers_df.columns if c not in ordered]
    all_cols = ordered + rest

    header = "".join(f"<th>{html_lib.escape(c)}</th>" for c in all_cols)
    rows = ""
    for _, row in outliers_df.iterrows():
        sev = str(row.get("severity", ""))
        css = "crit" if sev in ("critical", "extreme") else ""
        cells = _row_cells(row, all_cols)
        rows += f'<tr class="{css}">{cells}</tr>\n'

    return (
        f"<h2>All Outliers ({len(outliers_df):,})</h2>\n"
        f"<table><tr>{header}</tr>\n{rows}</table>"
    )


def _build_warnings_section(warnings: list[str]) -> str:
    if not warnings:
        return "<h2>Warnings</h2>\n" '<p class="no-data">No warnings.</p>'
    items = "".join(f"<li>{html_lib.escape(w)}</li>\n" for w in warnings)
    return (
        "<h2>Warnings</h2>\n"
        f'<div class="warn-box"><ul>\n{items}</ul></div>'
    )


def _build_chart_links(output_path: Path) -> str:
    chart_dir = output_path.parent
    charts = [
        ("Parity Plot", "parity_plot.png"),
        ("Delta Histogram", "delta_histogram.png"),
        ("Residuals Plot", "residuals_plot.png"),
    ]
    links = "".join(
        f'<a href="{html_lib.escape(fname)}">{html_lib.escape(label)}</a>'
        for label, fname in charts
        if (chart_dir / fname).exists()
    )
    if not links:
        return ""
    return "<h2>Charts</h2>\n" f'<div class="chart-links">{links}</div>'


def _row_cells(row: pd.Series[Any], cols: list[str]) -> str:
    cells = ""
    for c in cols:
        if c not in row.index:
            cells += "<td>-</td>"
            continue
        val = row[c]
        if isinstance(val, bool):
            cells += f"<td>{'Yes' if val else 'No'}</td>"
        elif isinstance(val, float):
            cells += f"<td>{val:.4f}</td>"
        elif isinstance(val, int | np.integer):
            cells += f"<td>{int(val)}</td>"
        else:
            cells += f"<td>{html_lib.escape(str(val))}</td>"
    return cells
