"""Tests for src/grimperium/ml/html_report.py."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from grimperium.ml.error_analysis import AnalysisResult, ErrorAnalysisConfig, analyze
from grimperium.ml.html_report import generate_html_report

# ── Helpers ───────────────────────────────────────────────────────────────────


def _minimal_result() -> AnalysisResult:
    """Build a real AnalysisResult from a tiny synthetic DataFrame."""
    df = pd.DataFrame(
        {
            "mol_id": ["m1", "m2", "m3", "m4", "m5"],
            "H298_cbs": [-100.0, -200.0, -300.0, -400.0, -500.0],
            "H298_predicted": [-101.0, -198.0, -303.0, -399.0, -550.0],
        }
    )
    return analyze(df, ErrorAnalysisConfig(threshold_kcalmol=5.0))


def _empty_outliers_result() -> AnalysisResult:
    """AnalysisResult with perfect predictions — outliers_df is guaranteed empty."""
    df = pd.DataFrame(
        {
            "H298_cbs": [-100.0, -200.0, -300.0],
            "H298_predicted": [-100.0, -200.0, -300.0],
        }
    )
    return analyze(df)


# ── Tests ─────────────────────────────────────────────────────────────────────


def test_html_file_created(tmp_path: Path) -> None:
    """generate_html_report creates the output file."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    assert out.exists()


def test_html_returns_path(tmp_path: Path) -> None:
    """generate_html_report returns the path where the file was saved."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    returned = generate_html_report(result, out)
    assert returned == out


def test_html_contains_title_tag(tmp_path: Path) -> None:
    """HTML output contains a <title> tag."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    content = out.read_text()
    assert "<title>" in content


def test_html_contains_h1(tmp_path: Path) -> None:
    """HTML output contains an <h1> heading."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    content = out.read_text()
    assert "<h1>" in content


def test_html_contains_outliers_section(tmp_path: Path) -> None:
    """HTML output contains the All Outliers section heading."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    content = out.read_text()
    assert "All Outliers" in content or "Outliers" in content


def test_html_contains_warnings_section(tmp_path: Path) -> None:
    """HTML output contains the Warnings section heading."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    content = out.read_text()
    assert "Warnings" in content


def test_html_contains_severity_section(tmp_path: Path) -> None:
    """HTML output contains the Severity Distribution section."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    content = out.read_text()
    assert "Severity" in content


def test_no_chart_files_no_chart_links(tmp_path: Path) -> None:
    """When chart PNGs are absent, no chart links appear in the HTML."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    content = out.read_text()
    # None of the chart filenames should appear as href links
    for fname in ("parity_plot.png", "delta_histogram.png", "residuals_plot.png"):
        assert f'href="{fname}"' not in content


def test_chart_files_present_creates_links(tmp_path: Path) -> None:
    """When chart PNGs exist in the same directory, links appear in the HTML."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    # Create fake chart files
    (tmp_path / "parity_plot.png").write_bytes(b"")
    (tmp_path / "delta_histogram.png").write_bytes(b"")
    generate_html_report(result, out)
    content = out.read_text()
    assert 'href="parity_plot.png"' in content
    assert 'href="delta_histogram.png"' in content


def test_empty_outliers_renders_fallback(tmp_path: Path) -> None:
    """When outliers_df is empty, the HTML shows a 'no outliers' message."""
    result = _empty_outliers_result()
    assert result.outliers_df.empty
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    content = out.read_text()
    assert "No outliers" in content or "no-data" in content


def test_output_directory_created_if_missing(tmp_path: Path) -> None:
    """generate_html_report creates parent directories as needed."""
    result = _minimal_result()
    nested = tmp_path / "a" / "b" / "report.html"
    generate_html_report(result, nested)
    assert nested.exists()


def test_html_is_utf8_encoded(tmp_path: Path) -> None:
    """Output file is valid UTF-8."""
    result = _minimal_result()
    out = tmp_path / "report.html"
    generate_html_report(result, out)
    content = out.read_bytes().decode("utf-8")  # raises on invalid UTF-8
    assert len(content) > 0
