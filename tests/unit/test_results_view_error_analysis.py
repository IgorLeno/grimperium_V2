"""Tests for ResultsView error-analysis menu actions.

Verifies that:
  - _handle_detailed_metrics delegates to analyze() with zero inline statistics
  - _handle_show_outliers renders an outlier table
  - _handle_export_outliers saves CSV to the expected path
  - _handle_html_report calls generate_html_report with the correct path
"""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from rich.console import Console

from grimperium.cli.views.results_view import ResultsView

# ── Setup helper ──────────────────────────────────────────────────────────────


def _create_view() -> tuple[ResultsView, io.StringIO]:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    controller = MagicMock()
    controller.console = console
    view = ResultsView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    return view, buf


def _make_csv(tmp_path: Path) -> Path:
    csv_path = tmp_path / "data.csv"
    df = pd.DataFrame(
        {
            "mol_id": ["m1", "m2", "m3"],
            "H298_cbs": [-100.0, -200.0, -300.0],
            "H298_predicted": [-101.0, -198.0, -303.0],
        }
    )
    df.to_csv(csv_path, index=False)
    return csv_path


# ── New menu options exist ────────────────────────────────────────────────────


def test_show_outliers_menu_option_exists() -> None:
    view, _ = _create_view()
    values = [o.value for o in view.get_menu_options()]
    assert "show_outliers" in values


def test_top_errors_menu_option_exists() -> None:
    view, _ = _create_view()
    values = [o.value for o in view.get_menu_options()]
    assert "top_errors" in values


def test_html_report_menu_option_exists() -> None:
    view, _ = _create_view()
    values = [o.value for o in view.get_menu_options()]
    assert "html_report" in values


def test_export_outliers_menu_option_exists() -> None:
    view, _ = _create_view()
    values = [o.value for o in view.get_menu_options()]
    assert "export_outliers" in values


# ── _handle_detailed_metrics calls analyze() ──────────────────────────────────


def test_handle_detailed_metrics_calls_analyze(tmp_path: Path) -> None:
    """_handle_detailed_metrics must call analyze() and read from its summary."""
    view, buf = _create_view()
    csv_path = _make_csv(tmp_path)

    fake_summary = {
        "mae": 2.0,
        "rmse": 2.5,
        "max_error": 5.0,
        "bias": -0.5,
        "r2": 0.9990,
        "pearson_r": 0.9995,
        "pct_within_1": 33.3,
        "pct_within_2": 66.7,
        "pct_within_5": 100.0,
        "n_molecules": 3,
    }
    fake_result = MagicMock()
    fake_result.summary = fake_summary

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch(
            "grimperium.ml.error_analysis.analyze", return_value=fake_result
        ) as mock_analyze,
    ):
        view._handle_detailed_metrics()

    mock_analyze.assert_called_once()
    output = buf.getvalue()
    # Values from the mock summary appear in the output
    assert "2.0000" in output  # mae
    assert "2.5000" in output  # rmse


def test_handle_detailed_metrics_no_csv_shows_error(tmp_path: Path) -> None:
    view, buf = _create_view()
    missing = tmp_path / "nope.csv"
    with patch.object(view, "_get_csv_path", return_value=missing):
        view._handle_detailed_metrics()
    assert "Error" in buf.getvalue() or "error" in buf.getvalue().lower()


# ── _handle_show_outliers ─────────────────────────────────────────────────────


def test_handle_show_outliers_renders_table(tmp_path: Path) -> None:
    """_handle_show_outliers renders an outlier table from analyze() output."""
    view, buf = _create_view()
    csv_path = _make_csv(tmp_path)

    outliers_df = pd.DataFrame(
        {
            "mol_id": ["m3"],
            "H298_cbs": [-300.0],
            "H298_predicted": [-303.0],
            "abs_error": [3.0],
            "severity": ["normal"],
            "outlier_score": [1],
        }
    )
    fake_result = MagicMock()
    fake_result.outliers_df = outliers_df

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch("grimperium.ml.error_analysis.analyze", return_value=fake_result),
    ):
        view._handle_show_outliers()

    output = buf.getvalue()
    assert "m3" in output


def test_handle_show_outliers_no_csv_shows_error(tmp_path: Path) -> None:
    view, buf = _create_view()
    missing = tmp_path / "nope.csv"
    with patch.object(view, "_get_csv_path", return_value=missing):
        view._handle_show_outliers()
    assert "Error" in buf.getvalue() or "error" in buf.getvalue().lower()


def test_handle_show_outliers_empty_shows_message(tmp_path: Path) -> None:
    """Empty outliers_df shows a 'no outliers' message."""
    view, buf = _create_view()
    csv_path = _make_csv(tmp_path)

    fake_result = MagicMock()
    fake_result.outliers_df = pd.DataFrame()

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch("grimperium.ml.error_analysis.analyze", return_value=fake_result),
    ):
        view._handle_show_outliers()

    assert "No outliers" in buf.getvalue() or "no outliers" in buf.getvalue().lower()


# ── _handle_export_outliers ───────────────────────────────────────────────────


def test_handle_export_outliers_saves_csv(tmp_path: Path) -> None:
    """_handle_export_outliers writes outliers.csv inside the charts directory."""
    view, buf = _create_view()
    csv_path = _make_csv(tmp_path)
    charts_dir = tmp_path / "charts"

    outliers_df = pd.DataFrame(
        {"mol_id": ["m1"], "abs_error": [5.0], "severity": ["elevated"]}
    )
    fake_result = MagicMock()
    fake_result.outliers_df = outliers_df

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_charts_dir", return_value=charts_dir),
        patch("grimperium.ml.error_analysis.analyze", return_value=fake_result),
    ):
        view._handle_export_outliers()

    expected_path = charts_dir / "outliers.csv"
    assert expected_path.exists()
    saved = pd.read_csv(expected_path)
    assert "mol_id" in saved.columns


def test_handle_export_outliers_no_csv_shows_error(tmp_path: Path) -> None:
    view, buf = _create_view()
    missing = tmp_path / "nope.csv"
    with patch.object(view, "_get_csv_path", return_value=missing):
        view._handle_export_outliers()
    assert "Error" in buf.getvalue() or "error" in buf.getvalue().lower()


# ── _handle_html_report ───────────────────────────────────────────────────────


def test_handle_html_report_calls_generate_with_correct_path(
    tmp_path: Path,
) -> None:
    """_handle_html_report calls generate_html_report with charts_dir/error_report.html."""
    view, buf = _create_view()
    csv_path = _make_csv(tmp_path)
    charts_dir = tmp_path / "charts"

    fake_result = MagicMock()
    fake_result.summary = {
        "mae": 1.0,
        "rmse": 1.0,
        "r2": 0.99,
        "n_molecules": 3,
        "n_outliers": 0,
    }

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_charts_dir", return_value=charts_dir),
        patch("grimperium.ml.error_analysis.analyze", return_value=fake_result),
        patch(
            "grimperium.ml.html_report.generate_html_report",
            return_value=charts_dir / "error_report.html",
        ) as mock_gen,
    ):
        view._handle_html_report()

    mock_gen.assert_called_once()
    call_args = mock_gen.call_args
    assert call_args[0][0] is fake_result
    assert str(call_args[0][1]).endswith("error_report.html")


def test_handle_html_report_no_csv_shows_error(tmp_path: Path) -> None:
    view, buf = _create_view()
    missing = tmp_path / "nope.csv"
    with patch.object(view, "_get_csv_path", return_value=missing):
        view._handle_html_report()
    assert "Error" in buf.getvalue() or "error" in buf.getvalue().lower()


# ── handle_action routing ─────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "action,handler",
    [
        ("show_outliers", "_handle_show_outliers"),
        ("top_errors", "_handle_top_errors"),
        ("html_report", "_handle_html_report"),
        ("export_outliers", "_handle_export_outliers"),
    ],
)
def test_handle_action_routes_to_new_handlers(action: str, handler: str) -> None:
    view, _ = _create_view()
    with patch.object(view, handler) as mock_handler:
        result = view.handle_action(action)
    mock_handler.assert_called_once()
    assert result is None
