"""Tests for Detailed Metrics feature in ResultsView."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from rich.console import Console

from grimperium.cli.views.results_view import ResultsView


def _create_results_view() -> tuple[ResultsView, io.StringIO]:
    """Create a ResultsView with a real Console capturing to StringIO."""
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    controller = MagicMock()
    controller.console = console
    view = ResultsView(controller)
    # Stub wait_for_enter to avoid blocking
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    return view, buf


# ── Fixture: 5 molecules with known errors ──
FIXTURE_CBS = [-100.0, -200.0, -300.0, -400.0, -500.0]
FIXTURE_PRED = [-101.0, -198.0, -303.0, -399.0, -505.0]


@pytest.fixture()
def known_csv(tmp_path: Path) -> Path:
    """Create CSV with known 5-molecule fixture."""
    df = pd.DataFrame({"H298_cbs": FIXTURE_CBS, "H298_predicted": FIXTURE_PRED})
    csv_path = tmp_path / "thermo_pm7.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


def test_detailed_metrics_menu_enabled() -> None:
    """MenuOption for 'Detailed Metrics' must NOT be disabled."""
    view, _ = _create_results_view()
    options = view.get_menu_options()
    detailed = next(o for o in options if o.value == "detailed")
    assert detailed.disabled is False


def test_handle_action_detailed_calls_handler() -> None:
    """handle_action('detailed') must call _handle_detailed_metrics."""
    view, _ = _create_results_view()
    with patch.object(view, "_handle_detailed_metrics") as mock_handler:
        result = view.handle_action("detailed")
    mock_handler.assert_called_once()
    assert result is None


def test_detailed_metrics_no_csv(tmp_path: Path) -> None:
    """Missing CSV → show_error is called."""
    view, buf = _create_results_view()
    missing_csv = tmp_path / "nonexistent.csv"
    with patch.object(view, "_get_csv_path", return_value=missing_csv):
        view._handle_detailed_metrics()
    output = buf.getvalue()
    assert "Error" in output


def test_detailed_metrics_no_predicted_column(tmp_path: Path) -> None:
    """CSV without H298_predicted column → show_error."""
    csv_path = tmp_path / "thermo_pm7.csv"
    df = pd.DataFrame({"H298_cbs": [-100.0], "other_col": [1.0]})
    df.to_csv(csv_path, index=False)

    view, buf = _create_results_view()
    with patch.object(view, "_get_csv_path", return_value=csv_path):
        view._handle_detailed_metrics()
    output = buf.getvalue()
    assert "Error" in output


def test_detailed_metrics_computation_perfect(tmp_path: Path) -> None:
    """When predicted == cbs, MAE=0 and R²=1."""
    csv_path = tmp_path / "thermo_pm7.csv"
    values = [-100.0, -200.0, -300.0]
    df = pd.DataFrame({"H298_cbs": values, "H298_predicted": values})
    df.to_csv(csv_path, index=False)

    view, buf = _create_results_view()
    with patch.object(view, "_get_csv_path", return_value=csv_path):
        view._handle_detailed_metrics()

    output = buf.getvalue()
    assert "0.0000" in output  # MAE = 0
    assert "1.0000" in output  # R² = 1


def test_detailed_metrics_computation_known_values(known_csv: Path) -> None:
    """Verify exact metric values with 5-molecule fixture."""
    view, buf = _create_results_view()
    with patch.object(view, "_get_csv_path", return_value=known_csv):
        view._handle_detailed_metrics()

    output = buf.getvalue()

    # MAE = 2.4
    assert "2.4000" in output
    # RMSE = sqrt(8) ≈ 2.8284
    assert "2.8284" in output
    # Max Error = 5.0
    assert "5.0000" in output
    # Bias = -1.2
    assert "-1.2000" in output
    # R² = 0.9996
    assert "0.9996" in output


def test_detailed_metrics_filters_nan_rows(tmp_path: Path) -> None:
    """Rows with NaN in predicted or cbs are excluded from computation."""
    csv_path = tmp_path / "thermo_pm7.csv"
    df = pd.DataFrame(
        {
            "H298_cbs": [-100.0, -200.0, np.nan, -400.0],
            "H298_predicted": [-100.0, np.nan, -300.0, -400.0],
        }
    )
    df.to_csv(csv_path, index=False)

    view, buf = _create_results_view()
    with patch.object(view, "_get_csv_path", return_value=csv_path):
        view._handle_detailed_metrics()

    output = buf.getvalue()
    # 2 valid rows (row 0 and 3) with perfect predictions → MAE = 0
    assert "0.0000" in output


def test_detailed_metrics_bias_direction(tmp_path: Path) -> None:
    """Positive bias when predictions consistently overestimate."""
    csv_path = tmp_path / "thermo_pm7.csv"
    # pred > cbs → positive error → positive bias
    df = pd.DataFrame(
        {
            "H298_cbs": [-100.0, -200.0, -300.0],
            "H298_predicted": [-98.0, -197.0, -296.0],
        }
    )
    df.to_csv(csv_path, index=False)

    view, buf = _create_results_view()
    with patch.object(view, "_get_csv_path", return_value=csv_path):
        view._handle_detailed_metrics()

    output = buf.getvalue()
    # Bias = mean([+2, +3, +4]) = +3.0
    assert "3.0000" in output
    assert "overestimates" in output


def test_detailed_metrics_percentages(tmp_path: Path) -> None:
    """Verify percentage within thresholds with custom dataset."""
    csv_path = tmp_path / "thermo_pm7.csv"
    # Errors: [0.5, 1.5, 2.5, 4.0, 6.0]
    # ±1: 1/5=20%, ±2: 2/5=40%, ±5: 4/5=80%
    df = pd.DataFrame(
        {
            "H298_cbs": [-100.0, -200.0, -300.0, -400.0, -500.0],
            "H298_predicted": [-100.5, -201.5, -302.5, -404.0, -506.0],
        }
    )
    df.to_csv(csv_path, index=False)

    view, buf = _create_results_view()
    with patch.object(view, "_get_csv_path", return_value=csv_path):
        view._handle_detailed_metrics()

    output = buf.getvalue()
    assert "20.0%" in output  # within ±1
    assert "40.0%" in output  # within ±2
    assert "80.0%" in output  # within ±5
