"""Tests for PM7 Baseline Analysis feature in DatabasesView."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from rich.console import Console

from grimperium.cli.views.databases_view import DatabasesView


def _create_view_with_registry(
    tmp_path: Path,
) -> tuple[DatabasesView, io.StringIO]:
    """Create a DatabasesView with real Console capturing to StringIO."""
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    controller = MagicMock()
    controller.console = console
    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        view = DatabasesView(controller)
    # Stub wait_for_enter to avoid blocking
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    return view, buf


# ── Fixture: 5 molecules with known PM7 vs CBS errors ──
FIXTURE_CBS = [-100.0, -200.0, -50.0, -400.0, -250.0]
FIXTURE_PM7 = [-105.0, -190.0, -52.0, -380.0, -245.0]


@pytest.fixture()
def pm7_csv(tmp_path: Path) -> Path:
    """Create CSV with PM7 baseline fixture data."""
    df = pd.DataFrame(
        {
            "mol_id": [f"mol_{i}" for i in range(5)],
            "H298_cbs": FIXTURE_CBS,
            "H298_pm7": FIXTURE_PM7,
            "status": ["OK"] * 5,
        }
    )
    csv_path = tmp_path / "thermo_pm7.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


def test_pm7_baseline_menu_appears_with_columns(tmp_path: Path) -> None:
    """PM7 Baseline option visible when H298_pm7 and H298_cbs columns exist."""
    df = pd.DataFrame({"H298_cbs": [-100.0], "H298_pm7": [-105.0], "status": ["OK"]})
    df.to_csv(tmp_path / "thermo_pm7.csv", index=False)

    view, _ = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    view.selected_db.csv_path = "thermo_pm7.csv"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        options = view.get_detail_menu_options()

    values = [o.value for o in options]
    assert "pm7_baseline" in values


def test_pm7_baseline_menu_hidden_without_columns(tmp_path: Path) -> None:
    """PM7 Baseline option NOT shown when columns are missing."""
    df = pd.DataFrame({"H298_cbs": [-100.0], "other": [1.0], "status": ["OK"]})
    df.to_csv(tmp_path / "thermo_pm7.csv", index=False)

    view, _ = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    view.selected_db.csv_path = "thermo_pm7.csv"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        options = view.get_detail_menu_options()

    values = [o.value for o in options]
    assert "pm7_baseline" not in values


def test_handle_action_pm7_baseline_calls_handler(tmp_path: Path) -> None:
    """handle_action('pm7_baseline') routes to _handle_pm7_baseline."""
    view, _ = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    with patch.object(view, "_handle_pm7_baseline") as mock_handler:
        result = view.handle_action("pm7_baseline")
    mock_handler.assert_called_once()
    assert result is None


def test_pm7_baseline_no_csv(tmp_path: Path) -> None:
    """Missing CSV → error message displayed."""
    view, buf = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    view.selected_db.csv_path = "nonexistent.csv"
    view.selected_db.alias = "PM7"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        view._handle_pm7_baseline()

    output = buf.getvalue()
    assert "not" in output.lower() or "Error" in output or "error" in output.lower()


def test_pm7_baseline_computation_known_values(
    tmp_path: Path, pm7_csv: Path
) -> None:
    """Verify exact metric values with 5-molecule fixture."""
    view, buf = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    view.selected_db.csv_path = "thermo_pm7.csv"
    view.selected_db.alias = "PM7"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", pm7_csv.parent):
        view._handle_pm7_baseline()

    output = buf.getvalue()

    # MARE = 8.4 kcal/mol
    assert "8.4000" in output
    # Bias = +5.6
    assert "5.6000" in output
    # MRE% = 4.2%
    assert "4.20" in output
    # MdRE% = 5.0%
    assert "5.00" in output


def test_pm7_baseline_excludes_zero_cbs(tmp_path: Path) -> None:
    """Rows with H298_cbs=0 are excluded from relative error calculation."""
    df = pd.DataFrame(
        {
            "H298_cbs": [0.0, -100.0, -200.0],
            "H298_pm7": [5.0, -105.0, -190.0],
        }
    )
    csv_path = tmp_path / "thermo_pm7.csv"
    df.to_csv(csv_path, index=False)

    view, buf = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    view.selected_db.csv_path = "thermo_pm7.csv"
    view.selected_db.alias = "PM7"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        view._handle_pm7_baseline()

    output = buf.getvalue()
    # RE% computed only for rows 1,2: RE%=[5.0, 5.0] → MRE%=5.0%
    assert "5.00" in output


def test_pm7_baseline_filters_nan_rows(tmp_path: Path) -> None:
    """Rows with NaN in H298_pm7 or H298_cbs are excluded."""
    df = pd.DataFrame(
        {
            "H298_cbs": [-100.0, np.nan, -300.0],
            "H298_pm7": [-105.0, -195.0, np.nan],
        }
    )
    csv_path = tmp_path / "thermo_pm7.csv"
    df.to_csv(csv_path, index=False)

    view, buf = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    view.selected_db.csv_path = "thermo_pm7.csv"
    view.selected_db.alias = "PM7"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        view._handle_pm7_baseline()

    output = buf.getvalue()
    # Only row 0 valid: cbs=-100, pm7=-105, RE%=5%
    assert "5.00" in output


def test_pm7_baseline_interpretation_text(
    tmp_path: Path, pm7_csv: Path
) -> None:
    """Interpretation panel contains formatted metric values."""
    view, buf = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    view.selected_db.csv_path = "thermo_pm7.csv"
    view.selected_db.alias = "PM7"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", pm7_csv.parent):
        view._handle_pm7_baseline()

    output = buf.getvalue()
    assert "PM7" in output
    assert "kcal/mol" in output


def test_pm7_baseline_percentage_thresholds(tmp_path: Path) -> None:
    """Verify RE% threshold percentages with custom dataset."""
    # RE%: [0.5, 3, 7, 12, 20]
    # <1%: 1/5=20%, <5%: 2/5=40%, <10%: 3/5=60%
    df = pd.DataFrame(
        {
            "H298_cbs": [-200.0, -100.0, -100.0, -100.0, -100.0],
            "H298_pm7": [-201.0, -103.0, -107.0, -112.0, -120.0],
        }
    )
    csv_path = tmp_path / "thermo_pm7.csv"
    df.to_csv(csv_path, index=False)

    view, buf = _create_view_with_registry(tmp_path)
    view.selected_db = MagicMock()
    view.selected_db.csv_path = "thermo_pm7.csv"
    view.selected_db.alias = "PM7"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        view._handle_pm7_baseline()

    output = buf.getvalue()
    assert "20.0%" in output  # RE% < 1%
    assert "40.0%" in output  # RE% < 5%
    assert "60.0%" in output  # RE% < 10%
