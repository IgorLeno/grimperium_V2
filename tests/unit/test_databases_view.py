"""Tests for databases view module."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from grimperium.cli.views.databases_view import DatabasesView


def _create_view_with_registry(tmp_path: Path) -> DatabasesView:
    """Create a DatabasesView pointing its registry at tmp_path."""
    controller = MagicMock()
    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        view = DatabasesView(controller)
    return view


def test_get_databases_returns_defaults(tmp_path: Path) -> None:
    """Registry auto-creates 3 default entries when no JSON exists."""
    view = _create_view_with_registry(tmp_path)
    databases = view.get_databases()

    assert len(databases) == 3
    aliases = {db.alias for db in databases}
    assert aliases == {"CBS", "PM7", "NIST"}


def test_get_databases_enriches_pm7_from_csv(tmp_path: Path) -> None:
    """PM7 entry gets molecule count from CSV (only OK rows)."""
    # Create CSV with mixed statuses
    df = pd.DataFrame(
        {
            "mol_id": ["m1", "m2", "m3"],
            "status": ["OK", "OK", "PENDING"],
            "smiles": ["C", "CC", "CCC"],
        }
    )
    df.to_csv(tmp_path / "thermo_pm7.csv", index=False)

    view = _create_view_with_registry(tmp_path)
    databases = view.get_databases()

    pm7 = next(db for db in databases if db.alias == "PM7")
    assert pm7.molecules == 2  # Only OK rows
    assert pm7.status == "ready"


def test_get_databases_pm7_empty_csv(tmp_path: Path) -> None:
    """PM7 with header-only CSV stays in_development."""
    df = pd.DataFrame(columns=["mol_id", "status", "smiles"])
    df.to_csv(tmp_path / "thermo_pm7.csv", index=False)

    view = _create_view_with_registry(tmp_path)
    databases = view.get_databases()

    pm7 = next(db for db in databases if db.alias == "PM7")
    assert pm7.molecules == 0
    assert pm7.status == "in_development"


def test_get_databases_no_csv(tmp_path: Path) -> None:
    """PM7 without any CSV file stays at 0 molecules."""
    view = _create_view_with_registry(tmp_path)
    databases = view.get_databases()

    pm7 = next(db for db in databases if db.alias == "PM7")
    assert pm7.molecules == 0


def test_menu_options_use_aliases(tmp_path: Path) -> None:
    """Menu options display alias-based labels."""
    view = _create_view_with_registry(tmp_path)
    options = view.get_menu_options()

    # First 3 options should be database entries
    db_options = [o for o in options if o.value.startswith("view_")]
    assert len(db_options) == 3
    labels = {o.label for o in db_options}
    assert labels == {"CBS", "PM7", "NIST"}
    values = {o.value for o in db_options}
    assert values == {"view_CBS", "view_PM7", "view_NIST"}


def test_detail_menu_has_run_config_analyze(tmp_path: Path) -> None:
    """Detail menu includes Run, Config, and Analyze options."""
    view = _create_view_with_registry(tmp_path)
    options = view.get_detail_menu_options()

    values = [o.value for o in options]
    assert "calculate_run" in values
    assert "calculate_config" in values
    assert "analyze" in values


def test_handle_action_view_alias(tmp_path: Path) -> None:
    """handle_action with view_{alias} selects correct database."""
    view = _create_view_with_registry(tmp_path)
    view.handle_action("view_PM7")

    assert view.selected_db is not None
    assert view.selected_db.alias == "PM7"


def test_refresh_databases_from_filesystem_finds_csvs(tmp_path: Path) -> None:
    """Refresh discovers CSV files in data/ directory."""
    controller = MagicMock()
    view = DatabasesView(controller)

    (tmp_path / "thermo_cbs_chon.csv").write_text(
        "smiles,H298_cbs\nCCO,-100.5\nCC,-50.2\n", encoding="utf-8"
    )
    (tmp_path / "thermo_pm7.csv").write_text(
        "smiles,H298_pm7\nCCO,-100.3\n", encoding="utf-8"
    )

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        count = view.refresh_databases_from_filesystem()

    assert count == 2


def test_refresh_databases_from_filesystem_missing_directory(tmp_path: Path) -> None:
    """Refresh handles missing data/ directory gracefully."""
    controller = MagicMock()
    view = DatabasesView(controller)

    nonexistent_dir = tmp_path / "does_not_exist"

    with patch("grimperium.cli.views.databases_view.DATA_DIR", nonexistent_dir):
        count = view.refresh_databases_from_filesystem()

    assert count == 0


def test_refresh_databases_from_filesystem_empty_directory(tmp_path: Path) -> None:
    """Refresh handles empty data/ directory."""
    controller = MagicMock()
    view = DatabasesView(controller)

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()

    with patch("grimperium.cli.views.databases_view.DATA_DIR", empty_dir):
        count = view.refresh_databases_from_filesystem()

    assert count == 0


def test_refresh_databases_from_filesystem_counts_molecules(tmp_path: Path) -> None:
    """Refresh counts actual CSV rows for each database."""
    controller = MagicMock()
    view = DatabasesView(controller)

    (tmp_path / "test_database.csv").write_text(
        "smiles,energy\nCCO,-100\nCC,-50\nC,-25\n", encoding="utf-8"
    )

    with patch("grimperium.cli.views.databases_view.DATA_DIR", tmp_path):
        count = view.refresh_databases_from_filesystem()

    assert count == 1
