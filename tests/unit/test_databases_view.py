"""Tests for databases view module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd

from grimperium.cli.views.databases_view import DatabasesView


def _create_view_with_registry(tmp_path: Path, monkeypatch) -> DatabasesView:
    """Create a DatabasesView pointing its registry at tmp_path."""
    data_dir = tmp_path / "data"
    data_dir.mkdir(exist_ok=True)
    monkeypatch.setenv("GRIMPERIUM_CONFIG_DIR", str(tmp_path / "config"))
    monkeypatch.setattr("grimperium.cli.views.databases_view.DATA_DIR", data_dir)
    controller = MagicMock()
    view = DatabasesView(controller)
    return view


def test_get_databases_returns_official_manifests(tmp_path: Path, monkeypatch) -> None:
    view = _create_view_with_registry(tmp_path, monkeypatch)
    databases = view.get_databases()

    assert len(databases) == 3
    assert {db.database_id for db in databases} == {
        "official.cbs_chon",
        "official.crest_pm7",
        "official.nist_experimental",
    }


def test_get_databases_enriches_pm7_from_csv(tmp_path: Path, monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "mol_id": ["m1", "m2", "m3"],
            "status": ["OK", "OK", "PENDING"],
            "smiles": ["C", "CC", "CCC"],
        }
    )
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    df.to_csv(data_dir / "thermo_pm7.csv", index=False)

    view = _create_view_with_registry(tmp_path, monkeypatch)
    databases = view.get_databases()

    pm7 = next(db for db in databases if db.database_id == "official.crest_pm7")
    assert pm7.molecules == 2
    assert pm7.status == "available"


def test_get_databases_pm7_empty_csv_is_available(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    pd.DataFrame(columns=["mol_id", "status", "smiles"]).to_csv(
        data_dir / "thermo_pm7.csv",
        index=False,
    )

    view = _create_view_with_registry(tmp_path, monkeypatch)
    pm7 = next(
        db for db in view.get_databases() if db.database_id == "official.crest_pm7"
    )

    assert pm7.molecules == 0
    assert pm7.status == "available"


def test_get_databases_no_csv_marks_missing(tmp_path: Path, monkeypatch) -> None:
    view = _create_view_with_registry(tmp_path, monkeypatch)
    pm7 = next(
        db for db in view.get_databases() if db.database_id == "official.crest_pm7"
    )

    assert pm7.molecules == 0
    assert pm7.status == "missing"


def test_menu_options_use_database_ids(tmp_path: Path, monkeypatch) -> None:
    view = _create_view_with_registry(tmp_path, monkeypatch)
    options = view.get_menu_options()

    db_options = [o for o in options if o.value.startswith("view_")]
    assert len(db_options) == 3
    assert {o.label for o in db_options} == {"CBS", "PM7", "NIST"}
    assert {o.value for o in db_options} == {
        "view_official.cbs_chon",
        "view_official.crest_pm7",
        "view_official.nist_experimental",
    }


def test_detail_menu_has_session_and_capability_gated_actions(
    tmp_path: Path,
    monkeypatch,
) -> None:
    view = _create_view_with_registry(tmp_path, monkeypatch)
    view.selected_db = view.registry.get_by_id("official.nist_experimental")

    options = {o.value: o for o in view.get_detail_menu_options()}

    assert "use_session_dataset" in options
    assert options["calculate_run"].disabled is True
    assert options["analyze"].disabled is True


def test_handle_action_view_id(tmp_path: Path, monkeypatch) -> None:
    view = _create_view_with_registry(tmp_path, monkeypatch)

    view.handle_action("view_official.crest_pm7")

    assert view.selected_db is not None
    assert view.selected_db.database_id == "official.crest_pm7"


def test_refresh_databases_from_filesystem_finds_known_csvs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "thermo_cbs_chon.csv").write_text(
        "smiles,H298_cbs\nCCO,-100.5\nCC,-50.2\n",
        encoding="utf-8",
    )
    (data_dir / "thermo_pm7.csv").write_text(
        "smiles,H298_pm7,status\nCCO,-100.3,OK\n",
        encoding="utf-8",
    )
    view = _create_view_with_registry(tmp_path, monkeypatch)

    count = view.refresh_databases_from_filesystem()

    assert count == 2


def test_refresh_databases_from_filesystem_missing_directory(
    tmp_path: Path,
    monkeypatch,
) -> None:
    view = _create_view_with_registry(tmp_path, monkeypatch)
    missing_dir = tmp_path / "missing"
    view.registry.data_dir = missing_dir

    count = view.refresh_databases_from_filesystem()

    assert count == 0


def test_refresh_databases_from_filesystem_counts_unregistered(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "custom.csv").write_text(
        "smiles,energy\nCCO,-100\nCC,-50\nC,-25\n",
        encoding="utf-8",
    )
    view = _create_view_with_registry(tmp_path, monkeypatch)

    count = view.refresh_databases_from_filesystem()

    assert count == 1
