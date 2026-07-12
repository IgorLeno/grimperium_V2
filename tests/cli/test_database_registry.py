"""Tests for DatabaseRegistry v2 manifests and overlay."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from grimperium.cli.database_registry import DatabaseRegistry


def _registry(tmp_path: Path, monkeypatch) -> DatabaseRegistry:
    config_dir = tmp_path / "config"
    monkeypatch.setenv("GRIMPERIUM_CONFIG_DIR", str(config_dir))
    return DatabaseRegistry(tmp_path / "data")


def test_loads_official_manifests_without_writing_data_registry(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    registry = _registry(tmp_path, monkeypatch)

    entries = registry.load()

    assert {entry.database_id for entry in entries} == {
        "official.cbs_chon",
        "official.crest_pm7",
        "official.nist_experimental",
    }
    assert not (data_dir / "databases_registry.json").exists()
    assert all(entry.status == "missing" for entry in entries)


def test_availability_is_computed_from_filesystem(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    pd.DataFrame(
        {
            "mol_id": ["m1", "m2", "m3"],
            "smiles": ["C", "CC", "CCC"],
            "status": ["OK", "PENDING", "OK"],
        }
    ).to_csv(data_dir / "thermo_pm7.csv", index=False)
    registry = _registry(tmp_path, monkeypatch)

    pm7 = registry.get_by_id("official.crest_pm7")

    assert pm7 is not None
    assert pm7.status == "available"
    assert pm7.molecules == 2


def test_invalid_schema_status_is_computed(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    pd.DataFrame({"energy": [-1.0]}).to_csv(data_dir / "thermo_pm7.csv", index=False)
    registry = _registry(tmp_path, monkeypatch)

    pm7 = registry.get_by_id("official.crest_pm7")

    assert pm7 is not None
    assert pm7.status == "invalid_schema"


def test_add_user_database_persists_overlay(
    tmp_path: Path,
    monkeypatch,
) -> None:
    registry = _registry(tmp_path, monkeypatch)
    csv_path = tmp_path / "custom.csv"
    csv_path.write_text("smiles,H298_cbs\nC,-10\n", encoding="utf-8")

    entry = registry.add_user_database(
        path=csv_path,
        name="Custom",
        alias="CUSTOM",
        description="Custom analysis data",
        role="analysis",
        capabilities={"readable", "analysis_input"},
    )

    assert entry.database_id.startswith("user.")
    overlay = json.loads(registry.overlay_path.read_text(encoding="utf-8"))
    assert overlay["schema_version"] == 2
    assert overlay["entries"][0]["database_id"] == entry.database_id
    assert registry.reload()[-1].alias == "CUSTOM"


def test_official_path_override_is_stored_in_overlay(
    tmp_path: Path,
    monkeypatch,
) -> None:
    registry = _registry(tmp_path, monkeypatch)
    override_path = tmp_path / "pm7_override.csv"
    override_path.write_text("smiles,status\nC,OK\n", encoding="utf-8")

    registry.update_entry("official.crest_pm7", path=override_path)
    pm7 = registry.reload()[1]

    assert pm7.database_id == "official.crest_pm7"
    assert pm7.path == override_path
    overlay = json.loads(registry.overlay_path.read_text(encoding="utf-8"))
    assert overlay["overrides"]["official.crest_pm7"]["path"] == str(override_path)


def test_migrates_legacy_data_registry_to_overlay(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    legacy_path = data_dir / "databases_registry.json"
    legacy_path.write_text(
        json.dumps(
            [
                {
                    "name": "Legacy PM7",
                    "alias": "PM7",
                    "csv_path": "legacy_pm7.csv",
                    "description": "legacy",
                    "pipeline": "crest_pm7",
                    "created_at": "2026-01-01",
                    "properties": ["smiles"],
                }
            ]
        ),
        encoding="utf-8",
    )
    registry = _registry(tmp_path, monkeypatch)

    pm7 = registry.get_by_id("official.crest_pm7")

    assert pm7 is not None
    assert pm7.name == "Legacy PM7"
    assert pm7.path == data_dir / "legacy_pm7.csv"
    assert not legacy_path.exists()
    assert (data_dir / "databases_registry.json.bak").exists()
