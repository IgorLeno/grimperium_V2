"""Tests for the DatabaseRegistry module."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from grimperium.cli.database_registry import DatabaseRegistry


# ── Test 1: auto_creates_defaults ──


def test_auto_creates_defaults(tmp_path: Path) -> None:
    """No JSON exists → creates with 3 entries."""
    registry = DatabaseRegistry(tmp_path)
    entries = registry.load()

    assert len(entries) == 3
    aliases = {e.alias for e in entries}
    assert aliases == {"CBS", "PM7", "NIST"}

    # File should now exist
    json_path = tmp_path / "databases_registry.json"
    assert json_path.exists()

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(data) == 3


# ── Test 2: load_existing ──


def test_load_existing(tmp_path: Path) -> None:
    """Read valid JSON."""
    registry_data = [
        {
            "name": "Test DB",
            "alias": "TEST",
            "csv_path": "",
            "description": "A test database",
            "pipeline": "reference",
            "created_at": "2026-01-01",
            "status": "ready",
            "properties": ["prop1"],
        }
    ]
    json_path = tmp_path / "databases_registry.json"
    json_path.write_text(json.dumps(registry_data), encoding="utf-8")

    registry = DatabaseRegistry(tmp_path)
    entries = registry.load()

    assert len(entries) == 1
    assert entries[0].alias == "TEST"
    assert entries[0].name == "Test DB"
    assert entries[0].properties == ["prop1"]


# ── Test 3: enrich_counts_ok_molecules ──


def test_enrich_counts_ok_molecules(tmp_path: Path) -> None:
    """PM7 CSV molecule count — only OK rows counted for crest_pm7 pipeline."""
    # Create a CSV with mixed statuses
    df = pd.DataFrame(
        {
            "mol_id": ["m1", "m2", "m3", "m4"],
            "status": ["OK", "OK", "PENDING", "RUNNING"],
            "smiles": ["C", "CC", "CCC", "CCCC"],
        }
    )
    df.to_csv(tmp_path / "thermo_pm7.csv", index=False)

    # Create registry with PM7 entry pointing to CSV
    registry_data = [
        {
            "name": "CREST PM7",
            "alias": "PM7",
            "csv_path": "thermo_pm7.csv",
            "description": "test",
            "pipeline": "crest_pm7",
            "created_at": "2026-01-01",
            "status": "in_development",
            "properties": [],
        }
    ]
    json_path = tmp_path / "databases_registry.json"
    json_path.write_text(json.dumps(registry_data), encoding="utf-8")

    registry = DatabaseRegistry(tmp_path)
    entries = registry.load()

    assert len(entries) == 1
    assert entries[0].molecules == 2  # Only OK rows
    assert entries[0].status == "ready"  # Updated because molecules > 0


# ── Test 4: get_by_alias ──


def test_get_by_alias(tmp_path: Path) -> None:
    """Alias lookup (case-insensitive)."""
    registry = DatabaseRegistry(tmp_path)
    registry.load()  # Creates defaults

    result = registry.get_by_alias("pm7")
    assert result is not None
    assert result.alias == "PM7"

    result = registry.get_by_alias("CBS")
    assert result is not None
    assert result.alias == "CBS"

    result = registry.get_by_alias("nonexistent")
    assert result is None


# ── Test 5: reload ──


def test_reload(tmp_path: Path) -> None:
    """Cache invalidation on reload."""
    registry = DatabaseRegistry(tmp_path)
    entries1 = registry.load()
    assert len(entries1) == 3

    # Modify the JSON file directly
    json_path = tmp_path / "databases_registry.json"
    data = json.loads(json_path.read_text(encoding="utf-8"))
    data.append(
        {
            "name": "New DB",
            "alias": "NEW",
            "csv_path": "",
            "description": "added",
            "pipeline": "reference",
            "created_at": "2026-03-01",
            "status": "in_development",
            "properties": [],
        }
    )
    json_path.write_text(json.dumps(data), encoding="utf-8")

    # Without reload, cache returns old data
    entries_cached = registry.load()
    assert len(entries_cached) == 3

    # After reload, new entry appears
    entries2 = registry.reload()
    assert len(entries2) == 4
    assert any(e.alias == "NEW" for e in entries2)
