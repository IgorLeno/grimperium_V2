"""Shared builders for the training tests."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager

CsvFactory = Callable[[list[dict[str, Any]]], Path]


@pytest.fixture
def make_csv(tmp_path: Path) -> CsvFactory:
    """Write a ``thermo_pm7.csv``-shaped file with the given rows."""

    def _make(rows: list[dict[str, Any]]) -> Path:
        csv_path = tmp_path / "thermo_pm7.csv"
        schema = BatchCSVManager(csv_path=None).get_schema()
        frame = pd.DataFrame(rows).reindex(columns=schema)
        frame.to_csv(csv_path, index=False)
        return csv_path

    return _make


def row(mol_id: str, status: str, **extra: Any) -> dict[str, Any]:
    """One CSV row with the identity columns the routine needs."""
    base: dict[str, Any] = {
        "mol_id": mol_id,
        "status": status,
        "smiles": "C",
        "multiplicity": 1,
        "charge": 0,
        "nheavy": 1,
        "H298_cbs": -17.8,
        "reruns": 0,
    }
    base.update(extra)
    return base
