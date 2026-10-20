"""Shared builders for the training tests."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

CbsFactory = Callable[[list[dict[str, Any]]], Path]


@pytest.fixture
def make_cbs(tmp_path: Path) -> CbsFactory:
    """Write a ``cbs_referencia.csv``-shaped file with the given rows."""

    def _make(rows: list[dict[str, Any]]) -> Path:
        csv_path = tmp_path / "cbs_referencia.csv"
        pd.DataFrame(rows).to_csv(csv_path, index=False)
        return csv_path

    return _make


def row(mol_id: str, **extra: Any) -> dict[str, Any]:
    """One CBS row with the columns the training queue reads."""
    base: dict[str, Any] = {
        "mol_id": mol_id,
        "smiles": f"smiles-{mol_id}",
        "multiplicity": 1,
        "charge": 0,
        "nheavy": 1,
        "H298_cbs": -17.8,
        "smiles_count": 1,
        "conformer_rank_by_h298": 1,
    }
    base.update(extra)
    return base
