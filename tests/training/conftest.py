"""Shared builders for the training tests."""

from __future__ import annotations

import sqlite3
from collections.abc import Callable
from contextlib import closing
from datetime import datetime, timezone
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


def set_running(db: Path, mol_id: str, cell_id: str | None = "central") -> None:
    """Put a queue row in flight without a logged attempt (legacy or fixture)."""
    with closing(sqlite3.connect(db)) as connection, connection:
        connection.execute(
            "UPDATE training_queue SET status = 'running', cell_id = ?,"
            " updated_at = ? WHERE mol_id = ?",
            (cell_id, datetime.now(timezone.utc).isoformat(), mol_id),
        )


def set_ok(
    db: Path,
    mol_id: str,
    *,
    hof: float,
    elapsed_s: float = 1.0,
    conformers: int = 1,
    batch_id: str = "b",
) -> None:
    """Store a result row as ``accept_success`` leaves it (no geometry)."""
    with closing(sqlite3.connect(db)) as connection, connection:
        connection.execute(
            "UPDATE training_queue SET status = 'ok', h298_pm7 = ?,"
            " delta_kcal_mol = h298_cbs - ?, conformers_evaluated = ?,"
            " elapsed_s = ?, batch_id = ?, updated_at = ? WHERE mol_id = ?",
            (
                hof,
                hof,
                conformers,
                elapsed_s,
                batch_id,
                datetime.now(timezone.utc).isoformat(),
                mol_id,
            ),
        )


def set_status(db: Path, mol_id: str, status: str, reruns: int = 0) -> None:
    """Force a queue row's status and rerun count."""
    with closing(sqlite3.connect(db)) as connection, connection:
        connection.execute(
            "UPDATE training_queue SET status = ?, reruns = ? WHERE mol_id = ?",
            (status, reruns, mol_id),
        )
