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


#: The admission references (``ADMISSION_REFERENCES``) as CBS rows.
REFERENCE_ROWS = (
    ("cbs_09698", "CO", 1, 2),
    ("cbs_00187", "[H][C]=C", 2, 2),
    ("cbs_24671", "C[N+](=O)[O-]", 1, 4),
    ("cbs_08982", "OCCO", 1, 4),
    ("cbs_01400", "CCC(=O)O", 1, 5),
)


def reference_rows() -> list[dict[str, Any]]:
    """CBS rows of the five admission references (seed them with the queue)."""
    return [
        row(mol_id, smiles=smiles, multiplicity=mult, nheavy=nheavy)
        for mol_id, smiles, mult, nheavy in REFERENCE_ROWS
    ]


def accept_references(db: Path, hof: float = -12.5) -> None:
    """Accept every reference at the central with ``hof``, as the queue would."""
    from semi_imperium.training.geometry_store import GeometryRecord, GeometryStore

    store = GeometryStore(db)
    for mol_id, smiles, _mult, _nheavy in REFERENCE_ROWS:
        set_ok(db, mol_id, hof=hof)
        store.save(
            GeometryRecord(
                mol_id=mol_id,
                smiles=smiles,
                hamiltonian="PM7",
                protocol="2",
                signature="sig",
                calculation_id=f"calc-{mol_id}",
                conformer_index=0,
                elements=("C",),
                coordinates=((0.0, 0.0, 0.0),),
                hof_kcal_mol=hof,
                lowest_frequency_cm1=100.0,
                minimum_confirmed=True,
            )
        )
