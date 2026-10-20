"""Final geometry of every protocol-2 molecule, one row per ``mol_id``.

Only the optimized geometry of the selected conformer is kept; the CSV never
holds coordinates. ``to_xyz`` exports one molecule for other programs.
"""

from __future__ import annotations

import json
import sqlite3
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from semi_imperium.training.db import connect


@dataclass(frozen=True)
class GeometryRecord:
    """One molecule's final geometry plus the numbers that justify it."""

    mol_id: str
    smiles: str
    hamiltonian: str
    protocol: str
    signature: str
    calculation_id: str
    conformer_index: int | None
    elements: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    hof_kcal_mol: float | None
    lowest_frequency_cm1: float | None
    minimum_confirmed: bool
    created_at: str = ""

    def to_xyz(self) -> str:
        """Render the geometry as one XYZ block, preserving atom order."""
        comment = (
            f"{self.mol_id} {self.hamiltonian} protocol={self.protocol} "
            f"hof_kcal_mol={self.hof_kcal_mol}"
        )
        lines = [str(len(self.elements)), comment]
        for element, (x, y, z) in zip(self.elements, self.coordinates, strict=True):
            lines.append(f"{element:<2} {x:16.10f} {y:16.10f} {z:16.10f}")
        return "\n".join(lines) + "\n"


def insert_geometry(connection: sqlite3.Connection, record: GeometryRecord) -> None:
    """Insert or replace ``record`` on an open connection.

    Lets the queue store the geometry and mark the molecule ``ok`` in one
    transaction, so an accepted molecule always has its geometry.
    """
    created_at = record.created_at or datetime.now(timezone.utc).isoformat()
    connection.execute(
        "INSERT OR REPLACE INTO final_geometry (mol_id, smiles, hamiltonian,"
        " protocol, signature, calculation_id, conformer_index, n_atoms,"
        " elements, coordinates, hof_kcal_mol, lowest_frequency_cm1,"
        " minimum_confirmed, created_at)"
        " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (
            record.mol_id,
            record.smiles,
            record.hamiltonian,
            record.protocol,
            record.signature,
            record.calculation_id,
            record.conformer_index,
            len(record.elements),
            json.dumps(list(record.elements)),
            json.dumps([list(point) for point in record.coordinates]),
            record.hof_kcal_mol,
            record.lowest_frequency_cm1,
            int(record.minimum_confirmed),
            created_at,
        ),
    )


class GeometryStore:
    """SQLite-backed table of final geometries (see :mod:`...training.db`)."""

    def __init__(self, path: Path | str) -> None:
        self._path = Path(path)

    def save(self, record: GeometryRecord) -> None:
        """Insert or replace the row of ``record.mol_id``."""
        with closing(connect(self._path)) as connection:
            insert_geometry(connection, record)

    def get(self, mol_id: str) -> GeometryRecord | None:
        """Return the stored geometry of ``mol_id``, or ``None``."""
        if not self._path.exists():
            return None
        with closing(connect(self._path)) as connection:
            row = connection.execute(
                "SELECT * FROM final_geometry WHERE mol_id = ?", (mol_id,)
            ).fetchone()
        if row is None:
            return None
        return GeometryRecord(
            mol_id=row["mol_id"],
            smiles=row["smiles"],
            hamiltonian=row["hamiltonian"],
            protocol=row["protocol"],
            signature=row["signature"],
            calculation_id=row["calculation_id"],
            conformer_index=row["conformer_index"],
            elements=tuple(json.loads(row["elements"])),
            coordinates=tuple(
                (float(p[0]), float(p[1]), float(p[2]))
                for p in json.loads(row["coordinates"])
            ),
            hof_kcal_mol=row["hof_kcal_mol"],
            lowest_frequency_cm1=row["lowest_frequency_cm1"],
            minimum_confirmed=bool(row["minimum_confirmed"]),
            created_at=row["created_at"],
        )

    def count(self) -> int:
        """Number of stored geometries (0 when the database does not exist)."""
        if not self._path.exists():
            return 0
        with closing(connect(self._path)) as connection:
            return int(
                connection.execute("SELECT COUNT(*) FROM final_geometry").fetchone()[0]
            )

    def mol_ids(self) -> set[str]:
        """Every ``mol_id`` that has a stored geometry."""
        if not self._path.exists():
            return set()
        with closing(connect(self._path)) as connection:
            return {
                str(row["mol_id"])
                for row in connection.execute("SELECT mol_id FROM final_geometry")
            }
