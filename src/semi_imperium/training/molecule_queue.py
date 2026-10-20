"""The training queue: which molecules to compute and how each one stands.

The queue is seeded from ``cbs_referencia.csv`` (read only) and lives in the
``training_queue`` table of the training database, so the CSV never changes.
One molecule per SMILES: among the rows that share a SMILES the reference
conformer wins (lowest ``conformer_rank_by_h298``, which is the lowest H298);
ties fall back to the lowest ``H298_cbs`` and then to ``mol_id``.

Only the coordinator writes here; workers never touch the database.
"""

from __future__ import annotations

import math
import sqlite3
from collections.abc import Iterator
from contextlib import closing, contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

from semi_imperium.training.db import connect

#: Columns read from the CBS CSV (``xyz`` and the heat-capacity fits are skipped).
CBS_COLUMNS = (
    "mol_id",
    "smiles",
    "charge",
    "multiplicity",
    "nheavy",
    "H298_cbs",
    "conformer_rank_by_h298",
)

PENDING = "pending"
RUNNING = "running"
OK = "ok"
RERUN = "rerun"
SKIP = "skip"

#: Rows the routine still has to compute. ``running`` is an orphan of an
#: interrupted run (the coordinator waits for in-flight work before exiting).
QUEUE_STATUSES = (PENDING, RERUN, RUNNING)


#: The cell that runs on the central machine itself (the local routine).
CENTRAL_CELL = "central"

#: How far ahead a claim or renewal pushes ``lease_until``. Informative only:
#: nothing returns a molecule when its lease lapses (decision of 2026-10-20);
#: the dashboard uses it to tell a quiet cell from a working one.
DEFAULT_LEASE_S = 300.0


class CbsReadError(ValueError):
    """The CBS CSV cannot be read as the queue expects."""


@dataclass(frozen=True)
class QueuedMolecule:
    """One molecule waiting for protocol 2."""

    mol_id: str
    smiles: str
    charge: int
    multiplicity: int
    nheavy: int
    h298_cbs: float


@dataclass(frozen=True)
class ClaimedMolecule(QueuedMolecule):
    """A molecule reserved to one cell; ``attempt_id`` names this try."""

    attempt_id: str


@dataclass(frozen=True)
class TrainingProgress:
    """Molecules per stage, with the denominator (``total``) explicit."""

    total: int
    ok: int
    pending: int
    rerun: int
    running: int
    skipped: int


def select_reference_rows(cbs_csv: Path) -> pd.DataFrame:
    """The queue's rows: one per SMILES, the reference conformer."""
    if not cbs_csv.is_file():
        raise FileNotFoundError(f"CBS CSV not found: {cbs_csv}")
    header = pd.read_csv(cbs_csv, nrows=0).columns
    missing = [name for name in CBS_COLUMNS if name not in header]
    if missing:
        raise CbsReadError(f"{cbs_csv}: missing columns {missing}")
    frame = pd.read_csv(cbs_csv, usecols=list(CBS_COLUMNS))
    if frame[list(CBS_COLUMNS)].isna().any().any():
        raise CbsReadError(f"{cbs_csv}: empty cells in {list(CBS_COLUMNS)}")
    ordered = frame.sort_values(
        ["smiles", "conformer_rank_by_h298", "H298_cbs", "mol_id"], kind="stable"
    )
    return ordered.drop_duplicates("smiles", keep="first").reset_index(drop=True)


# The CBS has ~27k rows and the dashboard polls progress, so the derived
# (not yet seeded) count is cached per file version.
_derived_total: dict[str, tuple[int, int]] = {}


def _cbs_total(cbs_csv: Path) -> int:
    version = cbs_csv.stat().st_mtime_ns if cbs_csv.is_file() else -1
    cached = _derived_total.get(str(cbs_csv))
    if cached and cached[0] == version:
        return cached[1]
    total = len(select_reference_rows(cbs_csv))
    _derived_total[str(cbs_csv)] = (version, total)
    return total


def read_progress(db_path: Path, cbs_csv: Path) -> TrainingProgress:
    """Count molecules by stage; never writes.

    Before the first run the queue table is empty, so every selected CBS
    molecule counts as pending (that is what seeding will insert).
    """
    if db_path.exists():
        with closing(connect(db_path)) as connection:
            counts = {
                str(r["status"]): int(r["n"])
                for r in connection.execute(
                    "SELECT status, COUNT(*) AS n FROM training_queue GROUP BY status"
                )
            }
        if counts:
            return TrainingProgress(
                total=sum(counts.values()),
                ok=counts.get(OK, 0),
                pending=counts.get(PENDING, 0),
                rerun=counts.get(RERUN, 0),
                running=counts.get(RUNNING, 0),
                skipped=counts.get(SKIP, 0),
            )
    total = _cbs_total(cbs_csv)
    return TrainingProgress(
        total=total, ok=0, pending=total, rerun=0, running=0, skipped=0
    )


@dataclass(frozen=True)
class InFlightMolecule:
    """A molecule flagged ``running``; ``started_at`` is when it was dispatched."""

    mol_id: str
    smiles: str
    nheavy: int
    started_at: str | None


@dataclass(frozen=True)
class CompletedMolecule:
    """A molecule accepted under protocol 2, as the progress feed shows it."""

    mol_id: str
    smiles: str
    h298_pm7: float
    delta_kcal_mol: float | None
    elapsed_s: float | None
    conformers_evaluated: int | None
    finished_at: str


@dataclass(frozen=True)
class TrainingActivity:
    """What is computing now and what finished last; never written by readers."""

    in_flight: tuple[InFlightMolecule, ...]
    recent: tuple[CompletedMolecule, ...]
    #: Mean worker time of the molecules accepted since ``since`` (``None``
    #: when none finished yet), so a fresh run does not inherit old timings.
    mean_elapsed_s: float | None


def read_activity(
    db_path: Path, *, since: str | None = None, recent: int = 8
) -> TrainingActivity:
    """In-flight molecules, latest successes and mean time since ``since``.

    ``updated_at`` is an ISO-8601 UTC string written by :func:`_now`, so the
    ``>=`` comparison is chronological.
    """
    if not db_path.exists():
        return TrainingActivity((), (), None)
    with closing(connect(db_path)) as connection:
        running = connection.execute(
            "SELECT mol_id, smiles, nheavy, updated_at FROM training_queue"
            " WHERE status = ? ORDER BY updated_at",
            (RUNNING,),
        ).fetchall()
        done = connection.execute(
            "SELECT mol_id, smiles, h298_pm7, delta_kcal_mol, elapsed_s,"
            " conformers_evaluated, updated_at FROM training_queue"
            " WHERE status = ? ORDER BY updated_at DESC LIMIT ?",
            (OK, recent),
        ).fetchall()
        mean = connection.execute(
            "SELECT AVG(elapsed_s) AS mean FROM training_queue"
            " WHERE status = ? AND elapsed_s IS NOT NULL AND updated_at >= ?",
            (OK, since or ""),
        ).fetchone()
    return TrainingActivity(
        in_flight=tuple(
            InFlightMolecule(
                mol_id=str(r["mol_id"]),
                smiles=str(r["smiles"]),
                nheavy=int(r["nheavy"]),
                started_at=r["updated_at"],
            )
            for r in running
        ),
        recent=tuple(
            CompletedMolecule(
                mol_id=str(r["mol_id"]),
                smiles=str(r["smiles"]),
                h298_pm7=float(r["h298_pm7"]),
                delta_kcal_mol=r["delta_kcal_mol"],
                elapsed_s=r["elapsed_s"],
                conformers_evaluated=r["conformers_evaluated"],
                finished_at=str(r["updated_at"]),
            )
            for r in done
        ),
        mean_elapsed_s=None if mean["mean"] is None else float(mean["mean"]),
    )


class TrainingQueue:
    """Seeding, queue order and state transitions of the training molecules."""

    def __init__(
        self, db_path: Path | str, cbs_csv: Path | str, *, max_reruns: int = 3
    ) -> None:
        self._db_path = Path(db_path)
        self._cbs_csv = Path(cbs_csv)
        self._max_reruns = max_reruns

    def ensure_seeded(self) -> int:
        """Insert the CBS molecules that are not in the queue yet.

        Existing rows are never touched, so a rerun of the routine keeps every
        status. Returns how many rows were inserted.
        """
        selected = select_reference_rows(self._cbs_csv)
        records = [
            (
                str(r.mol_id),
                str(r.smiles),
                int(r.charge),
                int(r.multiplicity),
                int(r.nheavy),
                float(r.H298_cbs),
            )
            for r in selected.itertuples(index=False)
        ]
        with closing(connect(self._db_path)) as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                before = connection.total_changes
                connection.executemany(
                    "INSERT OR IGNORE INTO training_queue"
                    " (mol_id, smiles, charge, multiplicity, nheavy, h298_cbs)"
                    " VALUES (?, ?, ?, ?, ?, ?)",
                    records,
                )
                inserted = connection.total_changes - before
                connection.execute("COMMIT")
            except BaseException:
                connection.execute("ROLLBACK")
                raise
        return inserted

    def queue(self, limit: int | None = None) -> list[QueuedMolecule]:
        """Queued molecules, smallest first (ties by ``mol_id``)."""
        with closing(connect(self._db_path)) as connection:
            rows = connection.execute(
                "SELECT mol_id, smiles, charge, multiplicity, nheavy, h298_cbs"
                " FROM training_queue WHERE status IN ('pending', 'rerun', 'running')"
                " ORDER BY nheavy, mol_id LIMIT ?",
                (-1 if limit is None else limit,),  # SQLite: negative = no limit
            ).fetchall()
        return [
            QueuedMolecule(
                mol_id=r["mol_id"],
                smiles=r["smiles"],
                charge=r["charge"],
                multiplicity=r["multiplicity"],
                nheavy=r["nheavy"],
                h298_cbs=r["h298_cbs"],
            )
            for r in rows
        ]

    def claim(
        self, cell_id: str, count: int, *, lease_s: float = DEFAULT_LEASE_S
    ) -> list[ClaimedMolecule]:
        """Reserve up to ``count`` queued molecules to ``cell_id``, smallest first.

        One ``UPDATE … RETURNING`` under the write lock, so two claimers (the
        central routine and a cell API request) never get the same molecule.
        Each reserved row gets a fresh ``attempt_id``, also logged in
        ``training_attempt``.
        """
        if not cell_id:
            raise ValueError("cell_id must not be empty")
        if count < 0:
            raise ValueError(f"count must be >= 0, got {count}")
        if count == 0:
            return []
        now = datetime.now(timezone.utc)
        lease_until = (now + timedelta(seconds=lease_s)).isoformat()
        with self._write() as connection:
            rows = connection.execute(
                "UPDATE training_queue SET status = 'running', cell_id = ?,"
                " attempt_id = lower(hex(randomblob(16))), lease_until = ?,"
                " updated_at = ? WHERE mol_id IN ("
                "  SELECT mol_id FROM training_queue"
                "  WHERE status IN ('pending', 'rerun')"
                "  ORDER BY nheavy, mol_id LIMIT ?)"
                " RETURNING mol_id, smiles, charge, multiplicity, nheavy, h298_cbs,"
                " attempt_id",
                (cell_id, lease_until, now.isoformat(), count),
            ).fetchall()
            connection.executemany(
                "INSERT INTO training_attempt (attempt_id, mol_id, cell_id, claimed_at)"
                " VALUES (?, ?, ?, ?)",
                [
                    (r["attempt_id"], r["mol_id"], cell_id, now.isoformat())
                    for r in rows
                ],
            )
        claimed = [
            ClaimedMolecule(
                mol_id=r["mol_id"],
                smiles=r["smiles"],
                charge=r["charge"],
                multiplicity=r["multiplicity"],
                nheavy=r["nheavy"],
                h298_cbs=r["h298_cbs"],
                attempt_id=r["attempt_id"],
            )
            for r in rows
        ]
        # RETURNING has no defined order; keep the queue's.
        return sorted(claimed, key=lambda m: (m.nheavy, m.mol_id))

    def mark_running(self, mol_id: str) -> None:
        """Flag ``mol_id`` as in flight (persisted immediately)."""
        self._expect_one(
            "UPDATE training_queue SET status = 'running', updated_at = ?"
            " WHERE mol_id = ?",
            (_now(), mol_id),
            mol_id,
        )

    def record_success(
        self,
        mol_id: str,
        *,
        hof_kcal_mol: float,
        conformers_evaluated: int,
        elapsed_s: float,
        batch_id: str,
    ) -> None:
        """Store a protocol-2 result: PM7 value and delta to CBS (CBS − PM7)."""
        if not math.isfinite(hof_kcal_mol):
            raise ValueError(f"non-finite PM7 value for {mol_id}: {hof_kcal_mol}")
        self._expect_one(
            "UPDATE training_queue SET status = 'ok', h298_pm7 = ?,"
            " delta_kcal_mol = h298_cbs - ?, conformers_evaluated = ?,"
            " elapsed_s = ?, batch_id = ?, updated_at = ? WHERE mol_id = ?",
            (
                hof_kcal_mol,
                hof_kcal_mol,
                conformers_evaluated,
                elapsed_s,
                batch_id,
                _now(),
                mol_id,
            ),
            mol_id,
        )

    def record_failure(self, mol_id: str, reason: str) -> str:
        """Count a failed attempt: ``rerun`` below the limit, else ``skip``.

        The reason is persisted by ``JobState.record_failure``.
        """
        del reason
        with closing(connect(self._db_path)) as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                found = connection.execute(
                    "SELECT reruns FROM training_queue WHERE mol_id = ?", (mol_id,)
                ).fetchone()
                if found is None:
                    raise KeyError(mol_id)
                reruns = int(found["reruns"]) + 1
                status = SKIP if reruns >= self._max_reruns else RERUN
                connection.execute(
                    "UPDATE training_queue SET status = ?, reruns = ?, updated_at = ?"
                    " WHERE mol_id = ?",
                    (status, reruns, _now(), mol_id),
                )
                connection.execute("COMMIT")
            except BaseException:
                connection.execute("ROLLBACK")
                raise
        return status

    @contextmanager
    def _write(self) -> Iterator[sqlite3.Connection]:
        """One write transaction, taken up front so readers never see half of it."""
        with closing(connect(self._db_path)) as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                yield connection
                connection.execute("COMMIT")
            except BaseException:
                connection.execute("ROLLBACK")
                raise

    def _expect_one(self, sql: str, params: tuple[object, ...], mol_id: str) -> None:
        with closing(connect(self._db_path)) as connection:
            changed = connection.execute(sql, params).rowcount
        if changed != 1:
            raise KeyError(mol_id)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()
