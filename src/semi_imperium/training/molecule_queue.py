"""The training queue: which molecules to compute and how each one stands.

The queue is seeded from ``cbs_referencia.csv`` (read only) and lives in the
``training_queue`` table of the training database, so the CSV never changes.
One molecule per SMILES: among the rows that share a SMILES the reference
conformer wins (lowest ``conformer_rank_by_h298``, which is the lowest H298);
ties fall back to the lowest ``H298_cbs`` and then to ``mol_id``.

Every consumer (the local routine as cell ``central``, later the remote
cells through the central API) takes molecules with ``claim`` and reports
them by ``attempt_id``; workers never touch the database.
"""

from __future__ import annotations

import json
import math
import sqlite3
from collections.abc import Collection, Iterator
from contextlib import closing, contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

from semi_imperium.training.db import connect
from semi_imperium.training.geometry_store import GeometryRecord, insert_geometry

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

#: Rows still to compute. ``running`` ones are reserved to a cell and only
#: return to the queue through ``TrainingQueue.release``.
QUEUE_STATUSES = (PENDING, RERUN, RUNNING)


#: The cell that runs on the central machine itself (the local routine).
CENTRAL_CELL = "central"

#: How far ahead a claim or renewal pushes ``lease_until``. Informative only:
#: nothing returns a molecule when its lease lapses (decision of 2026-10-20);
#: the dashboard uses it to tell a quiet cell from a working one.
DEFAULT_LEASE_S = 300.0


#: Outcomes of reporting a result for an attempt (besides ``rerun``/``skip``).
ACCEPTED = "accepted"
DUPLICATE = "duplicate"
DISCARDED = "discarded"
_FINAL_OUTCOMES = (ACCEPTED, "failed", DISCARDED)


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
    """A molecule flagged ``running``; ``started_at`` is when it was claimed."""

    mol_id: str
    smiles: str
    nheavy: int
    started_at: str | None
    #: Which cell holds it (``None`` only for rows from before the migration).
    cell_id: str | None


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
            "SELECT mol_id, smiles, nheavy, updated_at, cell_id FROM training_queue"
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
                cell_id=r["cell_id"],
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

    def claim(
        self,
        cell_id: str,
        count: int,
        *,
        lease_s: float = DEFAULT_LEASE_S,
        exclude: Collection[str] = (),
    ) -> list[ClaimedMolecule]:
        """Reserve up to ``count`` queued molecules to ``cell_id``, smallest first.

        One ``UPDATE … RETURNING`` under the write lock, so two claimers (the
        central routine and a cell API request) never get the same molecule.
        Each reserved row gets a fresh ``attempt_id``, also logged in
        ``training_attempt``. ``exclude`` keeps molecules out of this claim:
        the central run passes what it already tried, so a molecule that just
        failed (``rerun``, still smallest first) waits for the next run
        instead of being retried at once.
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
                "  AND mol_id NOT IN (SELECT value FROM json_each(?))"
                "  ORDER BY nheavy, mol_id LIMIT ?)"
                " RETURNING mol_id, smiles, charge, multiplicity, nheavy, h298_cbs,"
                " attempt_id",
                (
                    cell_id,
                    lease_until,
                    now.isoformat(),
                    json.dumps(sorted(exclude)),
                    count,
                ),
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

    def renew(
        self,
        cell_id: str,
        attempt_ids: list[str],
        *,
        lease_s: float = DEFAULT_LEASE_S,
    ) -> int:
        """Push ``lease_until`` of ``cell_id``'s running attempts; count them.

        Attempts of other cells, finished or released ones are ignored, so a
        late heartbeat cannot take a molecule back. ``updated_at`` is left
        alone: readers show it as the time the molecule was claimed.
        """
        if not attempt_ids:
            return 0
        lease_until = (
            datetime.now(timezone.utc) + timedelta(seconds=lease_s)
        ).isoformat()
        with self._write() as connection:
            return connection.execute(
                "UPDATE training_queue SET lease_until = ?"
                " WHERE status = 'running' AND cell_id = ?"
                " AND attempt_id IN (SELECT value FROM json_each(?))",
                (lease_until, cell_id, json.dumps(attempt_ids)),
            ).rowcount

    def release(self, cell_id: str, *, include_unowned: bool = False) -> list[str]:
        """Give ``cell_id``'s running molecules back to the queue; return them.

        This is the only way a reservation ends without a result (no lease
        expiry). No rerun is spent: a row returns to ``rerun`` if it had
        failed before, else to ``pending``. ``include_unowned`` also takes
        ``running`` rows with no cell, left by runs from before the cell
        columns existed. A result that still arrives for a released attempt
        is handled by the acceptance rules, not refused here.
        """
        now = _now()
        params = (cell_id, int(include_unowned))
        with self._write() as connection:
            # Attempts first: RETURNING would only show the cleared attempt_id.
            connection.execute(
                "UPDATE training_attempt SET outcome = 'released', released_at = ?"
                " WHERE outcome IS NULL AND attempt_id IN ("
                "  SELECT attempt_id FROM training_queue WHERE status = 'running'"
                "  AND (cell_id = ? OR (? AND cell_id IS NULL)))",
                (now, *params),
            )
            rows = connection.execute(
                "UPDATE training_queue SET"
                " status = CASE WHEN reruns > 0 THEN 'rerun' ELSE 'pending' END,"
                " cell_id = NULL, attempt_id = NULL, lease_until = NULL,"
                " updated_at = ? WHERE status = 'running'"
                " AND (cell_id = ? OR (? AND cell_id IS NULL)) RETURNING mol_id",
                (now, *params),
            ).fetchall()
        return sorted(str(r["mol_id"]) for r in rows)

    def accept_success(
        self,
        attempt_id: str,
        *,
        geometry: GeometryRecord,
        conformers_evaluated: int,
        elapsed_s: float,
        batch_id: str,
    ) -> str:
        """Record a protocol-2 result of ``attempt_id``; first valid result wins.

        ``accepted``: stored with its geometry in one transaction, also when
        the attempt had been released, as long as the molecule is not ``ok``.
        ``discarded``: another attempt already completed it; logged in
        ``training_attempt``. ``duplicate``: this attempt was already settled,
        nothing changes (a cell resending after a network drop).
        """
        hof = geometry.hof_kcal_mol
        if hof is None or not math.isfinite(hof):
            raise ValueError(f"non-finite PM7 value for {geometry.mol_id}: {hof}")
        with self._write() as connection:
            attempt = self._attempt(connection, attempt_id)
            if attempt["outcome"] in _FINAL_OUTCOMES:
                return DUPLICATE
            mol_id = str(attempt["mol_id"])
            if geometry.mol_id != mol_id:
                raise ValueError(
                    f"attempt {attempt_id} is for mol_id {mol_id},"
                    f" result is for {geometry.mol_id}"
                )
            current = connection.execute(
                "SELECT status, attempt_id FROM training_queue WHERE mol_id = ?",
                (mol_id,),
            ).fetchone()
            if current["status"] == OK:
                self._settle(
                    connection,
                    attempt_id,
                    DISCARDED,
                    f"já concluída pela tentativa {current['attempt_id']}",
                )
                return DISCARDED
            insert_geometry(connection, geometry)
            connection.execute(
                "UPDATE training_queue SET status = 'ok', h298_pm7 = ?,"
                " delta_kcal_mol = h298_cbs - ?, conformers_evaluated = ?,"
                " elapsed_s = ?, batch_id = ?, cell_id = ?, attempt_id = ?,"
                " lease_until = NULL, updated_at = ? WHERE mol_id = ?",
                (
                    hof,
                    hof,
                    conformers_evaluated,
                    elapsed_s,
                    batch_id,
                    attempt["cell_id"],
                    attempt_id,
                    _now(),
                    mol_id,
                ),
            )
            self._settle(connection, attempt_id, ACCEPTED, None)
        return ACCEPTED

    def accept_failure(self, attempt_id: str, reason: str) -> str:
        """Record a failed ``attempt_id``: ``rerun``/``skip``, or no effect.

        Only the attempt that currently holds the molecule spends a rerun. A
        failure of a released or superseded attempt is ``discarded`` (logged,
        the molecule untouched); a settled attempt gives ``duplicate``.
        """
        with self._write() as connection:
            attempt = self._attempt(connection, attempt_id)
            if attempt["outcome"] in _FINAL_OUTCOMES:
                return DUPLICATE
            current = connection.execute(
                "SELECT status, attempt_id, reruns FROM training_queue"
                " WHERE mol_id = ?",
                (attempt["mol_id"],),
            ).fetchone()
            if current["status"] != RUNNING or current["attempt_id"] != attempt_id:
                self._settle(
                    connection, attempt_id, DISCARDED, f"falha tardia: {reason}"
                )
                return DISCARDED
            reruns = int(current["reruns"]) + 1
            status = SKIP if reruns >= self._max_reruns else RERUN
            connection.execute(
                "UPDATE training_queue SET status = ?, reruns = ?, cell_id = NULL,"
                " attempt_id = NULL, lease_until = NULL, updated_at = ?"
                " WHERE mol_id = ?",
                (status, reruns, _now(), attempt["mol_id"]),
            )
            self._settle(connection, attempt_id, "failed", reason)
        return status

    @staticmethod
    def _attempt(connection: sqlite3.Connection, attempt_id: str) -> sqlite3.Row:
        found: sqlite3.Row | None = connection.execute(
            "SELECT mol_id, cell_id, outcome FROM training_attempt"
            " WHERE attempt_id = ?",
            (attempt_id,),
        ).fetchone()
        if found is None:
            raise KeyError(attempt_id)
        return found

    @staticmethod
    def _settle(
        connection: sqlite3.Connection,
        attempt_id: str,
        outcome: str,
        note: str | None,
    ) -> None:
        connection.execute(
            "UPDATE training_attempt SET outcome = ?, note = ?, finished_at = ?"
            " WHERE attempt_id = ?",
            (outcome, note, _now(), attempt_id),
        )

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


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()
