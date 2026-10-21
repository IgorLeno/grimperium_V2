"""The cells allowed to compute for the central, one row per machine.

A cell authenticates with a bearer token generated here and handed over once
(later inside the cell package). Only the token's sha256 is stored: the token
carries 256 random bits, so a fast hash is enough and a lookup by hash never
compares secrets byte by byte. Revoking clears the hash, so the old token
stops working at once while the row (and its history) stays.

Every time is the central's clock; a cell's own clock is never trusted.
"""

from __future__ import annotations

import hashlib
import json
import secrets
import sqlite3
from collections.abc import Iterator, Mapping
from contextlib import AbstractContextManager, closing, contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from semi_imperium.training.db import connect

ADMISSION_PENDING = "pending"
ADMISSION_PASSED = "passed"
ADMISSION_BLOCKED = "blocked"
ADMISSIONS = (ADMISSION_PENDING, ADMISSION_PASSED, ADMISSION_BLOCKED)

#: Heartbeat age limits of the dashboard states (decision of 2026-10-20).
ONLINE_S = 90.0
LATE_S = 300.0

#: The admission test (decision of 2026-10-20): a new cell computes these
#: five molecules, already accepted at the central with the pinned programs,
#: and each HOF must match the central's within the tolerance.
ADMISSION_REFERENCES = (
    "cbs_09698",  # CO
    "cbs_00187",  # [H][C]=C
    "cbs_24671",  # C[N+](=O)[O-]
    "cbs_08982",  # OCCO
    "cbs_01400",  # CCC(=O)O
)
ADMISSION_TOLERANCE_KCAL_MOL = 0.1
#: Admission attempts travel through the same claim/result routes as the
#: queue's; this prefix tells them apart (queue attempts are bare hex).
ADMISSION_PREFIX = "adm-"

# ``cell_admission_attempt.outcome``: NULL while open, then one of these.
_RELEASED = "released"  # the cell restarted; computed again unless it arrives
_WITHIN = "within"
_OUTSIDE = "outside"
_FAILED = "failed"
_DUPLICATE = "duplicate"  # another attempt of the same molecule came first
_SUPERSEDED = "superseded"  # the test ended or started over without it
_DONE = (_WITHIN, _OUTSIDE, _FAILED)


@dataclass(frozen=True)
class Cell:
    """One cell as the central knows it (never its token)."""

    cell_id: str
    name: str
    host: str | None
    cores: int | None
    workers: int | None
    crest_threads: int | None
    program_versions: dict[str, str]
    protocol_signature: str | None
    admission: str
    created_at: str
    #: Last accepted registration; ``None`` until one passes the version check.
    registered_at: str | None
    last_seen_at: str | None
    revoked_at: str | None
    #: Why the last registration was refused (programs differ from the central).
    refusal: str | None


def cell_state(cell: Cell, now: datetime) -> str:
    """Dashboard state of ``cell`` at ``now`` (the central's clock).

    Refusals outrank connectivity: a revoked, mismatched or blocked cell
    cannot receive work however recent its last signal.
    """
    if cell.revoked_at is not None:
        return "revoked"
    if cell.refusal is not None:
        return "version_mismatch"
    if cell.admission == ADMISSION_BLOCKED:
        return "blocked"
    if cell.registered_at is None or cell.last_seen_at is None:
        return "never_seen"
    age = (now - datetime.fromisoformat(cell.last_seen_at)).total_seconds()
    if age < ONLINE_S:
        return "online"
    return "late" if age <= LATE_S else "silent"


def _hash(token: str) -> str:
    return hashlib.sha256(token.encode()).hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class CellRegistry:
    """Create, authenticate and update the rows of the ``cells`` table."""

    def __init__(self, db_path: Path | str) -> None:
        self._db_path = Path(db_path)

    def create(self, name: str) -> tuple[Cell, str]:
        """Add a cell named ``name``; return it with its token (shown only now)."""
        name = _clean_name(name)
        token = secrets.token_urlsafe(32)
        cell_id = f"cell-{secrets.token_hex(4)}"
        with self._write() as connection:
            connection.execute(
                "INSERT INTO cells (cell_id, name, token_hash, created_at)"
                " VALUES (?, ?, ?, ?)",
                (cell_id, name, _hash(token), _now()),
            )
            return self._get(connection, cell_id), token

    def authenticate(self, token: str) -> Cell | None:
        """The cell holding ``token``; ``None`` if unknown or revoked."""
        if not token:
            return None
        with closing(connect(self._db_path)) as connection:
            found = connection.execute(
                "SELECT * FROM cells WHERE token_hash = ?",
                (_hash(token),),
            ).fetchone()
        return None if found is None else _cell(found)

    def record_registration(
        self,
        cell_id: str,
        *,
        host: str,
        cores: int,
        workers: int,
        crest_threads: int,
        program_versions: Mapping[str, str],
        protocol_signature: str,
        refusal: str | None,
    ) -> Cell:
        """Store what the cell reported; ``refusal`` says why it was turned away.

        A refused registration keeps the reported versions (the dashboard
        shows what differs) but leaves ``registered_at`` alone; an accepted
        one clears any earlier refusal.
        """
        now = _now()
        with self._write() as connection:
            self._require(connection, cell_id)
            connection.execute(
                "UPDATE cells SET host = ?, cores = ?, workers = ?,"
                " crest_threads = ?, program_versions = ?, protocol_signature = ?,"
                " refusal = ?, last_seen_at = ?,"
                " registered_at = CASE WHEN ? IS NULL THEN ? ELSE registered_at END"
                " WHERE cell_id = ?",
                (
                    host,
                    cores,
                    workers,
                    crest_threads,
                    json.dumps(dict(program_versions), sort_keys=True),
                    protocol_signature,
                    refusal,
                    now,
                    refusal,
                    now,
                    cell_id,
                ),
            )
            return self._get(connection, cell_id)

    def touch(self, cell_id: str) -> None:
        """Record a signal from ``cell_id`` now."""
        with self._write() as connection:
            self._require(connection, cell_id)
            connection.execute(
                "UPDATE cells SET last_seen_at = ? WHERE cell_id = ?",
                (_now(), cell_id),
            )

    def rename(self, cell_id: str, name: str) -> Cell:
        name = _clean_name(name)
        return self._update(
            cell_id, "UPDATE cells SET name = ? WHERE cell_id = ?", (name,)
        )

    def revoke(self, cell_id: str) -> Cell:
        """Invalidate the cell's token; its row and history stay."""
        return self._update(
            cell_id,
            "UPDATE cells SET token_hash = NULL,"
            " revoked_at = COALESCE(revoked_at, ?) WHERE cell_id = ?",
            (_now(),),
        )

    def rotate_token(self, cell_id: str) -> tuple[Cell, str]:
        """Give ``cell_id`` a new token (shown only now); the old one stops working.

        The cell package carries the token, and the central keeps only its
        hash, so every new package needs a new token. A revoked cell stays
        revoked: ``ValueError``.
        """
        token = secrets.token_urlsafe(32)
        with self._write() as connection:
            if self._get(connection, cell_id).revoked_at is not None:
                raise ValueError(f"cell {cell_id} is revoked")
            connection.execute(
                "UPDATE cells SET token_hash = ? WHERE cell_id = ?",
                (_hash(token), cell_id),
            )
            return self._get(connection, cell_id), token

    def set_admission(self, cell_id: str, admission: str) -> Cell:
        """Decide the admission by hand.

        Back to ``pending`` starts a new admission test: every earlier
        reference result is set aside. Any other decision only ends the
        test's open attempts, so the results that led to it stay visible.
        """
        if admission not in ADMISSIONS:
            raise ValueError(f"admission must be one of {ADMISSIONS}, got {admission}")
        with self._write() as connection:
            self._require(connection, cell_id)
            connection.execute(
                "UPDATE cells SET admission = ? WHERE cell_id = ?",
                (admission, cell_id),
            )
            _supersede(connection, cell_id, everything=admission == ADMISSION_PENDING)
            return self._get(connection, cell_id)

    def get(self, cell_id: str) -> Cell:
        with closing(connect(self._db_path)) as connection:
            return self._get(connection, cell_id)

    def list(self) -> list[Cell]:
        """Every cell, oldest first (the dashboard's row order)."""
        if not self._db_path.exists():
            return []
        with closing(connect(self._db_path)) as connection:
            rows = connection.execute(
                "SELECT * FROM cells ORDER BY created_at, cell_id"
            ).fetchall()
        return [_cell(r) for r in rows]

    def _update(self, cell_id: str, statement: str, params: tuple[str, ...]) -> Cell:
        """Run ``statement`` (its last parameter is ``cell_id``) on an existing cell."""
        with self._write() as connection:
            self._require(connection, cell_id)
            connection.execute(statement, (*params, cell_id))
            return self._get(connection, cell_id)

    @staticmethod
    def _require(connection: sqlite3.Connection, cell_id: str) -> None:
        if (
            connection.execute(
                "SELECT 1 FROM cells WHERE cell_id = ?", (cell_id,)
            ).fetchone()
            is None
        ):
            raise KeyError(cell_id)

    @staticmethod
    def _get(connection: sqlite3.Connection, cell_id: str) -> Cell:
        found = connection.execute(
            "SELECT * FROM cells WHERE cell_id = ?", (cell_id,)
        ).fetchone()
        if found is None:
            raise KeyError(cell_id)
        return _cell(found)

    def _write(self) -> AbstractContextManager[sqlite3.Connection]:
        return _write(self._db_path)


class AdmissionUnavailable(RuntimeError):
    """The central has no accepted HOF for some reference molecule."""


@dataclass(frozen=True)
class AdmissionMolecule:
    """One reference molecule handed to a cell under admission."""

    attempt_id: str
    mol_id: str
    smiles: str
    charge: int
    multiplicity: int
    nheavy: int


class AdmissionTest:
    """The five reference molecules a ``pending`` cell computes before work.

    The references never enter the queue: the attempts live in
    ``cell_admission_attempt`` and a result only becomes a verdict, never a
    geometry. Each HOF is compared with the central's ``final_geometry``
    when it arrives. One result outside the tolerance (or one failed
    calculation) blocks the cell; the five inside admit it.
    """

    def __init__(self, db_path: Path | str) -> None:
        self._db_path = Path(db_path)

    def claim(self, cell_id: str, count: int) -> list[AdmissionMolecule]:
        """Up to ``count`` references ``cell_id`` has neither open nor answered.

        ``AdmissionUnavailable`` if the central lacks a reference: without
        the accepted value there is nothing to compare with.
        """
        if count <= 0:
            return []
        now = _now()
        with _write(self._db_path) as connection:
            references = _references(connection)
            taken = {
                str(r["mol_id"])
                for r in connection.execute(
                    "SELECT mol_id FROM cell_admission_attempt WHERE cell_id = ?"
                    " AND (outcome IS NULL OR outcome IN (?, ?, ?))",
                    (cell_id, *_DONE),
                )
            }
            handed = [
                AdmissionMolecule(
                    attempt_id=f"{ADMISSION_PREFIX}{secrets.token_hex(16)}",
                    mol_id=mol_id,
                    smiles=str(references[mol_id]["smiles"]),
                    charge=int(references[mol_id]["charge"]),
                    multiplicity=int(references[mol_id]["multiplicity"]),
                    nheavy=int(references[mol_id]["nheavy"]),
                )
                for mol_id in ADMISSION_REFERENCES
                if mol_id not in taken
            ][:count]
            connection.executemany(
                "INSERT INTO cell_admission_attempt"
                " (attempt_id, cell_id, mol_id, claimed_at) VALUES (?, ?, ?, ?)",
                [(m.attempt_id, cell_id, m.mol_id, now) for m in handed],
            )
        return handed

    def held(self, cell_id: str) -> set[str]:
        """Open admission attempts of ``cell_id``."""
        if not self._db_path.exists():
            return set()
        with closing(connect(self._db_path)) as connection:
            rows = connection.execute(
                "SELECT attempt_id FROM cell_admission_attempt"
                " WHERE cell_id = ? AND outcome IS NULL",
                (cell_id,),
            ).fetchall()
        return {str(r["attempt_id"]) for r in rows}

    def release(self, cell_id: str) -> None:
        """The cell restarted: its open references will be handed out again."""
        with _write(self._db_path) as connection:
            connection.execute(
                "UPDATE cell_admission_attempt SET outcome = ?, released_at = ?"
                " WHERE cell_id = ? AND outcome IS NULL",
                (_RELEASED, _now(), cell_id),
            )

    def owner(self, attempt_id: str) -> tuple[str, str]:
        """``(mol_id, cell_id)`` of ``attempt_id``; ``KeyError`` if unknown."""
        with closing(connect(self._db_path)) as connection:
            found = connection.execute(
                "SELECT mol_id, cell_id FROM cell_admission_attempt"
                " WHERE attempt_id = ?",
                (attempt_id,),
            ).fetchone()
        if found is None:
            raise KeyError(attempt_id)
        return str(found["mol_id"]), str(found["cell_id"])

    def record(
        self, attempt_id: str, *, hof_kcal_mol: float | None, failure: str | None
    ) -> str:
        """Judge one reference result; return the outcome the cell is told.

        ``accepted`` when it counted, ``duplicate`` when this attempt (or
        another one of the same molecule) was already judged, ``discarded``
        when the test is over or started again. A released attempt still
        counts: the cell did compute it.
        """
        from semi_imperium.training.molecule_queue import (
            ACCEPTED,
            DISCARDED,
            DUPLICATE,
        )

        if (hof_kcal_mol is None) == (failure is None):
            raise ValueError("pass hof_kcal_mol or failure, exactly one")
        now = _now()
        with _write(self._db_path) as connection:
            attempt = connection.execute(
                "SELECT * FROM cell_admission_attempt WHERE attempt_id = ?",
                (attempt_id,),
            ).fetchone()
            if attempt is None:
                raise KeyError(attempt_id)
            cell_id, mol_id = str(attempt["cell_id"]), str(attempt["mol_id"])
            if attempt["outcome"] in (*_DONE, _DUPLICATE):
                return DUPLICATE
            admission = connection.execute(
                "SELECT admission FROM cells WHERE cell_id = ?", (cell_id,)
            ).fetchone()
            if attempt["outcome"] == _SUPERSEDED or (
                admission is None or admission["admission"] != ADMISSION_PENDING
            ):
                return DISCARDED
            if connection.execute(
                "SELECT 1 FROM cell_admission_attempt WHERE cell_id = ?"
                " AND mol_id = ? AND outcome IN (?, ?, ?)",
                (cell_id, mol_id, *_DONE),
            ).fetchone():
                connection.execute(
                    "UPDATE cell_admission_attempt SET outcome = ?, finished_at = ?"
                    " WHERE attempt_id = ?",
                    (_DUPLICATE, now, attempt_id),
                )
                return DUPLICATE
            reference = float(_references(connection)[mol_id]["hof_kcal_mol"])
            if failure is not None or hof_kcal_mol is None:
                outcome, note = _FAILED, failure
            # MOPAC reports 5 decimals: rounding drops the float noise, so a
            # difference of exactly the tolerance is inside it.
            elif round(abs(hof_kcal_mol - reference), 6) <= (
                ADMISSION_TOLERANCE_KCAL_MOL
            ):
                outcome, note = _WITHIN, None
            else:
                outcome, note = _OUTSIDE, None
            connection.execute(
                "UPDATE cell_admission_attempt SET outcome = ?, finished_at = ?,"
                " hof_kcal_mol = ?, reference_hof_kcal_mol = ?, note = ?"
                " WHERE attempt_id = ?",
                (outcome, now, hof_kcal_mol, reference, note, attempt_id),
            )
            within = connection.execute(
                "SELECT COUNT(DISTINCT mol_id) FROM cell_admission_attempt"
                " WHERE cell_id = ? AND outcome = ?",
                (cell_id, _WITHIN),
            ).fetchone()[0]
            verdict = None
            if outcome != _WITHIN:
                verdict = ADMISSION_BLOCKED
            elif within == len(ADMISSION_REFERENCES):
                verdict = ADMISSION_PASSED
            if verdict is not None:
                connection.execute(
                    "UPDATE cells SET admission = ? WHERE cell_id = ?",
                    (verdict, cell_id),
                )
                _supersede(connection, cell_id, everything=False)
        return ACCEPTED

    def note(self, cell_id: str, admission: str) -> str | None:
        """What the dashboard says about the test: progress, or why it blocked."""
        if admission == ADMISSION_PASSED or not self._db_path.exists():
            return None
        with closing(connect(self._db_path)) as connection:
            rows = connection.execute(
                "SELECT a.mol_id, a.outcome, a.hof_kcal_mol, a.reference_hof_kcal_mol,"
                " a.note, g.smiles FROM cell_admission_attempt a"
                " LEFT JOIN final_geometry g USING (mol_id)"
                " WHERE a.cell_id = ? AND a.outcome IN (?, ?, ?)"
                " ORDER BY a.finished_at, a.mol_id",
                (cell_id, *_DONE),
            ).fetchall()
        if admission == ADMISSION_PENDING:
            if not rows:
                return None
            return (
                f"Teste de admissão: {len(rows)} de {len(ADMISSION_REFERENCES)}"
                f" referências calculadas, todas dentro de"
                f" {ADMISSION_TOLERANCE_KCAL_MOL:g} kcal/mol."
            )
        reasons = [_reason(r) for r in rows if r["outcome"] != _WITHIN]
        return "; ".join(reasons) or None


def _reason(row: sqlite3.Row) -> str:
    name = f"{row['mol_id']} ({row['smiles']})" if row["smiles"] else row["mol_id"]
    if row["outcome"] == _FAILED:
        return f"{name} falhou: {row['note']}"
    hof, reference = float(row["hof_kcal_mol"]), float(row["reference_hof_kcal_mol"])
    return (
        f"{name}: HOF {hof:.5f} contra {reference:.5f} da central"
        f" (Δ {abs(hof - reference):.5f} > {ADMISSION_TOLERANCE_KCAL_MOL:g}"
        " kcal/mol)"
    )


def _references(connection: sqlite3.Connection) -> dict[str, sqlite3.Row]:
    """The reference molecules with the central's accepted HOF, by ``mol_id``."""
    rows = connection.execute(
        "SELECT g.mol_id, g.smiles, g.hof_kcal_mol, q.charge, q.multiplicity,"
        " q.nheavy FROM final_geometry g JOIN training_queue q USING (mol_id)"
        " WHERE g.hof_kcal_mol IS NOT NULL"
        " AND g.mol_id IN (SELECT value FROM json_each(?))",
        (json.dumps(ADMISSION_REFERENCES),),
    ).fetchall()
    found = {str(r["mol_id"]): r for r in rows}
    missing = [m for m in ADMISSION_REFERENCES if m not in found]
    if missing:
        raise AdmissionUnavailable(
            "a central não tem HOF aceito para " + ", ".join(missing)
        )
    return found


def _supersede(
    connection: sqlite3.Connection, cell_id: str, *, everything: bool
) -> None:
    """End the open admission attempts of ``cell_id`` (all of them if ``everything``).

    A result that still arrives for one is discarded.
    """
    statement = (
        "UPDATE cell_admission_attempt SET outcome = ?, released_at ="
        " COALESCE(released_at, ?) WHERE cell_id = ?"
    )
    if everything:
        statement += " AND outcome IS NOT ?"
    else:
        statement += " AND (outcome IS NULL OR outcome = ?)"
    connection.execute(
        statement,
        (_SUPERSEDED, _now(), cell_id, _SUPERSEDED if everything else _RELEASED),
    )


@contextmanager
def _write(db_path: Path) -> Iterator[sqlite3.Connection]:
    with closing(connect(db_path)) as connection:
        connection.execute("BEGIN IMMEDIATE")
        try:
            yield connection
            connection.execute("COMMIT")
        except BaseException:
            connection.execute("ROLLBACK")
            raise


def _clean_name(name: str) -> str:
    cleaned = name.strip()
    if not cleaned:
        raise ValueError("cell name must not be empty")
    return cleaned


def _cell(row: sqlite3.Row) -> Cell:
    return Cell(
        cell_id=str(row["cell_id"]),
        name=str(row["name"]),
        host=row["host"],
        cores=row["cores"],
        workers=row["workers"],
        crest_threads=row["crest_threads"],
        program_versions=json.loads(row["program_versions"] or "{}"),
        protocol_signature=row["protocol_signature"],
        admission=str(row["admission"]),
        created_at=str(row["created_at"]),
        registered_at=row["registered_at"],
        last_seen_at=row["last_seen_at"],
        revoked_at=row["revoked_at"],
        refusal=row["refusal"],
    )
