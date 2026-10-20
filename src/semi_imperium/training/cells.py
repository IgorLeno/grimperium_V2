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
from contextlib import closing, contextmanager
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

    def set_admission(self, cell_id: str, admission: str) -> Cell:
        if admission not in ADMISSIONS:
            raise ValueError(f"admission must be one of {ADMISSIONS}, got {admission}")
        return self._update(
            cell_id, "UPDATE cells SET admission = ? WHERE cell_id = ?", (admission,)
        )

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

    @contextmanager
    def _write(self) -> Iterator[sqlite3.Connection]:
        with closing(connect(self._db_path)) as connection:
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
