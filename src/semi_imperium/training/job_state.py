"""Single-row job state in SQLite: who runs, how far, and whether to stop.

The coordinator writes a heartbeat; readers (dashboard API, CLI) decide whether
a job that claims to run is really alive by checking the pid and the
heartbeat age. A job that is not alive is reported as ``stale`` and does not
block a new one.
"""

from __future__ import annotations

import os
import sqlite3
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from semi_imperium.training.db import connect

STALE_AFTER_SECONDS = 180.0
_ACTIVE = ("running", "stopping")


class JobAlreadyRunning(RuntimeError):
    """Another training job is alive."""


@dataclass(frozen=True)
class JobSnapshot:
    """What readers see of the job, with staleness already decided."""

    state: str
    pid: int | None
    started_at: str | None
    heartbeat_at: str | None
    finished_at: str | None
    workers: int
    limit: int | None
    done: int
    failed: int
    message: str
    stop_requested: bool

    @property
    def is_active(self) -> bool:
        """Whether a live process is running or winding down."""
        return self.state in _ACTIVE


_IDLE = JobSnapshot("idle", None, None, None, None, 0, None, 0, 0, "", False)


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _pid_alive(pid: int | None) -> bool:
    if pid is None:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OverflowError:
        return False
    return True


def _snapshot(row: sqlite3.Row) -> JobSnapshot:
    state = str(row["state"])
    if state in _ACTIVE:
        beat = row["heartbeat_at"]
        age = (
            (_now() - datetime.fromisoformat(beat)).total_seconds()
            if beat
            else float("inf")
        )
        if not _pid_alive(row["pid"]) or age > STALE_AFTER_SECONDS:
            state = "stale"
    return JobSnapshot(
        state=state,
        pid=row["pid"],
        started_at=row["started_at"],
        heartbeat_at=row["heartbeat_at"],
        finished_at=row["finished_at"],
        workers=int(row["workers"]),
        limit=row["limit_count"],
        done=int(row["done"]),
        failed=int(row["failed"]),
        message=str(row["message"]),
        stop_requested=bool(row["stop_requested"]),
    )


class JobState:
    """Reads and writes the ``training_job`` and ``training_failure`` tables."""

    def __init__(self, path: Path | str) -> None:
        self._path = Path(path)

    def snapshot(self) -> JobSnapshot:
        """Current job state (``idle`` if nothing ever ran)."""
        if not self._path.exists():
            return _IDLE
        with closing(connect(self._path)) as connection:
            row = connection.execute(
                "SELECT * FROM training_job WHERE id = 1"
            ).fetchone()
        return _IDLE if row is None else _snapshot(row)

    def begin(self, *, pid: int, workers: int, limit: int | None) -> None:
        """Claim the job slot, atomically refusing when a live job exists."""
        now = _now().isoformat()
        with closing(connect(self._path)) as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                row = connection.execute(
                    "SELECT * FROM training_job WHERE id = 1"
                ).fetchone()
                if row is not None and _snapshot(row).is_active:
                    raise JobAlreadyRunning(
                        f"training job already running (pid {row['pid']})"
                    )
                connection.execute(
                    "INSERT OR REPLACE INTO training_job (id, state, pid, started_at,"
                    " heartbeat_at, finished_at, workers, limit_count, done, failed,"
                    " message, stop_requested)"
                    " VALUES (1, 'running', ?, ?, ?, NULL, ?, ?, 0, 0, '', 0)",
                    (pid, now, now, workers, limit),
                )
                connection.execute("COMMIT")
            except BaseException:
                connection.execute("ROLLBACK")
                raise

    def heartbeat(self, *, done: int, failed: int) -> None:
        """Record progress and prove the coordinator is alive."""
        with closing(connect(self._path)) as connection:
            connection.execute(
                "UPDATE training_job SET heartbeat_at = ?, done = ?, failed = ?"
                " WHERE id = 1",
                (_now().isoformat(), done, failed),
            )

    def request_stop(self) -> bool:
        """Ask the running job to stop dispatching; ``False`` if none is active."""
        if not self.snapshot().is_active:
            return False
        with closing(connect(self._path)) as connection:
            connection.execute(
                "UPDATE training_job SET stop_requested = 1 WHERE id = 1"
            )
        return True

    def stop_requested(self) -> bool:
        """Whether a stop was requested for the current job."""
        return self.snapshot().stop_requested

    def mark_stopping(self) -> None:
        """Show that the coordinator saw the stop and is draining in-flight work."""
        with closing(connect(self._path)) as connection:
            connection.execute(
                "UPDATE training_job SET state = 'stopping' WHERE id = 1"
                " AND state = 'running'"
            )

    def finish(self, state: str, message: str) -> None:
        """Close the job as ``finished`` or ``failed``."""
        with closing(connect(self._path)) as connection:
            connection.execute(
                "UPDATE training_job SET state = ?, message = ?, finished_at = ?"
                " WHERE id = 1",
                (state, message, _now().isoformat()),
            )

    def record_failure(self, mol_id: str, reason: str) -> None:
        """Keep why one molecule failed (the CSV only keeps the status)."""
        with closing(connect(self._path)) as connection:
            connection.execute(
                "INSERT INTO training_failure (mol_id, reason, at) VALUES (?, ?, ?)",
                (mol_id, reason, _now().isoformat()),
            )

    def recent_failures(self, limit: int = 20) -> list[tuple[str, str, str]]:
        """``(mol_id, reason, at)`` of the latest failures, newest first."""
        if not self._path.exists():
            return []
        with closing(connect(self._path)) as connection:
            rows = connection.execute(
                "SELECT mol_id, reason, at FROM training_failure"
                " ORDER BY id DESC LIMIT ?",
                (limit,),
            ).fetchall()
        return [(str(r["mol_id"]), str(r["reason"]), str(r["at"])) for r in rows]
