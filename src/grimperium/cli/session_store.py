"""Persistent session state for distributed-mode runs.

Session state is stored in ~/.grimperium/session.json and survives process
restarts. It records which workers were active and what config they received,
allowing the CLI to re-attach to a running session.
"""

import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any


@dataclass
class WorkerSessionInfo:
    """Snapshot of a single worker's config at session start."""

    worker_id: str
    hostname: str
    batch_size: int = 10
    profile_name: str = "default"
    crest_timeout_minutes: int = 60
    mopac_timeout_minutes: int = 30

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "WorkerSessionInfo":
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in known})


@dataclass
class SessionState:
    """State of a distributed processing session."""

    started_at: str
    server_url: str = "http://localhost:8000"
    workers: list[WorkerSessionInfo] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "started_at": self.started_at,
            "server_url": self.server_url,
            "workers": [w.to_dict() for w in self.workers],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SessionState":
        workers = [
            WorkerSessionInfo.from_dict(w) for w in data.get("workers", [])
        ]
        return cls(
            started_at=data.get("started_at", ""),
            server_url=data.get("server_url", "http://localhost:8000"),
            workers=workers,
        )


def _session_path() -> Path:
    """Return ~/.grimperium/session.json, creating the directory if needed."""
    try:
        home = Path.home() / ".grimperium"
    except RuntimeError:
        raise RuntimeError(
            "Cannot determine home directory ($HOME not set). "
            "Set $HOME before running Grimperium."
        ) from None
    home.mkdir(parents=True, exist_ok=True)
    return home / "session.json"


def load_session() -> SessionState | None:
    """Load the persisted session, returning None if none exists or is corrupt."""
    path = _session_path()
    if not path.is_file():
        return None
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
        return SessionState.from_dict(raw)
    except (json.JSONDecodeError, OSError, TypeError, KeyError):
        return None


def save_session(state: SessionState) -> bool:
    """Persist session state to ~/.grimperium/session.json.

    Returns:
        True on success, False on I/O error.
    """
    path = _session_path()
    try:
        path.write_text(
            json.dumps(state.to_dict(), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        return True
    except OSError:
        return False


def delete_session() -> bool:
    """Remove the session file.

    Returns:
        True if the file was deleted, False if it did not exist.
    """
    path = _session_path()
    try:
        path.unlink()
        return True
    except FileNotFoundError:
        return False
