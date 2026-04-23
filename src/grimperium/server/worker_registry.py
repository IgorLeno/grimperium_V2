"""In-memory worker registry for the distributed processing server.

WorkerRegistry is the single source of truth for per-worker state during a
session. It replaces the two raw dicts (heartbeat_registry + running_molecules)
on app.state and adds per-worker metrics, config overrides, and shutdown flags.
"""

import threading
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any


@dataclass
class WorkerEntry:
    """All mutable state for a single registered worker."""

    worker_id: str
    hostname: str
    registered_at: datetime
    last_seen: datetime
    current_mol_id: str | None = None
    processed: int = 0
    successful: int = 0
    failed: int = 0
    skipped: int = 0
    config_override: dict[str, Any] | None = None
    shutdown_requested: bool = False


class WorkerRegistry:
    """Thread-safe registry of connected workers and their runtime state.

    Used by the FastAPI server to track all registered workers, their current
    molecule, processing metrics, and per-worker config overrides sent via
    POST /configure/{worker_id}.

    The registry also exposes a heartbeat_registry view compatible with the
    watchdog's HeartbeatRegistry type alias so the watchdog can still evict
    stale workers without knowing about the full registry.
    """

    def __init__(self) -> None:
        self._workers: dict[str, WorkerEntry] = {}
        self._lock = threading.Lock()

    # ── Registration ──────────────────────────────────────────────────────

    def register(self, worker_id: str, hostname: str) -> WorkerEntry:
        """Register or refresh a worker. Returns the (possibly new) entry."""
        now = datetime.now(timezone.utc)
        with self._lock:
            if worker_id in self._workers:
                entry = self._workers[worker_id]
                entry.last_seen = now
                entry.hostname = hostname
            else:
                entry = WorkerEntry(
                    worker_id=worker_id,
                    hostname=hostname,
                    registered_at=now,
                    last_seen=now,
                )
                self._workers[worker_id] = entry
        return entry

    def evict(self, worker_id: str) -> None:
        """Remove a worker from the registry (called by watchdog on timeout)."""
        with self._lock:
            self._workers.pop(worker_id, None)

    # ── Heartbeat ─────────────────────────────────────────────────────────

    def heartbeat(self, worker_id: str) -> bool:
        """Update last_seen for a registered worker.

        Returns:
            True if the worker was found and updated; False if unknown.
        """
        now = datetime.now(timezone.utc)
        with self._lock:
            if worker_id not in self._workers:
                return False
            self._workers[worker_id].last_seen = now
        return True

    # ── Molecule tracking ─────────────────────────────────────────────────

    def set_current_mol(self, worker_id: str, mol_id: str) -> None:
        """Record that a worker has claimed a molecule."""
        with self._lock:
            if worker_id in self._workers:
                self._workers[worker_id].current_mol_id = mol_id

    def clear_current_mol(self, worker_id: str) -> None:
        """Clear the current molecule after completion."""
        with self._lock:
            if worker_id in self._workers:
                self._workers[worker_id].current_mol_id = None

    def get_worker_for_mol(self, mol_id: str) -> str | None:
        """Return the worker_id that currently holds mol_id, or None."""
        with self._lock:
            for entry in self._workers.values():
                if entry.current_mol_id == mol_id:
                    return entry.worker_id
        return None

    # ── Metrics ───────────────────────────────────────────────────────────

    def record_success(self, worker_id: str) -> None:
        with self._lock:
            if worker_id in self._workers:
                e = self._workers[worker_id]
                e.processed += 1
                e.successful += 1
                e.current_mol_id = None

    def record_failure(self, worker_id: str) -> None:
        with self._lock:
            if worker_id in self._workers:
                e = self._workers[worker_id]
                e.processed += 1
                e.failed += 1
                e.current_mol_id = None

    def record_skip(self, worker_id: str) -> None:
        with self._lock:
            if worker_id in self._workers:
                e = self._workers[worker_id]
                e.processed += 1
                e.skipped += 1
                e.current_mol_id = None

    # ── Config override ───────────────────────────────────────────────────

    def set_config(self, worker_id: str, config: dict[str, Any]) -> bool:
        """Store a per-worker config override.

        Returns:
            True if the worker exists; False otherwise.
        """
        with self._lock:
            if worker_id not in self._workers:
                return False
            self._workers[worker_id].config_override = config
        return True

    def get_config(self, worker_id: str) -> dict[str, Any] | None:
        """Return the per-worker config override, or None if not set."""
        with self._lock:
            entry = self._workers.get(worker_id)
            return entry.config_override if entry else None

    # ── Shutdown signalling ────────────────────────────────────────────────

    def request_shutdown(self, worker_id: str) -> bool:
        """Signal a worker to shut down after finishing its current molecule.

        Returns:
            True if the worker exists; False if unknown.
        """
        with self._lock:
            if worker_id not in self._workers:
                return False
            self._workers[worker_id].shutdown_requested = True
        return True

    def is_shutdown_requested(self, worker_id: str) -> bool:
        with self._lock:
            entry = self._workers.get(worker_id)
            return entry.shutdown_requested if entry else False

    def request_shutdown_all(self) -> list[str]:
        """Signal all registered workers to shut down.

        Returns:
            List of worker_ids that were signalled.
        """
        with self._lock:
            ids = list(self._workers.keys())
            for entry in self._workers.values():
                entry.shutdown_requested = True
        return ids

    # ── Queries ───────────────────────────────────────────────────────────

    def all_workers(self) -> list[WorkerEntry]:
        """Snapshot of all registered workers (safe copy)."""
        with self._lock:
            return list(self._workers.values())

    def get_worker(self, worker_id: str) -> WorkerEntry | None:
        with self._lock:
            return self._workers.get(worker_id)

    def __len__(self) -> int:
        with self._lock:
            return len(self._workers)

    # ── Watchdog compatibility view ────────────────────────────────────────

    @property
    def heartbeat_registry(self) -> dict[str, tuple[str, datetime]]:
        """Live view of {worker_id: (hostname, last_seen)} for the watchdog.

        NOTE: Mutations (del registry[wid]) are NOT reflected back into the
        WorkerRegistry. The watchdog should call registry.evict() instead.
        This property is provided for gradual migration only.
        """
        with self._lock:
            return {
                wid: (e.hostname, e.last_seen)
                for wid, e in self._workers.items()
            }

    @property
    def running_molecules(self) -> dict[str, str]:
        """Live snapshot of {mol_id: worker_id} for legacy route compat."""
        with self._lock:
            return {
                e.current_mol_id: e.worker_id
                for e in self._workers.values()
                if e.current_mol_id is not None
            }


def make_worker_registry() -> WorkerRegistry:
    return WorkerRegistry()
