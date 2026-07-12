"""In-memory store for molecules claimed by a worker during a session."""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from grimperium.worker.offline_queue import new_result_id


@dataclass
class LocalRecord:
    """Tracks a single molecule from claim through completion."""

    mol_id: str
    smiles: str
    claimed_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    completed: bool = False
    success: bool | None = None
    result_update: dict[str, Any] = field(default_factory=dict)
    error: str | None = None
    completed_at: datetime | None = None
    result_id: str | None = None


class LocalStore:
    """In-memory store for molecules claimed by this worker.

    Thread-safety: not thread-safe. The runner calls this only from its main
    thread; the heartbeat thread does not touch it.
    """

    def __init__(self) -> None:
        self._records: dict[str, LocalRecord] = {}

    def add(self, mol_id: str, smiles: str) -> LocalRecord:
        if mol_id in self._records:
            raise ValueError(f"mol_id {mol_id!r} already in store")
        record = LocalRecord(mol_id=mol_id, smiles=smiles)
        self._records[mol_id] = record
        return record

    def get(self, mol_id: str) -> LocalRecord | None:
        return self._records.get(mol_id)

    def mark_success(self, mol_id: str, result_update: dict[str, Any]) -> None:
        record = self._records[mol_id]
        record.completed = True
        record.success = True
        record.result_update = result_update
        record.completed_at = datetime.now(timezone.utc)
        if not record.result_id:
            record.result_id = new_result_id()

    def mark_failure(self, mol_id: str, error: str) -> None:
        record = self._records[mol_id]
        record.completed = True
        record.success = False
        record.error = error
        record.completed_at = datetime.now(timezone.utc)
        if not record.result_id:
            record.result_id = new_result_id()

    def pending(self) -> list[LocalRecord]:
        return [r for r in self._records.values() if not r.completed]

    def completed(self) -> list[LocalRecord]:
        return [r for r in self._records.values() if r.completed]

    def clear_completed(self) -> int:
        done = [mid for mid, r in self._records.items() if r.completed]
        for mid in done:
            del self._records[mid]
        return len(done)

    def __len__(self) -> int:
        return len(self._records)
