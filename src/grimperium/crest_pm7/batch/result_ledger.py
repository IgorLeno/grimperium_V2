"""Append-only ledger with durable prepare/commit journal for sync idempotency."""

from __future__ import annotations

import hashlib
import json
import os
import threading
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1


class LedgerStatus(str, Enum):
    """Possible outcomes of a ledger idempotency check."""

    APPLIED = "applied"
    DUPLICATE = "duplicate"
    CONFLICT = "conflict"


class JournalTxnStatus(str, Enum):
    """Durable transaction states for crash recovery."""

    PREPARED = "prepared"
    COMMITTED = "committed"
    FAILED = "failed"


@dataclass(frozen=True)
class LedgerDecision:
    """Result of checking and optionally recording a synced result."""

    status: LedgerStatus
    result_id: str
    fingerprint: str

    @property
    def duplicate(self) -> bool:
        """Whether the result was already recorded with the same fingerprint."""
        return self.status is LedgerStatus.DUPLICATE

    @property
    def conflict(self) -> bool:
        """Whether the result ID was reused for different content."""
        return self.status is LedgerStatus.CONFLICT


@dataclass(frozen=True)
class JournalEntry:
    """One durable sync transaction record."""

    result_id: str
    fingerprint: str
    mol_id: str
    txn_status: JournalTxnStatus
    desired_success: bool
    previous_status: str | None = None
    previous_reruns: int | None = None
    final_status: str | None = None
    prepared_at: str | None = None
    committed_at: str | None = None
    error: str | None = None


class ResultLedger:
    """Idempotency ledger with in-memory index and durable journal."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.journal_path = self.path.with_name(
            self.path.stem + "_journal.jsonl"
            if self.path.suffix
            else self.path.name + "_journal.jsonl"
        )
        self._lock = threading.RLock()
        self._index: dict[str, str] = {}
        self._journal: dict[str, JournalEntry] = {}
        self._load_index()
        self._load_journal()

    def check(self, result_id: str, fingerprint: str) -> LedgerDecision:
        """Check a result ID without writing a ledger entry."""
        if not result_id:
            raise ValueError("result_id must not be blank")
        if not fingerprint:
            raise ValueError("fingerprint must not be blank")
        with self._lock:
            existing = self._index.get(result_id)
            if existing is None:
                journal = self._journal.get(result_id)
                if (
                    journal is not None
                    and journal.txn_status is JournalTxnStatus.COMMITTED
                ):
                    existing = journal.fingerprint
            if existing is None:
                return LedgerDecision(
                    status=LedgerStatus.APPLIED,
                    result_id=result_id,
                    fingerprint=fingerprint,
                )
            if existing == fingerprint:
                return LedgerDecision(
                    status=LedgerStatus.DUPLICATE,
                    result_id=result_id,
                    fingerprint=fingerprint,
                )
            return LedgerDecision(
                status=LedgerStatus.CONFLICT,
                result_id=result_id,
                fingerprint=fingerprint,
            )

    def check_and_record(
        self,
        result_id: str,
        mol_id: str,
        fingerprint: str,
    ) -> LedgerDecision:
        """Legacy helper: prepare+commit in one step after external apply."""
        decision = self.check(result_id, fingerprint)
        if decision.status is not LedgerStatus.APPLIED:
            return decision
        self.prepare(
            result_id=result_id,
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=True,
        )
        self.commit(result_id, final_status="applied")
        return LedgerDecision(
            status=LedgerStatus.APPLIED,
            result_id=result_id,
            fingerprint=fingerprint,
        )

    def prepare(
        self,
        *,
        result_id: str,
        mol_id: str,
        fingerprint: str,
        desired_success: bool,
        previous_status: str | None = None,
        previous_reruns: int | None = None,
    ) -> JournalEntry:
        """Durably mark a sync transaction as prepared before dual-write."""
        if not result_id or not fingerprint:
            raise ValueError("result_id and fingerprint must not be blank")
        now = datetime.now(timezone.utc).isoformat()
        entry = JournalEntry(
            result_id=result_id,
            fingerprint=fingerprint,
            mol_id=mol_id,
            txn_status=JournalTxnStatus.PREPARED,
            desired_success=desired_success,
            previous_status=previous_status,
            previous_reruns=previous_reruns,
            prepared_at=now,
        )
        with self._lock:
            self._append_journal(entry)
            self._journal[result_id] = entry
        return entry

    def commit(
        self, result_id: str, *, final_status: str | None = None
    ) -> JournalEntry:
        """Mark a prepared transaction committed and index the fingerprint."""
        with self._lock:
            previous = self._journal.get(result_id)
            if previous is None:
                raise ValueError(f"No prepared journal entry for {result_id}")
            if previous.txn_status is JournalTxnStatus.COMMITTED:
                return previous
            entry = JournalEntry(
                result_id=previous.result_id,
                fingerprint=previous.fingerprint,
                mol_id=previous.mol_id,
                txn_status=JournalTxnStatus.COMMITTED,
                desired_success=previous.desired_success,
                previous_status=previous.previous_status,
                previous_reruns=previous.previous_reruns,
                final_status=final_status or previous.final_status,
                prepared_at=previous.prepared_at,
                committed_at=datetime.now(timezone.utc).isoformat(),
                error=previous.error,
            )
            self._append_journal(entry)
            self._journal[result_id] = entry
            self._index[result_id] = entry.fingerprint
            self._append_committed_index(entry)
            return entry

    def mark_failed(self, result_id: str, *, error: str) -> JournalEntry:
        """Record that a prepared transaction failed before commit."""
        with self._lock:
            previous = self._journal.get(result_id)
            if previous is None:
                raise ValueError(f"No prepared journal entry for {result_id}")
            entry = JournalEntry(
                result_id=previous.result_id,
                fingerprint=previous.fingerprint,
                mol_id=previous.mol_id,
                txn_status=JournalTxnStatus.FAILED,
                desired_success=previous.desired_success,
                previous_status=previous.previous_status,
                previous_reruns=previous.previous_reruns,
                final_status=previous.final_status,
                prepared_at=previous.prepared_at,
                committed_at=None,
                error=error,
            )
            self._append_journal(entry)
            self._journal[result_id] = entry
            return entry

    def get_incomplete(self) -> list[JournalEntry]:
        """Return prepared (not committed/failed-final) transactions for recovery."""
        with self._lock:
            return [
                entry
                for entry in self._journal.values()
                if entry.txn_status is JournalTxnStatus.PREPARED
            ]

    def recover_incomplete(
        self,
        *,
        is_already_applied: Callable[[JournalEntry], bool],
    ) -> list[JournalEntry]:
        """Recover prepared transactions using an external applied-state probe.

        ``is_already_applied(entry)`` should return True when the dual-write
        effect is already visible (so we can commit without reapply).
        """
        recovered: list[JournalEntry] = []
        for entry in self.get_incomplete():
            if is_already_applied(entry):
                recovered.append(
                    self.commit(entry.result_id, final_status=entry.final_status)
                )
            else:
                # Leave prepared so the next sync can resume safely via check().
                recovered.append(entry)
        return recovered

    def _load_index(self) -> None:
        if not self.path.exists():
            self._index = {}
            return
        entries: dict[str, str] = {}
        with self.path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                payload = json.loads(line)
                result_id = str(payload.get("result_id", ""))
                fingerprint = str(payload.get("decision_fingerprint", ""))
                if result_id and fingerprint:
                    entries[result_id] = fingerprint
        self._index = entries

    def _load_journal(self) -> None:
        if not self.journal_path.exists():
            self._journal = {}
            return
        journal: dict[str, JournalEntry] = {}
        with self.journal_path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                payload = json.loads(line)
                entry = JournalEntry(
                    result_id=str(payload["result_id"]),
                    fingerprint=str(payload["fingerprint"]),
                    mol_id=str(payload["mol_id"]),
                    txn_status=JournalTxnStatus(str(payload["txn_status"])),
                    desired_success=bool(payload.get("desired_success", True)),
                    previous_status=(
                        str(payload["previous_status"])
                        if payload.get("previous_status") is not None
                        else None
                    ),
                    previous_reruns=(
                        int(payload["previous_reruns"])
                        if payload.get("previous_reruns") is not None
                        else None
                    ),
                    final_status=(
                        str(payload["final_status"])
                        if payload.get("final_status") is not None
                        else None
                    ),
                    prepared_at=(
                        str(payload["prepared_at"])
                        if payload.get("prepared_at") is not None
                        else None
                    ),
                    committed_at=(
                        str(payload["committed_at"])
                        if payload.get("committed_at") is not None
                        else None
                    ),
                    error=(
                        str(payload["error"])
                        if payload.get("error") is not None
                        else None
                    ),
                )
                journal[entry.result_id] = entry
                if entry.txn_status is JournalTxnStatus.COMMITTED:
                    self._index[entry.result_id] = entry.fingerprint
        self._journal = journal

    def _append_journal(self, entry: JournalEntry) -> None:
        self.journal_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": SCHEMA_VERSION,
            "result_id": entry.result_id,
            "fingerprint": entry.fingerprint,
            "mol_id": entry.mol_id,
            "txn_status": entry.txn_status.value,
            "desired_success": entry.desired_success,
            "previous_status": entry.previous_status,
            "previous_reruns": entry.previous_reruns,
            "final_status": entry.final_status,
            "prepared_at": entry.prepared_at,
            "committed_at": entry.committed_at,
            "error": entry.error,
        }
        _append_jsonl_fsync(self.journal_path, payload)

    def _append_committed_index(self, entry: JournalEntry) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        record = {
            "schema_version": SCHEMA_VERSION,
            "result_id": entry.result_id,
            "mol_id": entry.mol_id,
            "decision_fingerprint": entry.fingerprint,
            "applied_at": entry.committed_at or datetime.now(timezone.utc).isoformat(),
        }
        _append_jsonl_fsync(self.path, record)


def build_result_fingerprint(payload: dict[str, Any]) -> str:
    """Return a stable SHA-256 fingerprint for a synced result payload."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_legacy_result_id(payload: dict[str, Any]) -> str:
    """Deterministic result_id for legacy clients without an explicit ID.

    Uses only immutable payload fields — never mutable operational state like
    ``reruns``.
    """
    mol_id = str(payload.get("mol_id", ""))
    fingerprint = build_result_fingerprint(
        {
            "mol_id": payload.get("mol_id"),
            "success": payload.get("success"),
            "result_update": payload.get("result_update"),
            "error": payload.get("error"),
            "completed_at": payload.get("completed_at"),
        }
    )
    return f"legacy:{mol_id}:{fingerprint}"


def _append_jsonl_fsync(path: Path, payload: dict[str, Any]) -> None:
    line = json.dumps(payload, sort_keys=True) + "\n"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())
