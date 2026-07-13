"""Fila durável de resultados do worker com result_id estável."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def new_result_id() -> str:
    """Gerar um result_id explícito e estável para retries."""
    return str(uuid.uuid4())


@dataclass(frozen=True)
class OfflineResult:
    """Um resultado pendente de confirmação pelo servidor."""

    result_id: str
    mol_id: str
    success: bool
    result_update: dict[str, Any] | None
    error: str | None
    completed_at: str
    attempt_id: str | None = None

    def to_sync_dict(self) -> dict[str, Any]:
        """Payload compatível com SyncResult / POST /sync_results."""
        payload: dict[str, Any] = {
            "result_id": self.result_id,
            "mol_id": self.mol_id,
            "success": self.success,
            "result_update": self.result_update,
            "error": self.error,
            "completed_at": self.completed_at,
        }
        if self.attempt_id is not None:
            payload["attempt_id"] = self.attempt_id
        return payload

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> OfflineResult:
        return cls(
            result_id=str(payload["result_id"]),
            mol_id=str(payload["mol_id"]),
            success=bool(payload["success"]),
            result_update=(
                dict(payload["result_update"])
                if isinstance(payload.get("result_update"), dict)
                else None
            ),
            error=(str(payload["error"]) if payload.get("error") is not None else None),
            completed_at=str(payload["completed_at"]),
            attempt_id=(
                str(payload["attempt_id"])
                if payload.get("attempt_id") is not None
                else None
            ),
        )


class OfflineResultQueue:
    """JSONL durável: o mesmo result_id é reenviado até confirmação."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._entries: dict[str, OfflineResult] = {}
        self._load()

    def enqueue(
        self,
        *,
        mol_id: str,
        success: bool,
        result_update: dict[str, Any] | None = None,
        error: str | None = None,
        result_id: str | None = None,
        completed_at: str | None = None,
        attempt_id: str | None = None,
    ) -> OfflineResult:
        """Persistir um resultado com ID estável (gera um se omitido)."""
        entry = OfflineResult(
            result_id=result_id or new_result_id(),
            mol_id=mol_id,
            success=success,
            result_update=result_update,
            error=error,
            completed_at=completed_at
            or datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            attempt_id=attempt_id,
        )
        self._entries[entry.result_id] = entry
        self._persist()
        return entry

    def pending(self) -> list[OfflineResult]:
        return list(self._entries.values())

    def confirm(self, result_id: str) -> None:
        """Remover um resultado após confirmação do servidor."""
        if result_id in self._entries:
            del self._entries[result_id]
            self._persist()

    def _load(self) -> None:
        if not self.path.exists():
            return
        for line in self.path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            if isinstance(payload, dict) and payload.get("result_id"):
                entry = OfflineResult.from_dict(payload)
                self._entries[entry.result_id] = entry

    def _persist(self) -> None:
        lines = [
            json.dumps(entry.to_sync_dict(), ensure_ascii=False, sort_keys=True)
            for entry in self._entries.values()
        ]
        content = ("\n".join(lines) + "\n") if lines else ""
        _atomic_write_text(self.path, content)


def dead_letter_path_for(queue_path: Path) -> Path:
    """Derive the durable dead-letter sibling path from the offline queue path."""
    path = Path(queue_path)
    return path.with_name(f"{path.stem}_dead_letter{path.suffix}")


def pending_abort_path_for(queue_path: Path) -> Path:
    """Derive the durable pending-abort sibling path from the offline queue path."""
    path = Path(queue_path)
    return path.with_name(f"{path.stem}_pending_aborts{path.suffix}")


def _atomic_write_text(path: Path, content: str) -> None:
    """Persist ``content`` atomically via NamedTemporaryFile + fsync + replace."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temp_path = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink(missing_ok=True)


@dataclass(frozen=True)
class DeadLetterRecord:
    """Registro durável de rejeição terminal (conflict / stale_attempt)."""

    result_id: str
    mol_id: str
    attempt_id: str | None
    original_payload: dict[str, Any]
    returned_status: str
    detail: str | None
    worker_id: str
    rejected_at: str
    rejection_origin: str
    dead_letter_id: str | None = None

    def to_dict(self) -> dict[str, Any]:
        payload = {
            "result_id": self.result_id,
            "mol_id": self.mol_id,
            "attempt_id": self.attempt_id,
            "original_payload": self.original_payload,
            "returned_status": self.returned_status,
            "detail": self.detail,
            "worker_id": self.worker_id,
            "rejected_at": self.rejected_at,
            "rejection_origin": self.rejection_origin,
        }
        dead_letter_id = self.dead_letter_id or compute_dead_letter_id(
            result_id=self.result_id,
            returned_status=self.returned_status,
            original_payload=self.original_payload,
        )
        payload["dead_letter_id"] = dead_letter_id
        return payload

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> DeadLetterRecord:
        original = payload.get("original_payload")
        original_payload = dict(original) if isinstance(original, dict) else {}
        returned_status = str(payload["returned_status"])
        result_id = str(payload["result_id"])
        dead_letter_id = (
            str(payload["dead_letter_id"])
            if payload.get("dead_letter_id") is not None
            else compute_dead_letter_id(
                result_id=result_id,
                returned_status=returned_status,
                original_payload=original_payload,
            )
        )
        return cls(
            result_id=result_id,
            mol_id=str(payload["mol_id"]),
            attempt_id=(
                str(payload["attempt_id"])
                if payload.get("attempt_id") is not None
                else None
            ),
            original_payload=original_payload,
            returned_status=returned_status,
            detail=(
                str(payload["detail"]) if payload.get("detail") is not None else None
            ),
            worker_id=str(payload["worker_id"]),
            rejected_at=str(payload["rejected_at"]),
            rejection_origin=str(payload["rejection_origin"]),
            dead_letter_id=dead_letter_id,
        )


def compute_dead_letter_id(
    *,
    result_id: str,
    returned_status: str,
    original_payload: dict[str, Any],
) -> str:
    """Identidade estável de uma rejeição para deduplicar na janela de crash."""
    encoded = json.dumps(
        {
            "result_id": result_id,
            "returned_status": returned_status,
            "original_payload": original_payload,
        },
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


class DeadLetterQueue:
    """JSONL append-only para rejeições terminais do sync (auditoria local)."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._entries: list[DeadLetterRecord] = []
        self._index: set[str] = set()
        self._load()

    def append(self, record: DeadLetterRecord) -> None:
        """Gravar um registro de dead-letter de forma atômica e idempotente."""
        dead_letter_id = record.dead_letter_id or compute_dead_letter_id(
            result_id=record.result_id,
            returned_status=record.returned_status,
            original_payload=record.original_payload,
        )
        if dead_letter_id in self._index:
            return
        normalized = DeadLetterRecord(
            result_id=record.result_id,
            mol_id=record.mol_id,
            attempt_id=record.attempt_id,
            original_payload=record.original_payload,
            returned_status=record.returned_status,
            detail=record.detail,
            worker_id=record.worker_id,
            rejected_at=record.rejected_at,
            rejection_origin=record.rejection_origin,
            dead_letter_id=dead_letter_id,
        )
        new_entries = [*self._entries, normalized]
        self._persist_entries(new_entries)
        self._entries = new_entries
        self._index.add(dead_letter_id)

    def entries(self) -> list[DeadLetterRecord]:
        return list(self._entries)

    def _load(self) -> None:
        if not self.path.exists():
            return
        loaded: list[DeadLetterRecord] = []
        index: set[str] = set()
        for line in self.path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(payload, dict) and payload.get("result_id"):
                record = DeadLetterRecord.from_dict(payload)
                dead_letter_id = record.dead_letter_id or compute_dead_letter_id(
                    result_id=record.result_id,
                    returned_status=record.returned_status,
                    original_payload=record.original_payload,
                )
                if dead_letter_id in index:
                    continue
                index.add(dead_letter_id)
                loaded.append(
                    DeadLetterRecord(
                        result_id=record.result_id,
                        mol_id=record.mol_id,
                        attempt_id=record.attempt_id,
                        original_payload=record.original_payload,
                        returned_status=record.returned_status,
                        detail=record.detail,
                        worker_id=record.worker_id,
                        rejected_at=record.rejected_at,
                        rejection_origin=record.rejection_origin,
                        dead_letter_id=dead_letter_id,
                    )
                )
        self._entries = loaded
        self._index = index

    def _persist(self) -> None:
        self._persist_entries(self._entries)

    def _persist_entries(self, entries: list[DeadLetterRecord]) -> None:
        lines = [
            json.dumps(entry.to_dict(), ensure_ascii=False, sort_keys=True)
            for entry in entries
        ]
        content = ("\n".join(lines) + "\n") if lines else ""
        _atomic_write_text(self.path, content)


class PendingAbortQueue:
    """JSONL durável para aborts de lease-loss pendentes de confirmação no DL."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._entries: list[DeadLetterRecord] = []
        self._index: set[str] = set()
        self._load()

    def append(self, record: DeadLetterRecord) -> None:
        """Persistir abort pendente antes de tentar gravar no dead-letter."""
        dead_letter_id = record.dead_letter_id or compute_dead_letter_id(
            result_id=record.result_id,
            returned_status=record.returned_status,
            original_payload=record.original_payload,
        )
        if dead_letter_id in self._index:
            return
        normalized = DeadLetterRecord(
            result_id=record.result_id,
            mol_id=record.mol_id,
            attempt_id=record.attempt_id,
            original_payload=record.original_payload,
            returned_status=record.returned_status,
            detail=record.detail,
            worker_id=record.worker_id,
            rejected_at=record.rejected_at,
            rejection_origin=record.rejection_origin,
            dead_letter_id=dead_letter_id,
        )
        new_entries = [*self._entries, normalized]
        self._persist_entries(new_entries)
        self._entries = new_entries
        self._index.add(dead_letter_id)

    def remove(self, dead_letter_id: str) -> None:
        """Remover abort após confirmação durável no dead-letter."""
        if dead_letter_id not in self._index:
            return
        new_entries = [
            entry
            for entry in self._entries
            if (
                entry.dead_letter_id
                or compute_dead_letter_id(
                    result_id=entry.result_id,
                    returned_status=entry.returned_status,
                    original_payload=entry.original_payload,
                )
            )
            != dead_letter_id
        ]
        self._persist_entries(new_entries)
        self._entries = new_entries
        self._index.discard(dead_letter_id)

    def entries(self) -> list[DeadLetterRecord]:
        return list(self._entries)

    def _load(self) -> None:
        if not self.path.exists():
            return
        loaded: list[DeadLetterRecord] = []
        index: set[str] = set()
        for line in self.path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(payload, dict) and payload.get("result_id"):
                record = DeadLetterRecord.from_dict(payload)
                dead_letter_id = record.dead_letter_id or compute_dead_letter_id(
                    result_id=record.result_id,
                    returned_status=record.returned_status,
                    original_payload=record.original_payload,
                )
                if dead_letter_id in index:
                    continue
                index.add(dead_letter_id)
                loaded.append(
                    DeadLetterRecord(
                        result_id=record.result_id,
                        mol_id=record.mol_id,
                        attempt_id=record.attempt_id,
                        original_payload=record.original_payload,
                        returned_status=record.returned_status,
                        detail=record.detail,
                        worker_id=record.worker_id,
                        rejected_at=record.rejected_at,
                        rejection_origin=record.rejection_origin,
                        dead_letter_id=dead_letter_id,
                    )
                )
        self._entries = loaded
        self._index = index

    def _persist_entries(self, entries: list[DeadLetterRecord]) -> None:
        lines = [
            json.dumps(entry.to_dict(), ensure_ascii=False, sort_keys=True)
            for entry in entries
        ]
        content = ("\n".join(lines) + "\n") if lines else ""
        _atomic_write_text(self.path, content)
