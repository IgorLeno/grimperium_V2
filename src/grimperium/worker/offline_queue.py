"""Fila durável de resultados do worker com result_id estável."""

from __future__ import annotations

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

    def to_dict(self) -> dict[str, Any]:
        return {
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

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> DeadLetterRecord:
        original = payload.get("original_payload")
        return cls(
            result_id=str(payload["result_id"]),
            mol_id=str(payload["mol_id"]),
            attempt_id=(
                str(payload["attempt_id"])
                if payload.get("attempt_id") is not None
                else None
            ),
            original_payload=dict(original) if isinstance(original, dict) else {},
            returned_status=str(payload["returned_status"]),
            detail=(
                str(payload["detail"]) if payload.get("detail") is not None else None
            ),
            worker_id=str(payload["worker_id"]),
            rejected_at=str(payload["rejected_at"]),
            rejection_origin=str(payload["rejection_origin"]),
        )


class DeadLetterQueue:
    """JSONL append-only para rejeições terminais do sync (auditoria local)."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._entries: list[DeadLetterRecord] = []
        self._load()

    def append(self, record: DeadLetterRecord) -> None:
        """Gravar um registro de dead-letter de forma atômica."""
        self._entries.append(record)
        self._persist()

    def entries(self) -> list[DeadLetterRecord]:
        return list(self._entries)

    def _load(self) -> None:
        if not self.path.exists():
            return
        loaded: list[DeadLetterRecord] = []
        for line in self.path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            if isinstance(payload, dict) and payload.get("result_id"):
                loaded.append(DeadLetterRecord.from_dict(payload))
        self._entries = loaded

    def _persist(self) -> None:
        lines = [
            json.dumps(entry.to_dict(), ensure_ascii=False, sort_keys=True)
            for entry in self._entries
        ]
        content = ("\n".join(lines) + "\n") if lines else ""
        _atomic_write_text(self.path, content)
