"""Append-only ledger for idempotent distributed result sync."""

from __future__ import annotations

import hashlib
import json
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


class ResultLedger:
    """Append-only JSONL ledger keyed by result_id."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)

    def check_and_record(
        self,
        result_id: str,
        mol_id: str,
        fingerprint: str,
    ) -> LedgerDecision:
        """Record a result ID unless it is a duplicate or conflict."""
        if not result_id:
            raise ValueError("result_id must not be blank")
        if not fingerprint:
            raise ValueError("fingerprint must not be blank")

        decision = self.check(result_id, fingerprint)
        if decision.status is not LedgerStatus.APPLIED:
            return decision

        self.path.parent.mkdir(parents=True, exist_ok=True)
        record = {
            "schema_version": SCHEMA_VERSION,
            "result_id": result_id,
            "mol_id": mol_id,
            "decision_fingerprint": fingerprint,
            "applied_at": datetime.now(timezone.utc).isoformat(),
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")

        return LedgerDecision(
            status=LedgerStatus.APPLIED,
            result_id=result_id,
            fingerprint=fingerprint,
        )

    def check(self, result_id: str, fingerprint: str) -> LedgerDecision:
        """Check a result ID without writing a ledger entry."""
        if not result_id:
            raise ValueError("result_id must not be blank")
        if not fingerprint:
            raise ValueError("fingerprint must not be blank")
        existing = self._load_entries().get(result_id)
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

    def _load_entries(self) -> dict[str, str]:
        """Load result_id -> fingerprint from the append-only JSONL file."""
        if not self.path.exists():
            return {}

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
        return entries


def build_result_fingerprint(payload: dict[str, Any]) -> str:
    """Return a stable SHA-256 fingerprint for a synced result payload."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()
