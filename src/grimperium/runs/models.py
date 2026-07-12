"""Data models for persisted Grimperium runs."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any


class RunStatus(str, Enum):
    """Lifecycle states for a persisted run manifest."""

    CREATED = "created"
    RUNNING = "running"
    COMPLETED = "completed"
    PARTIAL = "partial"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INVALIDATED = "invalidated"


TERMINAL_STATUSES = frozenset(
    {
        RunStatus.COMPLETED,
        RunStatus.PARTIAL,
        RunStatus.FAILED,
        RunStatus.CANCELLED,
        RunStatus.INVALIDATED,
    }
)

# Explicit lifecycle transition matrix. Terminal states accept no further transitions.
ALLOWED_TRANSITIONS: dict[RunStatus, frozenset[RunStatus]] = {
    RunStatus.CREATED: frozenset(
        {RunStatus.RUNNING, RunStatus.CANCELLED, RunStatus.FAILED}
    ),
    RunStatus.RUNNING: frozenset(
        {
            RunStatus.COMPLETED,
            RunStatus.PARTIAL,
            RunStatus.FAILED,
            RunStatus.CANCELLED,
            RunStatus.INVALIDATED,
        }
    ),
    RunStatus.COMPLETED: frozenset(),
    RunStatus.PARTIAL: frozenset(),
    RunStatus.FAILED: frozenset(),
    RunStatus.CANCELLED: frozenset(),
    RunStatus.INVALIDATED: frozenset(),
}


@dataclass(frozen=True)
class RunManifest:
    """Persisted manifest describing one immutable scientific execution."""

    schema_version: str
    run_id: str
    status: RunStatus
    created_at: datetime
    started_at: datetime | None
    completed_at: datetime | None
    property_id: str
    method_id: str
    method_version: str
    method_snapshot: dict[str, Any]
    execution_overrides: dict[str, Any]
    dataset_ref: dict[str, Any] | None
    model_ref: dict[str, Any] | None
    molecule_count: int
    success_count: int
    failure_count: int
    output_paths: dict[str, Path]
    error: str | None

    def with_updates(self, **updates: Any) -> RunManifest:
        """Return a copy with selected fields updated."""
        return replace(self, **updates)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "schema_version": self.schema_version,
            "run_id": self.run_id,
            "status": self.status.value,
            "created_at": self.created_at.isoformat(),
            "started_at": (
                self.started_at.isoformat() if self.started_at is not None else None
            ),
            "completed_at": (
                self.completed_at.isoformat() if self.completed_at is not None else None
            ),
            "property_id": self.property_id,
            "method_id": self.method_id,
            "method_version": self.method_version,
            "method_snapshot": self.method_snapshot,
            "execution_overrides": self.execution_overrides,
            "dataset_ref": self.dataset_ref,
            "model_ref": self.model_ref,
            "molecule_count": self.molecule_count,
            "success_count": self.success_count,
            "failure_count": self.failure_count,
            "output_paths": {key: str(path) for key, path in self.output_paths.items()},
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> RunManifest:
        """Deserialize from JSON-compatible primitives."""
        return cls(
            schema_version=str(payload["schema_version"]),
            run_id=str(payload["run_id"]),
            status=RunStatus(str(payload["status"])),
            created_at=datetime.fromisoformat(str(payload["created_at"])),
            started_at=_optional_datetime(payload.get("started_at")),
            completed_at=_optional_datetime(payload.get("completed_at")),
            property_id=str(payload["property_id"]),
            method_id=str(payload["method_id"]),
            method_version=str(payload["method_version"]),
            method_snapshot=dict(payload.get("method_snapshot", {})),
            execution_overrides=dict(payload.get("execution_overrides", {})),
            dataset_ref=_optional_dict(payload.get("dataset_ref")),
            model_ref=_optional_dict(payload.get("model_ref")),
            molecule_count=int(payload.get("molecule_count", 0)),
            success_count=int(payload.get("success_count", 0)),
            failure_count=int(payload.get("failure_count", 0)),
            output_paths={
                str(key): Path(str(value))
                for key, value in dict(payload.get("output_paths", {})).items()
            },
            error=(str(payload["error"]) if payload.get("error") is not None else None),
        )


def _optional_datetime(value: object) -> datetime | None:
    if value is None:
        return None
    return datetime.fromisoformat(str(value))


def _optional_dict(value: object) -> dict[str, Any] | None:
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError("Run manifest reference fields must be objects or null")
    return dict(value)
