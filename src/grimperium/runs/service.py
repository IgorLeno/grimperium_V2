"""Service layer for run lifecycle management."""

from __future__ import annotations

import os
import uuid
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from grimperium.runs.models import TERMINAL_STATUSES, RunManifest, RunStatus
from grimperium.runs.persistence import RunManifestStore

DEFAULT_SCHEMA_VERSION = "1.0"
RUNS_DIR_ENV = "GRIMPERIUM_RUNS_DIR"


class RunService:
    """Coordinate run manifest creation and status transitions."""

    def __init__(self, runs_root: Path | str = Path("runs")) -> None:
        self.runs_root = Path(runs_root)
        self.store = RunManifestStore(self.runs_root)

    @classmethod
    def from_environment(cls) -> RunService:
        """Build a service from ``GRIMPERIUM_RUNS_DIR`` or cwd/runs."""
        return cls(Path(os.environ.get(RUNS_DIR_ENV, "runs")))

    def create_run(
        self,
        *,
        property_id: str,
        method_id: str,
        method_version: str,
        method_snapshot: dict[str, Any],
        execution_overrides: dict[str, Any],
        dataset_ref: dict[str, Any] | None,
        model_ref: dict[str, Any] | None,
        molecule_count: int,
    ) -> RunManifest:
        """Create and persist a new run manifest with a fresh run ID."""
        now = _utc_now()
        manifest = RunManifest(
            schema_version=DEFAULT_SCHEMA_VERSION,
            run_id=self._new_run_id(),
            status=RunStatus.CREATED,
            created_at=now,
            started_at=None,
            completed_at=None,
            property_id=property_id,
            method_id=method_id,
            method_version=method_version,
            method_snapshot=dict(method_snapshot),
            execution_overrides=dict(execution_overrides),
            dataset_ref=_copy_optional_dict(dataset_ref),
            model_ref=_copy_optional_dict(model_ref),
            molecule_count=molecule_count,
            success_count=0,
            failure_count=0,
            output_paths={},
            error=None,
        )
        self.store.write(manifest)
        return manifest

    def start_run(self, run_id: str) -> RunManifest:
        """Mark a created run as running."""
        manifest = self._load_mutable(run_id)
        started_at = manifest.started_at or _utc_now()
        return self._write(
            manifest.with_updates(
                status=RunStatus.RUNNING,
                started_at=started_at,
                error=None,
            )
        )

    def attach_output_paths(
        self,
        run_id: str,
        output_paths: Mapping[str, Path | str],
    ) -> RunManifest:
        """Attach or replace artifact paths before terminal finalization."""
        manifest = self._load_mutable(run_id)
        merged = dict(manifest.output_paths)
        merged.update({key: Path(value) for key, value in output_paths.items()})
        return self._write(manifest.with_updates(output_paths=merged))

    def complete_run(
        self,
        run_id: str,
        *,
        success_count: int,
        failure_count: int,
    ) -> RunManifest:
        """Finalize a run as completed or partial after validating outputs."""
        manifest = self._load_mutable(run_id)
        self._validate_completion_outputs(manifest)
        status = RunStatus.COMPLETED if failure_count == 0 else RunStatus.PARTIAL
        return self._write(
            manifest.with_updates(
                status=status,
                completed_at=_utc_now(),
                success_count=success_count,
                failure_count=failure_count,
                error=None,
            )
        )

    def fail_run(
        self,
        run_id: str,
        *,
        error: str,
        success_count: int = 0,
        failure_count: int = 0,
    ) -> RunManifest:
        """Finalize a run as failed."""
        manifest = self._load_mutable(run_id)
        return self._write(
            manifest.with_updates(
                status=RunStatus.FAILED,
                completed_at=_utc_now(),
                success_count=success_count,
                failure_count=failure_count,
                error=error,
            )
        )

    def invalidate_run(
        self,
        run_id: str,
        *,
        error: str,
        success_count: int = 0,
        failure_count: int = 0,
    ) -> RunManifest:
        """Finalize an ALL_OR_NOTHING run as invalidated without deleting artifacts."""
        manifest = self._load_mutable(run_id)
        return self._write(
            manifest.with_updates(
                status=RunStatus.INVALIDATED,
                completed_at=_utc_now(),
                success_count=success_count,
                failure_count=failure_count,
                error=error,
            )
        )

    def cancel_run(self, run_id: str, *, error: str | None = None) -> RunManifest:
        """Finalize a run as cancelled."""
        manifest = self._load_mutable(run_id)
        return self._write(
            manifest.with_updates(
                status=RunStatus.CANCELLED,
                completed_at=_utc_now(),
                error=error,
            )
        )

    def list_runs(self) -> list[RunManifest]:
        """List all persisted runs."""
        return self.store.list()

    def get_run(self, run_id: str) -> RunManifest:
        """Return one run manifest."""
        return self.store.read(run_id)

    def _load_mutable(self, run_id: str) -> RunManifest:
        manifest = self.store.read(run_id)
        if manifest.status in TERMINAL_STATUSES:
            raise ValueError(
                f"Run {run_id} is already terminal: {manifest.status.value}"
            )
        return manifest

    def _write(self, manifest: RunManifest) -> RunManifest:
        self.store.write(manifest)
        return manifest

    def _validate_completion_outputs(self, manifest: RunManifest) -> None:
        if not manifest.output_paths:
            raise ValueError("completed run requires output paths")
        missing = [
            str(path)
            for path in manifest.output_paths.values()
            if not Path(path).exists()
        ]
        if missing:
            raise ValueError(
                "completed run requires existing output paths: " + ", ".join(missing)
            )

    @staticmethod
    def _new_run_id() -> str:
        timestamp = _utc_now().strftime("%Y%m%dT%H%M%SZ")
        return f"run_{timestamp}_{uuid.uuid4().hex[:8]}"


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _copy_optional_dict(value: dict[str, Any] | None) -> dict[str, Any] | None:
    return dict(value) if value is not None else None
