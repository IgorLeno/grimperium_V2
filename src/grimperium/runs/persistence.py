"""Persistence helpers for run manifests."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

from grimperium.runs.models import RunManifest

MANIFEST_FILENAME = "manifest.json"


class RunManifestStore:
    """JSON store with atomic writes and run-root-relative artifact paths."""

    def __init__(self, runs_root: Path) -> None:
        self.runs_root = Path(runs_root)

    def manifest_path(self, run_id: str) -> Path:
        """Return the manifest path for a run ID."""
        return self.runs_root / run_id / MANIFEST_FILENAME

    def write(self, manifest: RunManifest) -> None:
        """Atomically write a manifest to disk."""
        run_dir = self.runs_root / manifest.run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        payload = self._manifest_payload(manifest)
        content = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        temp_path: Path | None = None

        try:
            with tempfile.NamedTemporaryFile(
                "w",
                encoding="utf-8",
                dir=run_dir,
                prefix=".manifest.",
                suffix=".tmp",
                delete=False,
            ) as handle:
                temp_path = Path(handle.name)
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
            temp_path.replace(self.manifest_path(manifest.run_id))
        finally:
            if temp_path is not None and temp_path.exists():
                temp_path.unlink()

    def read(self, run_id: str) -> RunManifest:
        """Read a single run manifest."""
        path = self.manifest_path(run_id)
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError(f"Run manifest is not a JSON object: {path}")
        manifest = RunManifest.from_dict(payload)
        return manifest.with_updates(
            output_paths={
                key: self._restore_path(value)
                for key, value in manifest.output_paths.items()
            }
        )

    def list(self) -> list[RunManifest]:
        """Return all readable manifests sorted by creation time descending."""
        if not self.runs_root.exists():
            return []
        manifests: list[RunManifest] = []
        for path in self.runs_root.glob(f"*/{MANIFEST_FILENAME}"):
            manifests.append(self.read(path.parent.name))
        return sorted(manifests, key=lambda item: item.created_at, reverse=True)

    def _manifest_payload(self, manifest: RunManifest) -> dict[str, Any]:
        payload = manifest.to_dict()
        payload["output_paths"] = {
            key: str(self._relative_to_runs_root(path))
            for key, path in manifest.output_paths.items()
        }
        return payload

    def _relative_to_runs_root(self, path: Path) -> Path:
        path = Path(path)
        try:
            return path.resolve(strict=False).relative_to(
                self.runs_root.resolve(strict=False)
            )
        except ValueError:
            return path

    def _restore_path(self, path: Path) -> Path:
        if path.is_absolute():
            return path
        return self.runs_root / path
