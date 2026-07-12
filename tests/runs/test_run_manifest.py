from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from grimperium.runs.models import RunManifest, RunStatus
from grimperium.runs.persistence import RunManifestStore


def _manifest(run_id: str, output_path: Path) -> RunManifest:
    return RunManifest(
        schema_version="1.0",
        run_id=run_id,
        status=RunStatus.COMPLETED,
        created_at=datetime.now(timezone.utc),
        started_at=datetime.now(timezone.utc),
        completed_at=datetime.now(timezone.utc),
        property_id="standard_enthalpy_of_formation",
        method_id="pm7_delta_learning",
        method_version="0.1.0",
        method_snapshot={"method_id": "pm7_delta_learning"},
        execution_overrides={"batch_size": 3},
        dataset_ref={"path": "data/input.csv"},
        model_ref={"path": "models/model.joblib"},
        molecule_count=1,
        success_count=1,
        failure_count=0,
        output_paths={"results_csv": output_path},
        error=None,
    )


def test_manifest_store_writes_relative_paths_atomically(tmp_path: Path) -> None:
    runs_root = tmp_path / "runs"
    output_path = runs_root / "run_abc" / "calculation_results.csv"
    output_path.parent.mkdir(parents=True)
    output_path.write_text("mol_id,H298_predicted\nm1,-10\n", encoding="utf-8")
    store = RunManifestStore(runs_root)

    store.write(_manifest("run_abc", output_path))

    manifest_path = runs_root / "run_abc" / "manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["output_paths"] == {"results_csv": "run_abc/calculation_results.csv"}
    assert list(runs_root.rglob("*.tmp")) == []


def test_manifest_store_round_trips_datetimes_and_paths(tmp_path: Path) -> None:
    runs_root = tmp_path / "runs"
    output_path = runs_root / "run_xyz" / "calculation_results.csv"
    output_path.parent.mkdir(parents=True)
    output_path.write_text("mol_id,H298_predicted\nm1,-10\n", encoding="utf-8")
    store = RunManifestStore(runs_root)
    manifest = _manifest("run_xyz", output_path)

    store.write(manifest)
    loaded = store.read("run_xyz")

    assert loaded.run_id == "run_xyz"
    assert loaded.status is RunStatus.COMPLETED
    assert loaded.created_at == manifest.created_at
    assert loaded.output_paths["results_csv"] == output_path
