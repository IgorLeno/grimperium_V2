from __future__ import annotations

from pathlib import Path

import pytest

from grimperium.runs.models import RunStatus
from grimperium.runs.service import RunService


def _service(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> RunService:
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    return RunService.from_environment()


def test_run_service_creates_starts_and_completes_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="pm7_delta_learning",
        method_version="0.1.0",
        method_snapshot={"method_id": "pm7_delta_learning"},
        execution_overrides={"n_conformers": 5},
        dataset_ref={"path": "data/input.csv"},
        model_ref={"path": "models/model.joblib"},
        molecule_count=1,
    )
    output_path = service.runs_root / manifest.run_id / "calculation_results.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("mol_id,H298_predicted\nm1,-10\n", encoding="utf-8")

    service.start_run(manifest.run_id)
    service.attach_output_paths(manifest.run_id, {"results_csv": output_path})
    completed = service.complete_run(
        manifest.run_id,
        success_count=1,
        failure_count=0,
    )

    assert completed.status is RunStatus.COMPLETED
    assert completed.success_count == 1
    assert completed.output_paths["results_csv"] == output_path
    assert service.get_run(manifest.run_id).status is RunStatus.COMPLETED
    assert [run.run_id for run in service.list_runs()] == [manifest.run_id]


def test_complete_run_requires_existing_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="pm7_delta_learning",
        method_version="0.1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    service.start_run(manifest.run_id)

    with pytest.raises(ValueError, match="completed run requires output paths"):
        service.complete_run(manifest.run_id, success_count=1, failure_count=0)


def test_invalidated_run_keeps_artifact_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0.0",
        method_snapshot={"method_id": "crest_pm7"},
        execution_overrides={"batch_size": 2},
        dataset_ref={"path": "data/thermo_pm7.csv"},
        model_ref=None,
        molecule_count=2,
    )
    output_path = service.runs_root / manifest.run_id / "calculation_results.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("mol_id,H298_pm7\nm1,-10\n", encoding="utf-8")

    service.start_run(manifest.run_id)
    service.attach_output_paths(manifest.run_id, {"results_csv": output_path})
    invalidated = service.invalidate_run(
        manifest.run_id,
        error="ALL_OR_NOTHING reset invalidated scientific completion",
        success_count=1,
        failure_count=1,
    )

    assert invalidated.status is RunStatus.INVALIDATED
    assert invalidated.output_paths["results_csv"] == output_path
    assert invalidated.error == "ALL_OR_NOTHING reset invalidated scientific completion"


def test_invalid_transition_created_to_completed_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    output_path = service.runs_root / manifest.run_id / "out.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("x\n", encoding="utf-8")
    service.attach_output_paths(manifest.run_id, {"results_csv": output_path})

    with pytest.raises(ValueError, match="Invalid run transition"):
        service.complete_run(manifest.run_id, success_count=1, failure_count=0)


def test_partial_requires_success_and_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=2,
    )
    output_path = service.runs_root / manifest.run_id / "out.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("x\n", encoding="utf-8")
    service.start_run(manifest.run_id)
    service.attach_output_paths(manifest.run_id, {"results_csv": output_path})

    with pytest.raises(ValueError, match="partial runs require"):
        service.complete_run(manifest.run_id, success_count=0, failure_count=2)


def test_portable_relative_output_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    output_path = service.runs_root / manifest.run_id / "calculation_results.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("mol_id\nm1\n", encoding="utf-8")
    service.start_run(manifest.run_id)
    service.attach_output_paths(
        manifest.run_id, {"calculation_results_csv": output_path}
    )
    service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    # Move runs root content conceptually by reloading with same relative layout.
    reloaded = RunService(service.runs_root).get_run(manifest.run_id)
    assert reloaded.output_paths["calculation_results_csv"].exists()
    assert "calculation_results.csv" in str(
        reloaded.output_paths["calculation_results_csv"]
    )


def test_invalid_transition_running_to_running_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    service.start_run(manifest.run_id)
    with pytest.raises(ValueError, match="Invalid run transition"):
        service.start_run(manifest.run_id)


def test_invalid_transition_created_to_invalidated_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    with pytest.raises(ValueError, match="Invalid run transition"):
        service.invalidate_run(manifest.run_id, error="nope")


def test_counts_exceeding_molecule_count_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    output_path = service.runs_root / manifest.run_id / "out.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("x\n", encoding="utf-8")
    service.start_run(manifest.run_id)
    service.attach_output_paths(manifest.run_id, {"results_csv": output_path})
    with pytest.raises(ValueError, match="molecule_count"):
        service.complete_run(manifest.run_id, success_count=2, failure_count=0)


def test_completed_run_requires_counts_match_molecule_count(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=2,
    )
    output_path = service.runs_root / manifest.run_id / "out.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("x\n", encoding="utf-8")
    service.start_run(manifest.run_id)
    service.attach_output_paths(manifest.run_id, {"results_csv": output_path})

    with pytest.raises(ValueError, match="must equal molecule_count"):
        service.complete_run(manifest.run_id, success_count=1, failure_count=0)


def test_partial_run_requires_counts_match_molecule_count(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=3,
    )
    output_path = service.runs_root / manifest.run_id / "out.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("x\n", encoding="utf-8")
    service.start_run(manifest.run_id)
    service.attach_output_paths(manifest.run_id, {"results_csv": output_path})

    with pytest.raises(ValueError, match="must equal molecule_count"):
        service.complete_run(manifest.run_id, success_count=1, failure_count=1)


def test_portable_paths_survive_runs_root_relocation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _service(tmp_path, monkeypatch)
    manifest = service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    output_path = service.runs_root / manifest.run_id / "calculation_results.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("mol_id\nm1\n", encoding="utf-8")
    service.start_run(manifest.run_id)
    service.attach_output_paths(
        manifest.run_id, {"calculation_results_csv": output_path}
    )
    service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    relocated = tmp_path / "moved_runs"
    service.runs_root.rename(relocated)
    reloaded = RunService(relocated).get_run(manifest.run_id)
    assert reloaded.output_paths["calculation_results_csv"].exists()
    assert reloaded.output_paths["calculation_results_csv"].is_relative_to(relocated)
