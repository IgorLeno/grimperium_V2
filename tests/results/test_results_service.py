from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from grimperium.results.service import ResultsService
from grimperium.runs.service import RunService


def test_analyze_dataset_returns_source_and_metrics(tmp_path: Path) -> None:
    csv_path = tmp_path / "predictions.csv"
    pd.DataFrame(
        {
            "mol_id": ["m1", "m2", "m3"],
            "H298_cbs": [-100.0, -200.0, -300.0],
            "H298_predicted": [-101.0, -198.0, -303.0],
        }
    ).to_csv(csv_path, index=False)

    report = ResultsService().analyze_dataset(csv_path)

    assert report.source_path == csv_path
    assert report.source_label == "predictions"
    assert report.summary["n_molecules"] == 3
    assert report.summary["mae"] == pytest.approx(2.0)
    assert sum(item.count for item in report.divergence_distribution) == 3


def test_analyze_run_uses_manifest_output_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    manifest = run_service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="pm7_delta_learning",
        method_version="0.1.0",
        method_snapshot={"display_name": "PM7 + Delta"},
        execution_overrides={},
        dataset_ref=None,
        model_ref={"name": "Model A", "path": "models/a.joblib"},
        molecule_count=2,
    )
    output_path = run_service.runs_root / manifest.run_id / "predictions.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "mol_id": ["m1", "m2"],
            "H298_cbs": [-100.0, -200.0],
            "H298_predicted": [-101.0, -198.0],
        }
    ).to_csv(output_path, index=False)
    run_service.start_run(manifest.run_id)
    run_service.attach_output_paths(manifest.run_id, {"results_csv": output_path})
    manifest = run_service.complete_run(
        manifest.run_id,
        success_count=2,
        failure_count=0,
    )

    report = ResultsService(run_service=run_service).analyze_run(manifest.run_id)

    assert report.run_id == manifest.run_id
    assert report.model_label == "Model A"
    assert report.method_label == "PM7 + Delta"
    assert report.summary["n_molecules"] == 2


def test_compare_runs_returns_basic_metric_rows(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    run_ids: list[str] = []
    for idx, err in enumerate((1.0, 3.0), start=1):
        manifest = run_service.create_run(
            property_id="standard_enthalpy_of_formation",
            method_id="pm7_delta_learning",
            method_version="0.1.0",
            method_snapshot={"display_name": f"Method {idx}"},
            execution_overrides={},
            dataset_ref=None,
            model_ref={"name": f"Model {idx}"},
            molecule_count=1,
        )
        output_path = run_service.runs_root / manifest.run_id / "predictions.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {"mol_id": ["m1"], "H298_cbs": [-100.0], "H298_predicted": [-100.0 - err]}
        ).to_csv(output_path, index=False)
        run_service.start_run(manifest.run_id)
        run_service.attach_output_paths(manifest.run_id, {"results_csv": output_path})
        run_service.complete_run(manifest.run_id, success_count=1, failure_count=0)
        run_ids.append(manifest.run_id)

    comparison = ResultsService(run_service=run_service).compare_runs(run_ids)

    assert [row.run_id for row in comparison.rows] == run_ids
    assert [row.mae for row in comparison.rows] == [1.0, 3.0]
    assert comparison.best_run_id == run_ids[0]
