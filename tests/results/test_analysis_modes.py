"""Unit tests for Results analysis capability modes."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from grimperium.results.models import ResultsAnalysisMode
from grimperium.results.service import ResultsService
from grimperium.runs.service import RunService


def _write_canonical_baseline(path: Path, *, mol_id: str, value: float) -> None:
    pd.DataFrame(
        [
            {
                "estimate_id": f"{mol_id}-baseline",
                "run_id": "run-x",
                "molecule_smiles": "CCO",
                "molecule_name": mol_id,
                "method_id": "crest_pm7",
                "method_version": "1.0",
                "property_id": "standard_enthalpy_of_formation",
                "role": "baseline",
                "hamiltonian": "PM7",
                "canonical_value": value,
                "canonical_unit": "kcal/mol",
                "value_kcal_mol": value,
                "value_kj_mol": "",
                "conformer_source_id": 1,
                "overall_status": "success",
                "execution_phase": "A",
                "started_at": "",
                "completed_at": "",
                "schema_version": 1,
            }
        ]
    ).to_csv(path, index=False)


def _complete_run_with_csv(
    run_service: RunService,
    *,
    csv_path: Path,
    dataset_ref: dict[str, object] | None = None,
    method_id: str = "crest_pm7",
) -> str:
    manifest = run_service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id=method_id,
        method_version="1.0",
        method_snapshot={"display_name": "CREST + MOPAC PM7"},
        execution_overrides={},
        dataset_ref=dataset_ref,
        model_ref=None,
        molecule_count=1,
    )
    run_service.start_run(manifest.run_id)
    run_service.attach_output_paths(
        manifest.run_id, {"calculation_results_csv": csv_path}
    )
    run_service.complete_run(manifest.run_id, success_count=1, failure_count=0)
    return manifest.run_id


def test_pm7_only_run_uses_scientific_summary_mode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    csv_path = tmp_path / "calculation_results.csv"
    _write_canonical_baseline(csv_path, mol_id="mol_ok", value=-55.0)
    run_id = _complete_run_with_csv(run_service, csv_path=csv_path)

    report = ResultsService(run_service=run_service).analyze_run(run_id)

    assert report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY
    assert report.analysis is None
    assert report.scientific_summary.n_molecules == 1
    assert "baseline" in report.scientific_summary.roles
    assert report.summary.get("mae") is None


def test_baseline_with_reference_from_dataset_ref(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    csv_path = tmp_path / "calculation_results.csv"
    _write_canonical_baseline(csv_path, mol_id="mol_ok", value=-55.0)
    ref_path = tmp_path / "reference.csv"
    pd.DataFrame({"mol_id": ["mol_ok"], "smiles": ["CCO"], "H298_cbs": [-50.0]}).to_csv(
        ref_path, index=False
    )
    original = csv_path.read_text(encoding="utf-8")
    run_id = _complete_run_with_csv(
        run_service,
        csv_path=csv_path,
        dataset_ref={"path": str(ref_path)},
    )

    report = ResultsService(run_service=run_service).analyze_run(run_id)

    assert report.analysis_mode is ResultsAnalysisMode.BASELINE_WITH_REFERENCE
    assert report.analysis is not None
    assert report.summary["mae"] == pytest.approx(5.0)
    assert report.scientific_summary.comparison_label == "PM7 baseline vs reference"
    assert csv_path.read_text(encoding="utf-8") == original
    assert "H298_predicted" not in original
    assert "final" not in original.lower()


def test_prediction_with_reference_legacy_wide(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    csv_path = tmp_path / "predictions.csv"
    pd.DataFrame(
        {
            "mol_id": ["m1"],
            "H298_cbs": [-100.0],
            "H298_predicted": [-101.0],
        }
    ).to_csv(csv_path, index=False)
    run_id = _complete_run_with_csv(
        run_service, csv_path=csv_path, method_id="pm7_delta_learning"
    )

    report = ResultsService(run_service=run_service).analyze_run(run_id)

    assert report.analysis_mode is ResultsAnalysisMode.PREDICTION_WITH_REFERENCE
    assert report.analysis is not None
    assert report.summary["mae"] == pytest.approx(1.0)
