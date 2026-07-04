"""Tests for the canonical PM7 + Delta Learning runner (Method B)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from grimperium.calculation.contracts.enums import OverallStatus, PropertyRole
from grimperium.calculation.contracts.models import MoleculeCalculationResult
from grimperium.calculation.runners.pm7_delta_runner import (
    PM7DeltaLearningRunner,
    PM7DeltaRunnerError,
)
from grimperium.cli.model_compatibility import ModelCompatibilityError
from grimperium.crest_pm7.config import CRESTStatus, MOPACStatus
from grimperium.crest_pm7.molecule_processor import ConformerData, PM7Result

H298_PM7 = -42.0
H298_FINAL = -50.0
EXPECTED_DELTA = H298_FINAL - H298_PM7  # -8.0


def _pm7_result(*, success: bool = True) -> PM7Result:
    selected = ConformerData(index=0, mol_id="mol-b", crest_rank=2)
    selected.crest_status = CRESTStatus.SUCCESS
    selected.crest_geometry_file = Path("work/mol-b/conf_2.xyz")
    selected.mopac_status = MOPACStatus.SUCCESS
    selected.mopac_output_file = Path("work/mol-b/conf_2.out")
    selected.energy_hof = H298_PM7
    selected.hof_extraction_successful = True

    return PM7Result(
        mol_id="mol-b",
        smiles="CCO",
        phase="B",
        nheavy=3,
        rdkit_descriptors={},
        crest_status=CRESTStatus.SUCCESS,
        crest_conformers_generated=3,
        crest_time=7.5,
        conformers=[selected],
        num_conformers_selected=1,
        k_selected_pm7=2,
        total_execution_time=31.0,
        success=success,
        error_message=None if success else "CREST failed",
    )


class _FakeLearner:
    def __init__(self, prediction: float = H298_FINAL, *, raises: bool = False) -> None:
        self._prediction = prediction
        self._raises = raises

    def predict(self, X: Any, y_pm7: Any) -> np.ndarray:
        if self._raises:
            raise RuntimeError("inference boom")
        return np.array([self._prediction])


class _FakePipeline:
    feature_cols = ["H298_pm7", "nheavy"]

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        return df


def _metadata() -> dict[str, Any]:
    return {
        "version": "1.2.3",
        "property_id": "standard_enthalpy_of_formation",
        "baseline_hamiltonian": "PM7",
        "feature_schema_id": "grimperium_dhf_v1",
        "feature_schema_hash": "deadbeef",
        "feature_columns": ["H298_pm7", "nheavy"],
    }


def _make_runner(
    *,
    tmp_path: Path,
    pm7_result: PM7Result | None = None,
    learner: _FakeLearner | None = None,
    feature_builder: Any = None,
    validator: Any = None,
) -> tuple[PM7DeltaLearningRunner, Path]:
    result = pm7_result if pm7_result is not None else _pm7_result()
    model_path = tmp_path / "model.pkl"
    model_path.write_text("stub", encoding="utf-8")

    runner = PM7DeltaLearningRunner(
        pm7_processor=lambda mol_id, smiles: result,
        model_loader=lambda path: (learner or _FakeLearner(), _FakePipeline()),
        metadata_loader=lambda path: _metadata(),
        feature_builder=feature_builder
        or (lambda smiles, res: pd.DataFrame([{"H298_pm7": H298_PM7, "nheavy": 3}])),
        model_validator=validator or (lambda path, method: None),
    )
    return runner, model_path


def _run(runner: PM7DeltaLearningRunner, model_path: Path) -> MoleculeCalculationResult:
    return runner.calculate_single_smiles(
        "CCO", molecule_id="mol-b", model_path=model_path
    )


def _estimate(result: MoleculeCalculationResult, role: PropertyRole) -> Any:
    matches = [e for e in result.estimates if e.role is role]
    assert len(matches) == 1, f"expected exactly one {role} estimate"
    return matches[0]


def test_runner_returns_canonical_success(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    result = _run(runner, model_path)
    assert isinstance(result, MoleculeCalculationResult)
    assert result.overall_status is OverallStatus.SUCCESS


def test_runner_method_id_and_version(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    result = _run(runner, model_path)
    assert result.run.method_ref.method_id == "pm7_delta_learning"
    assert result.run.method_ref.method_version == "0.1.0"
    final = _estimate(result, PropertyRole.FINAL)
    assert final.method_id == "pm7_delta_learning"
    assert final.method_version == "0.1.0"


def test_runner_baseline_estimate(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    baseline = _estimate(_run(runner, model_path), PropertyRole.BASELINE)
    assert baseline.hamiltonian == "PM7"
    assert baseline.value_kcal_mol == H298_PM7
    assert baseline.model_path is None


def test_runner_correction_estimate_is_delta(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    correction = _estimate(_run(runner, model_path), PropertyRole.CORRECTION)
    assert correction.value_kcal_mol == pytest.approx(EXPECTED_DELTA)
    assert correction.model_path is not None


def test_runner_final_estimate(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    final = _estimate(_run(runner, model_path), PropertyRole.FINAL)
    assert final.value_kcal_mol == pytest.approx(H298_FINAL)
    assert final.model_path is not None


def test_runner_delta_equals_final_minus_baseline(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    result = _run(runner, model_path)
    baseline = _estimate(result, PropertyRole.BASELINE).value_kcal_mol
    correction = _estimate(result, PropertyRole.CORRECTION).value_kcal_mol
    final = _estimate(result, PropertyRole.FINAL).value_kcal_mol
    assert correction == pytest.approx(final - baseline)


def test_runner_model_path_is_portable(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    correction = _estimate(_run(runner, model_path), PropertyRole.CORRECTION)
    assert not Path(correction.model_path).is_absolute()


def test_runner_records_stage_executions_with_times(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    result = _run(runner, model_path)
    stage_ids = {stage.stage_id for stage in result.stage_executions}
    assert "feature_assembly" in stage_ids
    assert "delta_learning_inference" in stage_ids
    inference = next(
        s for s in result.stage_executions if s.stage_id == "delta_learning_inference"
    )
    assert inference.execution_time_s is not None
    assert inference.settings.get("feature_schema_id") == "grimperium_dhf_v1"


def test_runner_preserves_selected_conformer(tmp_path: Path) -> None:
    runner, model_path = _make_runner(tmp_path=tmp_path)
    result = _run(runner, model_path)
    assert result.conformers  # conformer data preserved from the PM7 pipeline
    assert _estimate(result, PropertyRole.BASELINE).conformer_source_id == 2


def test_runner_raises_on_incompatible_model(tmp_path: Path) -> None:
    def _bad(path: Path, method: Any) -> None:
        raise ModelCompatibilityError("nope")

    runner, model_path = _make_runner(tmp_path=tmp_path, validator=_bad)
    with pytest.raises(ModelCompatibilityError):
        _run(runner, model_path)


def test_runner_raises_on_crest_failure(tmp_path: Path) -> None:
    runner, model_path = _make_runner(
        tmp_path=tmp_path, pm7_result=_pm7_result(success=False)
    )
    with pytest.raises(PM7DeltaRunnerError):
        _run(runner, model_path)


def test_runner_raises_on_feature_failure(tmp_path: Path) -> None:
    def _boom(smiles: str, res: PM7Result) -> pd.DataFrame:
        raise ValueError("bad features")

    runner, model_path = _make_runner(tmp_path=tmp_path, feature_builder=_boom)
    with pytest.raises(PM7DeltaRunnerError):
        _run(runner, model_path)


def test_runner_raises_on_inference_failure(tmp_path: Path) -> None:
    runner, model_path = _make_runner(
        tmp_path=tmp_path, learner=_FakeLearner(raises=True)
    )
    with pytest.raises(PM7DeltaRunnerError):
        _run(runner, model_path)
