"""End-to-end: real batch/calc outputs → Run → ResultsService.analyze_run."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from grimperium.calculation.contracts.enums import OverallStatus, PropertyRole
from grimperium.calculation.contracts.models import (
    CalculationMethodReference,
    MoleculeCalculationResult,
    MoleculeData,
    PropertyEstimate,
    RunMetadata,
)
from grimperium.calculation.contracts.quantity import Quantity
from grimperium.calculation.output.csv_writer import write_canonical_csv
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager
from grimperium.crest_pm7.batch.output_contracts import BatchOutputLayout
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import CRESTStatus, MOPACStatus, PM7Config
from grimperium.crest_pm7.csv_enhancements import CSVManagerExtensions
from grimperium.crest_pm7.molecule_processor import ConformerData, PM7Result
from grimperium.results.models import ResultsAnalysisMode
from grimperium.results.service import ResultsService
from grimperium.runs.service import RunService


class _FakeProcessor:
    def update_timeouts(self, **_kwargs: object) -> None:
        pass

    def process_with_fixed_timeout(
        self,
        mol_id: str,
        smiles: str,
        input_xyz: Path | None = None,
        progress_callback: object = None,
        charge: int = 0,
        multiplicity: int = 1,
    ) -> PM7Result:
        conformer = ConformerData(index=0, mol_id=mol_id, crest_rank=1)
        conformer.crest_status = CRESTStatus.SUCCESS
        conformer.mopac_status = MOPACStatus.SUCCESS
        conformer.energy_hof = -55.0
        conformer.hof_extraction_successful = True
        return PM7Result(
            mol_id=mol_id,
            smiles=smiles,
            phase="A",
            nheavy=3,
            rdkit_descriptors={},
            crest_status=CRESTStatus.SUCCESS,
            crest_conformers_generated=1,
            conformers=[conformer],
            num_conformers_selected=1,
            k_selected_pm7=1,
            total_execution_time=1.0,
            success=True,
        )


def _canonical_result(
    *,
    smiles: str,
    name: str,
    estimates: list[PropertyEstimate],
    method_id: str = "crest_pm7",
) -> MoleculeCalculationResult:
    return MoleculeCalculationResult(
        molecule=MoleculeData(smiles=smiles, name=name),
        run=RunMetadata(
            run_id="run-e2e",
            execution_phase="A",
            method_ref=CalculationMethodReference(
                method_id=method_id,
                method_version="1.0",
                property_id="standard_enthalpy_of_formation",
            ),
            started_at=datetime(2026, 7, 12, 10, 0, tzinfo=timezone.utc),
            completed_at=datetime(2026, 7, 12, 10, 1, tzinfo=timezone.utc),
            grimperium_version=None,
        ),
        overall_status=OverallStatus.SUCCESS,
        conformers=[],
        molecular_descriptors=None,
        estimates=estimates,
        artifacts=[],
        stage_executions=[],
    )


def test_batch_pm7_real_output_opens_in_results_as_baseline_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        [
            {
                "mol_id": "mol_ok",
                "smiles": "CCO",
                "nheavy": 3,
                "charge": 0,
                "multiplicity": 1,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            }
        ]
    ).to_csv(csv_path, index=False)
    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()
    state_manager = BatchStateManager(tmp_path / "batch_state.csv", PM7Config())
    state_manager.reconcile_molecules(csv_manager.state_seed_rows())
    layout = BatchOutputLayout(tmp_path / "batch_outputs")
    batch = csv_manager.select_batch(
        batch_id="batch_0001",
        batch_size=1,
        crest_timeout_minutes=30,
        mopac_timeout_minutes=60,
    )
    state_manager.mark_selected_from_batch(batch)
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=type(
            "DetailStub",
            (),
            {
                "pm7result_to_detail": lambda self, **kwargs: {},
                "save_detail": lambda self, detail: None,
            },
        )(),
        pm7_config=PM7Config(temp_dir=tmp_path),
        processor_adapter=_FakeProcessor(),
        output_layout=layout,
        write_xlsx=False,
    )
    with patch.object(
        CSVManagerExtensions, "update_molecule_with_mopac_results", return_value=True
    ):
        manager.execute_batch(batch)

    assert layout.calculation_results_csv.exists()
    canonical = layout.calculation_results_csv.read_text(encoding="utf-8")
    assert "baseline" in canonical.lower()
    assert "correction" not in {
        r.lower() for r in pd.read_csv(layout.calculation_results_csv)["role"]
    }
    assert "final" not in {
        r.lower() for r in pd.read_csv(layout.calculation_results_csv)["role"]
    }

    run_service = RunService.from_environment()
    manifest = run_service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={"display_name": "CREST + MOPAC PM7"},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    run_dir = run_service.runs_root / manifest.run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    run_csv = run_dir / "calculation_results.csv"
    run_csv.write_text(canonical, encoding="utf-8")
    run_service.start_run(manifest.run_id)
    run_service.attach_output_paths(
        manifest.run_id, {"calculation_results_csv": run_csv}
    )
    run_service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    report = ResultsService(run_service=run_service).analyze_run(manifest.run_id)

    assert report.run_id == manifest.run_id
    assert report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY
    assert report.analysis is None
    assert report.scientific_summary.n_molecules >= 1
    assert "baseline" in report.scientific_summary.roles
    assert report.summary.get("mae") is None


def test_pm7_with_reference_produces_baseline_metrics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    run_dir = tmp_path / "runs" / "placeholder"
    run_dir.mkdir(parents=True)
    csv_path = tmp_path / "calculation_results.csv"
    write_canonical_csv(
        [
            _canonical_result(
                smiles="CCO",
                name="mol_ok",
                estimates=[
                    PropertyEstimate(
                        estimate_id="e1",
                        property_id="standard_enthalpy_of_formation",
                        role=PropertyRole.BASELINE,
                        method_id="crest_pm7",
                        method_version="1.0",
                        hamiltonian="PM7",
                        value=Quantity(value=-55.0, unit="kcal/mol"),
                        value_kcal_mol=None,
                        value_kj_mol=None,
                        conformer_source_id=1,
                        uncertainty=None,
                        model_path=None,
                    )
                ],
            )
        ],
        csv_path,
    )
    original = csv_path.read_text(encoding="utf-8")
    ref_path = tmp_path / "ref.csv"
    pd.DataFrame({"mol_id": ["mol_ok"], "smiles": ["CCO"], "H298_cbs": [-50.0]}).to_csv(
        ref_path, index=False
    )

    manifest = run_service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={"display_name": "CREST + MOPAC PM7"},
        execution_overrides={},
        dataset_ref={"path": str(ref_path)},
        model_ref=None,
        molecule_count=1,
    )
    run_service.start_run(manifest.run_id)
    run_service.attach_output_paths(
        manifest.run_id, {"calculation_results_csv": csv_path}
    )
    run_service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    report = ResultsService(run_service=run_service).analyze_run(manifest.run_id)

    assert report.analysis_mode is ResultsAnalysisMode.BASELINE_WITH_REFERENCE
    assert report.analysis is not None
    assert report.summary["mae"] == pytest.approx(5.0)
    assert report.summary["rmse"] == pytest.approx(5.0)
    assert report.summary["bias"] is not None
    assert report.summary["r2"] is not None
    assert report.scientific_summary.comparison_label == "PM7 baseline vs reference"
    assert csv_path.read_text(encoding="utf-8") == original
    canonical = pd.read_csv(csv_path)
    assert "FINAL" not in set(canonical["role"].astype(str).str.upper())


def test_pm7_without_reference_scientific_summary_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    csv_path = tmp_path / "calculation_results.csv"
    write_canonical_csv(
        [
            _canonical_result(
                smiles="CCO",
                name="mol_ok",
                estimates=[
                    PropertyEstimate(
                        estimate_id="e1",
                        property_id="standard_enthalpy_of_formation",
                        role=PropertyRole.BASELINE,
                        method_id="crest_pm7",
                        method_version="1.0",
                        hamiltonian="PM7",
                        value=Quantity(value=-55.0, unit="kcal/mol"),
                        value_kcal_mol=None,
                        value_kj_mol=None,
                        conformer_source_id=1,
                        uncertainty=None,
                        model_path=None,
                    )
                ],
            )
        ],
        csv_path,
    )
    manifest = run_service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="crest_pm7",
        method_version="1.0",
        method_snapshot={"display_name": "CREST + MOPAC PM7"},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    run_service.start_run(manifest.run_id)
    run_service.attach_output_paths(
        manifest.run_id, {"calculation_results_csv": csv_path}
    )
    run_service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    report = ResultsService(run_service=run_service).analyze_run(manifest.run_id)

    assert report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY
    assert report.analysis is None
    assert report.scientific_summary.n_molecules == 1
    assert report.scientific_summary.n_estimates == 1
    assert report.summary.get("mae") is None
    assert report.summary.get("rmse") is None


def test_method_a_canonical_csv_opens_in_results(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    csv_path = tmp_path / "calculation_results.csv"
    write_canonical_csv(
        [
            _canonical_result(
                smiles="CCO",
                name="ethanol",
                method_id="semiempirical_am1_pm3_pm7",
                estimates=[
                    PropertyEstimate(
                        estimate_id="e1",
                        property_id="standard_enthalpy_of_formation",
                        role=PropertyRole.FINAL,
                        method_id="semiempirical_am1_pm3_pm7",
                        method_version="1.0",
                        hamiltonian="PM7",
                        value=Quantity(value=-56.0, unit="kcal/mol"),
                        value_kcal_mol=None,
                        value_kj_mol=None,
                        conformer_source_id=1,
                        uncertainty=None,
                        model_path=None,
                    )
                ],
            )
        ],
        csv_path,
    )
    manifest = run_service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="semiempirical_am1_pm3_pm7",
        method_version="1.0",
        method_snapshot={"display_name": "Semiempirical"},
        execution_overrides={},
        dataset_ref=None,
        model_ref=None,
        molecule_count=1,
    )
    run_service.start_run(manifest.run_id)
    compat = tmp_path / "single_result.json"
    compat.write_text(json.dumps({"smiles": "CCO"}), encoding="utf-8")
    run_service.attach_output_paths(
        manifest.run_id,
        {
            "calculation_results_csv": csv_path,
            "single_result": compat,
        },
    )
    run_service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    report = ResultsService(run_service=run_service).analyze_run(manifest.run_id)
    reloaded = run_service.get_run(manifest.run_id)

    assert "calculation_results_csv" in reloaded.output_paths
    assert report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY
    assert "PM7" in report.scientific_summary.hamiltonians


def test_method_b_preserves_baseline_correction_final(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    csv_path = tmp_path / "calculation_results.csv"
    write_canonical_csv(
        [
            _canonical_result(
                smiles="CCO",
                name="ethanol",
                method_id="pm7_delta_learning",
                estimates=[
                    PropertyEstimate(
                        estimate_id="b",
                        property_id="standard_enthalpy_of_formation",
                        role=PropertyRole.BASELINE,
                        method_id="pm7_delta_learning",
                        method_version="0.1.0",
                        hamiltonian="PM7",
                        value=Quantity(value=-60.0, unit="kcal/mol"),
                        value_kcal_mol=None,
                        value_kj_mol=None,
                        conformer_source_id=1,
                        uncertainty=None,
                        model_path=None,
                    ),
                    PropertyEstimate(
                        estimate_id="c",
                        property_id="standard_enthalpy_of_formation",
                        role=PropertyRole.CORRECTION,
                        method_id="pm7_delta_learning",
                        method_version="0.1.0",
                        hamiltonian=None,
                        value=Quantity(value=5.0, unit="kcal/mol"),
                        value_kcal_mol=None,
                        value_kj_mol=None,
                        conformer_source_id=None,
                        uncertainty=None,
                        model_path="models/a.joblib",
                    ),
                    PropertyEstimate(
                        estimate_id="f",
                        property_id="standard_enthalpy_of_formation",
                        role=PropertyRole.FINAL,
                        method_id="pm7_delta_learning",
                        method_version="0.1.0",
                        hamiltonian=None,
                        value=Quantity(value=-55.0, unit="kcal/mol"),
                        value_kcal_mol=None,
                        value_kj_mol=None,
                        conformer_source_id=None,
                        uncertainty=None,
                        model_path="models/a.joblib",
                    ),
                ],
            )
        ],
        csv_path,
    )
    roles = set(pd.read_csv(csv_path)["role"].str.lower())
    assert roles == {"baseline", "correction", "final"}

    manifest = run_service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id="pm7_delta_learning",
        method_version="0.1.0",
        method_snapshot={"display_name": "PM7 + Delta"},
        execution_overrides={},
        dataset_ref=None,
        model_ref={"name": "Model A"},
        molecule_count=1,
    )
    run_service.start_run(manifest.run_id)
    run_service.attach_output_paths(
        manifest.run_id, {"calculation_results_csv": csv_path}
    )
    run_service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    report = ResultsService(run_service=run_service).analyze_run(manifest.run_id)

    assert report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY
    assert set(report.scientific_summary.roles) >= {
        "baseline",
        "correction",
        "final",
    }
