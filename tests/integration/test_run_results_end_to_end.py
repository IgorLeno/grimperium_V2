"""End-to-end: real batch/calc outputs → Run → ResultsService.analyze_run."""

from __future__ import annotations

import io
import json
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from rich.console import Console

from grimperium.calculation.contracts.enums import OverallStatus, PropertyRole
from grimperium.calculation.contracts.models import (
    CalculationMethodReference,
    MoleculeCalculationResult,
    MoleculeData,
    PropertyEstimate,
    RunMetadata,
)
from grimperium.calculation.contracts.quantity import Quantity
from grimperium.calculation.methods import get_calculation_method
from grimperium.calculation.output.csv_writer import write_canonical_csv
from grimperium.cli.calc_pipeline import CalcPipelineResult
from grimperium.cli.session import SessionContext
from grimperium.cli.views.calc_view import CalcView
from grimperium.cli.views.databases_view import DatabasesView
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
from grimperium.runs.models import RunStatus
from grimperium.runs.service import RunService


class _FakeProcessor:
    def __init__(self, hofs: dict[str, float] | None = None) -> None:
        self._hofs = hofs or {}

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
        conformer.energy_hof = self._hofs.get(mol_id, -55.0)
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


def _controller(run_service: RunService) -> MagicMock:
    controller = MagicMock()
    controller.console = Console(file=io.StringIO(), highlight=False, width=140)
    controller.session = SessionContext()
    controller.run_service = run_service
    controller.current_model = "Model A"
    controller.current_model_path = None
    controller.current_method_definition = None
    controller.settings_manager = MagicMock()
    controller.session_summary.return_value = {
        "property": "Standard enthalpy of formation",
        "method": "Not selected",
        "dataset": "Not selected",
        "model": "Not required",
        "status": "Ready",
    }
    return controller


def _estimate(
    *,
    role: PropertyRole,
    value: float,
    method_id: str,
    method_version: str = "0.1.0",
    hamiltonian: str | None = None,
    model_path: str | None = None,
) -> PropertyEstimate:
    return PropertyEstimate(
        estimate_id=f"{role.value}-{hamiltonian or 'model'}",
        property_id="standard_enthalpy_of_formation",
        role=role,
        method_id=method_id,
        method_version=method_version,
        hamiltonian=hamiltonian,
        value=Quantity(value=value, unit="kcal/mol"),
        value_kcal_mol=value,
        value_kj_mol=None,
        conformer_source_id=0 if hamiltonian else None,
        uncertainty=None,
        model_path=model_path,
    )


def _method_a_canonical(
    *,
    smiles: str = "CCO",
    name: str = "ethanol",
    run_id: str = "run-method-a",
) -> MoleculeCalculationResult:
    result = _canonical_result(
        smiles=smiles,
        name=name,
        method_id="semiempirical_am1_pm3_pm7",
        estimates=[
            _estimate(
                role=PropertyRole.FINAL,
                value=-60.0,
                method_id="semiempirical_am1_pm3_pm7",
                hamiltonian="AM1",
            ),
            _estimate(
                role=PropertyRole.FINAL,
                value=-61.0,
                method_id="semiempirical_am1_pm3_pm7",
                hamiltonian="PM3",
            ),
            _estimate(
                role=PropertyRole.FINAL,
                value=-62.0,
                method_id="semiempirical_am1_pm3_pm7",
                hamiltonian="PM7",
            ),
        ],
    )
    return replace(
        result,
        run=RunMetadata(
            run_id=run_id,
            execution_phase="single_smiles",
            method_ref=CalculationMethodReference(
                method_id="semiempirical_am1_pm3_pm7",
                method_version="0.1.0",
                property_id="standard_enthalpy_of_formation",
            ),
            started_at=datetime(2026, 7, 12, 10, 0, tzinfo=timezone.utc),
            completed_at=datetime(2026, 7, 12, 10, 1, tzinfo=timezone.utc),
            grimperium_version=None,
        ),
    )


def _delta_canonical(run_id: str) -> MoleculeCalculationResult:
    result = _canonical_result(
        smiles="CCO",
        name="ethanol",
        method_id="pm7_delta_learning",
        estimates=[
            _estimate(
                role=PropertyRole.BASELINE,
                value=-65.0,
                method_id="pm7_delta_learning",
                hamiltonian="PM7",
            ),
            _estimate(
                role=PropertyRole.CORRECTION,
                value=-1.0,
                method_id="pm7_delta_learning",
                model_path="model.joblib",
            ),
            _estimate(
                role=PropertyRole.FINAL,
                value=-66.0,
                method_id="pm7_delta_learning",
                model_path="model.joblib",
            ),
        ],
    )
    return replace(
        result,
        run=RunMetadata(
            run_id=run_id,
            execution_phase="single_smiles",
            method_ref=CalculationMethodReference(
                method_id="pm7_delta_learning",
                method_version="0.1.0",
                property_id="standard_enthalpy_of_formation",
            ),
            started_at=datetime(2026, 7, 12, 10, 0, tzinfo=timezone.utc),
            completed_at=datetime(2026, 7, 12, 10, 1, tzinfo=timezone.utc),
            grimperium_version=None,
        ),
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
    run_service = RunService.from_environment()
    controller = _controller(run_service)
    view = DatabasesView(controller)
    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        [
            {
                "mol_id": mol_id,
                "smiles": "CCO",
                "nheavy": 3,
                "charge": 0,
                "multiplicity": 1,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            }
            for mol_id in ("mol_a", "mol_b")
        ]
    ).to_csv(csv_path, index=False)
    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()
    state_manager = BatchStateManager(tmp_path / "batch_state.csv", PM7Config())
    state_manager.reconcile_molecules(csv_manager.state_seed_rows())
    batch = csv_manager.select_batch(
        batch_id="batch_0001",
        batch_size=2,
        crest_timeout_minutes=30,
        mopac_timeout_minutes=60,
    )
    method_id, method_version, method_snapshot = view._method_run_fields()
    batch = batch.model_copy(
        update={
            "method_id": method_id,
            "method_version": method_version,
            "method_snapshot": method_snapshot,
        }
    )
    state_manager.mark_selected_from_batch(batch)
    manifest = run_service.create_run(
        property_id="standard_enthalpy_of_formation",
        method_id=method_id,
        method_version=method_version,
        method_snapshot=method_snapshot,
        execution_overrides=view._execution_overrides_snapshot(
            batch_size=batch.size,
            crest_timeout_minutes=30,
            mopac_timeout_minutes=60,
        ),
        dataset_ref=view._database_ref_snapshot(csv_path),
        model_ref=None,
        molecule_count=batch.size,
    )
    layout = BatchOutputLayout(view._pm7_batch_output_dir(manifest.run_id))
    manifest = run_service.start_run(manifest.run_id)
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
        processor_adapter=_FakeProcessor({"mol_a": -55.0, "mol_b": -57.0}),
        output_layout=layout,
        write_xlsx=False,
        canonical_run_id=manifest.run_id,
    )
    with patch.object(
        CSVManagerExtensions, "update_molecule_with_mopac_results", return_value=True
    ):
        batch_result = manager.execute_batch(batch)

    assert layout.calculation_results_csv.exists()
    canonical = pd.read_csv(layout.calculation_results_csv)
    assert set(canonical["run_id"]) == {manifest.run_id}
    assert set(canonical["method_id"]) == {"crest_pm7"}
    assert set(canonical["role"].str.lower()) == {"baseline"}
    assert "correction" not in set(canonical["role"].str.lower())
    assert "final" not in set(canonical["role"].str.lower())

    manifest = view._attach_existing_outputs(manifest, layout)
    terminal = run_service.complete_run(
        manifest.run_id,
        success_count=batch_result.success_count,
        failure_count=batch_result.failed_count
        + batch_result.rerun_count
        + batch_result.skip_count,
    )

    report = ResultsService(run_service=run_service).analyze_run(terminal.run_id)

    assert terminal.status is RunStatus.COMPLETED
    assert terminal.method_id == "crest_pm7"
    assert terminal.molecule_count == 2
    assert terminal.success_count == 2
    assert terminal.failure_count == 0
    assert "calculation_results_csv" in terminal.output_paths
    assert batch_result.total_count == terminal.molecule_count
    assert batch_result.success_count == terminal.success_count
    assert report.run_id == terminal.run_id
    assert report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY
    assert report.analysis is None
    assert report.scientific_summary.n_molecules == 2
    assert report.scientific_summary.n_estimates == 2
    assert report.scientific_summary.roles == ("baseline",)
    assert report.summary.get("mae") is None


def test_calc_view_method_a_do_prediction_writes_openable_authoritative_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    controller = _controller(run_service)
    method = get_calculation_method(
        "semiempirical_am1_pm3_pm7",
        property_id="standard_enthalpy_of_formation",
    )
    controller.current_method_definition = method
    controller.session.method_definition = method
    controller.session.property_id = method.property_id
    view = CalcView(controller)
    view.render = MagicMock()  # type: ignore[method-assign]
    view.render_method_a_result = MagicMock()  # type: ignore[method-assign]
    view._select_units = MagicMock(return_value="both")  # type: ignore[method-assign]

    class FakeMethodARunner:
        def __init__(self, **_kwargs: object) -> None:
            pass

        def calculate_single_smiles(self, smiles: str, **kwargs: object) -> object:
            return _method_a_canonical(
                smiles=smiles,
                name=str(kwargs["name"]),
                run_id=str(kwargs["run_id"]),
            )

    with (
        patch("grimperium.cli.views.calc_view.text_input", return_value="CCO"),
        patch(
            "grimperium.cli.views.calc_view.SemiempiricalFormationEnthalpyRunner",
            FakeMethodARunner,
        ),
    ):
        assert view.do_prediction() is None

    run_id = controller.session.run.run_id
    assert run_id is not None
    manifest = run_service.get_run(run_id)
    csv_path = manifest.output_paths["calculation_results_csv"]
    canonical = pd.read_csv(csv_path)
    report = ResultsService(run_service=run_service).analyze_run(run_id)

    assert manifest.status is RunStatus.COMPLETED
    assert manifest.method_id == "semiempirical_am1_pm3_pm7"
    assert manifest.molecule_count == 1
    assert manifest.success_count == 1
    assert manifest.failure_count == 0
    assert set(canonical["run_id"]) == {manifest.run_id}
    assert set(canonical["method_id"]) == {"semiempirical_am1_pm3_pm7"}
    assert set(canonical["role"].str.lower()) == {"final"}
    assert set(canonical["hamiltonian"]) == {"AM1", "PM3", "PM7"}
    assert csv_path.exists()
    assert "single_result" in manifest.output_paths
    assert report.run_id == manifest.run_id
    assert report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY
    assert report.scientific_summary.n_molecules == 1
    assert report.scientific_summary.n_estimates == 3
    assert "PM7" in report.scientific_summary.hamiltonians


def test_calc_view_method_b_without_canonical_fails_run_without_reusing_previous(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(tmp_path / "runs"))
    run_service = RunService.from_environment()
    controller = _controller(run_service)
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )
    model_path = tmp_path / "model.joblib"
    model_path.write_text("placeholder", encoding="utf-8")
    controller.current_method_definition = method
    controller.session.method_definition = method
    controller.session.property_id = method.property_id
    controller.current_model = "Delta model"
    controller.current_model_path = model_path
    controller.session.model.name = "Delta model"
    controller.session.model.path = model_path

    view = CalcView(controller)
    view.render = MagicMock()  # type: ignore[method-assign]
    view.render_result = MagicMock()  # type: ignore[method-assign]
    view._select_units = MagicMock(return_value="both")  # type: ignore[method-assign]
    view.last_calculation_result = _delta_canonical("stale-run")

    with (
        patch("grimperium.cli.views.calc_view.text_input", return_value="CCO"),
        patch("grimperium.cli.views.calc_view.validate_model_for_method"),
        patch(
            "grimperium.cli.views.calc_view.run_single_molecule_prediction",
            return_value=CalcPipelineResult(
                h298_pm7=-65.0,
                delta_correction=-1.0,
                h298_corrected=-66.0,
                n_conformers=3,
                execution_time=120.0,
                model_version="0.1.0",
                canonical=None,
            ),
        ),
    ):
        assert view.do_prediction() is None

    run_id = controller.session.run.run_id
    assert run_id is not None
    manifest = run_service.get_run(run_id)

    assert manifest.status is RunStatus.FAILED
    assert manifest.success_count == 0
    assert manifest.failure_count == 1
    assert manifest.output_paths == {}
    assert "MoleculeCalculationResult" in str(manifest.error)
    assert manifest.run_id != "stale-run"
    assert view.last_calculation_result is None
    assert not (run_service.runs_root / run_id / "calculation_results.csv").exists()


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
