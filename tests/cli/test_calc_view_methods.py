"""Tests for method-aware Calculate view behavior."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock

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
from grimperium.cli.calc_pipeline import CalcPipelineResult
from grimperium.cli.mock_data import PredictionResult
from grimperium.cli.views.calc_view import CalcView


@pytest.fixture
def controller() -> MagicMock:
    buffer = io.StringIO()
    console = Console(file=buffer, highlight=False, width=140)
    ctrl = MagicMock()
    ctrl.console = console
    ctrl.current_model = "test_model"
    ctrl.current_model_path = None
    ctrl.settings_manager = MagicMock()
    return ctrl


def test_render_available_methods_lists_standard_enthalpy_methods(
    controller: MagicMock,
) -> None:
    view = CalcView(controller)

    view.render_available_methods()

    output = controller.console.file.getvalue()
    assert "semiempirical_am1_pm3_pm7" in output
    assert "pm7_delta_learning" in output


def test_calculation_methods_action_is_available(controller: MagicMock) -> None:
    view = CalcView(controller)

    options = view.get_menu_options()

    assert any(option.value == "methods" for option in options)


def test_resolve_required_model_blocks_method_b_without_session_model(
    controller: MagicMock,
) -> None:
    view = CalcView(controller)
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )

    model_path = view._resolve_required_model(method)

    output = controller.console.file.getvalue()
    assert model_path is None
    assert "Select a compatible model" in output


def _method_a_result() -> MoleculeCalculationResult:
    return MoleculeCalculationResult(
        molecule=MoleculeData(smiles="CCO", name="calc_test"),
        run=RunMetadata(
            run_id="calc_test",
            execution_phase="single_smiles",
            method_ref=CalculationMethodReference(
                method_id="semiempirical_am1_pm3_pm7",
                method_version="0.1.0",
                property_id="standard_enthalpy_of_formation",
            ),
            started_at=None,
            completed_at=None,
            grimperium_version=None,
        ),
        overall_status=OverallStatus.SUCCESS,
        conformers=[MagicMock(), MagicMock()],
        molecular_descriptors=None,
        estimates=[
            PropertyEstimate(
                estimate_id="am1-final",
                property_id="standard_enthalpy_of_formation",
                role=PropertyRole.FINAL,
                method_id="semiempirical_am1_pm3_pm7",
                method_version="0.1.0",
                hamiltonian="AM1",
                value=Quantity(value=-60.0, unit="kcal/mol"),
                value_kcal_mol=-60.0,
                value_kj_mol=-251.04,
                conformer_source_id=0,
                uncertainty=None,
                model_path=None,
            ),
            PropertyEstimate(
                estimate_id="pm7-final",
                property_id="standard_enthalpy_of_formation",
                role=PropertyRole.FINAL,
                method_id="semiempirical_am1_pm3_pm7",
                method_version="0.1.0",
                hamiltonian="PM7",
                value=Quantity(value=-62.0, unit="kcal/mol"),
                value_kcal_mol=-62.0,
                value_kj_mol=-259.408,
                conformer_source_id=0,
                uncertainty=None,
                model_path=None,
            ),
        ],
        artifacts=[],
        stage_executions=[
            MagicMock(execution_time_s=1.5),
            MagicMock(execution_time_s=None),
            MagicMock(execution_time_s=2.5),
        ],
    )


def test_method_a_runs_without_session_model(
    monkeypatch: pytest.MonkeyPatch,
    controller: MagicMock,
) -> None:
    calls: dict[str, object] = {}

    class FakeRunner:
        def __init__(self, **kwargs: object) -> None:
            calls["xtb_enabled"] = kwargs["xtb_enabled"]

        def calculate_single_smiles(self, smiles: str, **kwargs: object) -> object:
            calls["smiles"] = smiles
            calls["molecule_id"] = kwargs["molecule_id"]
            return _method_a_result()

    from grimperium.cli.views import calc_view

    monkeypatch.setattr(calc_view, "SemiempiricalFormationEnthalpyRunner", FakeRunner)
    controller.current_model_path = None
    view = CalcView(controller)
    method = get_calculation_method(
        "semiempirical_am1_pm3_pm7",
        property_id="standard_enthalpy_of_formation",
    )

    result = view._run_method_a("CCO", method)

    assert result is True
    assert calls["xtb_enabled"] is True
    assert calls["smiles"] == "CCO"
    assert view.last_calculation_result is not None
    assert view.last_result is not None
    assert view.last_result.smiles == "CCO"
    assert view.last_result.h298_pm7 == -62.0
    assert view.last_result.delta_correction == 0.0
    assert view.last_result.h298_corrected == -62.0
    assert view.last_result.model_name == "Method A"
    assert view.last_result.model_version == method.version
    assert view.last_result.execution_time == 4.0
    assert view.last_result.n_conformers == 2


def test_method_b_validates_selected_model_and_uses_existing_pipeline(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    controller: MagicMock,
) -> None:
    model_path = tmp_path / "model.joblib"
    model_path.write_text("placeholder", encoding="utf-8")
    controller.current_model_path = model_path
    validations: list[tuple[Path, str]] = []
    predictions: list[tuple[str, Path]] = []

    def fake_validate(path: Path, method: object) -> None:
        validations.append((path, method.method_id))

    def fake_prediction(
        smiles: str,
        mol_id: str,
        path: Path,
        config: object,
        progress_callback: object,
    ) -> CalcPipelineResult:
        predictions.append((smiles, path))
        return CalcPipelineResult(
            h298_pm7=-65.0,
            delta_correction=-1.0,
            h298_corrected=-66.0,
            n_conformers=3,
            execution_time=120.0,
            model_version="1.0.0",
        )

    from grimperium.cli.views import calc_view

    monkeypatch.setattr(calc_view, "validate_model_for_method", fake_validate)
    monkeypatch.setattr(calc_view, "run_single_molecule_prediction", fake_prediction)
    view = CalcView(controller)
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )

    result = view._run_method_b("CCO", method)

    assert result is True
    assert validations == [(model_path, "pm7_delta_learning")]
    assert predictions == [("CCO", model_path)]
    assert view.history[-1].h298_corrected == -66.0


def test_method_a_result_can_render_both_units_without_mutating_canonical_value(
    controller: MagicMock,
) -> None:
    view = CalcView(controller)
    result = _method_a_result()

    view.render_method_a_result(result, units="both")

    output = controller.console.file.getvalue()
    assert "-60.00 kcal/mol" in output
    assert "-251.04 kJ/mol" in output
    assert result.estimates[0].value.unit == "kcal/mol"


def test_run_header_displays_method_a_last_result(
    monkeypatch: pytest.MonkeyPatch,
    controller: MagicMock,
) -> None:
    view = CalcView(controller)
    view.last_result = PredictionResult(
        smiles="CCO",
        h298_pm7=-62.0,
        delta_correction=0.0,
        h298_corrected=-62.0,
        model_name="Method A",
        model_version="0.1.0",
        execution_time=4.0,
        n_conformers=2,
    )
    calls = {"count": 0}

    def fake_show_back_menu(*args: object, **kwargs: object) -> str:
        calls["count"] += 1
        return "back"

    monkeypatch.setattr(
        "grimperium.cli.views.calc_view.show_back_menu", fake_show_back_menu
    )

    next_view = view.run()

    output = controller.console.file.getvalue()
    assert next_view == "main"
    assert calls["count"] == 1
    assert "Last prediction: CCO" in output
    assert "H298=-62.00 kcal/mol" in output
    assert "-259.41 kJ/mol" in output
