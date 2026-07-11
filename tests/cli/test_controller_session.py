"""Tests for CLI session controller (method/model state, no mock defaults)."""

from __future__ import annotations

from pathlib import Path

from grimperium.calculation.methods import get_calculation_method
from grimperium.cli.controller import CliController


def test_controller_starts_without_mock_model() -> None:
    ctrl = CliController()
    assert ctrl.current_model is None
    assert ctrl.current_model_path is None
    assert ctrl.current_method_definition is None
    assert ctrl.status == "No method selected"
    summary = ctrl.session_summary()
    assert summary["model"] == "No model selected"
    assert summary["method"] == "Not selected"
    assert "DeltaXGB" not in summary["model"]


def test_set_method_a_marks_model_not_required() -> None:
    ctrl = CliController()
    method = get_calculation_method(
        "semiempirical_am1_pm3_pm7",
        property_id="standard_enthalpy_of_formation",
    )
    ctrl.set_method(method)
    assert ctrl.current_method_id == "semiempirical_am1_pm3_pm7"
    assert ctrl.status == "Ready"
    summary = ctrl.session_summary()
    assert summary["model"] == "Not required"
    assert summary["property"] == method.property_name
    assert summary["method"] == method.display_name


def test_set_method_b_without_model_requires_model() -> None:
    ctrl = CliController()
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )
    ctrl.set_method(method)
    assert ctrl.current_method_id == "pm7_delta_learning"
    assert ctrl.status == "Model required"


def test_clear_method_resets_status() -> None:
    ctrl = CliController()
    method = get_calculation_method(
        "semiempirical_am1_pm3_pm7",
        property_id="standard_enthalpy_of_formation",
    )
    ctrl.set_method(method)
    ctrl.clear_method()
    assert ctrl.current_method_definition is None
    assert ctrl.status == "No method selected"


def test_set_model_refreshes_status_for_method_b(tmp_path: Path) -> None:
    ctrl = CliController()
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )
    ctrl.set_method(method)
    missing = tmp_path / "missing.joblib"
    ctrl.set_model("fake", model_path=missing)
    assert ctrl.status == "Model required"
