"""Provenance: PM7 batch must record crest_pm7 even if session is Delta Learning."""

from __future__ import annotations

from dataclasses import asdict
from unittest.mock import MagicMock

from grimperium.calculation.methods import get_calculation_method
from grimperium.cli.views.databases_view import DatabasesView


def test_method_run_fields_always_crest_pm7_ignores_session_delta() -> None:
    delta = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )
    controller = MagicMock()
    controller.current_method_definition = delta
    view = DatabasesView(controller)

    method_id, method_version, snapshot = view._method_run_fields()

    assert method_id == "crest_pm7"
    assert method_version == get_calculation_method("crest_pm7").version
    assert snapshot["method_id"] == "crest_pm7"
    assert asdict(delta)["method_id"] == "pm7_delta_learning"
