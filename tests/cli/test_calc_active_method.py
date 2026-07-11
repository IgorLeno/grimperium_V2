"""Tests for Calculate using the active session method."""

from __future__ import annotations

import io
from unittest.mock import MagicMock, patch

import pytest
from rich.console import Console

from grimperium.calculation.methods import get_calculation_method
from grimperium.cli.controller import CliController
from grimperium.cli.views.calc_view import CalcView


@pytest.fixture
def controller() -> CliController:
    ctrl = CliController()
    buf = io.StringIO()
    ctrl.console = Console(file=buf, highlight=False, width=140)
    return ctrl


def test_do_prediction_redirects_when_no_method(controller: CliController) -> None:
    view = CalcView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    view.show_error = MagicMock()  # type: ignore[method-assign]
    view.render = MagicMock()  # type: ignore[method-assign]

    assert view.do_prediction() == "methods"
    view.show_error.assert_called_once()


@patch("grimperium.cli.views.calc_view.text_input", return_value=None)
def test_do_prediction_uses_active_method_without_reselect(
    _mock_input: MagicMock,
    controller: CliController,
) -> None:
    method = get_calculation_method(
        "semiempirical_am1_pm3_pm7",
        property_id="standard_enthalpy_of_formation",
    )
    controller.set_method(method)
    view = CalcView(controller)
    view.render = MagicMock()  # type: ignore[method-assign]
    view._select_method = MagicMock(  # type: ignore[method-assign]
        side_effect=AssertionError("_select_method must not be called")
    )

    # User cancels SMILES → stay in calc; method already resolved from session.
    assert view.do_prediction() is None
    view._select_method.assert_not_called()


def test_handle_predict_without_method_navigates(
    controller: CliController,
) -> None:
    view = CalcView(controller)
    view.do_prediction = MagicMock(return_value="methods")  # type: ignore[method-assign]
    assert view.handle_action("predict") == "methods"
