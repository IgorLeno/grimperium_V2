"""Tests for MethodsView session method selection."""

from __future__ import annotations

import io
from unittest.mock import MagicMock

import pytest
from rich.console import Console

from grimperium.calculation.methods import get_calculation_method
from grimperium.cli.controller import CliController
from grimperium.cli.views.methods_view import MethodsView


@pytest.fixture
def controller() -> CliController:
    ctrl = CliController()
    buf = io.StringIO()
    ctrl.console = Console(file=buf, highlight=False, width=140)
    return ctrl


def test_methods_view_lists_registry_methods(controller: CliController) -> None:
    view = MethodsView(controller)
    view.render()
    output = controller.console.file.getvalue()  # type: ignore[attr-defined]
    assert "semiempirical_am1_pm3_pm7" in output
    assert "pm7_delta_learning" in output


def test_select_method_a_persists_on_controller(
    controller: CliController,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    view = MethodsView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    view.show_success = MagicMock()  # type: ignore[method-assign]
    monkeypatch.setattr(
        "grimperium.cli.views.methods_view.show_menu",
        lambda options, title="": "semiempirical_am1_pm3_pm7",
    )

    assert view.select_active_method() is True
    assert controller.current_method_id == "semiempirical_am1_pm3_pm7"
    assert controller.session_summary()["model"] == "Not required"
    assert controller.status == "Ready"


def test_select_method_b_requires_model(
    controller: CliController,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    view = MethodsView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    view.show_success = MagicMock()  # type: ignore[method-assign]
    monkeypatch.setattr(
        "grimperium.cli.views.methods_view.show_menu",
        lambda options, title="": "pm7_delta_learning",
    )

    assert view.select_active_method() is True
    assert controller.current_method_id == "pm7_delta_learning"
    assert controller.status == "Model required"


def test_cancel_selection_preserves_previous_method(
    controller: CliController,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    previous = get_calculation_method(
        "semiempirical_am1_pm3_pm7",
        property_id="standard_enthalpy_of_formation",
    )
    controller.set_method(previous)
    view = MethodsView(controller)
    monkeypatch.setattr(
        "grimperium.cli.views.methods_view.show_menu",
        lambda options, title="": None,
    )

    assert view.select_active_method() is False
    assert controller.current_method_id == "semiempirical_am1_pm3_pm7"
