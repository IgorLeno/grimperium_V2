"""Tests for redesigned main menu (session header + CALCULATION METHODS)."""

from __future__ import annotations

from typing import Any

from grimperium.cli import menu
from grimperium.cli.controller import CliController
from grimperium.cli.menu import format_session_header


class _FakeQuestion:
    def __init__(self, result: str) -> None:
        self._result = result

    def ask(self) -> str:
        return self._result


def test_main_menu_includes_calculation_methods(monkeypatch: Any) -> None:
    captured: dict[str, Any] = {}

    def fake_select(**kwargs: Any) -> _FakeQuestion:
        captured.update(kwargs)
        return _FakeQuestion("methods")

    monkeypatch.setattr(menu.questionary, "select", fake_select)

    result = menu.show_main_menu()

    titles = [
        choice.title for choice in captured["choices"] if hasattr(choice, "title")
    ]
    assert result == "methods"
    assert any("CALCULATION METHODS" in title for title in titles)
    assert any("CALCULATE" in title for title in titles)


def test_main_menu_header_has_no_mock_model(monkeypatch: Any) -> None:
    captured: dict[str, Any] = {}

    def fake_select(**kwargs: Any) -> _FakeQuestion:
        captured.update(kwargs)
        return _FakeQuestion("calc")

    monkeypatch.setattr(menu.questionary, "select", fake_select)

    menu.show_main_menu()

    message = captured.get("message", "")
    assert "DeltaXGB" not in message
    assert "Method:" in message
    assert "Property:" in message
    assert "Not selected" in message


def test_format_session_header_from_controller() -> None:
    ctrl = CliController()
    summary = ctrl.session_summary()
    header = format_session_header(
        property_label=summary["property"],
        method_label=summary["method"],
        model_label=summary["model"],
        status=summary["status"],
    )
    assert "Method: Not selected" in header
    assert "No model selected" in header
    assert "DeltaXGB" not in header


def test_main_menu_labels_calculate_without_changing_route(
    monkeypatch: Any,
) -> None:
    captured: dict[str, Any] = {}

    def fake_select(**kwargs: Any) -> _FakeQuestion:
        captured.update(kwargs)
        return _FakeQuestion("calc")

    monkeypatch.setattr(menu.questionary, "select", fake_select)

    result = menu.show_main_menu(model_label="model", status="Ready")

    titles = [
        choice.title for choice in captured["choices"] if hasattr(choice, "title")
    ]
    assert result == "calc"
    assert any("CALCULATE" in title for title in titles)
    assert not any(
        "CALC" in title.replace("CALCULATE", "").replace("CALCULATION METHODS", "")
        for title in titles
    )
