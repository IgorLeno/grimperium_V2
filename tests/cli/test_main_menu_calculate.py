"""Tests for the top-level Calculate menu label."""

from __future__ import annotations

from typing import Any

from grimperium.cli import menu


class _FakeQuestion:
    def __init__(self, result: str) -> None:
        self._result = result

    def ask(self) -> str:
        return self._result


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
    # Route value remains "calc"; label is CALCULATE (not a short "CALC" item).
    assert not any(
        "CALC" in title.replace("CALCULATE", "").replace("CALCULATION METHODS", "")
        for title in titles
    )
