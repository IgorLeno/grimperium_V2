"""Foundation contract for the independently launched Semi-Imperium shell."""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock, patch

from semi_imperium import menu
from semi_imperium.app import SemiImperiumCLI


def test_main_menu_exposes_only_three_focused_areas(monkeypatch: Any) -> None:
    captured: dict[str, Any] = {}

    def fake_show_menu_with_separator(**kwargs: Any) -> str:
        captured.update(kwargs)
        return "calc"

    monkeypatch.setattr(
        menu,
        "show_menu_with_separator",
        fake_show_menu_with_separator,
    )

    result = menu.show_main_menu()

    options = captured["options"]
    assert result == "calc"
    assert [(option.label, option.value) for option in options] == [
        ("CALCULATE", "calc"),
        ("DATABASE", "databases"),
        ("SETTINGS", "settings"),
    ]


def test_shell_registers_focused_areas_and_calculate_support() -> None:
    app = SemiImperiumCLI()

    assert set(app.controller._views) == {
        "calc",
        "methods",
        "databases",
        "settings",
    }
    assert app.controller.get_view("models") is None
    assert app.controller.get_view("results") is None
    assert app.controller.get_view("about") is None


def test_main_launches_semiimperium_without_preflight() -> None:
    with (
        patch("sys.argv", ["semi-imperium", "--skip-preflight"]),
        patch("grimperium.utils.logging.setup_logging"),
        patch("semi_imperium.app.SemiImperiumCLI") as cli_class,
    ):
        cli = MagicMock()
        cli.run.return_value = 0
        cli_class.return_value = cli

        from semi_imperium.app import main

        result = main()

    assert result == 0
    cli.run.assert_called_once_with(skip_preflight=True)
