"""Tests for grimperium_mini TUI application."""

from unittest.mock import patch

from grimperium_mini.app import MiniApp


def test_show_welcome_no_exception() -> None:
    app = MiniApp()
    app.console = app.console.__class__(force_terminal=False)
    app.show_welcome()


def test_main_menu_quit_returns_none() -> None:
    app = MiniApp()
    with patch("grimperium_mini.app.Prompt.ask", return_value="q"):
        result = app.display_main_menu()
    assert result is None


def test_main_menu_run_returns_run() -> None:
    app = MiniApp()
    with patch("grimperium_mini.app.Prompt.ask", return_value="1"):
        result = app.display_main_menu()
    assert result == "run"


def test_run_exits_cleanly_on_quit() -> None:
    app = MiniApp()
    with (
        patch.object(app, "show_welcome"),
        patch.object(app, "display_main_menu", return_value=None),
    ):
        code = app.run()
    assert code == 0
