"""Tests for ResultsView analysis source actions."""

from __future__ import annotations

import io
from unittest.mock import MagicMock

from rich.console import Console

from grimperium.cli.views.results_view import ResultsView


def _create_results_view() -> tuple[ResultsView, io.StringIO]:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    controller = MagicMock()
    controller.console = console
    controller.current_model_path = None
    controller.current_csv_path = None
    controller.session = MagicMock(analysis_path=None, run=None)
    view = ResultsView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    return view, buf


def test_analyze_source_menu_option_is_enabled() -> None:
    view, _ = _create_results_view()

    analysis = next(o for o in view.get_menu_options() if o.value == "analyze_source")

    assert analysis.label == "Analyze Active Source"
    assert analysis.disabled is False


def test_select_run_menu_option_is_enabled() -> None:
    view, _ = _create_results_view()

    select_run = next(o for o in view.get_menu_options() if o.value == "select_run")

    assert select_run.label == "Select Saved Run"
    assert select_run.disabled is False


def test_handle_action_routes_analyze_source() -> None:
    view, _ = _create_results_view()
    view._handle_analyze_source = MagicMock()  # type: ignore[method-assign]

    result = view.handle_action("analyze_source")

    view._handle_analyze_source.assert_called_once()
    assert result is None


def test_handle_action_routes_select_run() -> None:
    view, _ = _create_results_view()
    view._handle_select_run = MagicMock()  # type: ignore[method-assign]

    result = view.handle_action("select_run")

    view._handle_select_run.assert_called_once()
    assert result is None


def test_analyze_source_without_source_shows_error() -> None:
    view, buf = _create_results_view()
    view._load_analysis_report = MagicMock(return_value=None)  # type: ignore[method-assign]

    view._handle_analyze_source()

    view._load_analysis_report.assert_called_once_with(show_errors=True)
    view.wait_for_enter.assert_called_once()


def test_select_run_empty_history_does_not_call_set_run() -> None:
    view, buf = _create_results_view()
    view.controller.run_service = MagicMock()
    view.controller.run_service.list_runs.return_value = []
    view.controller.set_run = MagicMock()

    view._handle_select_run()

    view.controller.set_run.assert_not_called()
    assert "No saved runs found" in buf.getvalue()
