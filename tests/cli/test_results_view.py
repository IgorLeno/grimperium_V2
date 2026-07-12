"""PR7 boundary tests for ResultsView."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock

from rich.console import Console

from grimperium.cli.controller import CliController
from grimperium.cli.views.results_view import ResultsView


def _create_results_view() -> tuple[ResultsView, io.StringIO, MagicMock]:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    controller = MagicMock()
    controller.console = console
    controller.current_model_path = None
    controller.current_csv_path = None
    controller.session = MagicMock(analysis_path=None, run=None)
    view = ResultsView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    return view, buf, controller


def test_controller_tracks_current_csv_path() -> None:
    controller = CliController()
    csv_path = Path("data/custom.csv")

    assert controller.current_csv_path is None

    controller.set_csv(csv_path)
    assert controller.current_csv_path == csv_path

    controller.set_csv(None)
    assert controller.current_csv_path is None


def test_analysis_actions_are_enabled_in_results_menu() -> None:
    view, _, _ = _create_results_view()

    options = {option.value: option for option in view.get_menu_options()}

    assert options["analyze_source"].disabled is False
    assert options["select_run"].disabled is False
    assert "predict_batch" not in options
    assert "run_analysis" not in options


def test_pure_analysis_actions_remain_enabled() -> None:
    view, _, _ = _create_results_view()

    options = {option.value: option for option in view.get_menu_options()}

    for action in (
        "detailed",
        "charts",
        "show_outliers",
        "top_errors",
        "html_report",
        "export_outliers",
    ):
        assert options[action].disabled is False


def test_get_model_path_uses_controller_session_not_env(
    monkeypatch,
) -> None:
    view, buf, controller = _create_results_view()
    controller.current_model_path = Path("models/session.joblib")
    monkeypatch.setenv("GRIMPERIUM_MODEL_PATH", "models/env.joblib")

    assert view._get_model_path() == Path("models/session.joblib")
    assert "Model not selected" not in buf.getvalue()


def test_get_model_path_warns_and_returns_none_when_unselected() -> None:
    view, buf, controller = _create_results_view()
    controller.current_model_path = None

    assert view._get_model_path() is None
    assert "Model not selected" in buf.getvalue()


def test_get_csv_path_uses_controller_session_not_env(monkeypatch) -> None:
    view, buf, controller = _create_results_view()
    controller.current_csv_path = Path("data/session.csv")
    monkeypatch.setenv("GRIMPERIUM_DATA_PATH", "data/env.csv")

    assert view._get_csv_path() == Path("data/session.csv")
    assert "default data file" not in buf.getvalue()


def test_get_csv_path_returns_none_when_unselected() -> None:
    view, buf, controller = _create_results_view()
    controller.current_csv_path = None

    assert view._get_csv_path() is None
    assert "default data file" not in buf.getvalue()


def test_analyze_source_routes_without_predictor() -> None:
    view, _, _ = _create_results_view()
    view._load_analysis_report = MagicMock(return_value=None)  # type: ignore[method-assign]

    assert view.handle_action("analyze_source") is None
    view._load_analysis_report.assert_called_once_with(show_errors=True)


def test_select_run_with_empty_history_shows_empty_state() -> None:
    view, buf, controller = _create_results_view()
    controller.run_service = MagicMock()
    controller.run_service.list_runs.return_value = []

    assert view.handle_action("select_run") is None
    assert "No saved runs found" in buf.getvalue()
