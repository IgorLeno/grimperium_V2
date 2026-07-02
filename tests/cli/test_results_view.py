"""PR7 boundary tests for ResultsView."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock, patch

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


def test_calculation_actions_are_disabled_in_results_menu() -> None:
    view, _, _ = _create_results_view()

    options = {option.value: option for option in view.get_menu_options()}

    assert options["run_analysis"].disabled is True
    assert options["run_analysis"].disabled_reason == "Use Models view"
    assert options["predict_batch"].disabled is True
    assert options["predict_batch"].disabled_reason == "Use Models view"


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


def test_get_csv_path_falls_back_with_warning() -> None:
    view, buf, controller = _create_results_view()
    controller.current_csv_path = None

    assert view._get_csv_path() == Path("data/thermo_pm7.csv")
    assert "default data file" in buf.getvalue()


def test_predict_batch_redirect_does_not_call_predictor() -> None:
    view, buf, _ = _create_results_view()

    with patch("grimperium.ml.predictor.predict_batch") as mock_predict:
        view._handle_predict_batch()

    mock_predict.assert_not_called()
    output = buf.getvalue()
    assert "Batch prediction is handled in Models > Predict Batch" in output
    assert "Results is for analysis only" in output


def test_run_analysis_redirect_does_not_call_training_or_prediction() -> None:
    view, buf, _ = _create_results_view()

    with (
        patch("grimperium.ml.trainer.train") as mock_train,
        patch("grimperium.ml.predictor.predict_batch") as mock_predict,
    ):
        view._handle_run_analysis()

    mock_train.assert_not_called()
    mock_predict.assert_not_called()
    output = buf.getvalue()
    assert "Training and prediction are handled in Models" in output
    assert "Use Models > Train Model and Models > Predict Batch first" in output
