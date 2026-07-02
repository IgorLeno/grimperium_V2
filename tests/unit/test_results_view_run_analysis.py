"""Tests for redirected ResultsView calculation actions."""

from __future__ import annotations

import io
from unittest.mock import MagicMock, patch

from rich.console import Console

from grimperium.cli.views.results_view import ResultsView


def _create_results_view() -> tuple[ResultsView, io.StringIO]:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    controller = MagicMock()
    controller.console = console
    controller.current_model_path = None
    controller.current_csv_path = None
    view = ResultsView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    return view, buf


def test_run_analysis_menu_option_is_disabled() -> None:
    view, _ = _create_results_view()

    analysis = next(o for o in view.get_menu_options() if o.value == "run_analysis")

    assert analysis.label == "Run New Analysis"
    assert analysis.disabled is True
    assert analysis.disabled_reason == "Use Models view"


def test_predict_batch_menu_option_is_disabled() -> None:
    view, _ = _create_results_view()

    predict = next(o for o in view.get_menu_options() if o.value == "predict_batch")

    assert predict.label == "Predict Batch"
    assert predict.disabled is True
    assert predict.disabled_reason == "Use Models view"


def test_handle_action_routes_run_analysis_to_redirect_handler() -> None:
    view, _ = _create_results_view()

    with patch.object(view, "_handle_run_analysis") as mock_handler:
        result = view.handle_action("run_analysis")

    mock_handler.assert_called_once()
    assert result is None


def test_handle_action_routes_predict_batch_to_redirect_handler() -> None:
    view, _ = _create_results_view()

    with patch.object(view, "_handle_predict_batch") as mock_handler:
        result = view.handle_action("predict_batch")

    mock_handler.assert_called_once()
    assert result is None


def test_predict_batch_handler_redirects_without_prediction() -> None:
    view, buf = _create_results_view()

    with patch("grimperium.ml.predictor.predict_batch") as mock_predict:
        view._handle_predict_batch()

    mock_predict.assert_not_called()
    output = buf.getvalue()
    assert "Batch prediction is handled in Models > Predict Batch" in output
    assert "Results is for analysis only" in output


def test_run_analysis_handler_redirects_without_train_or_predict() -> None:
    view, buf = _create_results_view()

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
