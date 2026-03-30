"""Tests for Run New Analysis feature in ResultsView."""

from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from rich.console import Console

from grimperium.cli.views.results_view import ResultsView


def _create_results_view() -> tuple[ResultsView, io.StringIO]:
    """Create a ResultsView with a real Console capturing to StringIO."""
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    controller = MagicMock()
    controller.console = console
    view = ResultsView(controller)
    view.wait_for_enter = MagicMock()  # type: ignore[method-assign]
    return view, buf


def _make_csv(tmp_path: Path) -> Path:
    """Create a minimal CSV with eligible molecules."""
    csv_path = tmp_path / "thermo_pm7.csv"
    df = pd.DataFrame(
        {
            "mol_id": ["mol_1", "mol_2", "mol_3"],
            "status": ["OK", "OK", "Pending"],
            "cbs_quality_flag": ["OK", "OK", "OK"],
            "H298_cbs": [-100.0, -200.0, -300.0],
            "H298_pm7": [-101.0, -198.0, -303.0],
        }
    )
    df.to_csv(csv_path, index=False)
    return csv_path


# Fake train result for mocking
_FAKE_LEARNER = MagicMock()
_FAKE_TRAIN_M: dict[str, float] = {"rmse": 1.0, "mae": 0.8, "r2": 0.999}
_FAKE_TEST_M: dict[str, float] = {
    "rmse": 1.5,
    "mae": 1.0,
    "r2": 0.998,
    "max_error": 3.0,
    "mape": 0.5,
}
_FAKE_PIPELINE = MagicMock()
_FAKE_TRAIN_RESULT = (_FAKE_LEARNER, _FAKE_TRAIN_M, _FAKE_TEST_M, _FAKE_PIPELINE)

_FAKE_PRED_STATS = {
    "n_predicted": 2,
    "mean_delta": 1.5,
    "std_delta": 0.5,
    "min_predicted": -201.0,
    "max_predicted": -99.0,
}


@dataclass
class FakeChartResult:
    n_points: int = 2
    rmse: float = 1.5
    r2: float = 0.998
    parity_plot: str = "reports/charts/parity_plot.png"
    delta_histogram: str = "reports/charts/delta_histogram.png"
    residuals_plot: str = "reports/charts/residuals_plot.png"


def test_run_analysis_menu_option_exists() -> None:
    """MenuOption for 'Run New Analysis' must exist and not be disabled."""
    view, _ = _create_results_view()
    options = view.get_menu_options()
    analysis = next(o for o in options if o.value == "run_analysis")
    assert analysis.disabled is False
    assert analysis.label == "Run New Analysis"


def test_run_analysis_is_first_menu_option() -> None:
    """'Run New Analysis' must be the first menu option."""
    view, _ = _create_results_view()
    options = view.get_menu_options()
    assert options[0].value == "run_analysis"


def test_handle_action_routes_to_handler() -> None:
    """handle_action('run_analysis') must call _handle_run_analysis."""
    view, _ = _create_results_view()
    with patch.object(view, "_handle_run_analysis") as mock_handler:
        result = view.handle_action("run_analysis")
    mock_handler.assert_called_once()
    assert result is None


def test_no_csv_shows_error(tmp_path: Path) -> None:
    """Missing CSV → show_error called, train never called."""
    view, buf = _create_results_view()
    missing_csv = tmp_path / "nonexistent.csv"
    model_path = tmp_path / "model.joblib"

    with (
        patch.object(view, "_get_csv_path", return_value=missing_csv),
        patch.object(view, "_get_model_path", return_value=model_path),
        patch(
            "grimperium.cli.views.results_view.load_model_metadata",
            side_effect=FileNotFoundError,
        ),
        patch(
            "grimperium.ml.trainer.train",
        ) as mock_train,
    ):
        view._handle_run_analysis()

    output = buf.getvalue()
    assert "Error" in output
    mock_train.assert_not_called()


def test_confirmation_denied_aborts(tmp_path: Path) -> None:
    """User types 'no' → train not called."""
    view, buf = _create_results_view()
    csv_path = _make_csv(tmp_path)
    model_path = tmp_path / "model.joblib"

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_model_path", return_value=model_path),
        patch(
            "grimperium.cli.views.results_view.load_model_metadata",
            side_effect=FileNotFoundError,
        ),
        patch.object(view.console, "input", return_value="no"),
        patch(
            "grimperium.ml.trainer.train",
        ) as mock_train,
    ):
        view._handle_run_analysis()

    output = buf.getvalue()
    assert "cancelled" in output.lower()
    mock_train.assert_not_called()


def test_full_pipeline_success(tmp_path: Path) -> None:
    """All 3 steps complete successfully, summary panel displayed."""
    view, buf = _create_results_view()
    csv_path = _make_csv(tmp_path)
    model_path = tmp_path / "model.joblib"
    charts_dir = tmp_path / "charts"

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_model_path", return_value=model_path),
        patch.object(view, "_get_charts_dir", return_value=charts_dir),
        patch(
            "grimperium.cli.views.results_view.load_model_metadata",
            side_effect=FileNotFoundError,
        ),
        patch.object(view.console, "input", return_value="yes"),
        patch(
            "grimperium.cli.views.results_view.save_model",
        ) as mock_save,
        patch(
            "grimperium.cli.views.results_view.evaluate_gate",
            return_value=True,
        ),
        patch(
            "grimperium.ml.trainer.train",
            return_value=_FAKE_TRAIN_RESULT,
        ) as mock_train,
        patch(
            "grimperium.ml.predictor.predict_batch",
            return_value=(MagicMock(), _FAKE_PRED_STATS),
        ) as mock_predict,
        patch(
            "grimperium.ml.charts.generate_charts",
            return_value=FakeChartResult(),
        ) as mock_charts,
    ):
        view._handle_run_analysis()

    mock_train.assert_called_once()
    mock_save.assert_called_once()
    mock_predict.assert_called_once()
    mock_charts.assert_called_once()

    output = buf.getvalue()
    assert "Analysis Complete" in output
    assert "molecules predicted" in output
    assert "parity_plot.png" in output


def test_train_failure_stops_pipeline(tmp_path: Path) -> None:
    """train() raises → predict_batch and generate_charts not called."""
    view, buf = _create_results_view()
    csv_path = _make_csv(tmp_path)
    model_path = tmp_path / "model.joblib"

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_model_path", return_value=model_path),
        patch(
            "grimperium.cli.views.results_view.load_model_metadata",
            side_effect=FileNotFoundError,
        ),
        patch.object(view.console, "input", return_value="yes"),
        patch(
            "grimperium.ml.trainer.train",
            side_effect=RuntimeError("training boom"),
        ),
        patch(
            "grimperium.ml.predictor.predict_batch",
        ) as mock_predict,
        patch(
            "grimperium.ml.charts.generate_charts",
        ) as mock_charts,
    ):
        view._handle_run_analysis()

    output = buf.getvalue()
    assert "Training failed" in output
    mock_predict.assert_not_called()
    mock_charts.assert_not_called()


def test_predict_failure_stops_after_train(tmp_path: Path) -> None:
    """predict_batch raises → generate_charts not called."""
    view, buf = _create_results_view()
    csv_path = _make_csv(tmp_path)
    model_path = tmp_path / "model.joblib"

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_model_path", return_value=model_path),
        patch(
            "grimperium.cli.views.results_view.load_model_metadata",
            side_effect=FileNotFoundError,
        ),
        patch.object(view.console, "input", return_value="yes"),
        patch(
            "grimperium.cli.views.results_view.save_model",
        ),
        patch(
            "grimperium.cli.views.results_view.evaluate_gate",
            return_value=True,
        ),
        patch(
            "grimperium.ml.trainer.train",
            return_value=_FAKE_TRAIN_RESULT,
        ),
        patch(
            "grimperium.ml.predictor.predict_batch",
            side_effect=RuntimeError("predict boom"),
        ),
        patch(
            "grimperium.ml.charts.generate_charts",
        ) as mock_charts,
    ):
        view._handle_run_analysis()

    output = buf.getvalue()
    assert "Batch prediction failed" in output
    mock_charts.assert_not_called()


def test_charts_failure_non_fatal(tmp_path: Path) -> None:
    """generate_charts raises → warning shown, summary still displayed."""
    view, buf = _create_results_view()
    csv_path = _make_csv(tmp_path)
    model_path = tmp_path / "model.joblib"
    charts_dir = tmp_path / "charts"

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_model_path", return_value=model_path),
        patch.object(view, "_get_charts_dir", return_value=charts_dir),
        patch(
            "grimperium.cli.views.results_view.load_model_metadata",
            side_effect=FileNotFoundError,
        ),
        patch.object(view.console, "input", return_value="yes"),
        patch(
            "grimperium.cli.views.results_view.save_model",
        ),
        patch(
            "grimperium.cli.views.results_view.evaluate_gate",
            return_value=True,
        ),
        patch(
            "grimperium.ml.trainer.train",
            return_value=_FAKE_TRAIN_RESULT,
        ),
        patch(
            "grimperium.ml.predictor.predict_batch",
            return_value=(MagicMock(), _FAKE_PRED_STATS),
        ),
        patch(
            "grimperium.ml.charts.generate_charts",
            side_effect=RuntimeError("charts boom"),
        ),
    ):
        view._handle_run_analysis()

    output = buf.getvalue()
    assert "Chart generation failed" in output
    assert "Analysis Complete" in output


def test_before_after_comparison_shown(tmp_path: Path) -> None:
    """When old model exists, before/after comparison table is rendered."""
    view, buf = _create_results_view()
    csv_path = _make_csv(tmp_path)
    model_path = tmp_path / "model.joblib"
    model_path.touch()
    charts_dir = tmp_path / "charts"

    old_meta = {
        "metrics": {
            "test": {"mae": 2.0, "rmse": 3.0, "r2": 0.990, "max_error": 8.0}
        }
    }

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_model_path", return_value=model_path),
        patch.object(view, "_get_charts_dir", return_value=charts_dir),
        patch(
            "grimperium.cli.views.results_view.load_model_metadata",
            return_value=old_meta,
        ),
        patch.object(view.console, "input", return_value="yes"),
        patch(
            "grimperium.cli.views.results_view.save_model",
        ),
        patch(
            "grimperium.cli.views.results_view.evaluate_gate",
            return_value=True,
        ),
        patch(
            "grimperium.ml.trainer.train",
            return_value=_FAKE_TRAIN_RESULT,
        ),
        patch(
            "grimperium.ml.predictor.predict_batch",
            return_value=(MagicMock(), _FAKE_PRED_STATS),
        ),
        patch(
            "grimperium.ml.charts.generate_charts",
            return_value=FakeChartResult(),
        ),
    ):
        view._handle_run_analysis()

    output = buf.getvalue()
    assert "Before vs After" in output
    assert "2.0000" in output  # old MAE
    assert "1.0000" in output  # new MAE


def test_no_old_model_skips_comparison(tmp_path: Path) -> None:
    """No old model → analysis runs without comparison table."""
    view, buf = _create_results_view()
    csv_path = _make_csv(tmp_path)
    model_path = tmp_path / "model.joblib"
    charts_dir = tmp_path / "charts"

    with (
        patch.object(view, "_get_csv_path", return_value=csv_path),
        patch.object(view, "_get_model_path", return_value=model_path),
        patch.object(view, "_get_charts_dir", return_value=charts_dir),
        patch(
            "grimperium.cli.views.results_view.load_model_metadata",
            side_effect=FileNotFoundError,
        ),
        patch.object(view.console, "input", return_value="yes"),
        patch(
            "grimperium.cli.views.results_view.save_model",
        ),
        patch(
            "grimperium.cli.views.results_view.evaluate_gate",
            return_value=True,
        ),
        patch(
            "grimperium.ml.trainer.train",
            return_value=_FAKE_TRAIN_RESULT,
        ),
        patch(
            "grimperium.ml.predictor.predict_batch",
            return_value=(MagicMock(), _FAKE_PRED_STATS),
        ),
        patch(
            "grimperium.ml.charts.generate_charts",
            return_value=FakeChartResult(),
        ),
    ):
        view._handle_run_analysis()

    output = buf.getvalue()
    assert "Before vs After" not in output
    assert "Analysis Complete" in output
