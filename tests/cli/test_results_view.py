"""PR7 boundary tests for ResultsView."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
from rich.console import Console

from grimperium.cli.controller import CliController
from grimperium.cli.views.results_view import ResultsView
from grimperium.results.models import ResultsAnalysisMode


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


def test_handle_charts_uses_analysis_report_scored_dataframe(
    tmp_path: Path,
    monkeypatch,
) -> None:
    view, _, controller = _create_results_view()
    canonical_long_form = tmp_path / "calculation_results.csv"
    canonical_long_form.write_text(
        "estimate_id,run_id,role,canonical_value\nbaseline,run_1,baseline,-55\n",
        encoding="utf-8",
    )
    controller.current_csv_path = canonical_long_form
    report = MagicMock()
    report.analysis_mode = ResultsAnalysisMode.PREDICTION_WITH_REFERENCE
    report.scored_df = pd.DataFrame(
        {
            "mol_id": ["m1"],
            "H298_cbs": [-50.0],
            "H298_predicted": [-51.0],
        }
    )
    view._load_analysis_report = MagicMock(return_value=report)  # type: ignore[method-assign]
    view._get_charts_dir = MagicMock(return_value=tmp_path / "charts")  # type: ignore[method-assign]
    captured: dict[str, Path] = {}

    def fake_generate_charts(csv_path: Path, output_dir: Path):
        captured["csv_path"] = Path(csv_path)
        assert Path(csv_path) != canonical_long_form
        assert "H298_predicted" in pd.read_csv(csv_path).columns
        return MagicMock(
            n_points=1,
            rmse=1.0,
            r2=0.9,
            parity_plot=output_dir / "parity.png",
            delta_histogram=output_dir / "delta.png",
            residuals_plot=output_dir / "residuals.png",
        )

    monkeypatch.setattr(
        "grimperium.cli.views.results_view.charts_module.generate_charts",
        fake_generate_charts,
    )

    view._handle_charts()

    assert captured["csv_path"].parent == tmp_path / "charts"


def test_handle_charts_skips_scientific_summary_only() -> None:
    view, buf, _ = _create_results_view()
    report = MagicMock()
    report.analysis_mode = ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY
    view._load_analysis_report = MagicMock(return_value=report)  # type: ignore[method-assign]

    view._handle_charts()

    assert "comparative charts" in buf.getvalue().lower()
