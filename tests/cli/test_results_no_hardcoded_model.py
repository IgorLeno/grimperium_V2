from __future__ import annotations

import io
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from rich.console import Console

from grimperium.cli.session import ModelRef, ModelState, SessionContext
from grimperium.cli.views.results_view import ResultsView
from grimperium.results.models import ResultsAnalysisMode


def test_results_view_does_not_embed_legacy_model_names() -> None:
    source = Path("src/grimperium/cli/views/results_view.py").read_text(
        encoding="utf-8"
    )

    assert "DeltaLearner v1" not in source
    assert "KRR + XGBoost Ensemble" not in source
    assert "Models > Predict Batch" not in source
    assert "Predict Batch" not in source
    assert "data/thermo_pm7.csv" not in source
    assert "Analyze Active Source" in source
    assert "Select Saved Run" in source


def test_results_view_renders_session_model_metadata(monkeypatch) -> None:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    session = SessionContext(
        model=ModelRef(
            name="Custom Lab Model",
            path=Path("models/custom.joblib"),
            state=ModelState.READY,
        )
    )
    controller = SimpleNamespace(
        console=console,
        session=session,
        current_model_path=Path("models/custom.joblib"),
        current_csv_path=None,
        session_summary=lambda: {
            "property": "Standard enthalpy of formation",
            "method": "PM7 + Delta",
            "dataset": "Not selected",
            "model": "Custom Lab Model",
            "status": "Ready",
        },
    )
    monkeypatch.setattr(
        "grimperium.cli.views.results_view.ResultsService",
        lambda: MagicMock(
            analyze_dataset=MagicMock(return_value=None),
            model_metadata=MagicMock(
                return_value={
                    "model_label": "Custom Lab Model",
                    "algorithm": "metadata-driven",
                    "mae": 1.23,
                    "r2": 0.9876,
                    "status": "Ready",
                }
            ),
        ),
    )

    view = ResultsView(controller)  # type: ignore[arg-type]
    view._render_model_comparison()

    output = buf.getvalue()
    assert "Custom Lab Model" in output
    assert "metadata-driven" in output
    assert "DeltaLearner v1" not in output
    assert "KRR + XGBoost Ensemble" not in output


def test_results_view_prefers_active_run_model_label(monkeypatch) -> None:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    session = SessionContext(
        model=ModelRef(
            name="Session Model",
            path=Path("models/session.joblib"),
            state=ModelState.READY,
        )
    )
    controller = SimpleNamespace(
        console=console,
        session=session,
        current_model="Session Model",
        current_model_path=Path("models/session.joblib"),
        current_csv_path=None,
    )
    report = MagicMock(
        model_label="Run Manifest Model",
        analysis_mode=ResultsAnalysisMode.PREDICTION_WITH_REFERENCE,
        scientific_summary=MagicMock(method_id="pm7_delta_learning"),
    )

    view = ResultsView(controller)  # type: ignore[arg-type]
    view._load_analysis_report = MagicMock(return_value=report)  # type: ignore[method-assign]
    view._render_model_comparison()

    output = buf.getvalue()
    assert "Run Manifest Model" in output
    assert "Session Model" not in output


def test_results_view_shows_model_not_required_for_pm7_only_run() -> None:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=True, width=120)
    controller = SimpleNamespace(
        console=console,
        session=SessionContext(),
        current_model=None,
        current_model_path=None,
        current_csv_path=None,
    )
    report = MagicMock(
        model_label=None,
        analysis_mode=ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY,
        scientific_summary=MagicMock(method_id="crest_pm7"),
    )

    view = ResultsView(controller)  # type: ignore[arg-type]
    view._load_analysis_report = MagicMock(return_value=report)  # type: ignore[method-assign]
    view._render_model_comparison()

    assert "Not required" in buf.getvalue()
