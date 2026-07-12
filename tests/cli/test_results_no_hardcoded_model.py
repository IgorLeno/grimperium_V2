from __future__ import annotations

import io
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from rich.console import Console

from grimperium.cli.session import ModelRef, ModelState, SessionContext
from grimperium.cli.views.results_view import ResultsView


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
