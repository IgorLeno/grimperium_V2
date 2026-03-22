"""
Tests for unit display changes in CalcView:
- H298 shown in both kcal/mol and kJ/mol
- Execution time shown in minutes, not seconds
- History table header uses "Time (min)"
"""

import io
from unittest.mock import MagicMock

import pytest
from rich.console import Console

from grimperium.cli.mock_data import PredictionResult
from grimperium.cli.views.calc_view import KCAL_TO_KJ, CalcView


@pytest.fixture
def mock_controller() -> MagicMock:
    buf = io.StringIO()
    console = Console(file=buf, highlight=False, width=120)
    ctrl = MagicMock()
    ctrl.console = console
    ctrl.current_model = "test_model"
    return ctrl


@pytest.fixture
def sample_result() -> PredictionResult:
    return PredictionResult(
        smiles="CCO",
        h298_pm7=-65.0,
        delta_correction=-1.5,
        h298_corrected=-66.5,  # * 4.184 = -278.236 kJ/mol
        model_name="delta_model",
        model_version="1.0",
        execution_time=120.0,  # 2.00 min
        n_conformers=5,
    )


def test_kcal_to_kj_constant() -> None:
    """KCAL_TO_KJ must be exactly 4.184 per IUPAC definition."""
    assert KCAL_TO_KJ == 4.184


def test_kj_conversion_displayed(
    mock_controller: MagicMock, sample_result: PredictionResult
) -> None:
    """render_result must show H298 (Corrected) in kJ/mol alongside kcal/mol."""
    view = CalcView(mock_controller)
    view.render_result(sample_result)
    output = mock_controller.console.file.getvalue()
    expected_kj = f"{sample_result.h298_corrected * 4.184:.2f} kJ/mol"
    assert expected_kj in output


def test_execution_time_in_minutes(
    mock_controller: MagicMock, sample_result: PredictionResult
) -> None:
    """render_result must display execution time in minutes, not raw seconds."""
    view = CalcView(mock_controller)
    view.render_result(sample_result)
    output = mock_controller.console.file.getvalue()
    assert "2.00 min" in output
    # Raw seconds value must not appear (execution_time=120.0 → would show "120")
    assert "120" not in output


def test_history_column_header_minutes(
    mock_controller: MagicMock, sample_result: PredictionResult
) -> None:
    """render_history must use column header 'Time (min)', not 'Time (s)'."""
    view = CalcView(mock_controller)
    view.history.append(sample_result)
    view.render_history()
    output = mock_controller.console.file.getvalue()
    assert "Time (min)" in output
    assert "Time (s)" not in output
