"""Tests for calc view module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from grimperium.cli.calc_pipeline import (
    CalcPipelineError,
    CalcPipelineResult,
    run_single_molecule_prediction,
)
from grimperium.cli.controller import CliController
from grimperium.cli.mock_data import PredictionResult
from grimperium.cli.views.calc_view import CalcView
from grimperium.crest_pm7.config import PM7Config


@pytest.fixture
def mock_controller() -> MagicMock:
    """Create mock controller for CalcView."""
    controller = MagicMock()
    controller.current_model = "DeltaXGB_v1.0"
    controller.current_model_path = None
    return controller


@pytest.fixture
def calc_view(mock_controller: MagicMock) -> CalcView:
    """Create CalcView instance."""
    return CalcView(mock_controller)


# ──────────────────────────────────────────────────
# Tests: do_prediction pipeline integration
# ──────────────────────────────────────────────────


def _mock_pm7_delta_method() -> MagicMock:
    """Build a minimal method mock for do_prediction tests."""
    mock_method = MagicMock()
    mock_method.method_id = "pm7_delta_learning"
    mock_method.property_name = "Standard enthalpy of formation"
    mock_method.display_name = "PM7 + Delta Learning"
    mock_method.model_requirement.model_required = True
    mock_method.xtb.optional = True
    mock_method.xtb.enabled_by_default = False
    return mock_method


@patch("grimperium.cli.views.calc_view.text_input", return_value="CCO")
def test_do_prediction_sem_modelo(
    mock_input: MagicMock,
    calc_view: CalcView,
) -> None:
    """No model on disk → inline selection offers error, pipeline never executes.

    When no session model is set and no .joblib files exist,
    _select_model_inline shows 'No trained models found' via show_error
    and returns None so the pipeline is never reached.
    """
    mock_method = _mock_pm7_delta_method()
    calc_view._select_method = MagicMock(return_value=mock_method)
    calc_view._select_units = MagicMock(return_value="both")
    calc_view.render = MagicMock()
    calc_view.show_error = MagicMock()

    with (
        patch(
            "grimperium.cli.views.calc_view.discover_available_models",
            return_value=[],
        ),
        patch(
            "grimperium.cli.views.calc_view.run_single_molecule_prediction"
        ) as mock_pipeline,
    ):
        result = calc_view.do_prediction()

    calc_view.show_error.assert_called_once()
    error_msg = calc_view.show_error.call_args[0][0]
    assert "model" in error_msg.lower() or "Model" in error_msg
    mock_pipeline.assert_not_called()
    assert result is True


@patch("grimperium.cli.views.calc_view.text_input", return_value="CCO")
def test_do_prediction_crest_falha(
    mock_input: MagicMock,
    calc_view: CalcView,
) -> None:
    """CREST failure → CalcPipelineError → show_error."""
    mock_method = _mock_pm7_delta_method()
    calc_view._select_method = MagicMock(return_value=mock_method)
    calc_view._select_units = MagicMock(return_value="both")
    calc_view._resolve_required_model = MagicMock(  # type: ignore[method-assign]
        return_value=Path("/fake/model.joblib"),
    )
    calc_view.render = MagicMock()
    calc_view.show_error = MagicMock()

    with patch(
        "grimperium.cli.views.calc_view.run_single_molecule_prediction",
        side_effect=CalcPipelineError("Pipeline failed: CREST timeout"),
    ):
        result = calc_view.do_prediction()

    calc_view.show_error.assert_called_once()
    error_msg = calc_view.show_error.call_args[0][0]
    assert "CREST" in error_msg
    assert result is True
    assert len(calc_view.history) == 0


@patch("grimperium.cli.views.calc_view.text_input", return_value="CCO")
def test_do_prediction_mopac_falha(
    mock_input: MagicMock,
    calc_view: CalcView,
) -> None:
    """MOPAC failure → CalcPipelineError → show_error."""
    mock_method = _mock_pm7_delta_method()
    calc_view._select_method = MagicMock(return_value=mock_method)
    calc_view._select_units = MagicMock(return_value="both")
    calc_view._resolve_required_model = MagicMock(  # type: ignore[method-assign]
        return_value=Path("/fake/model.joblib"),
    )
    calc_view.render = MagicMock()
    calc_view.show_error = MagicMock()

    with patch(
        "grimperium.cli.views.calc_view.run_single_molecule_prediction",
        side_effect=CalcPipelineError("No successful MOPAC output found"),
    ):
        result = calc_view.do_prediction()

    calc_view.show_error.assert_called_once()
    error_msg = calc_view.show_error.call_args[0][0]
    assert "MOPAC" in error_msg
    assert result is True
    assert len(calc_view.history) == 0


@patch("grimperium.cli.views.calc_view.text_input", return_value="CCO")
def test_do_prediction_sucesso_completo(
    mock_input: MagicMock,
    calc_view: CalcView,
) -> None:
    """Full pipeline success → PredictionResult with all real fields."""
    mock_method = _mock_pm7_delta_method()
    calc_view._select_method = MagicMock(return_value=mock_method)
    calc_view._select_units = MagicMock(return_value="both")
    calc_view._resolve_required_model = MagicMock(  # type: ignore[method-assign]
        return_value=Path("/fake/model.joblib"),
    )
    calc_view.render = MagicMock()
    calc_view.render_result = MagicMock()

    fake_result = CalcPipelineResult(
        h298_pm7=-45.32,
        delta_correction=2.15,
        h298_corrected=-43.17,
        n_conformers=7,
        execution_time=120.5,
        model_version="v2.1",
    )

    with patch(
        "grimperium.cli.views.calc_view.run_single_molecule_prediction",
        return_value=fake_result,
    ):
        result = calc_view.do_prediction()

    assert result is True
    assert len(calc_view.history) == 1

    pred = calc_view.history[0]
    assert isinstance(pred, PredictionResult)
    assert pred.smiles == "CCO"
    assert pred.h298_pm7 == -45.32
    assert pred.delta_correction == 2.15
    assert isinstance(pred.delta_correction, float)
    assert pred.h298_corrected == -43.17
    assert pred.n_conformers == 7
    assert pred.execution_time == 120.5
    assert pred.model_version == "v2.1"
    assert pred.model_name == "DeltaXGB_v1.0"
    # molecule_name must NOT exist
    assert not hasattr(pred, "molecule_name")


def test_run_single_molecule_prediction_delta_correction_scalar() -> None:
    """Pipeline extracts pure delta from the final H298 prediction."""

    class FakeLearner:
        def __init__(self) -> None:
            self.y_pm7_seen: np.ndarray | None = None

        def predict(self, X: np.ndarray, y_pm7: np.ndarray) -> np.ndarray:
            self.y_pm7_seen = y_pm7
            delta_pred = np.array([3.24])
            return y_pm7 + delta_pred

    fake_learner = FakeLearner()
    fake_feature_pipeline = MagicMock()
    fake_feature_pipeline.transform.return_value = np.array([[1.0, 2.0, 3.0]])

    fake_conformer = MagicMock()
    fake_conformer.mopac_output_file = Path("/fake/mol.out")

    fake_pm7_result = MagicMock()
    fake_pm7_result.success = True
    fake_pm7_result.error_message = None
    fake_pm7_result.most_stable_hof = -72.76
    fake_pm7_result.nheavy = 2
    fake_pm7_result.crest_conformers_generated = 7
    fake_pm7_result.num_conformers_selected = 3
    fake_pm7_result.crest_time = 12.0
    fake_pm7_result.total_execution_time = 30.0
    fake_pm7_result.get_selected_conformer.return_value = fake_conformer

    with (
        patch("grimperium.cli.calc_pipeline.CRESTPM7Pipeline") as mock_pipeline_cls,
        patch(
            "grimperium.cli.calc_pipeline.load_model",
            return_value=(fake_learner, fake_feature_pipeline),
        ),
        patch(
            "grimperium.cli.calc_pipeline.load_model_metadata",
            return_value={"version": "v2.1"},
        ),
        patch(
            "grimperium.cli.calculation_features.extract_all_rdkit_descriptors",
            return_value={"rdkit_nrotbonds": 1.0},
        ),
        patch(
            "grimperium.cli.calculation_features.extract_mopac_descriptors",
            return_value={"mopac_homo_ev": -10.5},
        ),
    ):
        mock_pipeline_cls.return_value.process_molecule.return_value = fake_pm7_result

        result = run_single_molecule_prediction(
            smiles="CCO",
            mol_id="calc_123456",
            model_path=Path("/fake/model.joblib"),
            config=PM7Config(),
        )

    assert fake_learner.y_pm7_seen is not None
    assert fake_learner.y_pm7_seen.shape == (1,)
    assert isinstance(result.delta_correction, float)
    assert not isinstance(result.delta_correction, np.ndarray)
    assert isinstance(result.h298_corrected, float)
    assert result.h298_pm7 == pytest.approx(-72.76)
    assert result.delta_correction == pytest.approx(3.24)
    assert -10.0 < result.delta_correction < 10.0
    assert result.h298_corrected == pytest.approx(-69.52)
    assert result.h298_corrected != pytest.approx(result.h298_pm7 + result.h298_pm7)


def test_set_model_com_path() -> None:
    """set_model() with path updates current_model_path."""
    controller = CliController()
    model_path = Path("/models/delta_v2.joblib")

    controller.set_model("DeltaRF_v2.0", model_path)

    assert controller.current_model == "DeltaRF_v2.0"
    assert controller.current_model_path == model_path


def test_get_model_path_usa_controller() -> None:
    """_resolve_model_path() prioritizes controller path over env var."""
    controller = MagicMock()
    controller.current_model = "DeltaXGB_v1.0"
    controller_path = Path(__file__)  # use this test file (guaranteed to exist)
    controller.current_model_path = controller_path

    view = CalcView(controller)

    with patch.dict("os.environ", {"GRIMPERIUM_MODEL_PATH": "/env/model.joblib"}):
        resolved = view._resolve_model_path()

    assert resolved == controller_path


def test_validate_smiles_valid_molecules(calc_view: CalcView) -> None:
    """
    Bug #3: Test RDKit validation accepts valid SMILES.

    Valid SMILES should pass validation and return True.
    """
    valid_smiles = [
        "CCO",  # Ethanol
        "CC(=O)O",  # Acetic acid
        "c1ccccc1",  # Benzene
        "C",  # Methane
        "CC",  # Ethane
        "CC(C)C",  # Isobutane
    ]

    for smiles in valid_smiles:
        result = calc_view.validate_smiles(smiles)
        assert result is True, f"Valid SMILES '{smiles}' should pass validation"


def test_validate_smiles_invalid_molecules(calc_view: CalcView) -> None:
    """
    Bug #3: Test RDKit validation rejects invalid SMILES.

    Invalid SMILES should return error message string.
    """
    invalid_smiles = [
        "C(C",  # Unmatched parenthesis
        "CC==O",  # Invalid double bond
        "C1CCC",  # Unclosed ring
        "[Xx]",  # Invalid element
        "C(C)(C)(C)(C)(C)C",  # Carbon with too many bonds
    ]

    for smiles in invalid_smiles:
        result = calc_view.validate_smiles(smiles)
        assert isinstance(
            result, str
        ), f"Invalid SMILES '{smiles}' should return error message"
        assert len(result) > 0, "Error message should not be empty"


def test_validate_smiles_empty_string(calc_view: CalcView) -> None:
    """
    Bug #3: Test validation rejects empty strings.
    """
    result = calc_view.validate_smiles("")
    assert isinstance(result, str)
    assert "empty" in result.lower() or "valid" in result.lower()


def test_validate_smiles_whitespace_only(calc_view: CalcView) -> None:
    """
    Bug #3: Test validation rejects whitespace-only strings.
    """
    result = calc_view.validate_smiles("   ")
    assert isinstance(result, str)
    assert len(result) > 0


def test_validate_smiles_strips_whitespace(calc_view: CalcView) -> None:
    """
    Bug #3: Test validation handles whitespace around valid SMILES.

    Should strip whitespace and validate the inner SMILES.
    """
    result = calc_view.validate_smiles("  CCO  ")
    assert result is True, "Valid SMILES with whitespace should pass"


def test_validate_smiles_duplicate_detection() -> None:
    """
    Bug #3: Test that duplicate SMILES are tracked in prediction history.

    This is a future enhancement - just verify history is maintained.
    """
    controller = MagicMock()
    controller.current_model = "DeltaXGB_v1.0"
    view = CalcView(controller)

    # Should start with empty history
    assert len(view.history) == 0


def test_validate_molecules_method_exists(calc_view: CalcView) -> None:
    """
    Bug #3: Test that validate_molecules method exists for batch validation.

    This method should validate a list of molecule dictionaries.
    """
    # Check if method exists
    assert hasattr(
        calc_view, "validate_molecules"
    ), "CalcView should have validate_molecules method for batch operations"


def test_validate_molecules_filters_invalid_smiles(calc_view: CalcView) -> None:
    """
    Bug #3: Test validate_molecules filters out invalid SMILES.

    Should return only valid molecules with summary of rejections.
    """
    molecules = [
        {"smiles": "CCO", "name": "Ethanol"},
        {"smiles": "C(C", "name": "Invalid1"},  # Invalid
        {"smiles": "CC", "name": "Ethane"},
        {"smiles": "[Xx]", "name": "Invalid2"},  # Invalid
    ]

    valid, summary = calc_view.validate_molecules(molecules)

    assert len(valid) == 2, "Should return 2 valid molecules"
    assert valid[0]["smiles"] == "CCO"
    assert valid[1]["smiles"] == "CC"

    assert (
        "rejected" in summary.lower() or "2" in summary
    ), "Summary should mention rejected molecules"


def test_validate_molecules_removes_duplicates(calc_view: CalcView) -> None:
    """
    Bug #3: Test validate_molecules removes duplicate SMILES.
    """
    molecules = [
        {"smiles": "CCO", "name": "Ethanol"},
        {"smiles": "CC", "name": "Ethane"},
        {"smiles": "CCO", "name": "Ethanol_duplicate"},  # Duplicate
    ]

    valid, summary = calc_view.validate_molecules(molecules)

    assert len(valid) == 2, "Should remove duplicate SMILES"
    smiles_list = [m["smiles"] for m in valid]
    assert smiles_list.count("CCO") == 1, "CCO should appear only once"


def test_validate_molecules_empty_list(calc_view: CalcView) -> None:
    """
    Bug #3: Test validate_molecules handles empty molecule list.
    """
    valid, summary = calc_view.validate_molecules([])

    assert len(valid) == 0
    assert "no molecules" in summary.lower() or "empty" in summary.lower()
