"""Tests for calc view module."""

from unittest.mock import MagicMock

import pytest

from grimperium.cli.views.calc_view import CalcView


@pytest.fixture
def mock_controller() -> MagicMock:
    """Create mock controller for CalcView."""
    controller = MagicMock()
    controller.current_model = "DeltaXGB_v1.0"
    return controller


@pytest.fixture
def calc_view(mock_controller: MagicMock) -> CalcView:
    """Create CalcView instance."""
    return CalcView(mock_controller)


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
        assert isinstance(result, str), f"Invalid SMILES '{smiles}' should return error message"
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
    assert hasattr(calc_view, "validate_molecules"), (
        "CalcView should have validate_molecules method for batch operations"
    )


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

    assert "rejected" in summary.lower() or "2" in summary, (
        "Summary should mention rejected molecules"
    )


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
