"""Unit tests for CREST PM7 data models.

Tests ConformerData, PM7Result, and related classes.
"""

import pytest

# Skip all tests in this module if rdkit is not available
pytest.importorskip("rdkit", reason="rdkit not available")

from grimperium.crest_pm7.config import (
    CRESTStatus,
    HOFConfidence,
    MOPACStatus,
    QualityGrade,
)
from grimperium.crest_pm7.molecule_processor import (
    ConformerData,
    PM7Result,
    _collect_issues,
    _grade_from_issues,
)


class TestConformerData:
    """Tests for ConformerData dataclass."""

    def test_init_minimal(self) -> None:
        """Test minimal initialization."""
        conf = ConformerData(index=0, mol_id="TEST_001")
        assert conf.index == 0
        assert conf.mol_id == "TEST_001"
        assert conf.crest_status == CRESTStatus.NOT_ATTEMPTED
        assert conf.mopac_status == MOPACStatus.NOT_ATTEMPTED
        assert conf.energy_hof is None

    def test_is_successful_true(self) -> None:
        """Test is_successful when all conditions met."""
        conf = ConformerData(
            index=0,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=-50.5,
        )
        assert conf.is_successful is True

    def test_is_successful_false_mopac_failed(self) -> None:
        """Test is_successful when MOPAC failed."""
        conf = ConformerData(
            index=0,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.ERROR,  # Changed: FAILED -> ERROR
            hof_extraction_successful=True,
            energy_hof=-50.5,
        )
        assert conf.is_successful is False

    def test_is_successful_false_no_hof(self) -> None:
        """Test is_successful when HOF extraction failed."""
        conf = ConformerData(
            index=0,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=False,
            energy_hof=None,
        )
        assert conf.is_successful is False

    def test_to_dict(self) -> None:
        """Test serialization to dict."""
        conf = ConformerData(
            index=0,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.SUCCESS,
            energy_hof=-50.5,
            hof_confidence=HOFConfidence.HIGH,
        )
        result = conf.to_dict()

        assert result["index"] == 0
        assert result["mol_id"] == "TEST_001"
        assert result["energy_hof"] == -50.5
        assert "mopac_status" in result


class TestPM7Result:
    """Tests for PM7Result dataclass."""

    def test_init_empty(self) -> None:
        """Test initialization with no conformers."""
        result = PM7Result(
            mol_id="TEST_001",
            smiles="CCO",
            nheavy=3,
        )
        assert result.mol_id == "TEST_001"
        assert result.smiles == "CCO"
        assert result.nheavy == 3
        assert result.conformers == []
        assert (
            result.quality_grade == QualityGrade.FAILED
        )  # Changed: grade->quality_grade, C->FAILED (default)
        assert result.crest_time is None  # Changed: 0.0 -> None (default)

    def test_most_stable_hof_empty(self) -> None:
        """Test most_stable_hof with no conformers."""
        result = PM7Result(
            mol_id="TEST_001",
            smiles="CCO",
            nheavy=3,
        )
        assert result.most_stable_hof is None

    def test_most_stable_hof_with_conformers(self) -> None:
        """Test most_stable_hof with conformers."""
        conf1 = ConformerData(
            index=0,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=-50.0,
        )
        conf2 = ConformerData(
            index=1,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=-55.0,  # More negative = more stable
        )
        conf3 = ConformerData(
            index=2,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.ERROR,  # Changed: FAILED -> ERROR
            hof_extraction_successful=False,
            energy_hof=None,
        )

        result = PM7Result(
            mol_id="TEST_001",
            smiles="CCO",
            nheavy=3,
            conformers=[conf1, conf2, conf3],
        )

        # Most stable should be the most negative
        assert result.most_stable_hof == pytest.approx(-55.0)

    def test_successful_conformers(self) -> None:
        """Test successful_conformers property."""
        conf1 = ConformerData(
            index=0,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=-50.0,
        )
        conf2 = ConformerData(
            index=1,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.ERROR,  # Changed: FAILED -> ERROR
            hof_extraction_successful=False,
            energy_hof=None,
        )
        conf3 = ConformerData(
            index=2,
            mol_id="TEST_001",
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=-52.0,
        )

        result = PM7Result(
            mol_id="TEST_001",
            smiles="CCO",
            nheavy=3,
            conformers=[conf1, conf2, conf3],
        )

        successful = result.successful_conformers
        assert len(successful) == 2
        assert all(c.is_successful for c in successful)

    def test_to_dict(self) -> None:
        """Test serialization to dict."""
        result = PM7Result(
            mol_id="TEST_001",
            smiles="CCO",
            nheavy=3,
            quality_grade=QualityGrade.A,  # Changed: grade -> quality_grade
            crest_time=10.5,
            total_execution_time=45.2,  # Changed: total_mopac_time -> total_execution_time
        )
        data = result.to_dict()

        assert data["mol_id"] == "TEST_001"
        assert data["smiles"] == "CCO"
        assert data["nheavy"] == 3
        assert data["quality_grade"] == "A"  # Changed: grade -> quality_grade
        assert data["crest_time"] == pytest.approx(10.5)


class TestGradingFunctions:
    """Tests for grading helper functions."""

    def test_collect_issues_empty(self) -> None:
        """Test _collect_issues with perfect result - returns empty list for success."""
        conf = ConformerData(
            index=0,
            mol_id="TEST_001",
            crest_status=CRESTStatus.SUCCESS,
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=-50.0,
            hof_confidence=HOFConfidence.HIGH,
        )
        result = PM7Result(
            mol_id="TEST_001",
            smiles="CCO",
            nheavy=3,
            conformers=[conf],
            success=True,  # Required for _collect_issues to return detailed issues
        )

        issues = _collect_issues(result)
        assert isinstance(issues, list)  # Changed: dict -> list

    def test_grade_from_issues_A(self) -> None:
        """Test grade A assignment."""
        # No issues + success + conformers = Grade A
        issues: list[str] = []  # Changed: dict -> list
        grade = _grade_from_issues(issues, success=True, has_conformers=True)
        assert grade == QualityGrade.A

    def test_grade_from_issues_C(self) -> None:
        """Test grade C assignment for multiple issues."""
        # Multiple issues = Grade C
        issues = [
            "no_high_confidence_hof",
            "incomplete_conformer_coverage",
        ]  # 2+ issues
        grade = _grade_from_issues(issues, success=True, has_conformers=True)
        assert grade == QualityGrade.C


class TestQualityGradeEnum:
    """Tests for QualityGrade enum."""

    def test_grade_values(self) -> None:
        """Test grade enum values."""
        assert QualityGrade.A.value == "A"
        assert QualityGrade.B.value == "B"
        assert QualityGrade.C.value == "C"

    def test_grade_comparison(self) -> None:
        """Test grade comparison."""
        # Can compare string values
        assert QualityGrade.A.value < QualityGrade.B.value
        assert QualityGrade.B.value < QualityGrade.C.value


class TestStatusEnums:
    """Tests for status enums."""

    def test_crest_status_values(self) -> None:
        """Test CREST status enum values."""
        assert CRESTStatus.NOT_ATTEMPTED.value == "NOT_ATTEMPTED"
        assert CRESTStatus.SUCCESS.value == "SUCCESS"
        assert CRESTStatus.FAILED.value == "FAILED"
        # Note: CRESTStatus does not have TIMEOUT (removed in API update)

    def test_mopac_status_values(self) -> None:
        """Test MOPAC status enum values."""
        assert MOPACStatus.NOT_ATTEMPTED.value == "NOT_ATTEMPTED"
        assert MOPACStatus.SUCCESS.value == "SUCCESS"
        assert MOPACStatus.ERROR.value == "ERROR"  # Changed: FAILED -> ERROR
        assert MOPACStatus.TIMEOUT.value == "TIMEOUT"
        assert MOPACStatus.SCF_FAILED.value == "SCF_FAILED"
        assert MOPACStatus.GEOMETRY_ERROR.value == "GEOMETRY_ERROR"  # Added: new status

    def test_hof_confidence_values(self) -> None:
        """Test HOF confidence enum values."""
        assert HOFConfidence.HIGH.value == "HIGH"
        assert HOFConfidence.MEDIUM.value == "MEDIUM"
        assert HOFConfidence.LOW.value == "LOW"
