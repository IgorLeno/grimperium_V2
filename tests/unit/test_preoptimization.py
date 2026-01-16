"""Unit tests for xTBPreOptimizer.

Tests xTB pre-optimization for molecular structures.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

pytest.importorskip("rdkit", reason="rdkit not available")

from grimperium.crest_pm7.preoptimization import (
    ESTIMATED_TIME_PER_MOLECULE_SECONDS,
    xTBPreOptimizer,
    xTBPreOptResult,
)


class TestxTBPreOptResult:
    """Tests for xTBPreOptResult dataclass."""

    def test_init_success(self) -> None:
        """Test successful result initialization."""
        result = xTBPreOptResult(
            success=True,
            output_xyz=Path("/tmp/test.xyz"),
            error_message="",
            time_seconds=5.0,
        )
        assert result.success is True
        assert result.output_xyz == Path("/tmp/test.xyz")
        assert result.error_message == ""
        assert result.time_seconds == 5.0

    def test_init_failure(self) -> None:
        """Test failure result initialization."""
        result = xTBPreOptResult(
            success=False,
            output_xyz=None,
            error_message="xTB not found",
            time_seconds=0.0,
        )
        assert result.success is False
        assert result.output_xyz is None
        assert result.error_message == "xTB not found"


class TestxTBPreOptimizer:
    """Tests for xTBPreOptimizer class."""

    def test_disabled_returns_original_path(self, tmp_path: Path) -> None:
        """If disabled, return input unchanged."""
        preopt = xTBPreOptimizer(enabled=False)
        input_xyz = tmp_path / "input.xyz"
        input_xyz.write_text("test")
        work_dir = tmp_path / "work"

        result = preopt.preoptimize_structure("mol001", input_xyz, work_dir)

        assert result.success is True
        assert result.output_xyz == input_xyz
        assert result.error_message == ""
        assert result.time_seconds == 0.0

    def test_missing_input_returns_error(self, tmp_path: Path) -> None:
        """Handle missing XYZ gracefully."""
        preopt = xTBPreOptimizer(enabled=True)
        missing_xyz = tmp_path / "nonexistent.xyz"
        work_dir = tmp_path / "work"

        result = preopt.preoptimize_structure("mol001", missing_xyz, work_dir)

        assert result.success is False
        assert result.output_xyz is None
        assert "not found" in result.error_message
        assert result.time_seconds == 0.0

    def test_estimated_time_zero_when_disabled(self) -> None:
        """Zero hours if disabled."""
        preopt = xTBPreOptimizer(enabled=False)
        estimated = preopt.get_estimated_time_hours(100)
        assert estimated == 0.0

    def test_estimated_time_calculation(self) -> None:
        """~15s/molecule if enabled."""
        preopt = xTBPreOptimizer(enabled=True)
        n_molecules = 100
        expected_hours = (n_molecules * ESTIMATED_TIME_PER_MOLECULE_SECONDS) / 3600.0

        estimated = preopt.get_estimated_time_hours(n_molecules)

        assert estimated == expected_hours
        assert estimated == pytest.approx(0.4167, rel=0.01)

    def test_init_defaults(self) -> None:
        """Test default initialization."""
        preopt = xTBPreOptimizer()
        assert preopt.enabled is False
        assert preopt.timeout_seconds == 300

    def test_init_custom_values(self) -> None:
        """Test custom initialization."""
        preopt = xTBPreOptimizer(enabled=True, timeout_seconds=600)
        assert preopt.enabled is True
        assert preopt.timeout_seconds == 600


class TestxTBPreOptimizerExecution:
    """Tests for xTB execution (with mocking)."""

    @pytest.fixture
    def preopt(self) -> xTBPreOptimizer:
        """Create enabled pre-optimizer."""
        return xTBPreOptimizer(enabled=True, timeout_seconds=60)

    def test_xtb_not_installed(self, preopt: xTBPreOptimizer, tmp_path: Path) -> None:
        """Handle FileNotFoundError when xTB not installed."""
        input_xyz = tmp_path / "input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")
        work_dir = tmp_path / "work"

        preopt._xtb_command = ["nonexistent_xtb_binary"]
        result = preopt.preoptimize_structure("mol001", input_xyz, work_dir)

        assert result.success is False
        assert "not found" in result.error_message

    @patch("grimperium.crest_pm7.preoptimization.subprocess.run")
    def test_xtb_timeout(
        self, mock_run: MagicMock, preopt: xTBPreOptimizer, tmp_path: Path
    ) -> None:
        """Handle TimeoutExpired."""
        import subprocess

        input_xyz = tmp_path / "input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")
        work_dir = tmp_path / "work"

        mock_run.side_effect = subprocess.TimeoutExpired(cmd="xtb", timeout=60)
        result = preopt.preoptimize_structure("mol001", input_xyz, work_dir)

        assert result.success is False
        assert "timed out" in result.error_message

    @patch("grimperium.crest_pm7.preoptimization.subprocess.run")
    def test_xtb_success(
        self, mock_run: MagicMock, preopt: xTBPreOptimizer, tmp_path: Path
    ) -> None:
        """Test successful xTB execution."""
        input_xyz = tmp_path / "input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")
        work_dir = tmp_path / "work"
        work_dir.mkdir(parents=True, exist_ok=True)

        xtbopt_file = work_dir / "xtbopt.xyz"
        xtbopt_file.write_text("3\noptimized\nC 0 0 0\nH 1.1 0 0\nH 0 1.1 0\n")

        mock_run.return_value = MagicMock(returncode=0, stderr="")
        result = preopt.preoptimize_structure("mol001", input_xyz, work_dir)

        assert result.success is True
        assert result.output_xyz == work_dir / "mol001_preopt.xyz"
        assert result.output_xyz.exists()

    @patch("grimperium.crest_pm7.preoptimization.subprocess.run")
    def test_xtb_failure_nonzero_returncode(
        self, mock_run: MagicMock, preopt: xTBPreOptimizer, tmp_path: Path
    ) -> None:
        """Test xTB failure with non-zero return code."""
        input_xyz = tmp_path / "input.xyz"
        input_xyz.write_text("3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")
        work_dir = tmp_path / "work"
        work_dir.mkdir(parents=True, exist_ok=True)

        mock_run.return_value = MagicMock(returncode=1, stderr="SCF failed")
        result = preopt.preoptimize_structure("mol001", input_xyz, work_dir)

        assert result.success is False
        assert "xTB failed" in result.error_message
