"""
Tests for BatchOrchestrator.

Includes tests for v2.2 visibility features.
"""

from unittest.mock import MagicMock

import pytest

from grimperium.core.batch_orchestrator import (
    BatchOrchestrator,
    BatchOrchestratorError,
    BatchSummary,
    CalculationSettings,
)


@pytest.fixture
def valid_csv_content():
    """Minimal valid CSV content."""
    return (
        "mol_id,smiles,nheavy,status\n"
        "mol_001,C,1,Pending\n"
        "mol_002,CC,2,Pending\n"
        "mol_003,CCC,3,OK\n"
    )


@pytest.fixture
def orchestrator_with_csv(tmp_path, valid_csv_content):
    """Create orchestrator with valid CSV."""
    csv_file = tmp_path / "test.csv"
    csv_file.write_text(valid_csv_content)

    orch = BatchOrchestrator(csv_file, strict=False)
    orch.log_dir = tmp_path / ".grimperium"
    orch.validation_log = orch.log_dir / "validation_errors.log"
    orch.summary_log = orch.log_dir / "batch_summary.log"

    return orch


class TestCalculationSettings:
    """Tests for CalculationSettings."""

    def test_defaults(self):
        """Default settings are correct."""
        settings = CalculationSettings()

        assert settings.max_reruns == 3
        assert settings.crest_timeout_minutes == 30
        assert settings.mopac_timeout_minutes == 10

    def test_custom_values(self):
        """Custom settings are applied."""
        settings = CalculationSettings(
            max_reruns=5,
            crest_timeout_minutes=60,
        )

        assert settings.max_reruns == 5
        assert settings.crest_timeout_minutes == 60


class TestBatchSummary:
    """Tests for BatchSummary."""

    def test_to_dict(self):
        """to_dict returns all fields."""
        summary = BatchSummary(total=100, complete=80, errors=5)

        d = summary.to_dict()

        assert d["total"] == 100
        assert d["complete"] == 80
        assert d["errors"] == 5


class TestBatchOrchestratorInit:
    """Tests for BatchOrchestrator initialization."""

    def test_init_with_defaults(self, tmp_path, valid_csv_content):
        """Initialize with default settings."""
        csv_file = tmp_path / "test.csv"
        csv_file.write_text(valid_csv_content)

        orch = BatchOrchestrator(csv_file)

        assert orch.strict is False
        assert orch.settings.max_reruns == 3

    def test_init_with_custom_settings(self, tmp_path, valid_csv_content):
        """Initialize with custom settings."""
        csv_file = tmp_path / "test.csv"
        csv_file.write_text(valid_csv_content)

        settings = CalculationSettings(max_reruns=5)
        orch = BatchOrchestrator(csv_file, settings=settings)

        assert orch.settings.max_reruns == 5


class TestBatchOrchestratorRun:
    """Tests for BatchOrchestrator.run()."""

    def test_run_returns_summary(self, orchestrator_with_csv):
        """run() returns BatchSummary."""
        summary = orchestrator_with_csv.run(dry_run=True)

        assert isinstance(summary, BatchSummary)
        assert summary.total == 3

    def test_dry_run_does_not_process(self, orchestrator_with_csv):
        """dry_run=True doesn't actually process."""
        summary = orchestrator_with_csv.run(dry_run=True)

        # Pending molecules should still be pending
        assert summary.pending == 2

    def test_callback_is_called(self, orchestrator_with_csv):
        """Callback is called for each molecule."""
        callback = MagicMock()

        orchestrator_with_csv.run(callback=callback, dry_run=False)

        # Callback should be called for each pending molecule (2 pending)
        assert callback.call_count == 2

    def test_run_with_missing_file(self, tmp_path):
        """run() raises error for missing file."""
        csv_file = tmp_path / "nonexistent.csv"
        orch = BatchOrchestrator(csv_file)

        with pytest.raises(BatchOrchestratorError):
            orch.run()


class TestBatchOrchestratorVisibility:
    """Tests for v2.2 visibility features."""

    def test_creates_validation_log_permissive(
        self, tmp_path
    ):
        """v2.2: Creates validation log file in permissive mode."""
        # Create CSV with invalid rows
        csv_content = (
            "mol_id,smiles,nheavy,status\n"
            "mol_001,C,1,Pending\n"
            ",CC,2,Pending\n"  # Missing mol_id
            "mol_003,CCC,3,OK\n"
        )
        csv_file = tmp_path / "test.csv"
        csv_file.write_text(csv_content)

        orch = BatchOrchestrator(csv_file, strict=False)
        orch.log_dir = tmp_path / ".grimperium"
        orch.validation_log = orch.log_dir / "validation_errors.log"
        orch.summary_log = orch.log_dir / "batch_summary.log"

        orch.run(dry_run=True)

        # Check validation log was created
        assert orch.validation_log.exists()

    def test_saves_validation_log(self, tmp_path):
        """v2.2: Saves validation errors to file."""
        # Create CSV with invalid rows
        csv_content = (
            "mol_id,smiles,nheavy,status\n"
            "mol_001,C,1,Pending\n"
            ",CC,2,Pending\n"  # Missing mol_id
        )
        csv_file = tmp_path / "test.csv"
        csv_file.write_text(csv_content)

        orch = BatchOrchestrator(csv_file, strict=False)
        orch.log_dir = tmp_path / ".grimperium"
        orch.validation_log = orch.log_dir / "validation_errors.log"
        orch.summary_log = orch.log_dir / "batch_summary.log"

        orch.run(dry_run=True)

        # Check validation log exists and has content
        assert orch.validation_log.exists()
        content = orch.validation_log.read_text()
        assert "Validation Errors" in content

    def test_saves_batch_summary_log(self, orchestrator_with_csv):
        """v2.2: Saves batch summary to file (append)."""
        orchestrator_with_csv.run(dry_run=True)

        assert orchestrator_with_csv.summary_log.exists()
        content = orchestrator_with_csv.summary_log.read_text()
        assert "Total:" in content

    def test_batch_summary_log_appends(self, orchestrator_with_csv):
        """v2.2: Batch summary log appends (not overwrites)."""
        # Run twice
        orchestrator_with_csv.run(dry_run=True)
        orchestrator_with_csv.run(dry_run=True)

        content = orchestrator_with_csv.summary_log.read_text()

        # Should have two entries (count separators)
        assert content.count("=" * 70) >= 4  # 2 entries × 2 separators each

    def test_get_validation_report(self, tmp_path):
        """get_validation_report returns report from data manager."""
        # Create CSV with invalid rows
        csv_content = (
            "mol_id,smiles,nheavy,status\n"
            "mol_001,C,1,Pending\n"
            ",CC,2,Pending\n"  # Missing mol_id
        )
        csv_file = tmp_path / "test.csv"
        csv_file.write_text(csv_content)

        orch = BatchOrchestrator(csv_file, strict=False)
        orch.log_dir = tmp_path / ".grimperium"
        orch.validation_log = orch.log_dir / "validation_errors.log"
        orch.summary_log = orch.log_dir / "batch_summary.log"

        orch.run(dry_run=True)

        report = orch.get_validation_report()
        assert report.total_errors == 1


class TestBatchOrchestratorNoMolecules:
    """Tests for handling empty batches."""

    def test_no_molecules_to_process(self, tmp_path):
        """Handles case with no molecules to process."""
        csv_content = (
            "mol_id,smiles,nheavy,status\n"
            "mol_001,C,1,OK\n"  # Already complete
            "mol_002,CC,2,OK\n"
        )
        csv_file = tmp_path / "test.csv"
        csv_file.write_text(csv_content)

        orch = BatchOrchestrator(csv_file, strict=False)
        orch.log_dir = tmp_path / ".grimperium"
        orch.validation_log = orch.log_dir / "validation_errors.log"
        orch.summary_log = orch.log_dir / "batch_summary.log"

        summary = orch.run(dry_run=True)

        assert summary.total == 2
        assert summary.pending == 0  # No pending molecules
