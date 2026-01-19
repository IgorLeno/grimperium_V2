"""
Tests for BatchCSVManager.

Focus on FIX #1 (mol_id ordering), FIX #2 (batch naming), FIX #3 (timestamp format).
"""

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import BatchSortingStrategy, MoleculeStatus


@pytest.fixture
def csv_path(tmp_path: Path) -> Path:
    """Create a temporary CSV file for testing."""
    csv_file = tmp_path / "test_batch.csv"
    return csv_file


@pytest.fixture
def sample_csv_content() -> str:
    """Sample CSV content with 3 molecules."""
    return """mol_id,smiles,nheavy,nrotbonds,status,batch_id,batch_order,batch_failure_policy,retry_count,last_error_message,crest_status,crest_conformers_generated,crest_time,crest_error,num_conformers_selected,most_stable_hof,quality_grade,success,error_message,total_execution_time,actual_crest_timeout_used,actual_mopac_timeout_used,delta_e_12,delta_e_13,delta_e_15,timestamp,tpsa,aromatic_rings,has_heteroatoms,reference_hof,crest_v3,crest_quick,crest_nci,crest_gfnff,crest_ewin,crest_rthr,crest_optlev,crest_threads,crest_xtb_preopt,mopac_precise,mopac_scfcrt,mopac_itry,mopac_pulay,mopac_prtall,mopac_archive,max_retries
id_00001,C,1,0,PENDING,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,3
id_00002,CC,2,1,PENDING,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,3
id_00003,CCC,3,2,PENDING,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,3
"""


@pytest.fixture
def manager_with_csv(csv_path: Path, sample_csv_content: str) -> BatchCSVManager:
    """Create a BatchCSVManager with sample CSV loaded."""
    csv_path.write_text(sample_csv_content)
    manager = BatchCSVManager(csv_path)
    manager.load_csv()
    return manager


class TestMolIdOrdering:
    """Tests for FIX #1: mol_id ordering in RERUN_FIRST_THEN_EASY strategy."""

    def test_apply_sorting_strategy_respects_mol_id_order(
        self, manager_with_csv: BatchCSVManager
    ):
        """Verify that RERUN_FIRST_THEN_EASY sorts by mol_id, not nheavy."""
        # Create test DataFrame with mixed RERUN/PENDING status
        # mol_ids are out of order, nheavy is deliberately inverse to mol_id
        test_df = pd.DataFrame(
            {
                "mol_id": ["id_00003", "id_00001", "id_00002"],  # Out of order
                "status": ["PENDING", "RERUN", "PENDING"],
                "nheavy": [15, 10, 20],  # Deliberately inverse to mol_id
            }
        )

        result = manager_with_csv._apply_sorting_strategy(
            test_df, BatchSortingStrategy.RERUN_FIRST_THEN_EASY
        )

        # Verify order: RERUN first (id_00001), then PENDING (id_00002, id_00003)
        assert result["mol_id"].tolist() == ["id_00001", "id_00002", "id_00003"]
        assert result["status"].tolist() == ["RERUN", "PENDING", "PENDING"]

    def test_apply_sorting_strategy_no_nheavy_dependency(
        self, manager_with_csv: BatchCSVManager
    ):
        """Verify sorting is NOT based on nheavy column."""
        # Create test where nheavy order differs from mol_id order
        test_df = pd.DataFrame(
            {
                "mol_id": ["id_00001", "id_00002", "id_00003"],
                "status": ["PENDING", "PENDING", "PENDING"],
                "nheavy": [20, 10, 15],  # Different order than mol_id
            }
        )

        result = manager_with_csv._apply_sorting_strategy(
            test_df, BatchSortingStrategy.RERUN_FIRST_THEN_EASY
        )

        # Result must follow mol_id, not nheavy
        assert result["mol_id"].tolist() == ["id_00001", "id_00002", "id_00003"]


class TestSequentialBatchNaming:
    """Tests for FIX #2: Sequential batch numbering."""

    def test_generate_batch_id_sequential(
        self, manager_with_csv: BatchCSVManager
    ):
        """Verify batch IDs are sequential: batch_0001, batch_0002, ..."""
        # Clear existing batch IDs
        manager_with_csv.df["batch_id"] = None
        manager_with_csv.save_csv()

        # Generate first batch
        batch_id_1 = manager_with_csv.generate_batch_id()
        assert batch_id_1 == "batch_0001"

        # Manually add to DataFrame
        manager_with_csv.df.at[0, "batch_id"] = batch_id_1
        manager_with_csv.save_csv()

        # Generate second batch
        batch_id_2 = manager_with_csv.generate_batch_id()
        assert batch_id_2 == "batch_0002"

        # Manually add to DataFrame
        manager_with_csv.df.at[1, "batch_id"] = batch_id_2
        manager_with_csv.save_csv()

        # Generate third batch
        batch_id_3 = manager_with_csv.generate_batch_id()
        assert batch_id_3 == "batch_0003"

    def test_generate_batch_id_no_timestamp(
        self, manager_with_csv: BatchCSVManager
    ):
        """Verify batch ID format is batch_NNNN, not timestamp-based."""
        manager_with_csv.df["batch_id"] = None
        manager_with_csv.save_csv()

        batch_id = manager_with_csv.generate_batch_id()

        # Must NOT contain timestamp patterns
        assert batch_id.count("_") == 1  # Only one underscore (batch_NNNN)
        assert "T" not in batch_id  # No ISO format
        assert batch_id.startswith("batch_")
        assert batch_id[-4:].isdigit()  # Last 4 chars are digits
        assert len(batch_id) == len("batch_0001")  # Exact length


class TestTimestampFormat:
    """Tests for FIX #3: Simplified timestamp format."""

    def test_timestamp_format_dd_mm_hh_mm(
        self, manager_with_csv: BatchCSVManager
    ):
        """Verify timestamp format is dd/mm-HH:MM, not ISO."""

        # Create mock PM7Result with known timestamp
        class MockPM7Result:
            timestamp = datetime(2026, 1, 19, 14, 18, 20, 676511, tzinfo=timezone.utc)
            crest_status = None
            crest_conformers_generated = 0
            crest_time = None
            crest_error = None
            num_conformers_selected = 0
            most_stable_hof = None
            quality_grade = None
            success = False
            error_message = None
            total_execution_time = None
            delta_e_12 = None
            delta_e_13 = None
            delta_e_15 = None

        result = MockPM7Result()
        csv_update = manager_with_csv.pm7result_to_csv_update(
            mol_id="id_00001",
            result=result,
            batch_id="batch_0001",
            batch_order=1,
            crest_timeout_used=30.0,
            mopac_timeout_used=10.0,
        )

        # Verify timestamp format
        assert csv_update["timestamp"] == "19/01-14:18"
        # Must NOT be ISO format
        assert "T" not in csv_update["timestamp"]
        assert "+00:00" not in csv_update["timestamp"]
        assert "." not in csv_update["timestamp"]  # No fractional seconds
