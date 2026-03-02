"""Tests for retry/error handling in BatchCSVManager.

Focus on ensuring reruns counter is persisted correctly.
"""

from pathlib import Path

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus


@pytest.fixture
def csv_path(tmp_path: Path) -> Path:
    """Provide a temporary CSV path.

    Args:
        tmp_path: Temporary directory fixture provided by pytest.

    Returns:
        Path: Path to the temporary CSV file named "retry_test.csv".
    """
    return tmp_path / "retry_test.csv"


def _write_csv(csv_path: Path, rows: list[dict[str, object]]) -> None:
    """Write rows to a CSV file.

    Args:
        csv_path: Path to write CSV file.
        rows: List of row dictionaries.
    """
    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)


class TestRetryHandling:
    """Tests for retry/skip transitions using reruns counter."""

    def test_mark_rerun_first_attempt(self, csv_path: Path) -> None:
        """Ensure mark_rerun increments reruns counter."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_00014",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.RUNNING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                }
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        manager.mark_rerun("mol_00014", "CREST timeout after 3600s")

        df = pd.read_csv(csv_path)
        row = df.iloc[0]

        assert row["status"] == MoleculeStatus.RERUN.value
        assert int(row["reruns"]) == 1

    def test_mark_rerun_second_attempt(self, csv_path: Path) -> None:
        """Ensure mark_rerun increments reruns to 2 without skipping."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_00014",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.RUNNING.value,
                    "reruns": 1,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                }
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        manager.mark_rerun("mol_00014", "CREST timeout after 3600s")

        df = pd.read_csv(csv_path)
        row = df.iloc[0]

        assert row["status"] == MoleculeStatus.RERUN.value
        assert int(row["reruns"]) == 2

    def test_mark_rerun_third_attempt_skips(self, csv_path: Path) -> None:
        """Ensure third attempt transitions to Skip (reruns >= 3)."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_00014",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.RUNNING.value,
                    "reruns": 2,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                }
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        manager.mark_rerun("mol_00014", "CREST timeout after 3600s")

        df = pd.read_csv(csv_path)
        row = df.iloc[0]

        assert row["status"] == MoleculeStatus.SKIP.value
        assert int(row["reruns"]) == 3

    def test_mark_success_preserves_reruns(self, csv_path: Path) -> None:
        """Ensure mark_success does not alter reruns counter."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_00002",
                    "smiles": "CC",
                    "nheavy": 2,
                    "status": MoleculeStatus.RUNNING.value,
                    "reruns": 2,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                }
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        manager.mark_success("mol_00002", {})

        df = pd.read_csv(csv_path)
        row = df.iloc[0]

        assert row["status"] == MoleculeStatus.OK.value
        assert int(row["reruns"]) == 2


class TestResetStuckRunning:
    """Tests for reset_stuck_running() startup recovery."""

    def test_reset_stuck_running_resets_to_pending(self, csv_path: Path) -> None:
        """RUNNING molecules are reset to PENDING with reruns incremented."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_001",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.RUNNING.value,
                    "reruns": 0,
                },
                {
                    "mol_id": "mol_002",
                    "smiles": "CC",
                    "nheavy": 2,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        count = manager.reset_stuck_running()

        assert count == 1

        df = pd.read_csv(csv_path)
        assert df.loc[0, "status"] == MoleculeStatus.PENDING.value
        assert int(df.loc[0, "reruns"]) == 1
        # mol_002 should be unchanged
        assert df.loc[1, "status"] == MoleculeStatus.PENDING.value
        assert int(df.loc[1, "reruns"]) == 0

    def test_reset_stuck_running_no_running(self, csv_path: Path) -> None:
        """No RUNNING molecules returns 0."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_001",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        count = manager.reset_stuck_running()

        assert count == 0
