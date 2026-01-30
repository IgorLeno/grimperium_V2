"""Tests for retry/error handling in BatchCSVManager.

Focus on ensuring retry counters and error messages are persisted correctly.
"""

from pathlib import Path

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus


@pytest.fixture
def csv_path(tmp_path: Path) -> Path:
    """Provide a temporary CSV path."""
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
    """Tests for retry/skip transitions and error logging."""

    def test_mark_rerun_first_attempt_appends_suffix(self, csv_path: Path) -> None:
        """Ensure mark_rerun appends attempt 1 and increments counters."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_00014",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.RUNNING.value,
                    "retry_count": 0,
                    "max_retries": 3,
                    "reruns": 0,
                    "last_error_message": "",
                    "error_message": "",
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
        assert int(row["retry_count"]) == 1
        assert int(row["reruns"]) == 1
        assert row["last_error_message"] == "CREST timeout after 3600s (attempt 1)"
        assert row["error_message"] == "CREST timeout after 3600s (attempt 1)"

    def test_mark_rerun_second_attempt_appends_suffix(self, csv_path: Path) -> None:
        """Ensure mark_rerun appends attempt 2 without skipping."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_00014",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.RERUN.value,
                    "retry_count": 1,
                    "max_retries": 3,
                    "reruns": 1,
                    "last_error_message": "CREST timeout after 3600s (attempt 1)",
                    "error_message": "CREST timeout after 3600s (attempt 1)",
                }
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        manager.mark_rerun("mol_00014", "CREST timeout after 3600s")

        df = pd.read_csv(csv_path)
        row = df.iloc[0]

        assert row["status"] == MoleculeStatus.RERUN.value
        assert int(row["retry_count"]) == 2
        assert int(row["reruns"]) == 2
        assert row["last_error_message"] == "CREST timeout after 3600s (attempt 2)"
        assert row["error_message"] == "CREST timeout after 3600s (attempt 2)"

    def test_mark_rerun_third_attempt_skips(self, csv_path: Path) -> None:
        """Ensure third attempt transitions to Skip with attempt suffix."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_00014",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.RERUN.value,
                    "retry_count": 2,
                    "max_retries": 3,
                    "reruns": 2,
                    "last_error_message": "CREST timeout after 3600s (attempt 2)",
                    "error_message": "CREST timeout after 3600s (attempt 2)",
                }
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        manager.mark_rerun("mol_00014", "CREST timeout after 3600s")

        df = pd.read_csv(csv_path)
        row = df.iloc[0]

        assert row["status"] == MoleculeStatus.SKIP.value
        assert int(row["retry_count"]) == 3
        assert int(row["reruns"]) == 3
        assert row["last_error_message"] == "CREST timeout after 3600s (attempt 3)"
        assert row["error_message"] == "CREST timeout after 3600s (attempt 3)"

    def test_mark_success_preserves_retry_counts(self, csv_path: Path) -> None:
        """Ensure mark_success keeps reruns aligned with retry_count."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_00002",
                    "smiles": "CC",
                    "nheavy": 2,
                    "status": MoleculeStatus.RERUN.value,
                    "retry_count": 2,
                    "max_retries": 3,
                    "reruns": 2,
                    "last_error_message": "prev error",
                    "error_message": "prev error",
                }
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        manager.mark_success("mol_00002", {"error_message": None})

        df = pd.read_csv(csv_path)
        row = df.iloc[0]

        assert row["status"] == MoleculeStatus.OK.value
        assert int(row["retry_count"]) == 2
        assert int(row["reruns"]) == 2
