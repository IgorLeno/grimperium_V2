"""Tests for retry/error handling in BatchCSVManager.

Focus on ensuring reruns counter is persisted correctly.
"""

from pathlib import Path

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import (
    BatchFailurePolicy,
    BatchSortingStrategy,
    MoleculeStatus,
)


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
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_002",
                    "smiles": "CC",
                    "nheavy": 2,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
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
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        count = manager.reset_stuck_running()

        assert count == 0

    def test_reset_stuck_running_skips_at_max_reruns(self, csv_path: Path) -> None:
        """RUNNING molecule is skipped when recovery reaches max reruns."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_003",
                    "smiles": "CCC",
                    "nheavy": 3,
                    "status": MoleculeStatus.RUNNING.value,
                    "reruns": 2,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                }
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        count = manager.reset_stuck_running()

        assert count == 1

        df = pd.read_csv(csv_path)
        assert df.loc[0, "status"] == MoleculeStatus.SKIP.value
        assert int(df.loc[0, "reruns"]) == 3


class TestRecoverOrphanSelected:
    """Tests for recovering orphaned SELECTED molecules from interrupted batches."""

    def test_select_batch_recovers_orphan_selected(self, csv_path: Path) -> None:
        """Orphaned SELECTED molecule is recovered into new batch."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_001",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.SELECTED.value,
                    "batch_id": "batch_0001",
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_002",
                    "smiles": "CC",
                    "nheavy": 2,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_003",
                    "smiles": "CCC",
                    "nheavy": 3,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        batch = manager.select_batch(
            batch_id="batch_0002",
            batch_size=3,
            crest_timeout_minutes=60,
            mopac_timeout_minutes=30,
        )

        assert len(batch.molecules) == 3
        assert batch.molecules[0].mol_id == "mol_001"

        df = pd.read_csv(csv_path)
        orphan_row = df[df["mol_id"] == "mol_001"].iloc[0]
        assert orphan_row["batch_id"] == "batch_0002"

    def test_orphan_selected_does_not_increment_reruns(self, csv_path: Path) -> None:
        """Recovering orphan SELECTED must NOT increment reruns (never executed)."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_001",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.SELECTED.value,
                    "batch_id": "batch_0001",
                    "reruns": 1,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        batch = manager.select_batch(
            batch_id="batch_0002",
            batch_size=5,
            crest_timeout_minutes=60,
            mopac_timeout_minutes=30,
        )

        assert len(batch.molecules) == 1
        assert batch.molecules[0].reruns == 1

        df = pd.read_csv(csv_path)
        assert int(df.iloc[0]["reruns"]) == 1

    def test_orphan_selected_has_highest_priority(self, csv_path: Path) -> None:
        """Orphaned SELECTED molecules get priority over RERUN and PENDING."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_001",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_002",
                    "smiles": "CC",
                    "nheavy": 2,
                    "status": MoleculeStatus.RERUN.value,
                    "reruns": 1,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_003",
                    "smiles": "CCC",
                    "nheavy": 3,
                    "status": MoleculeStatus.SELECTED.value,
                    "batch_id": "batch_0001",
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        batch = manager.select_batch(
            batch_id="batch_0002",
            batch_size=3,
            crest_timeout_minutes=60,
            mopac_timeout_minutes=30,
            strategy=BatchSortingStrategy.RERUN_FIRST_THEN_EASY,
        )

        mol_ids = [m.mol_id for m in batch.molecules]
        assert mol_ids[0] == "mol_003"  # orphan first
        assert mol_ids[1] == "mol_002"  # rerun second
        assert mol_ids[2] == "mol_001"  # pending last

    def test_orphan_selected_fills_batch_before_pending(self, csv_path: Path) -> None:
        """When orphans >= batch_size, no PENDING molecules are included."""
        rows: list[dict[str, object]] = []
        for i in range(1, 4):
            rows.append(
                {
                    "mol_id": f"mol_orphan_{i:03d}",
                    "smiles": "C" * i,
                    "nheavy": i,
                    "status": MoleculeStatus.SELECTED.value,
                    "batch_id": "batch_0001",
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                }
            )
        for i in range(1, 6):
            rows.append(
                {
                    "mol_id": f"mol_pending_{i:03d}",
                    "smiles": "N" * i,
                    "nheavy": i,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                }
            )
        _write_csv(csv_path, rows)

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        batch = manager.select_batch(
            batch_id="batch_0002",
            batch_size=3,
            crest_timeout_minutes=60,
            mopac_timeout_minutes=30,
        )

        assert len(batch.molecules) == 3
        for mol in batch.molecules:
            assert mol.mol_id.startswith("mol_orphan_")

    def test_orphan_selected_partial_fill_then_pending(self, csv_path: Path) -> None:
        """1 orphan + remaining capacity filled from PENDING."""
        rows: list[dict[str, object]] = [
            {
                "mol_id": "mol_orphan_001",
                "smiles": "C",
                "nheavy": 1,
                "status": MoleculeStatus.SELECTED.value,
                "batch_id": "batch_0001",
                "reruns": 0,
                "crest_status": "NOT_ATTEMPTED",
                "mopac_status": "NOT_ATTEMPTED",
            },
        ]
        for i in range(1, 6):
            rows.append(
                {
                    "mol_id": f"mol_pending_{i:03d}",
                    "smiles": "N" * i,
                    "nheavy": i,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                }
            )
        _write_csv(csv_path, rows)

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        batch = manager.select_batch(
            batch_id="batch_0002",
            batch_size=3,
            crest_timeout_minutes=60,
            mopac_timeout_minutes=30,
        )

        assert len(batch.molecules) == 3
        assert batch.molecules[0].mol_id == "mol_orphan_001"
        assert batch.molecules[0].batch_order == 1
        assert batch.molecules[1].batch_order == 2
        assert batch.molecules[2].batch_order == 3

    def test_no_orphans_normal_behavior_unchanged(self, csv_path: Path) -> None:
        """Regression: no orphans means standard PENDING/RERUN selection."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_001",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_002",
                    "smiles": "CC",
                    "nheavy": 2,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_003",
                    "smiles": "CCC",
                    "nheavy": 3,
                    "status": MoleculeStatus.PENDING.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        batch = manager.select_batch(
            batch_id="batch_0002",
            batch_size=3,
            crest_timeout_minutes=60,
            mopac_timeout_minutes=30,
        )

        assert len(batch.molecules) == 3
        mol_ids = [m.mol_id for m in batch.molecules]
        assert mol_ids == ["mol_001", "mol_002", "mol_003"]

    def test_only_orphans_no_pending_available(self, csv_path: Path) -> None:
        """Batch from orphans alone when no PENDING/RERUN exist."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_001",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.SELECTED.value,
                    "batch_id": "batch_0001",
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_002",
                    "smiles": "CC",
                    "nheavy": 2,
                    "status": MoleculeStatus.SELECTED.value,
                    "batch_id": "batch_0001",
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
                {
                    "mol_id": "mol_ok",
                    "smiles": "CCC",
                    "nheavy": 3,
                    "status": MoleculeStatus.OK.value,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        batch = manager.select_batch(
            batch_id="batch_0002",
            batch_size=5,
            crest_timeout_minutes=60,
            mopac_timeout_minutes=30,
        )

        assert len(batch.molecules) == 2
        mol_ids = [m.mol_id for m in batch.molecules]
        assert "mol_001" in mol_ids
        assert "mol_002" in mol_ids

    def test_orphan_selected_updates_batch_metadata(self, csv_path: Path) -> None:
        """Orphan recovery updates batch_id, batch_order, policy, and timeouts."""
        _write_csv(
            csv_path,
            [
                {
                    "mol_id": "mol_001",
                    "smiles": "C",
                    "nheavy": 1,
                    "status": MoleculeStatus.SELECTED.value,
                    "batch_id": "batch_0001",
                    "batch_order": 5,
                    "batch_failure_policy": BatchFailurePolicy.ALL_OR_NOTHING.value,
                    "assigned_crest_timeout": 120,
                    "assigned_mopac_timeout": 60,
                    "reruns": 0,
                    "crest_status": "NOT_ATTEMPTED",
                    "mopac_status": "NOT_ATTEMPTED",
                },
            ],
        )

        manager = BatchCSVManager(csv_path)
        manager.load_csv()
        batch = manager.select_batch(
            batch_id="batch_0002",
            batch_size=3,
            crest_timeout_minutes=90,
            mopac_timeout_minutes=45,
            failure_policy=BatchFailurePolicy.PARTIAL_OK,
        )

        assert len(batch.molecules) == 1
        assert batch.molecules[0].batch_order == 1

        df = pd.read_csv(csv_path)
        row = df.iloc[0]
        assert row["batch_id"] == "batch_0002"
        assert int(row["batch_order"]) == 1
        assert row["batch_failure_policy"] == BatchFailurePolicy.PARTIAL_OK.value
        assert int(row["assigned_crest_timeout"]) == 90
        assert int(row["assigned_mopac_timeout"]) == 45
