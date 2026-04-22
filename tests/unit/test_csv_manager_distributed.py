"""Tests for distributed processing methods in BatchCSVManager.

Covers:
- distribute_molecules()
- mark_worker_offline()
- reassign_offline_molecules()
- reset_stuck_assigned()
"""

from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus, WorkerStatus

# 3 PENDING (mol_001, mol_002, mol_005), 1 RERUN (mol_003), 1 OK (mol_004)
# Available pool: mol_001, mol_002, mol_003, mol_005
# RERUN_FIRST order: mol_003, mol_001, mol_002, mol_005
_MINIMAL_CSV = """\
mol_id,smiles,nheavy,status
mol_001,C,1,Pending
mol_002,CC,2,Pending
mol_003,CCC,3,Rerun
mol_004,CCCC,4,OK
mol_005,CCCCC,5,Pending
"""


@pytest.fixture
def csv_path(tmp_path: Path) -> Path:
    return tmp_path / "distributed_test.csv"


@pytest.fixture
def manager(csv_path: Path) -> BatchCSVManager:
    csv_path.write_text(_MINIMAL_CSV)
    m = BatchCSVManager(csv_path)
    m.load_csv()
    return m


@pytest.fixture
def manager_max_reruns_1(csv_path: Path) -> BatchCSVManager:
    csv_path.write_text(_MINIMAL_CSV)
    m = BatchCSVManager(csv_path, max_reruns=1)
    m.load_csv()
    return m


# ── distribute_molecules ──────────────────────────────────────────────────────


class TestDistributeMolecules:
    def test_basic_allocation_returns_correct_mol_ids(
        self, manager: BatchCSVManager
    ) -> None:
        result = manager.distribute_molecules({"worker_A": 2}, 30, 60)
        assert len(result["worker_A"]) == 2

    def test_rerun_first_priority(self, manager: BatchCSVManager) -> None:
        result = manager.distribute_molecules({"worker_A": 2}, 30, 60)
        # mol_003 (RERUN) must be first
        assert result["worker_A"][0] == "mol_003"

    def test_assigned_status_written_to_csv(self, manager: BatchCSVManager) -> None:
        result = manager.distribute_molecules({"worker_A": 2}, 30, 60)
        df = manager._ensure_loaded()
        for mol_id in result["worker_A"]:
            row = df[df["mol_id"] == mol_id].iloc[0]
            assert row["status"] == MoleculeStatus.ASSIGNED.value
            assert row["assigned_worker"] == "worker_A"
            assert row["worker_status"] == WorkerStatus.ONLINE.value
            assert pd.notna(row["assigned_at"])

    def test_ok_molecules_never_assigned(self, manager: BatchCSVManager) -> None:
        result = manager.distribute_molecules({"worker_A": 10}, 30, 60)
        assert "mol_004" not in result["worker_A"]

    def test_multi_worker_no_overlap(self, manager: BatchCSVManager) -> None:
        result = manager.distribute_molecules({"worker_A": 2, "worker_B": 2}, 30, 60)
        assert not set(result["worker_A"]) & set(result["worker_B"])

    def test_pool_exhaustion_returns_available(self, manager: BatchCSVManager) -> None:
        # 4 available (3 PENDING + 1 RERUN), request 10
        result = manager.distribute_molecules({"worker_A": 10}, 30, 60)
        assert len(result["worker_A"]) == 4

    def test_second_worker_gets_remaining_pool(self, manager: BatchCSVManager) -> None:
        result = manager.distribute_molecules({"worker_A": 3, "worker_B": 10}, 30, 60)
        # 4 total available: worker_A gets 3, worker_B gets 1
        assert len(result["worker_A"]) == 3
        assert len(result["worker_B"]) == 1

    def test_existing_assigned_reused_in_allocation(
        self, manager: BatchCSVManager
    ) -> None:
        r1 = manager.distribute_molecules({"worker_A": 2}, 30, 60)
        # Second call with same count — should return same molecules, no new ones
        r2 = manager.distribute_molecules({"worker_A": 2}, 30, 60)
        assert set(r2["worker_A"]) == set(r1["worker_A"])

    def test_existing_assigned_refreshes_timeouts(
        self, manager: BatchCSVManager
    ) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        manager.distribute_molecules({"worker_A": 1}, 45, 90)
        df = manager._ensure_loaded()
        assigned = df[df["assigned_worker"] == "worker_A"]
        assert not assigned.empty
        assert int(assigned.iloc[0]["assigned_crest_timeout"]) == 45
        assert int(assigned.iloc[0]["assigned_mopac_timeout"]) == 90

    def test_timeouts_stored_in_csv(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 45, 90)
        df = manager._ensure_loaded()
        assigned = df[df["assigned_worker"] == "worker_A"]
        assert int(assigned.iloc[0]["assigned_crest_timeout"]) == 45
        assert int(assigned.iloc[0]["assigned_mopac_timeout"]) == 90

    def test_zero_allocation_returns_empty_list(self, manager: BatchCSVManager) -> None:
        result = manager.distribute_molecules({"worker_A": 0}, 30, 60)
        assert result["worker_A"] == []

    def test_empty_allocations_returns_empty_dict(
        self, manager: BatchCSVManager
    ) -> None:
        assert manager.distribute_molecules({}, 30, 60) == {}

    def test_negative_allocation_raises(self, manager: BatchCSVManager) -> None:
        with pytest.raises(ValueError, match="must be >= 0"):
            manager.distribute_molecules({"worker_A": -1}, 30, 60)

    def test_csv_persisted_after_distribution(
        self, manager: BatchCSVManager, csv_path: Path
    ) -> None:
        manager.distribute_molecules({"worker_A": 2}, 30, 60)
        # Reload from disk to verify persistence
        fresh = BatchCSVManager(csv_path)
        fresh.load_csv()
        assert fresh.get_status_counts().get(MoleculeStatus.ASSIGNED.value, 0) == 2


# ── mark_worker_offline ───────────────────────────────────────────────────────


class TestMarkWorkerOffline:
    def test_marks_assigned_molecules_offline(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 2}, 30, 60)
        count = manager.mark_worker_offline("worker_A")
        assert count == 2
        df = manager._ensure_loaded()
        offline = df[df["assigned_worker"] == "worker_A"]
        assert all(offline["worker_status"] == WorkerStatus.OFFLINE.value)

    def test_molecule_status_unchanged(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 2}, 30, 60)
        manager.mark_worker_offline("worker_A")
        df = manager._ensure_loaded()
        worker_a = df[df["assigned_worker"] == "worker_A"]
        # Status must still be ASSIGNED — only worker_status changed
        assert all(worker_a["status"] == MoleculeStatus.ASSIGNED.value)

    def test_only_affects_target_worker(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1, "worker_B": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        df = manager._ensure_loaded()
        worker_b = df[df["assigned_worker"] == "worker_B"]
        assert all(worker_b["worker_status"] == WorkerStatus.ONLINE.value)

    def test_returns_zero_for_unknown_worker(self, manager: BatchCSVManager) -> None:
        assert manager.mark_worker_offline("ghost_worker") == 0

    def test_includes_running_molecules(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        # Force RUNNING status
        df = manager._ensure_loaded()
        idx = df[df["assigned_worker"] == "worker_A"].index[0]
        df.at[idx, "status"] = MoleculeStatus.RUNNING.value
        manager.save_csv()

        count = manager.mark_worker_offline("worker_A")
        assert count == 1


# ── reassign_offline_molecules ────────────────────────────────────────────────


class TestReassignOfflineMolecules:
    def test_assigned_offline_back_to_pending_no_reruns(
        self, manager: BatchCSVManager
    ) -> None:
        manager.distribute_molecules({"worker_A": 2}, 30, 60)
        manager.mark_worker_offline("worker_A")
        count = manager.reassign_offline_molecules()
        assert count == 2
        df = manager._ensure_loaded()
        # Both molecules should now be PENDING with cleared worker fields
        for _, row in df[df["mol_id"].isin(["mol_003", "mol_001"])].iterrows():
            assert row["status"] == MoleculeStatus.PENDING.value
            assert pd.isna(row["assigned_worker"]) or row["assigned_worker"] is None
            assert row["worker_status"] == WorkerStatus.UNASSIGNED.value
            assert pd.isna(row["assigned_at"]) or row["assigned_at"] is None

    def test_assigned_offline_no_reruns_increment(
        self, manager: BatchCSVManager
    ) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        df = manager._ensure_loaded()
        idx = df[df["assigned_worker"] == "worker_A"].index[0]
        mol_id: str = str(df.at[idx, "mol_id"])

        manager.mark_worker_offline("worker_A")
        manager.reassign_offline_molecules()

        df2 = manager._ensure_loaded()
        row = df2[df2["mol_id"] == mol_id].iloc[0]
        # reruns column may not exist in minimal CSV — either absent or 0 is correct
        reruns_val = row["reruns"] if "reruns" in df2.columns else None
        assert manager._safe_int(reruns_val, 0) == 0

    def test_running_offline_increments_reruns(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        df = manager._ensure_loaded()
        idx = df[df["assigned_worker"] == "worker_A"].index[0]
        mol_id: str = str(df.at[idx, "mol_id"])
        df.at[idx, "status"] = MoleculeStatus.RUNNING.value
        df.at[idx, "reruns"] = 0
        manager.save_csv()

        manager.mark_worker_offline("worker_A")
        manager.reassign_offline_molecules()

        df2 = manager._ensure_loaded()
        row = df2[df2["mol_id"] == mol_id].iloc[0]
        assert row["status"] == MoleculeStatus.PENDING.value
        assert manager._safe_int(row["reruns"], 0) == 1

    def test_running_offline_max_reruns_gives_skip(
        self, manager_max_reruns_1: BatchCSVManager
    ) -> None:
        manager_max_reruns_1.distribute_molecules({"worker_A": 1}, 30, 60)
        df = manager_max_reruns_1._ensure_loaded()
        idx = df[df["assigned_worker"] == "worker_A"].index[0]
        mol_id: str = str(df.at[idx, "mol_id"])
        df.at[idx, "status"] = MoleculeStatus.RUNNING.value
        df.at[idx, "reruns"] = 0
        manager_max_reruns_1.save_csv()

        manager_max_reruns_1.mark_worker_offline("worker_A")
        manager_max_reruns_1.reassign_offline_molecules()

        df2 = manager_max_reruns_1._ensure_loaded()
        assert (
            df2[df2["mol_id"] == mol_id].iloc[0]["status"] == MoleculeStatus.SKIP.value
        )

    def test_filter_by_worker_id_leaves_other_offline(
        self, manager: BatchCSVManager
    ) -> None:
        manager.distribute_molecules({"worker_A": 1, "worker_B": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        manager.mark_worker_offline("worker_B")

        count = manager.reassign_offline_molecules(worker_id="worker_A")
        assert count == 1

        df = manager._ensure_loaded()
        worker_b = df[df["assigned_worker"] == "worker_B"]
        # worker_B molecules are still ASSIGNED+offline — not touched
        assert not worker_b.empty
        assert all(worker_b["worker_status"] == WorkerStatus.OFFLINE.value)
        assert all(worker_b["status"] == MoleculeStatus.ASSIGNED.value)

    def test_returns_zero_when_nothing_offline(self, manager: BatchCSVManager) -> None:
        assert manager.reassign_offline_molecules() == 0

    def test_reassigned_molecules_re_enter_pool(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 2}, 30, 60)
        manager.mark_worker_offline("worker_A")
        manager.reassign_offline_molecules()

        # Now distribute again — molecules should be available
        result = manager.distribute_molecules({"worker_B": 4}, 30, 60)
        assert len(result["worker_B"]) == 4  # All 4 available again


# ── reset_stuck_assigned ──────────────────────────────────────────────────────


class TestResetStuckAssigned:
    def _set_assigned_at(
        self,
        manager: BatchCSVManager,
        worker_id: str,
        assigned_at: str | None,
    ) -> None:
        """Helper to override assigned_at in the DataFrame directly."""
        df = manager._ensure_loaded()
        mask = df["assigned_worker"] == worker_id
        for idx in df[mask].index:
            df.at[idx, "assigned_at"] = assigned_at
        manager.save_csv()

    def test_no_offline_molecules_returns_zeros(self, manager: BatchCSVManager) -> None:
        assert manager.reset_stuck_assigned() == (0, 0)

    def test_old_assigned_auto_resets_to_pending(
        self, manager: BatchCSVManager
    ) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        old_ts = (datetime.now(timezone.utc) - timedelta(hours=3)).isoformat()
        self._set_assigned_at(manager, "worker_A", old_ts)

        auto_reset, warned = manager.reset_stuck_assigned(threshold_hours=2.0)
        assert auto_reset == 1
        assert warned == 0

        df = manager._ensure_loaded()
        # The previously assigned molecule should now be PENDING
        pending = df[df["status"] == MoleculeStatus.PENDING.value]
        assert not pending.empty

    def test_auto_reset_clears_worker_fields(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        df = manager._ensure_loaded()
        idx = df[df["assigned_worker"] == "worker_A"].index[0]
        mol_id: str = str(df.at[idx, "mol_id"])
        old_ts = (datetime.now(timezone.utc) - timedelta(hours=3)).isoformat()
        df.at[idx, "assigned_at"] = old_ts
        manager.save_csv()

        manager.reset_stuck_assigned(threshold_hours=2.0)

        df2 = manager._ensure_loaded()
        row = df2[df2["mol_id"] == mol_id].iloc[0]
        assert pd.isna(row["assigned_worker"]) or row["assigned_worker"] is None
        assert row["worker_status"] == WorkerStatus.UNASSIGNED.value
        assert pd.isna(row["assigned_at"]) or row["assigned_at"] is None

    def test_auto_reset_no_reruns_increment(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        df = manager._ensure_loaded()
        idx = df[df["assigned_worker"] == "worker_A"].index[0]
        mol_id: str = str(df.at[idx, "mol_id"])
        old_ts = (datetime.now(timezone.utc) - timedelta(hours=3)).isoformat()
        df.at[idx, "assigned_at"] = old_ts
        manager.save_csv()

        manager.reset_stuck_assigned(threshold_hours=2.0)

        df2 = manager._ensure_loaded()
        row = df2[df2["mol_id"] == mol_id].iloc[0]
        # reruns column may not exist in minimal CSV — either absent or 0 is correct
        reruns_val = row["reruns"] if "reruns" in df2.columns else None
        assert manager._safe_int(reruns_val, 0) == 0

    def test_recent_assigned_warns_only(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        # assigned_at is very recent (just set by distribute_molecules)

        auto_reset, warned = manager.reset_stuck_assigned(threshold_hours=2.0)
        assert auto_reset == 0
        assert warned == 1

    def test_recent_molecule_status_unchanged(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")

        manager.reset_stuck_assigned(threshold_hours=2.0)

        df = manager._ensure_loaded()
        worker_a = df[df["assigned_worker"] == "worker_A"]
        # Should still be ASSIGNED — only warned, not reset
        assert not worker_a.empty
        assert all(worker_a["status"] == MoleculeStatus.ASSIGNED.value)

    def test_missing_timestamp_auto_resets(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        self._set_assigned_at(manager, "worker_A", None)

        auto_reset, warned = manager.reset_stuck_assigned(threshold_hours=2.0)
        assert auto_reset == 1
        assert warned == 0

    def test_mixed_old_and_recent(self, manager: BatchCSVManager) -> None:
        # Assign 2 molecules to different workers
        manager.distribute_molecules({"worker_A": 1, "worker_B": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        manager.mark_worker_offline("worker_B")

        # worker_A: old timestamp → auto-reset
        old_ts = (datetime.now(timezone.utc) - timedelta(hours=5)).isoformat()
        self._set_assigned_at(manager, "worker_A", old_ts)
        # worker_B: recent timestamp (just now) → warn only

        auto_reset, warned = manager.reset_stuck_assigned(threshold_hours=2.0)
        assert auto_reset == 1
        assert warned == 1

    def test_configurable_threshold(self, manager: BatchCSVManager) -> None:
        manager.distribute_molecules({"worker_A": 1}, 30, 60)
        manager.mark_worker_offline("worker_A")
        # Set assigned_at to 30 minutes ago
        recent_ts = (datetime.now(timezone.utc) - timedelta(minutes=30)).isoformat()
        self._set_assigned_at(manager, "worker_A", recent_ts)

        # With 1-hour threshold: 30 min is recent → warn
        auto_reset, warned = manager.reset_stuck_assigned(threshold_hours=1.0)
        assert auto_reset == 0
        assert warned == 1

        # With 15-minute threshold: 30 min is old → reset
        manager.distribute_molecules({"worker_A": 0}, 30, 60)  # No-op re-call
        # Re-set the timestamp (reassign_offline was not called, so mol still ASSIGNED+offline)
        auto_reset2, warned2 = manager.reset_stuck_assigned(threshold_hours=0.25)
        assert auto_reset2 == 1
        assert warned2 == 0
