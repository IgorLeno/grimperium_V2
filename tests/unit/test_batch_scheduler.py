"""
Tests for BatchScheduler.

Includes tests for Ajuste #2 (explicit max_reruns parameter).
"""

import pytest

from grimperium.core.batch_scheduler import BatchScheduler
from grimperium.core.molecule import (
    Molecule,
    MoleculeIdentity,
    MoleculeMeta,
    MoleculeProperties,
    MoleculeStatus,
)


def create_molecule(
    mol_id: str,
    status: MoleculeStatus = MoleculeStatus.PENDING,
    reruns: int = 0,
) -> Molecule:
    """Helper to create test molecules."""
    return Molecule(
        identity=MoleculeIdentity(mol_id=mol_id, smiles="C"),
        properties=MoleculeProperties(nheavy=1),
        meta=MoleculeMeta(status=status, reruns=reruns),
    )


class TestScheduleAjuste2:
    """Tests for schedule() - specifically for Ajuste #2."""

    def test_schedule_explicit_max_reruns_includes(self):
        """AJUSTE #2: max_reruns is explicit parameter - includes eligible."""
        mol_failed = create_molecule(
            "mol_001",
            status=MoleculeStatus.RERUN,
            reruns=2,
        )

        # With max_reruns=3, should be included (2 < 3)
        scheduled = BatchScheduler.schedule([mol_failed], max_reruns=3)
        assert mol_failed in scheduled

    def test_schedule_explicit_max_reruns_excludes(self):
        """AJUSTE #2: max_reruns is explicit - excludes ineligible."""
        mol_failed = create_molecule(
            "mol_001",
            status=MoleculeStatus.RERUN,
            reruns=2,
        )

        # With max_reruns=2, should be excluded (2 >= 2)
        scheduled = BatchScheduler.schedule([mol_failed], max_reruns=2)
        assert mol_failed not in scheduled

    def test_schedule_testable_no_global_dependency(self):
        """AJUSTE #2: No global settings dependency - fully testable."""
        mol = create_molecule("mol_001", status=MoleculeStatus.PENDING)

        # Can test without mocking any settings
        scheduled = BatchScheduler.schedule([mol], max_reruns=5)

        assert len(scheduled) == 1
        assert mol in scheduled


class TestSchedulePriority:
    """Tests for schedule() priority order."""

    def test_failed_before_pending(self):
        """Failed molecules are scheduled before pending."""
        mol_pending = create_molecule("mol_pending", status=MoleculeStatus.PENDING)
        mol_failed = create_molecule(
            "mol_failed",
            status=MoleculeStatus.RERUN,
            reruns=0,
        )

        # Pass in opposite order
        scheduled = BatchScheduler.schedule(
            [mol_pending, mol_failed], max_reruns=3
        )

        # Failed should come first
        assert scheduled[0] is mol_failed
        assert scheduled[1] is mol_pending

    def test_preserves_order_within_group(self):
        """Order is preserved within each priority group."""
        mol_p1 = create_molecule("pending_1", status=MoleculeStatus.PENDING)
        mol_p2 = create_molecule("pending_2", status=MoleculeStatus.PENDING)
        mol_p3 = create_molecule("pending_3", status=MoleculeStatus.PENDING)

        scheduled = BatchScheduler.schedule(
            [mol_p1, mol_p2, mol_p3], max_reruns=3
        )

        # Same order as input
        assert scheduled == [mol_p1, mol_p2, mol_p3]

    def test_complete_not_scheduled(self):
        """Complete molecules are not scheduled."""
        mol_complete = create_molecule("mol_001", status=MoleculeStatus.OK)

        scheduled = BatchScheduler.schedule([mol_complete], max_reruns=3)

        assert mol_complete not in scheduled
        assert len(scheduled) == 0

    def test_skipped_not_scheduled(self):
        """Skipped molecules are not scheduled."""
        mol_skipped = create_molecule("mol_001", status=MoleculeStatus.SKIP)

        scheduled = BatchScheduler.schedule([mol_skipped], max_reruns=3)

        assert mol_skipped not in scheduled

    def test_running_not_scheduled(self):
        """Running molecules are not scheduled."""
        mol_running = create_molecule("mol_001", status=MoleculeStatus.RUNNING)

        scheduled = BatchScheduler.schedule([mol_running], max_reruns=3)

        assert mol_running not in scheduled


class TestGetSkipReasons:
    """Tests for get_skip_reasons()."""

    def test_counts_max_reruns_exceeded(self):
        """Counts molecules that exceeded max reruns."""
        mol = create_molecule("mol_001", status=MoleculeStatus.RERUN, reruns=3)

        reasons = BatchScheduler.get_skip_reasons([mol], max_reruns=3)

        assert reasons["max_reruns_exceeded"] == 1

    def test_counts_complete(self):
        """Counts complete molecules."""
        mol = create_molecule("mol_001", status=MoleculeStatus.OK)

        reasons = BatchScheduler.get_skip_reasons([mol], max_reruns=3)

        assert reasons["complete"] == 1

    def test_counts_skipped(self):
        """Counts skipped molecules."""
        mol = create_molecule("mol_001", status=MoleculeStatus.SKIP)

        reasons = BatchScheduler.get_skip_reasons([mol], max_reruns=3)

        assert reasons["skipped"] == 1


class TestPartitionByStatus:
    """Tests for partition_by_status()."""

    def test_partitions_correctly(self):
        """Molecules are partitioned by status."""
        mol_pending = create_molecule("pending", status=MoleculeStatus.PENDING)
        mol_ok = create_molecule("ok", status=MoleculeStatus.OK)
        mol_rerun = create_molecule("rerun", status=MoleculeStatus.RERUN)

        partitions = BatchScheduler.partition_by_status(
            [mol_pending, mol_ok, mol_rerun]
        )

        assert len(partitions["pending"]) == 1
        assert len(partitions["ok"]) == 1
        assert len(partitions["rerun"]) == 1


class TestCountSchedulable:
    """Tests for count_schedulable()."""

    def test_counts_correctly(self):
        """Counts schedulable molecules correctly."""
        molecules = [
            create_molecule("p1", status=MoleculeStatus.PENDING),
            create_molecule("p2", status=MoleculeStatus.PENDING),
            create_molecule("f1", status=MoleculeStatus.RERUN, reruns=0),
            create_molecule("ok", status=MoleculeStatus.OK),
            create_molecule("skip", status=MoleculeStatus.SKIP),
        ]

        counts = BatchScheduler.count_schedulable(molecules, max_reruns=3)

        assert counts["total"] == 5
        assert counts["schedulable"] == 3  # 2 pending + 1 failed
        assert counts["pending"] == 2
        assert counts["failed_rerunnable"] == 1
        assert counts["complete"] == 1
        assert counts["skipped"] == 1
