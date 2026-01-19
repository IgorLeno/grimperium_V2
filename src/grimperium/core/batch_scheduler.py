"""
Batch scheduler: prioritize molecules for processing.

This module determines which molecules to process and in what order.
Priority is given to failed molecules (eligible for rerun) before pending ones.

**AJUSTE v2.2:** max_reruns passed explicitly (not from global)
- Improves testability
- Removes implicit dependency
- Clear function signature
"""

from __future__ import annotations

import logging

from grimperium.core.molecule import Molecule

logger = logging.getLogger(__name__)


class BatchScheduler:
    """
    Determine processing order for molecules.

    **AJUSTE v2.2:** max_reruns is explicit parameter

    Priority order:
    1. FAILED molecules (reruns < max_reruns) - retry these first
    2. PENDING molecules - process new ones after retries

    Molecules with status COMPLETE, SKIP, or RUNNING are not scheduled.

    Usage:
        >>> molecules = load_molecules()
        >>> scheduled = BatchScheduler.schedule(molecules, max_reruns=3)
        >>> for mol in scheduled:
        ...     process(mol)
    """

    @staticmethod
    def schedule(
        molecules: list[Molecule],
        max_reruns: int,
    ) -> list[Molecule]:
        """
        Schedule molecules for processing.

        Args:
            molecules: List of molecules
            max_reruns: Maximum rerun attempts (explicit parameter)

        Returns:
            List[Molecule] in processing order
            (preserves CSV order within priority groups)

        Priority:
        1. FAILED molecules (reruns < max_reruns)
        2. PENDING molecules
        3. Everything else: SKIPPED

        Complexity: O(n)

        Examples:
            >>> from grimperium.core.molecule import Molecule, MoleculeStatus
            >>> mol_pending = Molecule(...)
            >>> mol_pending.meta.status = MoleculeStatus.PENDING
            >>> mol_failed = Molecule(...)
            >>> mol_failed.meta.status = MoleculeStatus.RERUN
            >>> mol_failed.meta.reruns = 1
            >>> scheduled = BatchScheduler.schedule(
            ...     [mol_pending, mol_failed], max_reruns=3
            ... )
            >>> # Failed first, then pending
            >>> scheduled[0] is mol_failed
            True
        """
        # Priority 1: FAILED molecules (eligible for rerun)
        failed_rerunnable = [
            m for m in molecules if m.is_failed and m.can_rerun(max_reruns)
        ]

        # Priority 2: PENDING molecules
        pending = [m for m in molecules if m.is_pending]

        # Combine (preserves order within each group)
        scheduled = failed_rerunnable + pending

        logger.info(
            f"Scheduled {len(scheduled)} molecules "
            f"({len(failed_rerunnable)} failed, {len(pending)} pending) "
            f"with max_reruns={max_reruns}"
        )

        return scheduled

    @staticmethod
    def get_skip_reasons(
        molecules: list[Molecule],
        max_reruns: int,
    ) -> dict[str, int]:
        """
        Get reasons why molecules are skipped.

        Args:
            molecules: List of molecules
            max_reruns: Maximum rerun attempts

        Returns:
            {
                'max_reruns_exceeded': count,
                'complete': count,
                'running': count,
                'skipped': count,
            }

        Examples:
            >>> reasons = BatchScheduler.get_skip_reasons(molecules, max_reruns=3)
            >>> if reasons['max_reruns_exceeded'] > 0:
            ...     print(f"Warning: {reasons['max_reruns_exceeded']} molecules exhausted retries")
        """
        reasons = {
            "max_reruns_exceeded": 0,
            "complete": 0,
            "running": 0,
            "skipped": 0,
        }

        for mol in molecules:
            if mol.is_failed and not mol.can_rerun(max_reruns):
                reasons["max_reruns_exceeded"] += 1
            elif mol.is_complete:
                reasons["complete"] += 1
            elif mol.is_running:
                reasons["running"] += 1
            elif mol.is_skipped:
                reasons["skipped"] += 1

        return reasons

    @staticmethod
    def partition_by_status(
        molecules: list[Molecule],
    ) -> dict[str, list[Molecule]]:
        """
        Partition molecules by their current status.

        Args:
            molecules: List of molecules

        Returns:
            Dictionary with status as key and list of molecules as value.
            Keys are lowercase status values.

        Examples:
            >>> partitioned = BatchScheduler.partition_by_status(molecules)
            >>> print(f"Pending: {len(partitioned['pending'])}")
            >>> print(f"Complete: {len(partitioned['ok'])}")
        """
        partitions: dict[str, list[Molecule]] = {
            "pending": [],
            "selected": [],
            "running": [],
            "ok": [],
            "rerun": [],
            "skip": [],
        }

        for mol in molecules:
            key = mol.status.value.lower()
            if key in partitions:
                partitions[key].append(mol)

        return partitions

    @staticmethod
    def count_schedulable(
        molecules: list[Molecule],
        max_reruns: int,
    ) -> dict[str, int]:
        """
        Count how many molecules are schedulable.

        Args:
            molecules: List of molecules
            max_reruns: Maximum rerun attempts

        Returns:
            {
                'total': total molecules,
                'schedulable': molecules that will be scheduled,
                'failed_rerunnable': failed but eligible for rerun,
                'pending': new molecules,
                'complete': already finished,
                'skipped': permanently skipped,
            }
        """
        failed_rerunnable = pending = complete = skipped = running = 0
        for m in molecules:
            if m.is_failed and m.can_rerun(max_reruns):
                failed_rerunnable += 1
            if m.is_pending:
                pending += 1
            if m.is_complete:
                complete += 1
            if m.is_skipped:
                skipped += 1
            if m.is_running:
                running += 1

        return {
            "total": len(molecules),
            "schedulable": failed_rerunnable + pending,
            "failed_rerunnable": failed_rerunnable,
            "pending": pending,
            "complete": complete,
            "skipped": skipped,
            "running": running,
        }
