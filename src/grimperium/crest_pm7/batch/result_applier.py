"""Shared dual-write service for batch result application."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.batch.state_manager import BatchStateManager


@dataclass(frozen=True)
class ResultApplicationDecision:
    """Structured outcome of applying one molecule result."""

    mol_id: str
    success: bool
    final_status: str
    reruns: int
    error: str | None = None


class BatchResultApplier:
    """Apply result decisions to operational state and mirror them to CSV."""

    def __init__(
        self,
        *,
        state_manager: BatchStateManager,
        csv_manager: BatchCSVManager,
    ) -> None:
        self.state_manager = state_manager
        self.csv_manager = csv_manager

    def apply_success(
        self,
        mol_id: str,
        result_update: dict[str, Any] | None = None,
    ) -> ResultApplicationDecision:
        """Mark a molecule OK in operational state and mirror scientific data."""
        state_snapshot = self.state_manager.snapshot_row(mol_id)
        self.state_manager.mark_success(mol_id)
        try:
            self.csv_manager.mark_success(mol_id, result_update or {})
        except Exception:
            self.state_manager.restore_row(mol_id, state_snapshot)
            raise

        return ResultApplicationDecision(
            mol_id=mol_id,
            success=True,
            final_status=MoleculeStatus.OK.value,
            reruns=self.state_manager.get_reruns(mol_id),
        )

    def apply_failure(
        self,
        mol_id: str,
        error: str,
        *,
        result_update: dict[str, Any] | None = None,
        force_skip: bool = False,
    ) -> ResultApplicationDecision:
        """Record a failure in state and mirror the decided status to CSV."""
        state_snapshot = self.state_manager.snapshot_row(mol_id)
        if force_skip:
            final_status = MoleculeStatus.SKIP.value
            self.state_manager.update_molecule_status(mol_id, final_status)
        else:
            final_status = self.state_manager.record_rerun_or_skip(mol_id, error)
        reruns = self.state_manager.get_reruns(mol_id)

        try:
            self.csv_manager.apply_operational_status(
                mol_id,
                final_status,
                reruns=reruns,
                result_update=result_update,
            )
        except Exception:
            self.state_manager.restore_row(mol_id, state_snapshot)
            raise

        return ResultApplicationDecision(
            mol_id=mol_id,
            success=False,
            final_status=final_status,
            reruns=reruns,
            error=error,
        )
