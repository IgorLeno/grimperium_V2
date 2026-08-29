"""Integrated conformer-to-MOPAC scientific calculation workflow.

This is the production composition point between Semi-Imperium's conformer
preparation and the concrete/adapter-backed minimum workflow.  A caller hands
over one molecular request and one effective configuration; the exact finite
selection produced by the first stage is the only selected-conformer set the
AM1, PM3 and PM7 calculations may optimize.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from grimperium.crest_pm7.config import PM7Config
from semi_imperium.conformers import (
    ConformerPreparation,
    ConformerRequest,
    ConformerWorkflow,
    MoleculeTopology,
)
from semi_imperium.domain import EffectiveConfiguration
from semi_imperium.mopac import (
    SUPPORTED_HAMILTONIANS,
    CommandRunner,
    JsonWorkflowJournal,
    MinimumWorkflowResult,
    MopacExecutableBackend,
    MopacMinimumBackend,
    MopacMinimumWorkflow,
    WorkflowJournal,
)


@dataclass(frozen=True)
class SemiImperiumCalculationResult:
    """Traceable output of conformer preparation plus minimum verification."""

    conformers: ConformerPreparation
    minima: MinimumWorkflowResult

    def to_dict(self) -> dict[str, Any]:
        """Serialize selection provenance and every MOPAC attempt together."""
        return {
            "conformers": self.conformers.to_dict(),
            "minima": self.minima.to_dict(),
        }


class SemiImperiumCalculationWorkflow:
    """Prepare conformers, then execute independent verified MOPAC surfaces."""

    def __init__(
        self,
        *,
        conformer_workflow: ConformerWorkflow,
        mopac_backend: MopacMinimumBackend,
        journal: WorkflowJournal,
    ) -> None:
        self.conformer_workflow = conformer_workflow
        self.mopac_backend = mopac_backend
        self.journal = journal

    @classmethod
    def from_pm7_config(
        cls,
        *,
        conformer_workflow: ConformerWorkflow,
        config: PM7Config,
        calculation_id: str,
        journal_path: Path,
        work_dir: Path | None = None,
        charge: int = 0,
        multiplicity: int = 1,
        runner: CommandRunner | None = None,
    ) -> SemiImperiumCalculationWorkflow:
        """Compose the real executable backend and durable workflow journal.

        ``calculation_id`` is the safe, unique namespace of one persisted
        calculation record.  It keeps identically numbered AM1/PM3/PM7 attempts
        from different molecules or runs in separate artifact trees.
        """
        backend = MopacExecutableBackend.from_pm7_config(
            config,
            calculation_id=calculation_id,
            work_dir=work_dir,
            charge=charge,
            multiplicity=multiplicity,
            runner=runner,
        )
        return cls(
            conformer_workflow=conformer_workflow,
            mopac_backend=backend,
            journal=JsonWorkflowJournal(journal_path),
        )

    def run(
        self,
        request: ConformerRequest,
        configuration: EffectiveConfiguration,
        *,
        hamiltonians: tuple[str, ...] = SUPPORTED_HAMILTONIANS,
        topology: MoleculeTopology | None = None,
    ) -> SemiImperiumCalculationResult:
        """Run the configured finite selection and persist all MOPAC transitions."""
        prepared = self.conformer_workflow.prepare(
            request,
            search_settings=configuration.conformer_search,
            selection_settings=configuration.conformer_selection,
            topology=topology,
        )
        minima = MopacMinimumWorkflow(
            self.mopac_backend,
            verification=configuration.verification,
            journal=self.journal,
        ).run(prepared.selection, hamiltonians=hamiltonians)
        return SemiImperiumCalculationResult(conformers=prepared, minima=minima)


__all__ = ["SemiImperiumCalculationResult", "SemiImperiumCalculationWorkflow"]
