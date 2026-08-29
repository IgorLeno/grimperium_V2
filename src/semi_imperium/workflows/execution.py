"""The production executor behind the Calculate area.

This is the one place where the Calculate area meets the real scientific
pipeline: it turns one pending calculation into a
:class:`~semi_imperium.calculation.SemiImperiumCalculationWorkflow` run —
conformer preparation, then an independent MOPAC optimization and
minimum verification for the single requested Hamiltonian — and projects
the result back onto the persisted lifecycle vocabulary.

Nothing is inferred when a program is missing. If a molecule asks for a
CREST search and no CREST execution backend is configured, the
calculation fails with that reason instead of quietly falling back to a
single embedded structure, because the two are different evidence.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from semi_imperium.domain import (
    CalculationResultData,
    CalculationState,
    VerificationOutcome,
    VerificationPolicy,
)
from semi_imperium.mopac.models import CandidateState, HamiltonianResult
from semi_imperium.settings import SemiImperiumSettings
from semi_imperium.workflows.calculation import ExecutionOutcome, ExecutionRequest

if TYPE_CHECKING:  # pragma: no cover - typing only
    from semi_imperium.conformers import ConformerEnsemble, ConformerRequest
    from semi_imperium.domain import ConformerSearchSettings

#: Subdirectory of the store that holds MOPAC inputs, outputs and journals.
ARTIFACTS_DIRNAME = "artifacts"


class ScientificCalculationExecutor:
    """Run one pending calculation through the real conformer/MOPAC pipeline."""

    def __init__(
        self,
        settings: SemiImperiumSettings,
        *,
        crest_runner: Any | None = None,
        command_runner: Any | None = None,
        confpass_backend: Any | None = None,
    ) -> None:
        self.settings = settings
        self.crest_runner = crest_runner
        self.command_runner = command_runner
        self.confpass_backend = confpass_backend

    def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
        """Compute one molecule under one Hamiltonian, or say why it could not."""
        configuration = request.configuration
        if configuration.conformer_search.enabled and self.crest_runner is None:
            return ExecutionOutcome(
                state=CalculationState.FAILED,
                verification=VerificationOutcome.FAILED,
                error_message=(
                    "This molecule asks for a CREST conformer search but no CREST "
                    "execution backend is configured; disable CREST for it in the "
                    "Calculate table or configure a CREST backend"
                ),
            )

        from semi_imperium.calculation import SemiImperiumCalculationWorkflow
        from semi_imperium.conformers import ConformerRequest, ConformerWorkflow
        from semi_imperium.conformers.crest import CrestConformerSearch
        from semi_imperium.conformers.initial_structure import RDKitInitialStructure

        identity = request.identity
        artifacts = (
            self.settings.runtime.store_root
            / ARTIFACTS_DIRNAME
            / request.calculation_id
        )
        conformer_workflow = ConformerWorkflow(
            search_backend=(
                CrestConformerSearch(runner=self.crest_runner)
                if self.crest_runner is not None
                else _UnavailableSearch()
            ),
            initial_structure_backend=RDKitInitialStructure(),
            confpass_backend=self.confpass_backend,
        )
        workflow = SemiImperiumCalculationWorkflow.from_pm7_config(
            conformer_workflow=conformer_workflow,
            config=self._pm7_config(),
            calculation_id=request.calculation_id,
            journal_path=artifacts / "workflow.json",
            work_dir=artifacts / "mopac",
            charge=identity.charge,
            multiplicity=identity.multiplicity,
            runner=self.command_runner,
        )
        result = workflow.run(
            ConformerRequest(
                molecule_id=identity.molecule_id,
                smiles=identity.canonical_smiles,
                charge=identity.charge,
                multiplicity=identity.multiplicity,
                run_id=request.run_id,
            ),
            configuration,
            hamiltonians=(request.hamiltonian,),
        )
        return outcome_from_hamiltonian_result(
            result.minima.for_hamiltonian(request.hamiltonian),
            policy=configuration.verification.policy,
            store_root=self.settings.runtime.store_root,
        )

    def _pm7_config(self) -> Any:
        """Project the runtime settings onto Grimperium's execution config."""
        from grimperium.crest_pm7.config import PM7Config

        runtime = self.settings.runtime
        return PM7Config(
            crest_executable=runtime.crest_executable,
            mopac_executable=runtime.mopac_executable,
            temp_dir=runtime.work_dir,
            crest_threads=runtime.crest_threads,
            crest_timeout=runtime.crest_timeout_seconds,
            mopac_timeout_base=runtime.mopac_timeout_seconds,
            mopac_timeout_margin=1.0,
            mopac_scf_threshold=self.settings.semiempirical.scf_convergence,
            mopac_precise_scf="PRECISE" in self.settings.semiempirical.keywords,
        )


def outcome_from_hamiltonian_result(
    result: HamiltonianResult,
    *,
    policy: VerificationPolicy,
    store_root: Path | None = None,
) -> ExecutionOutcome:
    """Project one independent Hamiltonian result onto the stored lifecycle.

    The mapping keeps the distinction the pipeline draws: a heat of
    formation that has not been confirmed as a minimum is reported as
    ``unverified`` and never promoted to ``verified``.
    """
    attempts = {attempt.attempt_id: attempt for attempt in result.attempts}
    chosen = attempts.get(result.verified_attempt_id or "") or attempts.get(
        result.provisional_lowest_attempt_id or ""
    )

    if result.state is CandidateState.OPTIMIZATION_FAILED or chosen is None:
        errors = [
            attempt.optimization_error
            for attempt in result.attempts
            if attempt.optimization_error
        ]
        return ExecutionOutcome(
            state=CalculationState.FAILED,
            verification=VerificationOutcome.FAILED,
            error_message=(
                "; ".join(errors)
                or f"{result.hamiltonian} produced no usable optimized geometry"
            ),
        )

    diagnostics = chosen.diagnostics
    data = CalculationResultData(
        energy_hof_kcal_mol=chosen.provisional_heat_of_formation_kcal_mol,
        conformer_index=chosen.source_conformer_index,
        conformers_evaluated=len(result.selected_conformer_indices),
        lowest_frequency_cm1=(
            None if diagnostics is None else diagnostics.lowest_frequency_cm1
        ),
        artifact_paths=_relative_artifacts(
            (chosen.optimization_output_path, chosen.verification_output_path),
            store_root,
        ),
    )

    if result.state is CandidateState.MINIMUM_VERIFIED:
        return ExecutionOutcome(
            state=CalculationState.VERIFIED,
            verification=VerificationOutcome.MINIMUM_CONFIRMED,
            result=data,
        )
    if result.state is CandidateState.SADDLE_DETECTED:
        return ExecutionOutcome(
            state=CalculationState.SADDLE,
            verification=VerificationOutcome.SADDLE_POINT,
            result=data,
        )
    if result.state is CandidateState.VERIFICATION_FAILED:
        return ExecutionOutcome(
            state=CalculationState.UNVERIFIED,
            verification=VerificationOutcome.FAILED,
            result=data,
        )
    return ExecutionOutcome(
        state=CalculationState.UNVERIFIED,
        verification=(
            VerificationOutcome.NOT_REQUESTED
            if policy is VerificationPolicy.NONE
            else VerificationOutcome.INCONCLUSIVE
        ),
        result=data,
    )


def _relative_artifacts(
    paths: tuple[str | None, ...], store_root: Path | None
) -> tuple[str, ...]:
    """Keep only artifact paths that can be expressed inside the store."""
    if store_root is None:
        return ()
    relative: list[str] = []
    for path in paths:
        if not path:
            continue
        try:
            relative.append(str(Path(path).relative_to(store_root)))
        except ValueError:
            # Artifacts written outside the store are execution detail; the
            # record refuses absolute paths rather than storing a broken one.
            continue
    return tuple(relative)


class _UnavailableSearch:
    """Search backend used when CREST is off; calling it is a programming error."""

    def search(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        """Refuse loudly: the workflow must have taken the initial-3D route."""
        from semi_imperium.conformers import ConformerBackendError

        del settings
        raise ConformerBackendError(
            "No CREST execution backend is configured, so "
            f"{request.molecule_id!r} cannot run a conformer search",
            code="crest_unavailable",
        )


__all__ = [
    "ARTIFACTS_DIRNAME",
    "ScientificCalculationExecutor",
    "outcome_from_hamiltonian_result",
]
