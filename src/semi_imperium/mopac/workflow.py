"""Independent AM1/PM3/PM7 optimization and bounded minimum recovery."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any, Protocol

from semi_imperium.conformers import Conformer, ConformerGeometry, SelectionResult
from semi_imperium.domain import VerificationPolicy, VerificationSettings
from semi_imperium.mopac.force_parser import classify_force_output
from semi_imperium.mopac.models import (
    SUPPORTED_HAMILTONIANS,
    AttemptOrigin,
    CandidateAttempt,
    CandidateState,
    DisplacementLineage,
    ForceRun,
    HamiltonianResult,
    MinimumWorkflowResult,
    OptimizationRun,
    SelectionLineage,
)


class MopacMinimumBackend(Protocol):
    """Injected MOPAC execution boundary used by the workflow."""

    def optimize(
        self,
        *,
        hamiltonian: str,
        geometry: ConformerGeometry,
        source_conformer_index: int,
        attempt_id: str,
        displacement: DisplacementLineage | None,
    ) -> OptimizationRun:
        """Optimize exactly one geometry under exactly one Hamiltonian."""
        ...

    def verify_force(
        self,
        *,
        hamiltonian: str,
        optimization: OptimizationRun,
        attempt_id: str,
    ) -> ForceRun:
        """Run MOPAC FORCE on the optimized geometry without changing it."""
        ...


class WorkflowJournal(Protocol):
    """Durable sink for attempt transitions and terminal outcomes."""

    def record_attempt(self, attempt: CandidateAttempt) -> None:
        """Persist one attempt snapshot before execution continues."""
        ...

    def record_hamiltonian_result(self, result: HamiltonianResult) -> None:
        """Persist one Hamiltonian's terminal result."""
        ...

    def record_workflow_result(self, result: MinimumWorkflowResult) -> None:
        """Persist the complete terminal result."""
        ...


class NullWorkflowJournal:
    """No-op journal for callers that persist returned models elsewhere."""

    def record_attempt(self, attempt: CandidateAttempt) -> None:
        """Accept an attempt without external storage."""

    def record_hamiltonian_result(self, result: HamiltonianResult) -> None:
        """Accept a terminal Hamiltonian result without external storage."""

    def record_workflow_result(self, result: MinimumWorkflowResult) -> None:
        """Accept a workflow result without external storage."""


class JsonWorkflowJournal:
    """Atomic JSON event journal preserving every attempt state transition."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self._payload: dict[str, Any] = {
            "schema_version": 1,
            "attempt_events": [],
            "hamiltonian_results": {},
            "workflow_result": None,
        }

    def record_attempt(self, attempt: CandidateAttempt) -> None:
        """Append an immutable snapshot, including OPTIMIZED_UNVERIFIED."""
        events = self._payload["attempt_events"]
        if not isinstance(events, list):  # pragma: no cover - internal invariant
            raise TypeError("JsonWorkflowJournal.attempt_events is not a list")
        events.append(attempt.to_dict())
        self._flush()

    def record_hamiltonian_result(self, result: HamiltonianResult) -> None:
        """Store an independent terminal result keyed by Hamiltonian."""
        results = self._payload["hamiltonian_results"]
        if not isinstance(results, dict):  # pragma: no cover - internal invariant
            raise TypeError("JsonWorkflowJournal.hamiltonian_results is not a dict")
        results[result.hamiltonian] = result.to_dict()
        self._flush()

    def record_workflow_result(self, result: MinimumWorkflowResult) -> None:
        """Store the complete workflow after all requested methods terminate."""
        self._payload["workflow_result"] = result.to_dict()
        self._flush()

    def _flush(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_name(f".{self.path.name}.tmp")
        temporary.write_text(
            json.dumps(self._payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(self.path)


class MopacMinimumWorkflow:
    """Optimize, verify and recover a finite candidate set per Hamiltonian.

    Every requested Hamiltonian starts from the same selected conformers but
    receives its own optimization calls, energy ordering, FORCE analyses and
    final selection.  There is no shared "best conformer" after CREST.
    """

    def __init__(
        self,
        backend: MopacMinimumBackend,
        *,
        verification: VerificationSettings | None = None,
        journal: WorkflowJournal | None = None,
    ) -> None:
        self.backend = backend
        self.verification = verification or VerificationSettings(
            policy=VerificationPolicy.REQUIRE_MINIMUM
        )
        if not self.verification.requires_minimum:
            raise ValueError(
                "MopacMinimumWorkflow requires the 'require_minimum' "
                "verification policy"
            )
        self.journal = journal or NullWorkflowJournal()

    def run(
        self,
        selection: SelectionResult,
        *,
        hamiltonians: tuple[str, ...] = SUPPORTED_HAMILTONIANS,
    ) -> MinimumWorkflowResult:
        """Run each requested Hamiltonian independently and in request order."""
        normalized = tuple(method.strip().upper() for method in hamiltonians)
        if not normalized:
            raise ValueError("At least one Hamiltonian must be requested")
        if len(set(normalized)) != len(normalized):
            raise ValueError("Hamiltonian requests must be unique")
        unsupported = tuple(
            method for method in normalized if method not in SUPPORTED_HAMILTONIANS
        )
        if unsupported:
            raise ValueError(
                f"Unsupported Hamiltonian(s): {', '.join(unsupported)}; "
                f"expected {', '.join(SUPPORTED_HAMILTONIANS)}"
            )

        outcomes = tuple(
            self._run_hamiltonian(method, selection.selected) for method in normalized
        )
        result = MinimumWorkflowResult(
            selection_lineage=SelectionLineage.from_selection(selection),
            hamiltonian_results=outcomes,
        )
        self.journal.record_workflow_result(result)
        return result

    def _run_hamiltonian(
        self, hamiltonian: str, conformers: tuple[Conformer, ...]
    ) -> HamiltonianResult:
        attempts: list[CandidateAttempt] = []
        optimizations: dict[str, OptimizationRun] = {}
        force_runs: dict[str, ForceRun] = {}

        for conformer in conformers:
            sequence = len(attempts)
            attempt_id = _attempt_id(hamiltonian, sequence)
            optimization = self._optimize(
                hamiltonian=hamiltonian,
                geometry=conformer.geometry,
                source_conformer_index=conformer.index,
                attempt_id=attempt_id,
                displacement=None,
            )
            attempt = _attempt_from_optimization(
                attempt_id=attempt_id,
                hamiltonian=hamiltonian,
                sequence=sequence,
                origin=AttemptOrigin.SELECTED_CONFORMER,
                source_conformer_index=conformer.index,
                optimization=optimization,
                displacement=None,
            )
            attempts.append(attempt)
            optimizations[attempt_id] = optimization
            self.journal.record_attempt(attempt)

        provisional = _lowest_optimized(attempts, selected_only=True)
        for position in _optimized_positions_by_energy(attempts):
            verified, force_run = self._verify(
                attempts[position], optimizations[attempts[position].attempt_id]
            )
            attempts[position] = verified
            force_runs[verified.attempt_id] = force_run
            self.journal.record_attempt(verified)

        selected_verified = _lowest_verified(attempts)
        if selected_verified is not None:
            return self._finish(
                hamiltonian,
                conformers,
                attempts,
                provisional,
                selected_verified,
                recovery_used=0,
            )

        recovery_used = 0
        while recovery_used < self.verification.max_displacement_reoptimizations:
            source = _recovery_source(attempts, force_runs)
            if source is None:
                break
            parent, force_run, mode = source
            direction = 1 if recovery_used % 2 == 0 else -1
            # Budgets above two must not silently repeat the same two Cartesian
            # geometries.  Each later +/- pair uses a smaller displacement,
            # still bounded by the configured maximum step.
            amplitude = self.verification.displacement_step_angstrom / (
                1 + recovery_used // 2
            )
            displacement = DisplacementLineage(
                parent_attempt_id=parent.attempt_id,
                normal_mode=mode,
                direction=direction,
                amplitude_angstrom=amplitude,
            )
            vector = force_run.normal_modes[mode]
            if parent.optimized_geometry is None:  # pragma: no cover - model invariant
                break
            geometry = displace_geometry(
                parent.optimized_geometry,
                vector,
                direction=direction,
                amplitude_angstrom=amplitude,
            )
            sequence = len(attempts)
            attempt_id = _attempt_id(hamiltonian, sequence)
            optimization = self._optimize(
                hamiltonian=hamiltonian,
                geometry=geometry,
                source_conformer_index=parent.source_conformer_index,
                attempt_id=attempt_id,
                displacement=displacement,
            )
            attempt = _attempt_from_optimization(
                attempt_id=attempt_id,
                hamiltonian=hamiltonian,
                sequence=sequence,
                origin=AttemptOrigin.DISPLACEMENT_REOPTIMIZATION,
                source_conformer_index=parent.source_conformer_index,
                optimization=optimization,
                displacement=displacement,
            )
            recovery_used += 1
            attempts.append(attempt)
            optimizations[attempt_id] = optimization
            self.journal.record_attempt(attempt)
            if attempt.state is CandidateState.OPTIMIZATION_FAILED:
                continue
            verified, new_force_run = self._verify(attempt, optimization)
            attempts[-1] = verified
            force_runs[attempt_id] = new_force_run
            self.journal.record_attempt(verified)
            if verified.is_verified_minimum:
                return self._finish(
                    hamiltonian,
                    conformers,
                    attempts,
                    provisional,
                    verified,
                    recovery_used=recovery_used,
                )

        return self._finish(
            hamiltonian,
            conformers,
            attempts,
            provisional,
            None,
            recovery_used=recovery_used,
        )

    def _optimize(
        self,
        *,
        hamiltonian: str,
        geometry: ConformerGeometry,
        source_conformer_index: int,
        attempt_id: str,
        displacement: DisplacementLineage | None,
    ) -> OptimizationRun:
        try:
            return self.backend.optimize(
                hamiltonian=hamiltonian,
                geometry=geometry,
                source_conformer_index=source_conformer_index,
                attempt_id=attempt_id,
                displacement=displacement,
            )
        except Exception as exc:  # backend boundary: preserve the failed candidate
            return OptimizationRun(
                converged=False,
                error_message=f"Optimization backend failed: {exc}",
            )

    def _verify(
        self, attempt: CandidateAttempt, optimization: OptimizationRun
    ) -> tuple[CandidateAttempt, ForceRun]:
        try:
            force_run = self.backend.verify_force(
                hamiltonian=attempt.hamiltonian,
                optimization=optimization,
                attempt_id=attempt.attempt_id,
            )
        except Exception as exc:  # backend boundary: never promote on an exception
            force_run = ForceRun(
                output="",
                execution_error=f"FORCE backend failed: {exc}",
            )
        atom_count = (
            None
            if attempt.optimized_geometry is None
            else attempt.optimized_geometry.atom_count
        )
        classification = classify_force_output(
            force_run.output,
            self.verification,
            atom_count=atom_count,
            execution_error=force_run.execution_error,
        )
        return (
            replace(
                attempt,
                state=classification.state,
                state_history=(*attempt.state_history, classification.state),
                verification_output_path=force_run.output_path,
                diagnostics=classification.diagnostics,
            ),
            force_run,
        )

    def _finish(
        self,
        hamiltonian: str,
        conformers: tuple[Conformer, ...],
        attempts: list[CandidateAttempt],
        provisional: CandidateAttempt | None,
        verified: CandidateAttempt | None,
        *,
        recovery_used: int,
    ) -> HamiltonianResult:
        if verified is not None:
            state = CandidateState.MINIMUM_VERIFIED
        elif any(item.state is CandidateState.SADDLE_DETECTED for item in attempts):
            state = CandidateState.SADDLE_DETECTED
        elif any(item.state is CandidateState.VERIFICATION_FAILED for item in attempts):
            state = CandidateState.VERIFICATION_FAILED
        else:
            state = CandidateState.OPTIMIZATION_FAILED
        result = HamiltonianResult(
            hamiltonian=hamiltonian,
            state=state,
            selected_conformer_indices=tuple(item.index for item in conformers),
            attempts=tuple(attempts),
            provisional_lowest_attempt_id=(
                None if provisional is None else provisional.attempt_id
            ),
            verified_attempt_id=None if verified is None else verified.attempt_id,
            recovery_attempts_used=recovery_used,
            recovery_attempt_limit=(self.verification.max_displacement_reoptimizations),
        )
        self.journal.record_hamiltonian_result(result)
        return result


def displace_geometry(
    geometry: ConformerGeometry,
    normal_mode: tuple[tuple[float, float, float], ...],
    *,
    direction: int,
    amplitude_angstrom: float,
) -> ConformerGeometry:
    """Displace along a normal mode with a bounded maximum atom movement."""
    if direction not in {-1, 1}:
        raise ValueError("direction must be -1 or 1")
    if amplitude_angstrom <= 0:
        raise ValueError("amplitude_angstrom must be > 0")
    if len(normal_mode) != geometry.atom_count:
        raise ValueError(
            "Normal mode needs one vector per atom: "
            f"{len(normal_mode)} vectors for {geometry.atom_count} atoms"
        )
    maximum_component = max(abs(value) for vector in normal_mode for value in vector)
    if maximum_component == 0:
        raise ValueError("Normal mode vector must not be all zero")
    scale = direction * amplitude_angstrom / maximum_component
    coordinates = tuple(
        (
            position[0] + scale * vector[0],
            position[1] + scale * vector[1],
            position[2] + scale * vector[2],
        )
        for position, vector in zip(geometry.coordinates, normal_mode)
    )
    return ConformerGeometry(
        elements=geometry.elements,
        coordinates=coordinates,
    )


def _attempt_from_optimization(
    *,
    attempt_id: str,
    hamiltonian: str,
    sequence: int,
    origin: AttemptOrigin,
    source_conformer_index: int,
    optimization: OptimizationRun,
    displacement: DisplacementLineage | None,
) -> CandidateAttempt:
    state = (
        CandidateState.OPTIMIZED_UNVERIFIED
        if optimization.converged
        else CandidateState.OPTIMIZATION_FAILED
    )
    return CandidateAttempt(
        attempt_id=attempt_id,
        hamiltonian=hamiltonian,
        sequence=sequence,
        origin=origin,
        source_conformer_index=source_conformer_index,
        state=state,
        state_history=(state,),
        provisional_heat_of_formation_kcal_mol=(
            optimization.heat_of_formation_kcal_mol
        ),
        optimized_geometry=optimization.geometry,
        optimization_output_path=optimization.output_path,
        optimization_error=optimization.error_message,
        displacement=displacement,
    )


def _attempt_id(hamiltonian: str, sequence: int) -> str:
    return f"{hamiltonian.lower()}-attempt-{sequence:03d}"


def _optimized_positions_by_energy(
    attempts: list[CandidateAttempt],
) -> tuple[int, ...]:
    positions = [
        position
        for position, attempt in enumerate(attempts)
        if attempt.state is CandidateState.OPTIMIZED_UNVERIFIED
    ]
    return tuple(
        sorted(
            positions,
            key=lambda position: (
                _required_provisional_energy(attempts[position]),
                attempts[position].sequence,
            ),
        )
    )


def _lowest_optimized(
    attempts: list[CandidateAttempt], *, selected_only: bool
) -> CandidateAttempt | None:
    candidates = [
        attempt
        for attempt in attempts
        if attempt.provisional_heat_of_formation_kcal_mol is not None
        and (not selected_only or attempt.origin is AttemptOrigin.SELECTED_CONFORMER)
    ]
    if not candidates:
        return None
    return min(
        candidates, key=lambda item: (_required_provisional_energy(item), item.sequence)
    )


def _lowest_verified(attempts: list[CandidateAttempt]) -> CandidateAttempt | None:
    candidates = [attempt for attempt in attempts if attempt.is_verified_minimum]
    if not candidates:
        return None
    return min(
        candidates, key=lambda item: (_required_provisional_energy(item), item.sequence)
    )


def _required_provisional_energy(attempt: CandidateAttempt) -> float:
    value = attempt.provisional_heat_of_formation_kcal_mol
    if value is None:  # pragma: no cover - callers filter optimized attempts
        raise ValueError(f"Attempt {attempt.attempt_id!r} has no provisional energy")
    return value


def _recovery_source(
    attempts: list[CandidateAttempt], force_runs: dict[str, ForceRun]
) -> tuple[CandidateAttempt, ForceRun, int] | None:
    saddles = sorted(
        (
            attempt
            for attempt in attempts
            if attempt.state is CandidateState.SADDLE_DETECTED
        ),
        key=lambda item: (_required_provisional_energy(item), -item.sequence),
    )
    for attempt in saddles:
        diagnostics = attempt.diagnostics
        force_run = force_runs.get(attempt.attempt_id)
        if diagnostics is None or force_run is None:
            continue
        for mode in diagnostics.imaginary_mode_numbers:
            if mode in force_run.normal_modes:
                return attempt, force_run, mode
    return None


__all__ = [
    "JsonWorkflowJournal",
    "MopacMinimumBackend",
    "MopacMinimumWorkflow",
    "NullWorkflowJournal",
    "WorkflowJournal",
    "displace_geometry",
]
