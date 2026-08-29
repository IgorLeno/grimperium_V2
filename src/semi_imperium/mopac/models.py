"""Auditable data contracts for MOPAC optimization and minimum verification.

The heat of formation produced by geometry optimization is deliberately named
``provisional``.  It becomes reportable only when a subsequent MOPAC ``FORCE``
calculation classifies that exact optimized geometry as a minimum.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from semi_imperium.conformers import ConformerGeometry, SelectionResult

SUPPORTED_HAMILTONIANS = ("AM1", "PM3", "PM7")


class CandidateState(str, Enum):
    """Scientific state of one independently optimized candidate."""

    OPTIMIZATION_FAILED = "optimization_failed"
    OPTIMIZED_UNVERIFIED = "optimized_unverified"
    MINIMUM_VERIFIED = "minimum_verified"
    SADDLE_DETECTED = "saddle_detected"
    VERIFICATION_FAILED = "verification_failed"


@dataclass(frozen=True)
class SelectionLineage:
    """Persistable account of why this finite conformer set was used."""

    strategy: str
    considered: int
    ranking_basis: str
    selected_indices: tuple[int, ...]
    evidence: tuple[str, ...]
    experimental: bool

    @classmethod
    def from_selection(cls, selection: SelectionResult) -> SelectionLineage:
        """Capture the immutable scientific selection fields."""
        return cls(
            strategy=selection.strategy.value,
            considered=selection.considered,
            ranking_basis=selection.ranking_basis,
            selected_indices=selection.selected_indices,
            evidence=selection.evidence,
            experimental=selection.is_experimental,
        )

    def to_dict(self) -> dict[str, Any]:
        """Serialize the selection account next to MOPAC attempts."""
        return {
            "strategy": self.strategy,
            "considered": self.considered,
            "ranking_basis": self.ranking_basis,
            "selected_indices": list(self.selected_indices),
            "evidence": list(self.evidence),
            "experimental": self.experimental,
        }


class AttemptOrigin(str, Enum):
    """Why an optimization entered the finite candidate set."""

    SELECTED_CONFORMER = "selected_conformer"
    DISPLACEMENT_REOPTIMIZATION = "displacement_reoptimization"


@dataclass(frozen=True)
class DisplacementLineage:
    """Explicit lineage of a recovery geometry."""

    parent_attempt_id: str
    normal_mode: int
    direction: int
    amplitude_angstrom: float

    def __post_init__(self) -> None:
        if not self.parent_attempt_id:
            raise ValueError("DisplacementLineage.parent_attempt_id must not be empty")
        if self.normal_mode < 1:
            raise ValueError("DisplacementLineage.normal_mode must be >= 1")
        if self.direction not in {-1, 1}:
            raise ValueError("DisplacementLineage.direction must be -1 or 1")
        if self.amplitude_angstrom <= 0:
            raise ValueError("DisplacementLineage.amplitude_angstrom must be > 0")

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "parent_attempt_id": self.parent_attempt_id,
            "normal_mode": self.normal_mode,
            "direction": self.direction,
            "amplitude_angstrom": self.amplitude_angstrom,
        }


@dataclass(frozen=True)
class OptimizationRun:
    """Raw result returned by an optimization backend."""

    converged: bool
    geometry: ConformerGeometry | None = None
    heat_of_formation_kcal_mol: float | None = None
    output_path: str | None = None
    error_message: str | None = None

    def __post_init__(self) -> None:
        if self.converged and (
            self.geometry is None or self.heat_of_formation_kcal_mol is None
        ):
            raise ValueError(
                "A converged OptimizationRun needs geometry and heat of formation"
            )
        if not self.converged and not self.error_message:
            raise ValueError("A failed OptimizationRun needs an error_message")


@dataclass(frozen=True)
class ForceRun:
    """MOPAC ``FORCE`` output and optional vectors needed for recovery.

    ``normal_modes`` maps the one-based MOPAC vibration number to one Cartesian
    vector per atom.  The parser classifies the text; the vectors remain typed
    backend data because parsing the several version-dependent matrix layouts
    is an adapter concern.
    """

    output: str
    output_path: str | None = None
    normal_modes: dict[int, tuple[tuple[float, float, float], ...]] = field(
        default_factory=dict
    )
    execution_error: str | None = None


@dataclass(frozen=True)
class FrequencyDiagnostics:
    """Parsed evidence behind a minimum/saddle/failure classification."""

    force_calculation_detected: bool
    force_calculation_completed: bool
    stationary_point_confirmed: bool
    frequencies_cm1: tuple[float, ...]
    nontrivial_frequencies_cm1: tuple[float, ...]
    trivial_frequencies_cm1: tuple[float, ...]
    numerical_low_frequencies_cm1: tuple[float, ...]
    imaginary_frequencies_cm1: tuple[float, ...]
    imaginary_mode_numbers: tuple[int, ...]
    expected_trivial_modes: int | None
    frequency_source: str
    gradient_norm: float | None = None
    failure_reason: str | None = None
    notes: tuple[str, ...] = ()

    @property
    def lowest_frequency_cm1(self) -> float | None:
        """Lowest parsed non-trivial frequency, when one was available."""
        if not self.nontrivial_frequencies_cm1:
            return None
        return min(self.nontrivial_frequencies_cm1)

    def to_dict(self) -> dict[str, Any]:
        """Serialize all classification evidence without dropping diagnostics."""
        return {
            "force_calculation_detected": self.force_calculation_detected,
            "force_calculation_completed": self.force_calculation_completed,
            "stationary_point_confirmed": self.stationary_point_confirmed,
            "frequencies_cm1": list(self.frequencies_cm1),
            "nontrivial_frequencies_cm1": list(self.nontrivial_frequencies_cm1),
            "trivial_frequencies_cm1": list(self.trivial_frequencies_cm1),
            "numerical_low_frequencies_cm1": list(self.numerical_low_frequencies_cm1),
            "imaginary_frequencies_cm1": list(self.imaginary_frequencies_cm1),
            "imaginary_mode_numbers": list(self.imaginary_mode_numbers),
            "expected_trivial_modes": self.expected_trivial_modes,
            "frequency_source": self.frequency_source,
            "gradient_norm": self.gradient_norm,
            "failure_reason": self.failure_reason,
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class ForceClassification:
    """Verdict plus diagnostics for one completed verification request."""

    state: CandidateState
    diagnostics: FrequencyDiagnostics

    def __post_init__(self) -> None:
        if self.state not in {
            CandidateState.MINIMUM_VERIFIED,
            CandidateState.SADDLE_DETECTED,
            CandidateState.VERIFICATION_FAILED,
        }:
            raise ValueError("ForceClassification needs a verification state")


@dataclass(frozen=True)
class CandidateAttempt:
    """Final snapshot of one attempt, including every state it passed through."""

    attempt_id: str
    hamiltonian: str
    sequence: int
    origin: AttemptOrigin
    source_conformer_index: int
    state: CandidateState
    state_history: tuple[CandidateState, ...]
    provisional_heat_of_formation_kcal_mol: float | None = None
    optimized_geometry: ConformerGeometry | None = None
    optimization_output_path: str | None = None
    verification_output_path: str | None = None
    optimization_error: str | None = None
    diagnostics: FrequencyDiagnostics | None = None
    displacement: DisplacementLineage | None = None

    def __post_init__(self) -> None:
        if not self.attempt_id:
            raise ValueError("CandidateAttempt.attempt_id must not be empty")
        if self.hamiltonian not in SUPPORTED_HAMILTONIANS:
            raise ValueError(f"Unsupported Hamiltonian {self.hamiltonian!r}")
        if self.sequence < 0:
            raise ValueError("CandidateAttempt.sequence must be >= 0")
        if self.source_conformer_index < 0:
            raise ValueError("CandidateAttempt.source_conformer_index must be >= 0")
        if not self.state_history or self.state_history[-1] is not self.state:
            raise ValueError("CandidateAttempt.state_history must end in state")
        if self.origin is AttemptOrigin.DISPLACEMENT_REOPTIMIZATION:
            if self.displacement is None:
                raise ValueError("A displacement attempt needs displacement lineage")
        elif self.displacement is not None:
            raise ValueError("A selected-conformer attempt cannot have displacement")
        if self.state is not CandidateState.OPTIMIZATION_FAILED and (
            self.provisional_heat_of_formation_kcal_mol is None
            or self.optimized_geometry is None
        ):
            raise ValueError(
                "An optimized CandidateAttempt needs provisional energy and geometry"
            )

    @property
    def is_verified_minimum(self) -> bool:
        """Whether this exact candidate may supply a final heat of formation."""
        return self.state is CandidateState.MINIMUM_VERIFIED

    def to_dict(self) -> dict[str, Any]:
        """Serialize the complete attempt lineage and diagnostics."""
        geometry = self.optimized_geometry
        return {
            "attempt_id": self.attempt_id,
            "hamiltonian": self.hamiltonian,
            "sequence": self.sequence,
            "origin": self.origin.value,
            "source_conformer_index": self.source_conformer_index,
            "state": self.state.value,
            "state_history": [state.value for state in self.state_history],
            "provisional_heat_of_formation_kcal_mol": (
                self.provisional_heat_of_formation_kcal_mol
            ),
            "optimized_geometry": (
                None
                if geometry is None
                else {
                    "elements": list(geometry.elements),
                    "coordinates": [list(point) for point in geometry.coordinates],
                }
            ),
            "optimization_output_path": self.optimization_output_path,
            "verification_output_path": self.verification_output_path,
            "optimization_error": self.optimization_error,
            "diagnostics": (
                None if self.diagnostics is None else self.diagnostics.to_dict()
            ),
            "displacement": (
                None if self.displacement is None else self.displacement.to_dict()
            ),
        }


@dataclass(frozen=True)
class HamiltonianResult:
    """Independent outcome for one of AM1, PM3 or PM7."""

    hamiltonian: str
    state: CandidateState
    selected_conformer_indices: tuple[int, ...]
    attempts: tuple[CandidateAttempt, ...]
    provisional_lowest_attempt_id: str | None
    verified_attempt_id: str | None
    recovery_attempts_used: int
    recovery_attempt_limit: int

    def __post_init__(self) -> None:
        if self.hamiltonian not in SUPPORTED_HAMILTONIANS:
            raise ValueError(f"Unsupported Hamiltonian {self.hamiltonian!r}")
        if self.recovery_attempts_used > self.recovery_attempt_limit:
            raise ValueError("Recovery attempts exceeded their configured limit")
        attempt_ids = {attempt.attempt_id for attempt in self.attempts}
        if self.provisional_lowest_attempt_id not in attempt_ids | {None}:
            raise ValueError("Unknown provisional_lowest_attempt_id")
        if self.verified_attempt_id not in attempt_ids | {None}:
            raise ValueError("Unknown verified_attempt_id")
        if self.state is CandidateState.MINIMUM_VERIFIED:
            if self.verified_attempt_id is None:
                raise ValueError("A verified result needs verified_attempt_id")
        elif self.verified_attempt_id is not None:
            raise ValueError("A non-minimum result cannot name a verified attempt")

    @property
    def verified_heat_of_formation_kcal_mol(self) -> float | None:
        """Return a final value only for a ``MINIMUM_VERIFIED`` candidate."""
        if self.verified_attempt_id is None:
            return None
        for attempt in self.attempts:
            if attempt.attempt_id == self.verified_attempt_id:
                if not attempt.is_verified_minimum:
                    return None
                return attempt.provisional_heat_of_formation_kcal_mol
        return None

    def to_dict(self) -> dict[str, Any]:
        """Serialize the independent result and its bounded candidate set."""
        return {
            "hamiltonian": self.hamiltonian,
            "state": self.state.value,
            "selected_conformer_indices": list(self.selected_conformer_indices),
            "attempts": [attempt.to_dict() for attempt in self.attempts],
            "provisional_lowest_attempt_id": self.provisional_lowest_attempt_id,
            "verified_attempt_id": self.verified_attempt_id,
            "verified_heat_of_formation_kcal_mol": (
                self.verified_heat_of_formation_kcal_mol
            ),
            "recovery_attempts_used": self.recovery_attempts_used,
            "recovery_attempt_limit": self.recovery_attempt_limit,
        }


@dataclass(frozen=True)
class MinimumWorkflowResult:
    """Outcomes for all requested Hamiltonians, kept separate by construction."""

    selection_lineage: SelectionLineage
    hamiltonian_results: tuple[HamiltonianResult, ...]

    def __post_init__(self) -> None:
        methods = [result.hamiltonian for result in self.hamiltonian_results]
        if len(set(methods)) != len(methods):
            raise ValueError("MinimumWorkflowResult has duplicate Hamiltonians")

    def for_hamiltonian(self, hamiltonian: str) -> HamiltonianResult:
        """Return one independent result by normalized Hamiltonian name."""
        normalized = hamiltonian.strip().upper()
        for result in self.hamiltonian_results:
            if result.hamiltonian == normalized:
                return result
        raise KeyError(normalized)

    def to_dict(self) -> dict[str, Any]:
        """Serialize all independently selected results."""
        return {
            "selection_lineage": self.selection_lineage.to_dict(),
            "hamiltonian_results": [
                result.to_dict() for result in self.hamiltonian_results
            ],
        }


__all__ = [
    "SUPPORTED_HAMILTONIANS",
    "AttemptOrigin",
    "CandidateAttempt",
    "CandidateState",
    "DisplacementLineage",
    "ForceRun",
    "ForceClassification",
    "FrequencyDiagnostics",
    "HamiltonianResult",
    "MinimumWorkflowResult",
    "OptimizationRun",
    "SelectionLineage",
]
