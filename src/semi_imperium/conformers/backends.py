"""Adapter boundary for every external conformer program.

CREST, the RDKit initial-structure route and CONFPASS are reached only
through the protocols declared here. The orchestration depends on these
protocols and never on a binary, so unit tests drive the whole workflow
with plain in-memory doubles and nothing external is ever spawned.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol

from semi_imperium.conformers.ensemble import ConformerEnsemble
from semi_imperium.domain.configuration import ConformerSearchSettings


class ConformerBackendError(RuntimeError):
    """An external conformer program failed, is absent or returned junk.

    Carries a stable ``code`` so callers can react to the reason without
    matching on message text.
    """

    def __init__(self, message: str, *, code: str) -> None:
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class ConformerRequest:
    """What an external conformer program needs about one molecule."""

    molecule_id: str
    smiles: str
    charge: int = 0
    multiplicity: int = 1
    run_id: str | None = None

    def __post_init__(self) -> None:
        if not self.molecule_id.strip():
            raise ValueError("ConformerRequest.molecule_id must not be empty")
        if not self.smiles.strip():
            raise ValueError("ConformerRequest.smiles must not be empty")
        if self.multiplicity < 1:
            raise ValueError(
                f"ConformerRequest.multiplicity must be >= 1, got {self.multiplicity}"
            )


@dataclass(frozen=True)
class ConfPassCandidate:
    """One conformer handed to CONFPASS, already adapted to SDF."""

    index: int
    sd_record: str
    energy_kcal_mol: float | None = None


@dataclass(frozen=True)
class ConfPassRanking:
    """CONFPASS's opinion about one candidate.

    ``pas_completeness_class`` is the PAS classifier's label. It is
    carried here so it can be *recorded*, never so it can be used as
    evidence about the science; the selection layer enforces that.
    """

    index: int
    priority: int
    pas_completeness_class: str | None = None

    def __post_init__(self) -> None:
        if self.priority < 0:
            raise ValueError(
                f"ConfPassRanking.priority must be >= 0, got {self.priority}"
            )


class ConformerSearchBackend(Protocol):
    """An external conformer search program, e.g. CREST."""

    def search(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        """Return the complete sampled ensemble, provenance included."""
        ...


class InitialStructureBackend(Protocol):
    """The route that builds one starting structure without a search."""

    def build(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        """Return a single-conformer ensemble usable by MOPAC."""
        ...


class ConfPassBackend(Protocol):
    """The experimental CONFPASS prioritizer."""

    def prioritize(
        self,
        candidates: Sequence[ConfPassCandidate],
    ) -> Sequence[ConfPassRanking]:
        """Return one ranking per candidate it was given."""
        ...


__all__ = [
    "ConfPassBackend",
    "ConfPassCandidate",
    "ConfPassRanking",
    "ConformerBackendError",
    "ConformerRequest",
    "ConformerSearchBackend",
    "InitialStructureBackend",
]
