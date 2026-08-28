"""Conformer geometries, ensembles and their search provenance.

A conformer here is inert data: the element order it belongs to, the
coordinates, and the energy the producing program reported for it.
Nothing in this module talks to CREST, RDKit or CONFPASS, which is what
lets the orchestration around it be exercised without any of them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any

from semi_imperium.domain.configuration import ConformerSearchSettings
from semi_imperium.domain.enums import ConformerSource
from semi_imperium.domain.records import utc_now

#: Version string recorded when a program does not report one. It is
#: spelled out on purpose: a blank version would read as "not recorded".
UNKNOWN_PROGRAM_VERSION = "unknown"


@dataclass(frozen=True)
class ConformerGeometry:
    """Element symbols and Cartesian coordinates, in one fixed order.

    The order is meaningful: connectivity, MOPAC input files and any
    SDF adaptation all address atoms by position, so a geometry that
    reorders its atoms silently invalidates everything downstream.
    """

    elements: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]

    def __post_init__(self) -> None:
        if not self.elements:
            raise ValueError("ConformerGeometry.elements must not be empty")
        if len(self.elements) != len(self.coordinates):
            raise ValueError(
                "ConformerGeometry needs one coordinate triple per element: "
                f"{len(self.elements)} elements, {len(self.coordinates)} triples"
            )
        for position in self.coordinates:
            if len(position) != 3:
                raise ValueError(
                    "ConformerGeometry coordinates must be (x, y, z) triples, "
                    f"got {position!r}"
                )

    @property
    def atom_count(self) -> int:
        """Number of atoms in this geometry."""
        return len(self.elements)


@dataclass(frozen=True)
class Conformer:
    """One sampled structure with the energy its producer reported.

    ``energy_kcal_mol`` stays ``None`` when the producing route has no
    energy to report — the single embedded structure of the initial-3D
    route, for instance. That is not a missing value to be guessed at;
    it is why energy-based ranking refuses to run on such an ensemble.
    """

    index: int
    geometry: ConformerGeometry
    energy_kcal_mol: float | None = None
    label: str = ""

    def __post_init__(self) -> None:
        if self.index < 0:
            raise ValueError(f"Conformer.index must be >= 0, got {self.index}")
        if not self.label:
            object.__setattr__(self, "label", f"conf{self.index:03d}")

    @property
    def atom_count(self) -> int:
        """Number of atoms in this conformer."""
        return self.geometry.atom_count

    def require_energy(self) -> float:
        """Return the reported energy, or fail loudly when there is none."""
        if self.energy_kcal_mol is None:
            raise ValueError(
                f"Conformer {self.label!r} carries no energy; an energy-based "
                "ranking cannot be applied to it"
            )
        return self.energy_kcal_mol


@dataclass(frozen=True)
class ConformerSearchProvenance:
    """How an ensemble was produced: settings, version and run.

    This is what makes a stored geometry defensible later: the effective
    search settings, the version of the executable that produced it and
    the run it belongs to, recorded next to the structures themselves.
    """

    source: ConformerSource
    program: str
    program_version: str
    settings: ConformerSearchSettings
    run_id: str | None = None
    command: tuple[str, ...] = ()
    generated_at: datetime = field(default_factory=utc_now)

    def __post_init__(self) -> None:
        if not self.program.strip():
            raise ValueError("ConformerSearchProvenance.program must not be empty")
        if not self.program_version.strip():
            raise ValueError(
                "ConformerSearchProvenance.program_version must not be empty; "
                f"record {UNKNOWN_PROGRAM_VERSION!r} explicitly instead"
            )

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "source": self.source.value,
            "program": self.program,
            "program_version": self.program_version,
            "settings": self.settings.to_dict(),
            "run_id": self.run_id,
            "command": list(self.command),
            "generated_at": self.generated_at.isoformat(),
        }


@dataclass(frozen=True)
class ConformerEnsemble:
    """Every conformer one search produced, before any selection.

    The ensemble is deliberately the *complete* output of the search:
    strategies receive it whole and record how many conformers they
    considered, so a narrow selection can never be mistaken for a narrow
    search.
    """

    conformers: tuple[Conformer, ...]
    provenance: ConformerSearchProvenance

    def __post_init__(self) -> None:
        if not self.conformers:
            raise ValueError("ConformerEnsemble must contain at least one conformer")
        indices = [conformer.index for conformer in self.conformers]
        if len(set(indices)) != len(indices):
            raise ValueError(
                f"ConformerEnsemble indices must be unique, got {sorted(indices)}"
            )
        atom_counts = {conformer.atom_count for conformer in self.conformers}
        if len(atom_counts) != 1:
            raise ValueError(
                "ConformerEnsemble conformers must all describe the same molecule, "
                f"got atom counts {sorted(atom_counts)}"
            )

    @property
    def size(self) -> int:
        """How many conformers the search produced."""
        return len(self.conformers)

    @property
    def source(self) -> ConformerSource:
        """Which route produced these structures."""
        return self.provenance.source

    @property
    def has_energies(self) -> bool:
        """Whether every conformer carries an energy to rank by."""
        return all(
            conformer.energy_kcal_mol is not None for conformer in self.conformers
        )

    def ranked_by_energy(self) -> tuple[Conformer, ...]:
        """Return the conformers ordered by the energy the search reported.

        Ties keep the producing program's own order, so the ranking is
        reproducible for ensembles that contain degenerate structures.

        Raises:
            ValueError: If any conformer carries no energy.
        """

        def sort_key(conformer: Conformer) -> tuple[float, int]:
            return (conformer.require_energy(), conformer.index)

        return tuple(sorted(self.conformers, key=sort_key))


__all__ = [
    "UNKNOWN_PROGRAM_VERSION",
    "Conformer",
    "ConformerEnsemble",
    "ConformerGeometry",
    "ConformerSearchProvenance",
]
