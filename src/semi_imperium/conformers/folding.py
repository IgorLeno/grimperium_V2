"""Optional pre-selection filter for extremely folded conformers.

A CREST ensemble can contain structures whose chain ends fold back onto
each other. This stage measures that directly from geometry and
connectivity, drops such conformers before any selection strategy runs,
and records what it did. It never empties an ensemble: when every
conformer would be dropped it keeps them all and says so in the
evidence (``folding_filter_bypassed_all_folded``).

The descriptor is the one calibrated by ``scripts/calibrate_folding_filter.py``
(see ``reports/folding_calibration/report.md``), recomputed here from
:class:`MoleculeTopology` so the conformer stage stays free of RDKit.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass
from typing import Any

from semi_imperium.conformers.confpass import MoleculeTopology
from semi_imperium.conformers.ensemble import (
    Conformer,
    ConformerEnsemble,
    ConformerGeometry,
)
from semi_imperium.domain.configuration import FoldingFilterSettings

#: Van der Waals radii (Å) from the RDKit periodic table (2026.03.3), the
#: same source the calibration used. The calibration itself only covered
#: C, N and O; the other entries apply the same contact criterion.
VDW_RADII_ANGSTROM: dict[str, float] = {
    "B": 1.8,
    "C": 1.7,
    "N": 1.6,
    "O": 1.55,
    "F": 1.5,
    "Si": 2.1,
    "P": 1.95,
    "S": 1.8,
    "Cl": 1.8,
    "Se": 1.9,
    "Br": 1.9,
    "I": 2.1,
}

#: Elements that act as both H-bond donor (when carrying H) and acceptor.
HBOND_ELEMENTS = frozenset({"N", "O"})

BYPASSED_ALL_FOLDED = "folding_filter_bypassed_all_folded"


@dataclass(frozen=True)
class FoldingMeasurement:
    """Descriptor values for one conformer."""

    index: int
    fold_contacts: int
    rg_ratio: float
    folded: bool

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "index": self.index,
            "fold_contacts": self.fold_contacts,
            "rg_ratio": self.rg_ratio,
            "folded": self.folded,
        }


@dataclass(frozen=True)
class FoldingFilterOutcome:
    """The ensemble the selection strategy receives, and why."""

    settings: FoldingFilterSettings
    kept: ConformerEnsemble
    measurements: tuple[FoldingMeasurement, ...]
    bypassed: bool

    @property
    def discarded_indices(self) -> tuple[int, ...]:
        """Ensemble indices removed by the filter (empty when bypassed)."""
        if self.bypassed:
            return ()
        return tuple(item.index for item in self.measurements if item.folded)

    @property
    def evidence(self) -> tuple[str, ...]:
        """Short, replayable statements about what the filter did."""
        folded = sum(item.folded for item in self.measurements)
        entries = [
            f"folding_filter_folded={folded}/{len(self.measurements)}",
        ]
        if self.bypassed:
            entries.append(BYPASSED_ALL_FOLDED)
        return tuple(entries)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "settings": self.settings.to_dict(),
            "bypassed": self.bypassed,
            "discarded_indices": list(self.discarded_indices),
            "measurements": [item.to_dict() for item in self.measurements],
            "evidence": list(self.evidence),
        }


def apply_folding_filter(
    ensemble: ConformerEnsemble,
    topology: MoleculeTopology,
    settings: FoldingFilterSettings,
) -> FoldingFilterOutcome:
    """Drop extremely folded conformers, never emptying the ensemble.

    Raises:
        ValueError: If the filter is disabled, the topology does not match
            the ensemble's atom count, or an element has no known radius.
    """
    if not settings.enabled:
        raise ValueError("apply_folding_filter was called with a disabled filter")
    elements = ensemble.conformers[0].geometry.elements
    if topology.atom_count != len(elements):
        raise ValueError(
            f"Folding filter topology has {topology.atom_count} atoms but the "
            f"ensemble has {len(elements)}"
        )

    neighbors = _neighbors(topology)
    distances = _topological_distances(neighbors)
    heavy = [i for i, element in enumerate(elements) if element != "H"]
    radii = {i: _vdw_radius(elements[i]) for i in heavy}
    candidate_pairs: list[tuple[int, int]] = []
    for position, i in enumerate(heavy):
        for j in heavy[position + 1 :]:
            bonds = distances[i][j]
            if bonds is not None and bonds >= settings.min_topological_distance:
                candidate_pairs.append((i, j))

    gyration = {
        conformer.index: _radius_of_gyration(conformer.geometry, heavy)
        for conformer in ensemble.conformers
    }
    largest = max(gyration.values())

    measurements = []
    for conformer in ensemble.conformers:
        contacts = _fold_contacts(
            conformer, candidate_pairs, radii, neighbors, settings
        )
        rg_ratio = gyration[conformer.index] / largest if largest > 0 else 1.0
        folded = contacts >= settings.min_contacts and (
            settings.max_rg_ratio is None or rg_ratio < settings.max_rg_ratio
        )
        measurements.append(
            FoldingMeasurement(
                index=conformer.index,
                fold_contacts=contacts,
                rg_ratio=rg_ratio,
                folded=folded,
            )
        )

    folded_indices = {item.index for item in measurements if item.folded}
    kept = tuple(c for c in ensemble.conformers if c.index not in folded_indices)
    bypassed = not kept
    return FoldingFilterOutcome(
        settings=settings,
        kept=(
            ensemble
            if bypassed or not folded_indices
            else ConformerEnsemble(conformers=kept, provenance=ensemble.provenance)
        ),
        measurements=tuple(measurements),
        bypassed=bypassed,
    )


def _fold_contacts(
    conformer: Conformer,
    pairs: list[tuple[int, int]],
    radii: dict[int, float],
    neighbors: list[list[int]],
    settings: FoldingFilterSettings,
) -> int:
    """Count heavy-atom pairs in extreme contact, H-bonded pairs excluded."""
    geometry = conformer.geometry
    hbonded = _hbond_pairs(geometry, neighbors, settings.hbond_max_h_acceptor_angstrom)
    return sum(
        1
        for i, j in pairs
        if (i, j) not in hbonded
        and _distance(geometry, i, j) < settings.vdw_scale * (radii[i] + radii[j])
    )


def _hbond_pairs(
    geometry: ConformerGeometry,
    neighbors: list[list[int]],
    cutoff: float,
) -> set[tuple[int, int]]:
    """Heavy-atom (donor, acceptor) pairs joined by a classic H-bond."""
    elements = geometry.elements
    polar = [i for i, element in enumerate(elements) if element in HBOND_ELEMENTS]
    pairs: set[tuple[int, int]] = set()
    for donor in polar:
        hydrogens = [n for n in neighbors[donor] if elements[n] == "H"]
        for hydrogen in hydrogens:
            for acceptor in polar:
                if acceptor == donor:
                    continue
                if _distance(geometry, hydrogen, acceptor) < cutoff:
                    pairs.add((min(donor, acceptor), max(donor, acceptor)))
    return pairs


def _neighbors(topology: MoleculeTopology) -> list[list[int]]:
    adjacency: list[list[int]] = [[] for _ in range(topology.atom_count)]
    for first, second, _order in topology.bonds:
        adjacency[first].append(second)
        adjacency[second].append(first)
    return adjacency


def _topological_distances(neighbors: list[list[int]]) -> list[list[int | None]]:
    """Bond counts between every atom pair; ``None`` across fragments.

    Pairs in different fragments are never counted as folding contacts:
    two fragments lying close together is not a chain folding back.
    """
    size = len(neighbors)
    table: list[list[int | None]] = []
    for start in range(size):
        row: list[int | None] = [None] * size
        row[start] = 0
        queue = deque([(start, 0)])
        while queue:
            atom, step = queue.popleft()
            for neighbor in neighbors[atom]:
                if row[neighbor] is None:
                    row[neighbor] = step + 1
                    queue.append((neighbor, step + 1))
        table.append(row)
    return table


def _radius_of_gyration(geometry: ConformerGeometry, atoms: list[int]) -> float:
    """Unweighted radius of gyration over ``atoms`` (heavy atoms)."""
    if not atoms:
        return 0.0
    points = [geometry.coordinates[i] for i in atoms]
    centre = tuple(sum(axis) / len(points) for axis in zip(*points, strict=True))
    squared = sum(
        sum((p - c) ** 2 for p, c in zip(point, centre, strict=True))
        for point in points
    )
    return math.sqrt(squared / len(points))


def _distance(geometry: ConformerGeometry, first: int, second: int) -> float:
    return math.dist(geometry.coordinates[first], geometry.coordinates[second])


def _vdw_radius(element: str) -> float:
    try:
        return VDW_RADII_ANGSTROM[element]
    except KeyError as exc:
        known = ", ".join(sorted(VDW_RADII_ANGSTROM))
        raise ValueError(
            f"Folding filter has no van der Waals radius for {element!r}; "
            f"known elements: {known}"
        ) from exc


__all__ = [
    "BYPASSED_ALL_FOLDED",
    "HBOND_ELEMENTS",
    "VDW_RADII_ANGSTROM",
    "FoldingFilterOutcome",
    "FoldingMeasurement",
    "apply_folding_filter",
]
