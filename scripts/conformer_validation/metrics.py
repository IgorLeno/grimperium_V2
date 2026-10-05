"""Recovery and dispersion metrics of the validation experiment.

RMSD: both structures are reduced to a heavy-atom connectivity graph
with every bond single and no charges, so symmetry-equivalent mappings
come from the graph alone and bond-order perception never decides
whether two structures can be compared. The reference graph is
perceived from its XYZ (``DetermineConnectivity``); the probe graph is
the SMILES topology in RDKit ``AddHs`` order, which is the order of
every CREST/MOPAC geometry. ``rdMolAlign.GetBestRMS`` then aligns over
all graph automorphisms.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolAlign
from rdkit.Geometry import Point3D

#: A selected geometry "recovers" the reference below this RMSD (angstrom).
RECOVERY_THRESHOLD_ANGSTROM = 0.5


def _skeleton(
    atomic_numbers: Sequence[int],
    bonds: Sequence[tuple[int, int]],
    coordinates: Sequence[tuple[float, float, float]],
) -> Any:
    """Heavy-atom molecule with only single bonds and the given coordinates."""
    editable = Chem.RWMol()
    for number in atomic_numbers:
        atom = Chem.Atom(int(number))
        atom.SetNoImplicit(True)
        editable.AddAtom(atom)
    for first, second in bonds:
        editable.AddBond(int(first), int(second), Chem.BondType.SINGLE)
    molecule = editable.GetMol()
    molecule.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(molecule)
    conformer = Chem.Conformer(len(atomic_numbers))
    for index, (x, y, z) in enumerate(coordinates):
        conformer.SetAtomPosition(index, Point3D(x, y, z))
    molecule.AddConformer(conformer, assignId=True)
    return molecule


def _heavy_subgraph(
    elements: Sequence[str],
    bonds: Sequence[tuple[int, int]],
    coordinates: Sequence[tuple[float, float, float]],
) -> Any:
    """Drop hydrogens, renumber heavy atoms, keep heavy-heavy bonds."""
    table = Chem.GetPeriodicTable()
    heavy = [i for i, element in enumerate(elements) if element != "H"]
    position = {atom: new for new, atom in enumerate(heavy)}
    return _skeleton(
        [table.GetAtomicNumber(elements[i]) for i in heavy],
        [
            (position[a], position[b])
            for a, b in bonds
            if a in position and b in position
        ],
        [coordinates[i] for i in heavy],
    )


def reference_skeleton(xyz_block: str) -> Any:
    """Heavy-atom skeleton of the reference geometry, bonds perceived."""
    perceived = Chem.MolFromXYZBlock(xyz_block)
    if perceived is None:
        raise ValueError("Reference XYZ block could not be parsed")
    rdDetermineBonds.DetermineConnectivity(perceived)
    positions = perceived.GetConformer().GetPositions()
    return _heavy_subgraph(
        [atom.GetSymbol() for atom in perceived.GetAtoms()],
        [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in perceived.GetBonds()],
        [tuple(float(v) for v in row) for row in positions],
    )


class ProbeTemplate:
    """SMILES topology in ``AddHs`` order, ready to receive coordinates."""

    def __init__(self, smiles: str) -> None:
        parsed = Chem.MolFromSmiles(smiles)
        if parsed is None:
            raise ValueError(f"RDKit could not parse SMILES {smiles!r}")
        molecule = Chem.AddHs(parsed)
        self.elements = tuple(atom.GetSymbol() for atom in molecule.GetAtoms())
        self.bonds = tuple(
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in molecule.GetBonds()
        )

    def skeleton(
        self,
        elements: Sequence[str],
        coordinates: Sequence[tuple[float, float, float]],
    ) -> Any:
        """Heavy-atom skeleton of a geometry in this template's atom order.

        Raises:
            ValueError: If the geometry's elements are not the template's.
        """
        if tuple(elements) != self.elements:
            raise ValueError(
                "Geometry elements do not follow the SMILES AddHs atom order"
            )
        return _heavy_subgraph(self.elements, self.bonds, coordinates)


def best_rmsd(probe: Any, reference: Any) -> float:
    """Symmetry-aware heavy-atom RMSD; ``probe`` is aligned, not changed."""
    return float(rdMolAlign.GetBestRMS(Chem.Mol(probe), reference))


def sample_std(values: Sequence[float]) -> float | None:
    """Sample standard deviation (n - 1), or ``None`` below two values."""
    if len(values) < 2:
        return None
    mean = sum(values) / len(values)
    return math.sqrt(sum((v - mean) ** 2 for v in values) / (len(values) - 1))


def median(values: Sequence[float]) -> float | None:
    """Median, or ``None`` for no values."""
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def median_absolute_deviation(values: Sequence[float]) -> float | None:
    """Unscaled median absolute deviation from the median."""
    centre = median(values)
    if centre is None:
        return None
    return median([abs(v - centre) for v in values])


__all__ = [
    "RECOVERY_THRESHOLD_ANGSTROM",
    "ProbeTemplate",
    "best_rmsd",
    "median",
    "median_absolute_deviation",
    "reference_skeleton",
    "sample_std",
]
