"""Recovery RMSD and dispersion statistics of the validation experiment."""

from __future__ import annotations

import math

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from tests.unit.conformer_validation_support import load

metrics = load("metrics")

SMILES = "CCCC(=O)O"


def _embedded(smiles: str = SMILES) -> Chem.Mol:
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(molecule, randomSeed=11) == 0
    AllChem.MMFFOptimizeMolecule(molecule)
    return molecule


def _elements_coords(
    molecule: Chem.Mol,
) -> tuple[list[str], list[tuple[float, float, float]]]:
    positions = molecule.GetConformer().GetPositions()
    return (
        [atom.GetSymbol() for atom in molecule.GetAtoms()],
        [tuple(float(v) for v in row) for row in positions],
    )


def _xyz(elements: list[str], coords: list[tuple[float, float, float]]) -> str:
    lines = [str(len(elements)), "ref"]
    lines += [f"{e} {x:.6f} {y:.6f} {z:.6f}" for e, (x, y, z) in zip(elements, coords)]
    return "\n".join(lines) + "\n"


def test_reference_in_another_atom_order_and_frame_gives_zero_rmsd() -> None:
    molecule = _embedded()
    elements, coords = _elements_coords(molecule)
    # Reference: atoms reversed, rotated and translated.
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    moved = [tuple(rotation @ np.array(c) + 3.0) for c in coords]
    reference = metrics.reference_skeleton(_xyz(elements[::-1], moved[::-1]))
    probe = metrics.ProbeTemplate(SMILES).skeleton(elements, coords)

    assert metrics.best_rmsd(probe, reference) == pytest.approx(0.0, abs=1e-4)


def test_carboxylic_oxygen_swap_is_absorbed_by_symmetry() -> None:
    # Skeleton bonds are all single, so C(=O)O oxygens are equivalent.
    molecule = _embedded()
    elements, coords = _elements_coords(molecule)
    swapped = list(coords)
    swapped[4], swapped[5] = coords[5], coords[4]
    # The probe's carbonyl O now sits where the reference has the OH O.
    reference = metrics.reference_skeleton(_xyz(elements, coords))
    probe = metrics.ProbeTemplate(SMILES).skeleton(elements, swapped)

    assert metrics.best_rmsd(probe, reference) == pytest.approx(0.0, abs=1e-4)


def test_rotated_dihedral_gives_a_positive_rmsd() -> None:
    molecule = _embedded()
    elements, coords = _elements_coords(molecule)
    reference = metrics.reference_skeleton(_xyz(elements, coords))
    twisted = Chem.Mol(molecule)
    conformer = twisted.GetConformer()
    angle = rdMolTransforms.GetDihedralDeg(conformer, 0, 1, 2, 3)
    rdMolTransforms.SetDihedralDeg(conformer, 0, 1, 2, 3, angle + 120.0)
    probe = metrics.ProbeTemplate(SMILES).skeleton(*_elements_coords(twisted))

    assert metrics.best_rmsd(probe, reference) > 0.3


def test_probe_refuses_geometry_out_of_smiles_order() -> None:
    elements, coords = _elements_coords(_embedded())
    with pytest.raises(ValueError, match="AddHs atom order"):
        metrics.ProbeTemplate(SMILES).skeleton(elements[::-1], coords[::-1])


def test_dispersion_statistics() -> None:
    values = [1.0, 2.0, 3.0, 4.0, 100.0]
    assert metrics.sample_std(values) == pytest.approx(
        math.sqrt(sum((v - 22.0) ** 2 for v in values) / 4)
    )
    assert metrics.median(values) == 3.0
    assert metrics.median([1.0, 2.0]) == 1.5
    assert metrics.median_absolute_deviation(values) == 1.0
    assert metrics.sample_std([1.0]) is None
    assert metrics.median([]) is None
    assert metrics.median_absolute_deviation([]) is None
