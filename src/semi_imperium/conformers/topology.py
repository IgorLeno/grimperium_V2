"""Molecule topology from the SMILES, checked against the real ensemble.

CONFPASS and the folding filter both need bonds addressed by atom
position, and a topology is only correct if it uses the ensemble's own
atom order. The structures CREST samples start from an RDKit embedding
of the SMILES with explicit hydrogens, and CREST keeps its input order,
so the same RDKit atom order is the topology's natural source.

That chain is an assumption about whoever produced the ensemble, so it
is never trusted on its own: every derived topology is checked against
the first conformer — element by element, and every bond by length —
and a mismatch is an error with a stable code, never a silent pass.

RDKit is imported through :func:`importlib.import_module` because its
generated bindings carry no type information; this is the same access
pattern used by :mod:`semi_imperium.conformers.initial_structure`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from importlib import import_module
from typing import Any

from semi_imperium.conformers.backends import (
    ConformerBackendError,
    ConformerRequest,
)
from semi_imperium.conformers.confpass import MoleculeTopology
from semi_imperium.conformers.ensemble import ConformerEnsemble

Chem: Any = import_module("rdkit.Chem")

#: Error codes raised by this module.
TOPOLOGY_PARSE_FAILED = "topology_parse_failed"
TOPOLOGY_BOND_ORDER_UNSUPPORTED = "topology_bond_order_unsupported"
TOPOLOGY_ATOM_ORDER_MISMATCH = "topology_atom_order_mismatch"

#: A bond is accepted up to this multiple of the summed covalent radii.
#: Loose enough for any bonded distance a conformer search produces, and
#: far below the distance between atoms a reordering would pair up.
DEFAULT_BOND_TOLERANCE = 1.3


@dataclass(frozen=True)
class SmilesTopology:
    """Topology provider for :class:`~semi_imperium.conformers.ConformerWorkflow`.

    Derives the topology from ``request.smiles`` and refuses it unless
    it matches the ensemble it will be applied to.
    """

    bond_tolerance: float = DEFAULT_BOND_TOLERANCE

    def __call__(
        self,
        request: ConformerRequest,
        ensemble: ConformerEnsemble,
    ) -> MoleculeTopology:
        """Return the topology of ``request`` in ``ensemble``'s atom order.

        Raises:
            ConformerBackendError: If the SMILES cannot be read, carries a
                bond SDF cannot express, or does not match the ensemble.
        """
        topology, elements = topology_from_smiles(request.smiles)
        require_matching_order(
            topology,
            elements,
            ensemble,
            molecule_id=request.molecule_id,
            bond_tolerance=self.bond_tolerance,
        )
        return topology


def topology_from_smiles(smiles: str) -> tuple[MoleculeTopology, tuple[str, ...]]:
    """Build the explicit-hydrogen topology and element order of ``smiles``.

    The hydrogens are added exactly as the initial-3D route and the CREST
    input generation add them, so the atom order is the one they embed.
    Aromatic bonds are kekulized because SDF bond orders are integers.

    Raises:
        ConformerBackendError: If RDKit cannot parse or kekulize the
            SMILES, or a bond order is not an integer after kekulization.
    """
    parsed = Chem.MolFromSmiles(smiles)
    if parsed is None:
        raise ConformerBackendError(
            f"RDKit could not parse SMILES {smiles!r} to build its topology",
            code=TOPOLOGY_PARSE_FAILED,
        )
    molecule = Chem.AddHs(parsed)
    try:
        Chem.Kekulize(molecule, clearAromaticFlags=True)
    except Exception as exc:  # RDKit raises its own KekulizeException
        raise ConformerBackendError(
            f"RDKit could not kekulize SMILES {smiles!r}: {exc}",
            code=TOPOLOGY_PARSE_FAILED,
        ) from exc

    bonds: list[tuple[int, int, int]] = []
    for bond in molecule.GetBonds():
        order = float(bond.GetBondTypeAsDouble())
        if order < 1 or order != int(order):
            raise ConformerBackendError(
                f"SMILES {smiles!r} has a bond of order {order} between atoms "
                f"{bond.GetBeginAtomIdx()} and {bond.GetEndAtomIdx()}, which "
                "an SDF topology cannot express",
                code=TOPOLOGY_BOND_ORDER_UNSUPPORTED,
            )
        bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), int(order)))
    elements = tuple(str(atom.GetSymbol()) for atom in molecule.GetAtoms())
    topology = MoleculeTopology(atom_count=len(elements), bonds=tuple(bonds))
    return topology, elements


def require_matching_order(
    topology: MoleculeTopology,
    elements: tuple[str, ...],
    ensemble: ConformerEnsemble,
    *,
    molecule_id: str,
    bond_tolerance: float = DEFAULT_BOND_TOLERANCE,
) -> None:
    """Refuse ``topology`` unless it describes ``ensemble`` atom by atom.

    Checked on the first conformer: the same element at every position,
    and every bond no longer than ``bond_tolerance`` times the summed
    covalent radii of its atoms.

    Raises:
        ConformerBackendError: With code ``topology_atom_order_mismatch``
            when any of the checks fails.
    """
    geometry = ensemble.conformers[0].geometry
    observed = tuple(symbol.capitalize() for symbol in geometry.elements)
    if observed != elements:
        position = next(
            (
                index
                for index, (seen, expected) in enumerate(zip(observed, elements))
                if seen != expected
            ),
            min(len(observed), len(elements)),
        )
        raise ConformerBackendError(
            f"The ensemble of {molecule_id!r} does not follow its SMILES atom "
            f"order: {len(observed)} atoms in the geometry, {len(elements)} "
            f"from the SMILES, first difference at atom {position}",
            code=TOPOLOGY_ATOM_ORDER_MISMATCH,
        )

    table = Chem.GetPeriodicTable()
    coordinates = geometry.coordinates
    for first, second, _ in topology.bonds:
        limit = bond_tolerance * (
            float(table.GetRcovalent(elements[first]))
            + float(table.GetRcovalent(elements[second]))
        )
        distance = math.dist(coordinates[first], coordinates[second])
        if distance > limit:
            raise ConformerBackendError(
                f"The ensemble of {molecule_id!r} does not follow its SMILES atom "
                f"order: bond {first}-{second} ({elements[first]}-"
                f"{elements[second]}) is {distance:.2f} Å apart, above the "
                f"{limit:.2f} Å limit",
                code=TOPOLOGY_ATOM_ORDER_MISMATCH,
            )


__all__ = [
    "DEFAULT_BOND_TOLERANCE",
    "TOPOLOGY_ATOM_ORDER_MISMATCH",
    "TOPOLOGY_BOND_ORDER_UNSUPPORTED",
    "TOPOLOGY_PARSE_FAILED",
    "SmilesTopology",
    "require_matching_order",
    "topology_from_smiles",
]
