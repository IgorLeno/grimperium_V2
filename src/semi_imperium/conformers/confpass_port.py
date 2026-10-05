"""EXPERIMENTAL port of CONFPASS PART 1: torsion-diversity prioritization.

Ported from CONFPASS (https://github.com/Goodman-lab/CONFPASS), commit
``1b5efb69585ea1f51bedccfed1d9d07133c18a53``, modules
``isolate_key_dihedral_v5``, ``cal_dihedral_v2``, ``dihedral_parameter_v2``,
``correcting_dihedral_v1``, ``clustering_dih_v7`` and ``GetPriority_v3``.
Only PART 1 is ported; the PART 2 PAS classifier is not, so no ranking
carries a completeness class.

    MIT License

    Copyright (c) 2022 Goodman Lab

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The port is held to golden files produced by the original
(``tests/fixtures/confpass_golden``). Its rules are kept as the original
wrote them, quirks included, because the golden files are the contract:
atoms are named ``"<symbol> <1-based index>"`` and several tests read
only the first character(s) of that name. Deliberate differences:

* the candidates are read in the order given (CREST's order); like the
  original, energies are never read, so "first member of a cluster" means
  "first in the file";
* where the original iterates a ``set`` (dihedral end atom fallback,
  surviving dihedral columns), the port uses atom/bond order, so results
  no longer depend on ``PYTHONHASHSEED``;
* the ward tree is built once and cut at every cluster count instead of
  refitting per count; sklearn builds the same full scipy ward tree for
  unstructured data either way;
* an ensemble with no variable dihedral, which makes every original
  method crash, raises :class:`ConformerBackendError` with code
  ``confpass_no_variable_dihedral`` so the caller can fall back openly.

RDKit and scikit-learn are imported through
:func:`importlib.import_module` because they carry no type information,
the same pattern as :mod:`semi_imperium.conformers.initial_structure`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from importlib import import_module
from typing import Any

import numpy as np

from semi_imperium.conformers.backends import (
    ConformerBackendError,
    ConfPassCandidate,
    ConfPassRanking,
)
from semi_imperium.conformers.confpass import CONFPASS_NO_VARIABLE_DIHEDRAL

Chem: Any = import_module("rdkit.Chem")
_cluster: Any = import_module("sklearn.cluster")

#: Priority methods of ``GetPriority_v3`` that this port reproduces.
METHODS = ("pipe_x_as", "pipe_x", "pipe_as")
DEFAULT_METHOD = "pipe_x_as"
DEFAULT_X = 0.8
DEFAULT_X_AS = 0.2

NO_VARIABLE_DIHEDRAL = CONFPASS_NO_VARIABLE_DIHEDRAL

_END_TAG = "M  END"
_HALOGEN_PREFIXES = {"F ", "Cl", "I ", "Br"}
_TRIHALO_NEIGHBOURS = (
    ["F ", "F ", "F "],
    ["Cl", "Cl", "Cl"],
    ["I ", "I ", "I "],
    ["Br", "Br", "Br"],
)

Bond = list[Any]
"""``[first_name, second_name, bond_order]`` as the original keeps it."""


@dataclass(frozen=True)
class ConfPassAnalysis:
    """Every intermediate the priority lists are derived from.

    Positions are 0-based in candidate order. ``clusterings[k - 1]`` holds
    the partition into ``k`` clusters, each cluster sorted, clusters in
    order of their first member.
    """

    descriptor_columns: tuple[str, ...]
    clusterings: tuple[tuple[tuple[int, ...], ...], ...]

    @property
    def conformer_count(self) -> int:
        """How many conformers were clustered."""
        return len(self.clusterings)

    def clusters_at(self, fraction: float) -> tuple[tuple[int, ...], ...]:
        """The partition into ``round(fraction * conformers)`` clusters."""
        count = round(self.conformer_count * fraction)
        if not 1 <= count <= self.conformer_count:
            raise ConformerBackendError(
                f"CONFPASS cannot cut {self.conformer_count} conformers into "
                f"{count} clusters (fraction {fraction})",
                code="confpass_invalid_parameters",
            )
        return self.clusterings[count - 1]

    def priority(
        self,
        method: str = DEFAULT_METHOD,
        *,
        x: float = DEFAULT_X,
        x_as: float = DEFAULT_X_AS,
    ) -> tuple[int, ...]:
        """The priority list of ``method``, as ``GetPriority.priority_df``."""
        if method == "pipe_as":
            return self._pipe_as()
        if method == "pipe_x":
            return self._pipe_x(x)
        if method == "pipe_x_as":
            head = list(self._pipe_x(x)[: round(self.conformer_count * x_as)])
            for position in self._pipe_as():
                if position not in head:
                    head.append(position)
            return tuple(head)
        raise ConformerBackendError(
            f"Unknown CONFPASS method {method!r}; expected one of {METHODS}",
            code="confpass_invalid_parameters",
        )

    def _pipe_as(self) -> tuple[int, ...]:
        """First member of every cluster, from 1 cluster up to one each."""
        ordered: list[int] = []
        for partition in self.clusterings:
            for cluster in partition:
                if cluster[0] not in ordered:
                    ordered.append(cluster[0])
        return tuple(ordered)

    def _pipe_x(self, x: float) -> tuple[int, ...]:
        """First member of each cluster at ``x``, then the other members."""
        partition = self.clusters_at(x)
        firsts = [cluster[0] for cluster in partition]
        rest = [member for cluster in partition for member in cluster[1:]]
        return tuple(firsts + rest)


@dataclass(frozen=True)
class PortedConfPass:
    """EXPERIMENTAL CONFPASS PART 1 backend; never reads energies.

    Defaults are CONFPASS's own: ``pipe_x_as`` with ``x=0.8`` and
    ``x_as=0.2``.
    """

    method: str = DEFAULT_METHOD
    x: float = DEFAULT_X
    x_as: float = DEFAULT_X_AS

    def __post_init__(self) -> None:
        if self.method not in METHODS:
            raise ValueError(
                f"Unknown CONFPASS method {self.method!r}; expected one of {METHODS}"
            )
        for name, value in (("x", self.x), ("x_as", self.x_as)):
            if not 0.0 < value <= 1.0:
                raise ValueError(
                    f"PortedConfPass.{name} must be in (0, 1], got {value}"
                )

    def prioritize(
        self,
        candidates: Sequence[ConfPassCandidate],
    ) -> Sequence[ConfPassRanking]:
        """Rank every candidate; priority 0 is the first to optimize.

        Raises:
            ConformerBackendError: ``confpass_too_few_conformers`` below two
                candidates, ``confpass_no_variable_dihedral`` when no
                rotatable dihedral varies across the ensemble, or another
                code when the input or the parameters are unusable.
        """
        analysis = analyse([candidate.sd_record for candidate in candidates])
        order = analysis.priority(self.method, x=self.x, x_as=self.x_as)
        return tuple(
            ConfPassRanking(index=candidates[position].index, priority=priority)
            for priority, position in enumerate(order)
        )


def analyse(sd_records: Sequence[str]) -> ConfPassAnalysis:
    """Run ``clustering_dih_v7.get_cluster_df`` on SD records, in order."""
    columns, descriptor = dihedral_descriptor(sd_records)
    return ConfPassAnalysis(
        descriptor_columns=columns,
        clusterings=_ward_partitions(descriptor),
    )


def dihedral_descriptor(
    sd_records: Sequence[str],
) -> tuple[tuple[str, ...], np.ndarray[Any, Any]]:
    """The corrected dihedral matrix CONFPASS clusters, one row per record."""
    if len(sd_records) < 2:
        raise ConformerBackendError(
            f"CONFPASS needs at least two conformers to cluster, got "
            f"{len(sd_records)}",
            code="confpass_too_few_conformers",
        )
    mol, mol_h = _read_molecule(sd_records[0])
    dihedrals, bonds = isolate_dihedrals(mol, mol_h)
    columns = [f"{first}_{second}" for first, second in bonds]
    angles = np.array(
        [_dihedral_values(dihedrals, _coordinates(record)) for record in sd_records],
        dtype=float,
    ).reshape(len(sd_records), len(columns))

    fixed = set(_fixed_columns(angles, columns))
    kept = [column for column in range(len(columns)) if columns[column] not in fixed]
    if not kept:
        raise ConformerBackendError(
            "CONFPASS found no rotatable dihedral that varies across the "
            f"{len(sd_records)} conformers, so it has nothing to cluster on",
            code=NO_VARIABLE_DIHEDRAL,
        )
    descriptor = np.column_stack([_correct_by_gap(angles[:, i]) for i in kept])
    if not np.isfinite(descriptor).all():
        raise ConformerBackendError(
            "CONFPASS computed a non-finite dihedral angle; the geometry has "
            "collinear atoms in a selected dihedral",
            code="confpass_invalid_dihedral",
        )
    return tuple(columns[i] for i in kept), descriptor


# ---------------------------------------------------------------------------
# isolate_key_dihedral_v5
# ---------------------------------------------------------------------------


def isolate_dihedrals(mol: Any, mol_h: Any) -> tuple[list[list[str]], list[list[str]]]:
    """``isolate_dihedral_mol``: the dihedrals to measure and their bonds."""
    names_h = [f"{atom.GetSymbol()} {atom.GetIdx() + 1}" for atom in mol_h.GetAtoms()]
    names = [name for name in names_h if not name.startswith("H ")]
    neighbours = _neighbour_map(mol, names)
    neighbours_h = _neighbour_map(mol_h, names_h)
    bond_list: list[Bond] = [
        [
            names[bond.GetBeginAtomIdx()],
            names[bond.GetEndAtomIdx()],
            bond.GetBondTypeAsDouble(),
        ]
        for bond in mol.GetBonds()
    ]

    candidates = _remove_cx(_remove_fixed_bond(mol, names, bond_list))
    inner = [bond for bond in candidates if not _is_terminal(bond, neighbours)]
    terminal = [bond for bond in candidates if _is_terminal(bond, neighbours)]
    terminal = [bond for bond in terminal if not _is_terminal(bond, neighbours_h)]
    terminal = [bond for bond in terminal if not _is_methyl(bond, neighbours_h)]
    inner = [bond for bond in inner if not _is_trihalomethyl(bond, neighbours)]
    inner = _remove_three_membered_ring(mol, names, inner)

    inner_bonds = [bond[0:2] for bond in inner]
    terminal_bonds = [bond[0:2] for bond in terminal]
    bonds = inner_bonds + terminal_bonds
    dihedrals = [
        _compose_dihedral(mol, names, neighbours, bonds, bond) for bond in inner_bonds
    ] + [
        _compose_dihedral(mol_h, names_h, neighbours_h, bonds, bond)
        for bond in terminal_bonds
    ]
    return dihedrals, bonds


def _neighbour_map(mol: Any, names: list[str]) -> dict[str, list[str]]:
    """Bonded atoms of every atom, in atom order (a bond-matrix column)."""
    neighbours: dict[str, list[str]] = {name: [] for name in names}
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 0.0:
            continue
        first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        neighbours[names[first]].append(names[second])
        neighbours[names[second]].append(names[first])
    order = {name: position for position, name in enumerate(names)}
    return {
        name: sorted(atoms, key=order.__getitem__) for name, atoms in neighbours.items()
    }


def _dedupe(bonds: list[Bond]) -> list[Bond]:
    unique: list[Bond] = []
    for bond in bonds:
        if bond not in unique:
            unique.append(bond)
    return unique


def _ring_names(mol: Any, names: list[str]) -> list[str]:
    """Ring atoms, repeated once per ring, as ``remove_FixedBond`` counts them."""
    return [names[atom] for ring in mol.GetRingInfo().AtomRings() for atom in ring]


def _remove_fixed_bond(mol: Any, names: list[str], bond_list: list[Bond]) -> list[Bond]:
    """``remove_FixedBond``: keep bonds that can rotate."""
    single = [bond for bond in bond_list if bond[2] <= 1.0]
    aromatic = [bond for bond in bond_list if bond[2] == 1.5]
    aromatic_hetero = [
        bond for bond in aromatic if bond[0][0] != "C" or bond[1][0] != "C"
    ]
    double = [bond for bond in bond_list if bond[2] == 2.0]
    double_hetero = [bond for bond in double if bond[0][0] != "C" and bond[1][0] != "C"]

    ring = _ring_names(mol, names)
    ring_n = [name for name in ring if name[0] == "N"]
    ring_c = [name for name in ring if name[0] == "C"]
    imine = [bond for bond in double if bond[0][0] == "N" and bond[1][0] == "C"] + [
        bond for bond in double if bond[0][0] == "C" and bond[1][0] == "N"
    ]
    imine_open = [
        bond
        for bond in imine
        if sum(name in bond for name in ring_n) + sum(name in bond for name in ring_c)
        <= 1
    ]

    carbonyl_like = [
        bond for bond in double if bond[0][0] in {"O", "N"} or bond[1][0] in {"O", "N"}
    ]
    carbonyl_carbons = [
        name for bond in carbonyl_like for name in bond[0:2] if name[0] == "C"
    ]
    aromatic_next_to_carbonyl = [
        bond for bond in aromatic if sum(name in bond for name in carbonyl_carbons) >= 1
    ]

    return _dedupe(
        single
        + aromatic_hetero
        + imine_open
        + aromatic_next_to_carbonyl
        + double_hetero
    )


def _remove_cx(bonds: list[Bond]) -> list[Bond]:
    """``remove_CX``: compares one character, so only F and I ever match."""
    blocked = {"F", "Cl", "Br", "I"}
    return [
        bond
        for bond in bonds
        if bond[0][0] not in blocked and bond[1][0] not in blocked
    ]


def _is_terminal(bond: Bond, neighbours: dict[str, list[str]]) -> bool:
    """Either atom has a single neighbour in this bond matrix."""
    return len(neighbours[bond[0]]) == 1 or len(neighbours[bond[1]]) == 1


def _is_methyl(bond: Bond, neighbours_h: dict[str, list[str]]) -> bool:
    """``remove_CH3``: an end with three H on a C (or Cl) or Si atom."""
    for end in (bond[0], bond[1]):
        hydrogens = [name for name in neighbours_h[end] if name[0:2] == "H "]
        if len(hydrogens) == 3 and (end[0] == "C" or end[0:2] == "Si"):
            return True
    return False


def _is_trihalomethyl(bond: Bond, neighbours: dict[str, list[str]]) -> bool:
    """``remove_terminalCCX3``: an end bonded to exactly three like halogens."""
    for end in (bond[0], bond[1]):
        halogens = [
            name[0:2] for name in neighbours[end] if name[0:2] in _HALOGEN_PREFIXES
        ]
        if halogens in _TRIHALO_NEIGHBOURS:
            return True
    return False


def _remove_three_membered_ring(
    mol: Any, names: list[str], bonds: list[Bond]
) -> list[Bond]:
    """``remove_3MemRing``: drop bonds with both atoms in a 3-membered ring."""
    in_ring = {
        names[atom.GetIdx()]: 1 if atom.IsInRingSize(3) else -1
        for atom in mol.GetAtoms()
    }
    return [bond for bond in bonds if in_ring[bond[0]] + in_ring[bond[1]] < 1]


def _same_ring_points(mol: Any, names: list[str], chosen: str) -> list[int]:
    """``atom_within_same_ring``: +2 per ring an atom shares with ``chosen``."""
    rings = [[names[atom] for atom in ring] for ring in mol.GetRingInfo().AtomRings()]
    shared = [ring for ring in rings if chosen in ring]
    return [2 * sum(name in ring for ring in shared) for name in names]


def _compose_dihedral(
    mol: Any,
    names: list[str],
    neighbours: dict[str, list[str]],
    bonds: list[list[str]],
    bond: list[str],
) -> list[str]:
    """``compose_dihedral``: extend a bond to the best-scoring neighbours.

    Neighbours that lie on another selected bond are preferred; otherwise
    any neighbour is allowed. Heteroatoms score 2, carbon 1, hydrogen 0,
    plus 2 per ring shared with the bond atom; the first highest wins.
    """
    bonded = {name for pair in bonds for name in pair}
    atom_points = [
        2 if atom.GetAtomicNum() not in {1, 6} else int(atom.GetAtomicNum() == 6)
        for atom in mol.GetAtoms()
    ]

    ends: list[str] = []
    for this, other in ((bond[0], bond[1]), (bond[1], bond[0])):
        options = [n for n in neighbours[this] if n in bonded and n != other]
        if not options:
            options = [n for n in neighbours[this] if n != other]
        ring_points = _same_ring_points(mol, names, this)
        points = {
            name: atom_points[position] + ring_points[position]
            for position, name in enumerate(names)
        }
        scores = [points[name] for name in options]
        ends.append(options[scores.index(max(scores))])
    return [ends[0], bond[0], bond[1], ends[1]]


# ---------------------------------------------------------------------------
# cal_dihedral_v2
# ---------------------------------------------------------------------------


def _read_molecule(sd_record: str) -> tuple[Any, Any]:
    """The first record, without and with hydrogens, as SDMolSupplier reads it."""
    lines = sd_record.splitlines()
    if _END_TAG not in lines:
        raise ConformerBackendError(
            "CONFPASS candidate is not an SD record: no 'M  END' line",
            code="sdf_parse_failed",
        )
    block = "\n".join(lines[: lines.index(_END_TAG) + 1]) + "\n"
    mol = Chem.MolFromMolBlock(block)
    mol_h = Chem.MolFromMolBlock(block, removeHs=False)
    if mol is None or mol_h is None:
        raise ConformerBackendError(
            "RDKit could not read the first CONFPASS candidate's molblock",
            code="sdf_parse_failed",
        )
    # Newer RDKit keeps a hydrogen that defines double-bond stereo perceived
    # from the 3D geometry (an imine N-H); RDKit 2020.09.5, which CONFPASS
    # was validated on, drops it. The heavy-atom naming assumes every
    # hydrogen is gone, so drop those too.
    if any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        params = Chem.RemoveHsParameters()
        params.removeDefiningBondStereo = True
        mol = Chem.RemoveHs(mol, params)
    return mol, mol_h


def _coordinates(sd_record: str) -> list[list[float]]:
    """Atom-block coordinates read the way ``dihedral_df`` reads them."""
    lines = sd_record.splitlines()
    try:
        atom_count = int(lines[3][:3])
        return [
            [float(value) for value in line[0:35].split()[0:3]]
            for line in lines[4 : 4 + atom_count]
        ]
    except (IndexError, ValueError) as exc:
        raise ConformerBackendError(
            "CONFPASS candidate holds an unreadable atom block",
            code="sdf_parse_failed",
        ) from exc


def _dihedral_angle(
    a1: np.ndarray[Any, Any],
    a2: np.ndarray[Any, Any],
    b1: np.ndarray[Any, Any],
    b2: np.ndarray[Any, Any],
) -> float:
    """``dihedral_angle`` in degrees, (-180, 180]."""
    vec_a = a2 - a1
    vec_b = b1 - a2
    vec_c = b2 - b1
    m = np.cross(vec_a, vec_b)
    n = np.cross(vec_b, vec_c)
    psi = np.sign(np.dot(vec_a, n)) * np.arccos(
        np.dot(n, m) / (np.sqrt(np.dot(n, n)) * np.sqrt(np.dot(m, m)))
    )
    return float(psi * 180 / np.pi)


def _dihedral_values(
    dihedrals: list[list[str]], coordinates: list[list[float]]
) -> list[float]:
    """``dihedral_descriptor`` for one conformer."""
    xyz = np.array(coordinates, dtype=float)
    values = []
    for dihedral in dihedrals:
        rows = [int(name.split(" ")[1]) - 1 for name in dihedral]
        values.append(_dihedral_angle(*(xyz[row] for row in rows)))
    return values


# ---------------------------------------------------------------------------
# dihedral_parameter_v2 and correcting_dihedral_v1
# ---------------------------------------------------------------------------


def _fixed_columns(angles: np.ndarray[Any, Any], columns: list[str]) -> list[str]:
    """``list_fixed_bond``: dihedrals that barely move across the ensemble."""
    fixed = []
    for position, column in enumerate(columns):
        values = angles[:, position]
        absolute = np.absolute(values)
        spread = np.amax(values) - np.amin(values)
        spread_abs = np.amax(absolute) - np.amin(absolute)
        still = np.std(values) < 2.4 and spread < 8.7
        # Flipping between about -180 and +180 is still one position.
        still_across_wrap = (
            spread > 359.2 and np.std(absolute) < 1.5 and spread_abs < 6.6
        )
        if still or still_across_wrap:
            fixed.append(column)
    return fixed


def _correct_by_gap(values: np.ndarray[Any, Any]) -> list[float]:
    """``correct_dih`` then ``shift_axis_df`` for one dihedral column.

    Angles are moved to 0-360; if the low end (below 100) has a gap of at
    least 15 degrees, everything left of it wraps by 360. Values are then
    made relative to the first conformer and shifted to start at zero.
    """
    shifted = [float(value) + 180 for value in values]
    low = sorted([value for value in shifted if value < 100] + [100])
    gaps = [after - before for before, after in zip(low[:-1], low[1:])]
    if gaps and max(gaps) >= 15:
        edge = low[gaps.index(max(gaps))]
        shifted = [value + 360 if value <= edge else value for value in shifted]
    relative = [value - shifted[0] for value in shifted]
    lowest = min(relative)
    return [value - lowest for value in relative]


# ---------------------------------------------------------------------------
# clustering_dih_v7
# ---------------------------------------------------------------------------


def _ward_partitions(
    descriptor: np.ndarray[Any, Any],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Ward partitions into 1..N clusters, as N separate sklearn fits give.

    ``AgglomerativeClustering(n_clusters=k)`` without connectivity always
    builds the full ward tree and undoes its last ``k - 1`` merges, so
    replaying the merges of one full tree yields every partition.
    """
    count = len(descriptor)
    model = _cluster.AgglomerativeClustering(
        n_clusters=None, distance_threshold=0.0, compute_full_tree=True
    ).fit(descriptor)
    members: dict[int, list[int]] = {leaf: [leaf] for leaf in range(count)}
    partitions = [_partition(members.values())]
    for step, (left, right) in enumerate(model.children_):
        members[count + step] = members.pop(int(left)) + members.pop(int(right))
        partitions.append(_partition(members.values()))
    partitions.reverse()
    return tuple(partitions)


def _partition(clusters: Any) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted(tuple(sorted(cluster)) for cluster in clusters))


__all__ = [
    "DEFAULT_METHOD",
    "DEFAULT_X",
    "DEFAULT_X_AS",
    "METHODS",
    "NO_VARIABLE_DIHEDRAL",
    "ConfPassAnalysis",
    "PortedConfPass",
    "analyse",
    "dihedral_descriptor",
    "isolate_dihedrals",
]
