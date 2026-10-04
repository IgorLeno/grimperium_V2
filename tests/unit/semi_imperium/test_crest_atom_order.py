"""CREST keeps the atom order of its input, so the SMILES topology holds.

``SmilesTopology`` derives bonds from the SMILES in RDKit ``AddHs`` order
and checks them only against the first conformer. That is enough only if
CREST writes every conformer in the order of its input structure. These
tests pin that on stored CREST runs (``tests/fixtures/crest_atom_order``)
and show that a canonical SMILES, re-parsed, can break the order.
"""

from __future__ import annotations

import json
from importlib import import_module
from pathlib import Path
from typing import Any

import pytest

from semi_imperium.conformers import (
    ConformerBackendError,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerRequest,
    ConformerSearchProvenance,
    parse_crest_ensemble,
)
from semi_imperium.conformers.topology import (
    TOPOLOGY_ATOM_ORDER_MISMATCH,
    SmilesTopology,
    require_matching_order,
    topology_from_smiles,
)
from semi_imperium.domain import ConformerSearchSettings, ConformerSource

Chem: Any = import_module("rdkit.Chem")
rdDetermineBonds: Any = import_module("rdkit.Chem.rdDetermineBonds")

FIXTURES = Path(__file__).resolve().parents[2] / "fixtures" / "crest_atom_order"
CASES: dict[str, dict[str, Any]] = {
    case["run"]: case
    for case in json.loads((FIXTURES / "cases.json").read_text())["cases"]
}
REORDERED_BY_CANONICAL = [
    run for run, case in CASES.items() if not case["canonical_keeps_order"]
]


def ensemble_of(run: str) -> ConformerEnsemble:
    conformers = parse_crest_ensemble(
        (FIXTURES / run / "crest_conformers.xyz").read_text()
    )
    return ConformerEnsemble(
        conformers=conformers,
        provenance=ConformerSearchProvenance(
            source=ConformerSource.CREST,
            program="crest",
            program_version="unknown",
            settings=ConformerSearchSettings(),
            run_id=run,
        ),
    )


def single(ensemble: ConformerEnsemble, position: int) -> ConformerEnsemble:
    return ConformerEnsemble(
        conformers=(ensemble.conformers[position],), provenance=ensemble.provenance
    )


def xyz_elements(path: Path) -> tuple[str, ...]:
    lines = path.read_text().splitlines()
    count = int(lines[0].split()[0])
    return tuple(line.split()[0].capitalize() for line in lines[2 : 2 + count])


def connectivity(geometry: ConformerGeometry) -> set[tuple[int, int]]:
    rows = [str(geometry.atom_count), ""] + [
        f"{symbol} {x:.8f} {y:.8f} {z:.8f}"
        for symbol, (x, y, z) in zip(geometry.elements, geometry.coordinates)
    ]
    molecule = Chem.MolFromXYZBlock("\n".join(rows))
    rdDetermineBonds.DetermineConnectivity(molecule)
    return {
        (
            min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
        )
        for bond in molecule.GetBonds()
    }


def input_connectivity(run: str) -> set[tuple[int, int]]:
    lines = (FIXTURES / run / "input.xyz").read_text().splitlines()
    count = int(lines[0].split()[0])
    rows = [line.split() for line in lines[2 : 2 + count]]
    geometry = ConformerGeometry(
        elements=tuple(row[0].capitalize() for row in rows),
        coordinates=tuple((float(r[1]), float(r[2]), float(r[3])) for r in rows),
    )
    return connectivity(geometry)


@pytest.mark.parametrize("run", sorted(CASES))
def test_every_crest_conformer_keeps_the_input_element_order(run: str) -> None:
    ensemble = ensemble_of(run)
    initial = xyz_elements(FIXTURES / run / "rdkit_initial.xyz")

    assert ensemble.size == CASES[run]["conformers"]
    assert xyz_elements(FIXTURES / run / "input.xyz") == initial
    for conformer in ensemble.conformers:
        observed = tuple(s.capitalize() for s in conformer.geometry.elements)
        assert observed == initial


@pytest.mark.parametrize("run", sorted(CASES))
def test_every_crest_conformer_keeps_the_input_connectivity(run: str) -> None:
    reference = input_connectivity(run)

    for conformer in ensemble_of(run).conformers:
        assert connectivity(conformer.geometry) == reference


@pytest.mark.parametrize("run", sorted(CASES))
def test_smiles_topology_accepts_the_real_crest_ensemble(run: str) -> None:
    smiles = CASES[run]["source_order_smiles"]
    ensemble = ensemble_of(run)
    request = ConformerRequest(molecule_id=run, smiles=smiles, run_id=run)

    topology = SmilesTopology()(request, ensemble)

    _, elements = topology_from_smiles(smiles)
    assert topology.atom_count == CASES[run]["atoms"]
    # The provider only checks conformer 0; the order holds for every one.
    for position in range(ensemble.size):
        require_matching_order(
            topology, elements, single(ensemble, position), molecule_id=run
        )


@pytest.mark.parametrize("run", REORDERED_BY_CANONICAL)
def test_a_canonical_smiles_does_not_keep_the_source_order(run: str) -> None:
    smiles = CASES[run]["canonical_smiles"]
    request = ConformerRequest(molecule_id=run, smiles=smiles, run_id=run)

    with pytest.raises(ConformerBackendError) as caught:
        SmilesTopology()(request, ensemble_of(run))

    assert caught.value.code == TOPOLOGY_ATOM_ORDER_MISMATCH


def test_the_fixtures_cover_both_canonical_outcomes() -> None:
    assert REORDERED_BY_CANONICAL
    assert len(REORDERED_BY_CANONICAL) < len(CASES)
