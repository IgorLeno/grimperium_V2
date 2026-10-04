"""Build the SDF inputs for the CONFPASS golden files from stored CREST runs.

Each case is a ``crest_conformers.xyz`` from ``runs/`` read with the
production parser (:func:`parse_crest_ensemble`) and written with the
production adapter (:func:`build_confpass_candidates`), so the golden
files pin CONFPASS's behaviour on exactly the SDF the Semi-Imperium
backend will hand it. Topology is perceived once from the first
conformer (``rdDetermineBonds``) and checked against every other one.

The SDF keeps CREST's own order, including the small energy inversions
CREST sometimes writes; they are counted in ``inputs.json`` because
CONFPASS PART 1 never reads energies and assumes file order is energy
order.

Usage (project environment)::

    poetry run python scripts/confpass_golden/build_inputs.py \\
        --runs runs --output tests/fixtures/confpass_golden
"""

from __future__ import annotations

import argparse
import hashlib
import json
from importlib import import_module
from pathlib import Path
from typing import Any

from semi_imperium.conformers import (
    ConformerEnsemble,
    ConformerSearchProvenance,
    MoleculeTopology,
    build_confpass_candidates,
    parse_crest_ensemble,
)
from semi_imperium.domain import ConformerSearchSettings, ConformerSource

Chem: Any = import_module("rdkit.Chem")
rdDetermineBonds: Any = import_module("rdkit.Chem.rdDetermineBonds")
rdBase: Any = import_module("rdkit.rdBase")

#: (case id, run directory). Chosen for 2-121 conformers and varied
#: chemistry: carboxylic acids, alcohols, methyl esters, a polyol and a
#: branched triester; the two smallest are edge cases.
CASES: tuple[tuple[str, str], ...] = (
    ("propanoic_acid", "cbt_002"),
    ("methyl_acetate", "cbt_020"),
    ("propanol", "cbt_014"),
    ("butanoic_acid", "cbt_003"),
    ("pentanoic_acid", "cbt_004"),
    ("butanol", "cbt_015"),
    ("methyl_propanoate", "cbt_021"),
    ("methyl_butanoate", "cbt_022"),
    ("hexanoic_acid", "cbt_005"),
    ("methyl_pentanoate", "cbt_023"),
    ("glycerol", "cbt_013"),
    ("triacetin", "cbt_030"),
    ("heptanoic_acid", "cbt_006"),
)

ENSEMBLE_FILE = "crest_conformers.xyz"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def xyz_block(ensemble: ConformerEnsemble, position: int) -> str:
    geometry = ensemble.conformers[position].geometry
    lines = [str(geometry.atom_count), ""]
    for symbol, (x, y, z) in zip(geometry.elements, geometry.coordinates):
        lines.append(f"{symbol} {x:.8f} {y:.8f} {z:.8f}")
    return "\n".join(lines)


def perceive(ensemble: ConformerEnsemble, position: int, orders: bool) -> Any:
    mol = Chem.MolFromXYZBlock(xyz_block(ensemble, position))
    if orders:
        rdDetermineBonds.DetermineBonds(mol, charge=0)
    else:
        rdDetermineBonds.DetermineConnectivity(mol)
    return mol


def bond_pairs(mol: Any) -> set[tuple[int, int]]:
    return {
        (
            min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
        )
        for b in mol.GetBonds()
    }


def topology_of(ensemble: ConformerEnsemble) -> tuple[MoleculeTopology, str]:
    """Bond orders from the first conformer, connectivity checked on all."""
    mol = perceive(ensemble, 0, orders=True)
    bonds: list[tuple[int, int, int]] = []
    for bond in mol.GetBonds():
        order = bond.GetBondTypeAsDouble()
        if order != int(order):
            raise ValueError(f"non-integer bond order {order} cannot go to SDF")
        bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), int(order)))
    reference = bond_pairs(mol)
    for position in range(1, ensemble.size):
        if bond_pairs(perceive(ensemble, position, orders=False)) != reference:
            raise ValueError(
                f"conformer {position} has a different connectivity than conformer 0"
            )
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    return MoleculeTopology(atom_count=mol.GetNumAtoms(), bonds=tuple(bonds)), smiles


def energy_inversions(ensemble: ConformerEnsemble) -> tuple[int, float]:
    energies = [conformer.require_energy() for conformer in ensemble.conformers]
    drops = [a - b for a, b in zip(energies, energies[1:]) if b < a]
    return len(drops), max(drops, default=0.0)


def build_case(case_id: str, run_dir: Path, output: Path) -> dict[str, Any]:
    source = (run_dir / ENSEMBLE_FILE).read_text()
    conformers = parse_crest_ensemble(source)
    ensemble = ConformerEnsemble(
        conformers=conformers,
        provenance=ConformerSearchProvenance(
            source=ConformerSource.CREST,
            program="crest",
            program_version="unknown",
            settings=ConformerSearchSettings(),
            run_id=run_dir.name,
        ),
    )
    topology, smiles = topology_of(ensemble)
    candidates = build_confpass_candidates(ensemble, topology, molecule_id=case_id)
    sdf = "\n".join(candidate.sd_record for candidate in candidates) + "\n"
    (output / f"{case_id}.sdf").write_text(sdf)
    inversions, largest = energy_inversions(ensemble)
    return {
        "case": case_id,
        "run": run_dir.name,
        "smiles": smiles,
        "atoms": topology.atom_count,
        "conformers": ensemble.size,
        "energy_inversions": inversions,
        "largest_inversion_kcal_mol": round(largest, 6),
        "source_sha256": sha256_text(source),
        "sdf": f"{case_id}.sdf",
        "sdf_sha256": sha256_text(sdf),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--runs", type=Path, default=Path("runs"))
    parser.add_argument(
        "--output", type=Path, default=Path("tests/fixtures/confpass_golden")
    )
    args = parser.parse_args(argv)

    rdBase.DisableLog("rdApp.*")
    output: Path = args.output / "inputs"
    output.mkdir(parents=True, exist_ok=True)
    cases = [build_case(case, args.runs / run, output) for case, run in CASES]
    manifest = {
        "generator": "scripts/confpass_golden/build_inputs.py",
        "rdkit": rdBase.rdkitVersion,
        "cases": cases,
    }
    (args.output / "inputs.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for case in cases:
        print(
            f"{case['case']:<18} {case['conformers']:>4} conformers "
            f"{case['energy_inversions']:>3} inversions  {case['smiles']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
