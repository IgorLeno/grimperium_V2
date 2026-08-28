"""The structure route that stays available when CREST is disabled.

Turning the conformer search off must never leave MOPAC without a
geometry, so this module embeds one starting structure with RDKit and
returns it as a regular single-conformer ensemble. Downstream code sees
the same types either way; only the recorded provenance differs, and it
says plainly that no conformational sampling happened.

RDKit is imported through :func:`importlib.import_module` because its
generated bindings carry no type information; this is the same access
pattern used by :mod:`semi_imperium.molecules.validation`.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from typing import Any

from semi_imperium.conformers.backends import (
    ConformerBackendError,
    ConformerRequest,
)
from semi_imperium.conformers.ensemble import (
    UNKNOWN_PROGRAM_VERSION,
    Conformer,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerSearchProvenance,
)
from semi_imperium.domain.configuration import ConformerSearchSettings
from semi_imperium.domain.enums import ConformerSource

rdkit: Any = import_module("rdkit")
Chem: Any = import_module("rdkit.Chem")
AllChem: Any = import_module("rdkit.Chem.AllChem")


@dataclass(frozen=True)
class RDKitInitialStructure:
    """Embeds one 3D structure to feed MOPAC without a conformer search."""

    program: str = "rdkit"
    random_seed: int = 0xF00D
    optimize: bool = True
    max_optimization_steps: int = 200

    def build(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        """Return a single-conformer ensemble for ``request``.

        Args:
            request: The molecule to embed.
            settings: The effective search settings, recorded as-is so
                the provenance shows the search was configured off.

        Raises:
            ConformerBackendError: If RDKit cannot parse the SMILES or
                cannot produce a 3D structure for it.
        """
        molecule = _embed(request.smiles, seed=self.random_seed)
        steps: tuple[str, ...] = ("rdkit.EmbedMolecule",)
        if self.optimize and _uff_optimize(molecule, self.max_optimization_steps):
            steps = (*steps, "rdkit.UFFOptimizeMolecule")
        conformer = Conformer(
            index=0,
            geometry=_geometry_of(molecule),
            energy_kcal_mol=None,
            label="initial",
        )
        provenance = ConformerSearchProvenance(
            source=ConformerSource.RDKIT_INITIAL_3D,
            program=self.program,
            program_version=_rdkit_version(),
            settings=settings,
            run_id=request.run_id,
            command=steps,
        )
        return ConformerEnsemble(conformers=(conformer,), provenance=provenance)


def _embed(smiles: str, *, seed: int) -> Any:
    """Parse ``smiles``, add hydrogens and embed one 3D conformer."""
    parsed = Chem.MolFromSmiles(smiles)
    if parsed is None:
        raise ConformerBackendError(
            f"RDKit could not parse SMILES {smiles!r} for the initial-3D route",
            code="initial_structure_parse_failed",
        )
    molecule = Chem.AddHs(parsed)
    status = AllChem.EmbedMolecule(molecule, randomSeed=seed)
    if status != 0:
        status = AllChem.EmbedMolecule(
            molecule,
            useRandomCoords=True,
            randomSeed=seed,
        )
    if status != 0:
        raise ConformerBackendError(
            f"RDKit could not embed an initial 3D structure for {smiles!r}",
            code="initial_structure_embed_failed",
        )
    return molecule


def _uff_optimize(molecule: Any, max_steps: int) -> bool:
    """Relax the embedded structure, reporting whether UFF could run.

    A molecule outside UFF's parameter set is not an error here: the raw
    embedded geometry is still a valid MOPAC starting point. The caller
    records which steps actually ran, so the skip is never invisible.
    """
    if not AllChem.UFFHasAllMoleculeParams(molecule):
        return False
    AllChem.UFFOptimizeMolecule(molecule, maxIters=max_steps)
    return True


def _geometry_of(molecule: Any) -> ConformerGeometry:
    """Read the embedded coordinates, preserving RDKit's atom order."""
    conformer = molecule.GetConformer()
    elements: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        elements.append(str(atom.GetSymbol()))
        coordinates.append((float(position.x), float(position.y), float(position.z)))
    return ConformerGeometry(
        elements=tuple(elements),
        coordinates=tuple(coordinates),
    )


def _rdkit_version() -> str:
    """Report the RDKit version, or say explicitly that it is unknown."""
    version = str(getattr(rdkit, "__version__", "")).strip()
    return version or UNKNOWN_PROGRAM_VERSION


__all__ = ["RDKitInitialStructure"]
