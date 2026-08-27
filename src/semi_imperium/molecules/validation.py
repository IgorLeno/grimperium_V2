"""Early local validation for Semi-Imperium computational inputs."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from typing import Any

Chem: Any = import_module("rdkit.Chem")
AllChem: Any = import_module("rdkit.Chem.AllChem")

SUPPORTED_HAMILTONIANS = frozenset({"AM1", "PM3", "PM7"})
SUPPORTED_MULTIPLICITIES = frozenset(range(1, 7))


class MoleculeValidationError(ValueError):
    """An actionable failure detected before CREST or MOPAC execution."""

    def __init__(self, message: str, *, code: str) -> None:
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class ComputationalPath:
    """Relevant near-term requirements of the requested calculation route."""

    hamiltonians: tuple[str, ...] = ("AM1", "PM3", "PM7")
    crest_enabled: bool = True
    require_initial_3d: bool = True

    def __post_init__(self) -> None:
        normalized = tuple(
            dict.fromkeys(item.strip().upper() for item in self.hamiltonians)
        )
        object.__setattr__(self, "hamiltonians", normalized)
        if not normalized:
            raise MoleculeValidationError(
                "Select at least one MOPAC Hamiltonian (AM1, PM3 or PM7)",
                code="no_hamiltonian",
            )
        unsupported = sorted(set(normalized) - SUPPORTED_HAMILTONIANS)
        if unsupported:
            raise MoleculeValidationError(
                f"Unsupported MOPAC Hamiltonian(s): {', '.join(unsupported)}. "
                "Choose AM1, PM3 or PM7.",
                code="unsupported_hamiltonian",
            )


@dataclass(frozen=True)
class ValidatedStructure:
    """Locally parsed structure safe to pass to the configured runner path."""

    canonical_smiles: str
    formal_charge: int
    atom_count: int
    inchi: str | None
    inchikey: str | None


class MoleculeValidator:
    """Validate parsing and known RDKit-to-CREST/MOPAC preconditions."""

    def canonicalize_smiles(self, smiles: str) -> str:
        """Canonicalize for structural comparison without route filtering.

        This deliberately performs only parsing and sanitization. Candidate
        disambiguation must not silently discard a distinct structure merely
        because that structure later fails the selected computational path.
        """
        mol = self._parse_smiles(smiles)
        return str(Chem.MolToSmiles(mol, canonical=True))

    def validate_smiles(
        self,
        smiles: str,
        *,
        path: ComputationalPath | None = None,
        expected_charge: int | None = None,
        multiplicity: int = 1,
    ) -> ValidatedStructure:
        """Parse, sanitize, canonicalize and preflight the requested path."""
        path = path or ComputationalPath()
        if multiplicity not in SUPPORTED_MULTIPLICITIES:
            raise MoleculeValidationError(
                f"Spin multiplicity {multiplicity} is unsupported by the current "
                "MOPAC path; choose a value from 1 through 6",
                code="unsupported_multiplicity",
            )

        mol = self._parse_smiles(smiles)
        if mol.GetNumAtoms() == 0:
            raise MoleculeValidationError(
                "SMILES does not contain any atoms",
                code="empty_structure",
            )
        if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
            raise MoleculeValidationError(
                "Wildcard/query atoms are not supported by the CREST/MOPAC path; "
                "enter an explicit molecular structure",
                code="query_atom",
            )
        if len(Chem.GetMolFrags(mol)) != 1:
            raise MoleculeValidationError(
                "Disconnected fragments (for example salts or mixtures) are not "
                "supported as one enthalpy calculation; enter one connected species",
                code="disconnected_structure",
            )

        formal_charge = int(Chem.GetFormalCharge(mol))
        if expected_charge is not None and expected_charge != formal_charge:
            raise MoleculeValidationError(
                f"Requested charge {expected_charge} does not match the SMILES "
                f"formal charge {formal_charge}; correct the charge or structure",
                code="charge_mismatch",
            )

        with_hydrogens = Chem.AddHs(mol)
        electron_count = (
            sum(atom.GetAtomicNum() for atom in with_hydrogens.GetAtoms())
            - formal_charge
        )
        if electron_count % 2 != (multiplicity - 1) % 2:
            suggested_parity = "even" if electron_count % 2 else "odd"
            raise MoleculeValidationError(
                f"Spin multiplicity {multiplicity} is inconsistent with the "
                f"structure's electron count; choose an {suggested_parity} "
                "multiplicity for this charge state",
                code="multiplicity_electron_mismatch",
            )

        if path.require_initial_3d:
            self._validate_initial_geometry(with_hydrogens, path)

        canonical_smiles = str(Chem.MolToSmiles(mol, canonical=True))
        inchi: str | None = None
        inchikey: str | None = None
        try:
            inchi = Chem.MolToInchi(mol) or None
            inchikey = Chem.MolToInchiKey(mol) or None
        except Exception:  # pragma: no cover - depends on RDKit InChI build
            inchi = None
            inchikey = None
        return ValidatedStructure(
            canonical_smiles=canonical_smiles,
            formal_charge=formal_charge,
            atom_count=mol.GetNumAtoms(),
            inchi=inchi,
            inchikey=inchikey,
        )

    @staticmethod
    def _parse_smiles(smiles: str) -> Any:
        if not isinstance(smiles, str) or not smiles.strip():
            raise MoleculeValidationError(
                "SMILES is empty; enter a molecular SMILES string",
                code="empty_smiles",
            )
        clean_smiles = smiles.strip()
        try:
            mol = Chem.MolFromSmiles(clean_smiles, sanitize=True)
        except Exception as exc:
            raise MoleculeValidationError(
                f"Malformed SMILES {clean_smiles!r}: RDKit could not sanitize it",
                code="malformed_smiles",
            ) from exc
        if mol is None:
            raise MoleculeValidationError(
                f"Malformed SMILES {clean_smiles!r}: RDKit could not parse it",
                code="malformed_smiles",
            )
        return mol

    @staticmethod
    def _validate_initial_geometry(mol: Any, path: ComputationalPath) -> None:
        """Exercise the same RDKit 3D/UFF prerequisite used by the runner."""
        if not AllChem.UFFHasAllMoleculeParams(mol):
            route = "CREST/MOPAC" if path.crest_enabled else "MOPAC"
            raise MoleculeValidationError(
                f"The {route} route cannot create its required initial 3D geometry "
                "because RDKit UFF lacks parameters for one or more atoms",
                code="unsupported_3d_parameters",
            )
        embedded = Chem.Mol(mol)
        result = AllChem.EmbedMolecule(embedded, randomSeed=0xF00D)
        if result != 0:
            result = AllChem.EmbedMolecule(
                embedded,
                useRandomCoords=True,
                randomSeed=0xF00D,
            )
        if result != 0:
            raise MoleculeValidationError(
                "RDKit could not generate the initial 3D structure required by the "
                "CREST/MOPAC route; provide a supported structure",
                code="initial_3d_failed",
            )


__all__ = [
    "ComputationalPath",
    "MoleculeValidationError",
    "MoleculeValidator",
    "SUPPORTED_HAMILTONIANS",
    "ValidatedStructure",
]
