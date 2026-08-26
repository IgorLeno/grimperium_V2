"""Molecular identity for Semi-Imperium.

Two calculations are about "the same molecule" only when their canonical
SMILES, charge and multiplicity all match. Display names are deliberately
excluded from the identity: renaming a molecule must never invalidate a
stored result, and two datasets that spell the same species differently
must still collapse onto one identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from typing import Any

from semi_imperium.domain.hashing import stable_digest

MOLECULE_ID_PREFIX = "mol"
MOLECULE_ID_LENGTH = 16


@dataclass(frozen=True)
class MolecularIdentity:
    """Canonical identity of one chemical species in a given charge state.

    Attributes:
        canonical_smiles: RDKit-canonical SMILES. Callers that start from
            arbitrary user input should build the identity through
            :meth:`from_smiles` so canonicalization is not skipped.
        charge: Net molecular charge.
        multiplicity: Spin multiplicity (``2S + 1``), always ``>= 1``.
        inchikey: Optional cross-reference key. It is provenance, not
            identity: it never participates in ``molecule_id`` so that an
            environment without InChI support produces the same ids.
        display_name: Optional human label. Never part of the identity.
    """

    canonical_smiles: str
    charge: int = 0
    multiplicity: int = 1
    inchikey: str | None = None
    display_name: str | None = None

    def __post_init__(self) -> None:
        if not self.canonical_smiles.strip():
            raise ValueError("MolecularIdentity.canonical_smiles must not be empty")
        if self.canonical_smiles != self.canonical_smiles.strip():
            raise ValueError(
                "MolecularIdentity.canonical_smiles must not have surrounding "
                f"whitespace: {self.canonical_smiles!r}"
            )
        if self.multiplicity < 1:
            raise ValueError(
                f"MolecularIdentity.multiplicity must be >= 1, got {self.multiplicity}"
            )
        if self.inchikey is not None and not self.inchikey.strip():
            raise ValueError(
                "MolecularIdentity.inchikey must be a non-empty string or None"
            )

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        *,
        charge: int = 0,
        multiplicity: int = 1,
        display_name: str | None = None,
        with_inchikey: bool = True,
    ) -> MolecularIdentity:
        """Build an identity by canonicalizing ``smiles`` with RDKit.

        Raises:
            RuntimeError: If RDKit is not importable.
            ValueError: If RDKit cannot parse ``smiles``.
        """
        try:
            chem: Any = import_module("rdkit.Chem")
        except ImportError as exc:  # pragma: no cover - optional environment
            raise RuntimeError(
                "RDKit is required to canonicalize SMILES into a MolecularIdentity"
            ) from exc

        mol = chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"RDKit cannot parse SMILES: {smiles}")

        inchikey: str | None = None
        if with_inchikey:
            # The InChIKey is a cross-reference, not part of the identity, and
            # RDKit builds without InChI support are legal. A missing key must
            # therefore never block identity creation.
            try:
                inchikey = str(chem.MolToInchiKey(mol)) or None
            except Exception:  # pragma: no cover - depends on RDKit build
                inchikey = None

        return cls(
            canonical_smiles=str(chem.MolToSmiles(mol)),
            charge=charge,
            multiplicity=multiplicity,
            inchikey=inchikey,
            display_name=display_name,
        )

    @property
    def identity_payload(self) -> dict[str, Any]:
        """Return exactly the fields that define this molecular identity."""
        return {
            "canonical_smiles": self.canonical_smiles,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }

    @property
    def molecule_id(self) -> str:
        """Deterministic identifier derived from the identity fields only."""
        digest = stable_digest(self.identity_payload)
        return f"{MOLECULE_ID_PREFIX}-{digest[:MOLECULE_ID_LENGTH]}"

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives, including the derived id."""
        return {
            "molecule_id": self.molecule_id,
            "canonical_smiles": self.canonical_smiles,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "inchikey": self.inchikey,
            "display_name": self.display_name,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> MolecularIdentity:
        """Deserialize from JSON-compatible primitives.

        ``molecule_id`` is recomputed rather than trusted, so a hand-edited
        or corrupted file cannot silently redirect reuse lookups.
        """
        identity = cls(
            canonical_smiles=str(payload["canonical_smiles"]),
            charge=int(payload["charge"]),
            multiplicity=int(payload["multiplicity"]),
            inchikey=_optional_str(payload.get("inchikey")),
            display_name=_optional_str(payload.get("display_name")),
        )
        stored_id = payload.get("molecule_id")
        if stored_id is not None and str(stored_id) != identity.molecule_id:
            raise ValueError(
                "Stored molecule_id does not match the identity fields: "
                f"{stored_id!r} != {identity.molecule_id!r}"
            )
        return identity


def _optional_str(value: Any) -> str | None:
    """Return ``value`` as a string, mapping absent values to ``None``."""
    if value is None:
        return None
    return str(value)


__all__ = ["MOLECULE_ID_PREFIX", "MolecularIdentity"]
