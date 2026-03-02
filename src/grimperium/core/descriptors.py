"""RDKit descriptor extraction for delta-learning features."""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any

logger = logging.getLogger(__name__)


def _safe_descriptor(
    func: Callable[[Any], float] | None,
    mol: Any,
) -> float | None:
    """Safely evaluate an RDKit descriptor callable.

    Args:
        func: RDKit descriptor function or None
        mol: RDKit molecule

    Returns:
        Descriptor value or None if unavailable

    Example:
        >>> from rdkit import Chem
        >>> from rdkit.Chem import Descriptors
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> result = _safe_descriptor(Descriptors.MolWt, mol)
        >>> 46.0 < result < 47.0  # Ethanol MW ~46.07
        True
    """
    if func is None:
        return None
    try:
        return float(func(mol))
    except Exception as exc:
        logger.debug("Descriptor failed: %s", exc)
        return None


def _count_atoms_by_symbol(mol: Any, symbol: str) -> int:
    """Count atoms of a given element in a molecule (including implicit H).

    Args:
        mol: RDKit molecule (should have explicit Hs added for H counts)
        symbol: Atomic symbol (e.g. "C", "H", "O", "N")

    Returns:
        Count of atoms with the given symbol

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.AddHs(Chem.MolFromSmiles("O"))
        >>> _count_atoms_by_symbol(mol, "H")
        2
    """
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == symbol)


def _count_bonds_by_type(mol: Any) -> dict[str, int]:
    """Count bonds by type in a molecule.

    Args:
        mol: RDKit molecule

    Returns:
        Dict with keys: single, double, triple, aromatic

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles("C=O")
        >>> counts = _count_bonds_by_type(mol)
        >>> counts == {"single": 0, "double": 1, "triple": 0, "aromatic": 0}
        True
    """
    from rdkit.Chem import rdchem

    counts = {"single": 0, "double": 0, "triple": 0, "aromatic": 0}
    bond_type_map = {
        rdchem.BondType.SINGLE: "single",
        rdchem.BondType.DOUBLE: "double",
        rdchem.BondType.TRIPLE: "triple",
        rdchem.BondType.AROMATIC: "aromatic",
    }
    for bond in mol.GetBonds():
        btype = bond_type_map.get(bond.GetBondType())
        if btype:
            counts[btype] += 1
    return counts


def extract_all_rdkit_descriptors(smiles: str) -> dict[str, float | int]:
    """Extract RDKit descriptors matching the CSV schema (15 rdkit_* columns).

    Args:
        smiles: SMILES string

    Returns:
        Dictionary with 15 rdkit_* descriptor values (empty if SMILES invalid)

    Example:
        >>> descriptors = extract_all_rdkit_descriptors("CCO")
        >>> "rdkit_mol_weight" in descriptors
        True
        >>> len(descriptors) == 15
        True
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        logger.warning("RDKit not available, returning empty descriptor set")
        return {}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Invalid SMILES for RDKit descriptors: %s", smiles)
        return {}

    # Molecule with explicit hydrogens for atom counting
    mol_h = Chem.AddHs(mol)

    # Warn if molecule contains atoms outside {C, H, O, N} — counts remain CHON-only
    _CHON = {"C", "H", "O", "N"}
    non_chon = {atom.GetSymbol() for atom in mol_h.GetAtoms()} - _CHON
    if non_chon:
        logger.warning(
            "Molecule contains non-CHON atoms: %s; rdkit_nC/nH/nO/nN count only C/H/O/N",
            sorted(non_chon),
        )

    # Bond counts (on molecule without explicit H, to match aromatic perception)
    bond_counts = _count_bonds_by_type(mol)

    descriptors: dict[str, float | int | None] = {
        "rdkit_nrotbonds": _safe_descriptor(
            getattr(Descriptors, "NumRotatableBonds", None), mol
        ),
        "rdkit_tpsa": _safe_descriptor(getattr(Descriptors, "TPSA", None), mol),
        "rdkit_num_rings": _safe_descriptor(
            getattr(Descriptors, "NumAromaticRings", None), mol
        ),
        "rdkit_fsp3": _safe_descriptor(getattr(Descriptors, "FractionCSP3", None), mol),
        "rdkit_mol_weight": _safe_descriptor(getattr(Descriptors, "MolWt", None), mol),
        "rdkit_hbond_donors": _safe_descriptor(
            getattr(Descriptors, "NumHDonors", None), mol
        ),
        "rdkit_hbond_acceptors": _safe_descriptor(
            getattr(Descriptors, "NumHAcceptors", None), mol
        ),
        # Atom counts (on mol with explicit H); only C, H, O, N are tallied
        "rdkit_nC": _count_atoms_by_symbol(mol_h, "C"),
        "rdkit_nH": _count_atoms_by_symbol(mol_h, "H"),
        "rdkit_nO": _count_atoms_by_symbol(mol_h, "O"),
        "rdkit_nN": _count_atoms_by_symbol(mol_h, "N"),
        # Bond counts
        "rdkit_bonds_single": bond_counts["single"],
        "rdkit_bonds_double": bond_counts["double"],
        "rdkit_bonds_triple": bond_counts["triple"],
        "rdkit_bonds_aromatic": bond_counts["aromatic"],
    }

    return {key: value for key, value in descriptors.items() if value is not None}
