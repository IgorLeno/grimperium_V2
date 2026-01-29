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
    """
    if func is None:
        return None
    try:
        return float(func(mol))
    except Exception as exc:
        logger.debug("Descriptor failed: %s", exc)
        return None


def extract_all_rdkit_descriptors(smiles: str) -> dict[str, float]:
    """Extract RDKit descriptors useful for delta-learning.

    Args:
        smiles: SMILES string

    Returns:
        Dictionary of descriptor values (empty if SMILES invalid)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors
    except ImportError:
        logger.warning("RDKit not available, returning empty descriptor set")
        return {}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Invalid SMILES for RDKit descriptors: %s", smiles)
        return {}

    descriptors: dict[str, float | None] = {
        "mol_wt": _safe_descriptor(getattr(Descriptors, "MolWt", None), mol),
        "hbd": _safe_descriptor(getattr(Descriptors, "NumHDonors", None), mol),
        "hba": _safe_descriptor(getattr(Descriptors, "NumHAcceptors", None), mol),
        "logp": _safe_descriptor(getattr(Crippen, "MolLogP", None), mol),
        "ring_count": _safe_descriptor(getattr(Descriptors, "RingCount", None), mol),
        "sat_ring_count": _safe_descriptor(
            getattr(Descriptors, "NumSaturatedRings", None), mol
        ),
        "arom_ring_count": _safe_descriptor(
            getattr(Descriptors, "NumAromaticRings", None), mol
        ),
        "num_heteroatoms": _safe_descriptor(
            getattr(Descriptors, "NumHeteroatoms", None), mol
        ),
        "fraction_csp3": _safe_descriptor(
            getattr(Descriptors, "FractionCSP3", None), mol
        ),
        "labute_asa": _safe_descriptor(getattr(Descriptors, "LabuteASA", None), mol),
        "bertz_ct": _safe_descriptor(getattr(Descriptors, "BertzCT", None), mol),
        "tpsa": _safe_descriptor(getattr(Descriptors, "TPSA", None), mol),
        "nrotbonds": _safe_descriptor(
            getattr(Descriptors, "NumRotatableBonds", None), mol
        ),
    }

    return {key: value for key, value in descriptors.items() if value is not None}
