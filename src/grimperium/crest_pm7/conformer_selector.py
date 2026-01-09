"""Conformer selection logic for CREST-PM7 Pipeline.

Determines optimal number of conformers based on molecular flexibility.
"""

import logging
from typing import Optional

from .config import PM7Config

LOG = logging.getLogger("grimperium.crest_pm7.conformer_selector")


def classify_flexibility(nrotbonds: int, config: PM7Config) -> str:
    """Classify molecule flexibility based on rotatable bonds.

    Args:
        nrotbonds: Number of rotatable bonds
        config: Pipeline configuration

    Returns:
        Flexibility class: "rigid", "medium", or "flexible"
    """
    if nrotbonds <= config.nrotbonds_threshold_rigid_to_medium:
        return "rigid"
    elif nrotbonds <= config.nrotbonds_threshold_medium_to_flexible:
        return "medium"
    else:
        return "flexible"


def get_num_conformers(
    nrotbonds: int,
    config: PM7Config,
    max_available: Optional[int] = None,
) -> tuple[int, str]:
    """Determine optimal number of conformers to process.

    Based on molecular flexibility:
    - Rigid (0-1 rotbonds): 1 conformer
    - Medium (2-4 rotbonds): 3 conformers
    - Flexible (5+ rotbonds): 5 conformers

    Args:
        nrotbonds: Number of rotatable bonds
        config: Pipeline configuration
        max_available: Maximum conformers available from CREST

    Returns:
        Tuple of (num_conformers, decision_reason)
    """
    flexibility = classify_flexibility(nrotbonds, config)

    if flexibility == "rigid":
        target = 1
        reason = f"rigid molecule (nrotbonds={nrotbonds}<=1)"
    elif flexibility == "medium":
        target = 3
        reason = f"medium flexibility (nrotbonds={nrotbonds}, 2-4)"
    else:
        target = 5
        reason = f"flexible molecule (nrotbonds={nrotbonds}>=5)"

    # Respect max_conformers from config
    target = min(target, config.max_conformers)

    # Respect available conformers
    if max_available is not None and max_available < target:
        old_target = target
        target = max_available
        reason += f", limited by available ({max_available} < {old_target})"

    LOG.debug(f"Selected {target} conformers: {reason}")
    return target, reason


def calculate_delta_e(
    energies: list[float],
    reference_index: int = 0,
) -> dict[str, Optional[float]]:
    """Calculate energy differences between conformers.

    Args:
        energies: List of conformer energies (kcal/mol)
        reference_index: Index of reference conformer (usually lowest energy)

    Returns:
        Dictionary with delta_e_12, delta_e_13, delta_e_15
    """
    if not energies:
        return {"delta_e_12": None, "delta_e_13": None, "delta_e_15": None}

    # Sort energies (lowest first)
    sorted_e = sorted(energies)
    n = len(sorted_e)

    result = {
        "delta_e_12": None,
        "delta_e_13": None,
        "delta_e_15": None,
    }

    e1 = sorted_e[0]

    if n >= 2:
        result["delta_e_12"] = sorted_e[1] - e1
    if n >= 3:
        result["delta_e_13"] = sorted_e[2] - e1
    if n >= 5:
        result["delta_e_15"] = sorted_e[4] - e1

    return result


def analyze_conformer_distribution(
    energies: list[float],
    energy_window: float = 6.0,
) -> dict:
    """Analyze conformer energy distribution.

    Args:
        energies: List of conformer energies
        energy_window: Energy window for counting (kcal/mol)

    Returns:
        Dictionary with distribution statistics
    """
    if not energies:
        return {
            "n_conformers": 0,
            "n_within_window": 0,
            "energy_spread": None,
            "mean_energy": None,
        }

    sorted_e = sorted(energies)
    e_min = sorted_e[0]

    within_window = sum(1 for e in sorted_e if (e - e_min) <= energy_window)

    return {
        "n_conformers": len(energies),
        "n_within_window": within_window,
        "energy_spread": sorted_e[-1] - sorted_e[0] if len(sorted_e) > 1 else 0.0,
        "mean_energy": sum(energies) / len(energies),
        "min_energy": e_min,
        "max_energy": sorted_e[-1],
    }
