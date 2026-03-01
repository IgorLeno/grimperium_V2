"""MOPAC .aux file descriptor parser.

Extracts electronic descriptors from MOPAC2016 .aux output files.
The .aux format uses Fortran D notation for floats (e.g., +0.577997D+02).
"""

import logging
import re
from pathlib import Path
from typing import Any

import numpy as np

LOG = logging.getLogger("grimperium.crest_pm7.mopac_descriptors")

# Descriptor keys with their default types
_DESCRIPTOR_KEYS: list[str] = [
    "mopac_dipole_debye",
    "mopac_ionization_potential_ev",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
    "mopac_cosmo_area_a2",
    "mopac_cosmo_volume_a3",
    "mopac_gradient_norm",
    "mopac_num_scf_cycles",
    "mopac_point_group",
    "mopac_time_s",
]


def parse_fortran_float(s: str) -> float:
    """Convert Fortran D notation to Python float.

    Args:
        s: String like '+0.577997D+02' or '-0.123456D-03' where D/d
            is the Fortran exponent separator (equivalent to E notation).

    Returns:
        Python float value after replacing the D/d exponent marker with E.

    Example:
        >>> parse_fortran_float('+0.577997D+02')
        57.7997
    """
    return float(s.replace("D", "E").replace("d", "e"))


def _empty_descriptors() -> dict[str, Any]:
    """Return dict with all descriptor keys set to NaN/None defaults.

    Returns:
        Dictionary with all keys from _DESCRIPTOR_KEYS initialized to np.nan,
        except 'mopac_point_group' which is set to None.
    """
    result: dict[str, Any] = {}
    for key in _DESCRIPTOR_KEYS:
        if key == "mopac_point_group":
            result[key] = None
        else:
            result[key] = np.nan
    return result


def _parse_aux_file(aux_content: str) -> dict[str, Any]:
    """Parse MOPAC .aux file content and extract descriptors.

    Args:
        aux_content: Full text content of the .aux file

    Returns:
        Dictionary with extracted descriptor values
    """
    result = _empty_descriptors()

    # Fortran D notation patterns
    fortran_re = r"([+-]?\d+\.\d+[Dd][+-]\d+)"

    # Simple Fortran D extractions
    patterns: list[tuple[str, str]] = [
        (rf"DIPOLE:DEBYE={fortran_re}", "mopac_dipole_debye"),
        (rf"IONIZATION_POTENTIAL:EV={fortran_re}", "mopac_ionization_potential_ev"),
        (rf"AREA:SQUARE ANGSTROMS={fortran_re}", "mopac_cosmo_area_a2"),
        (rf"VOLUME:CUBIC ANGSTROMS={fortran_re}", "mopac_cosmo_volume_a3"),
        (rf"GRADIENT_NORM:KCAL/MOL/ANGSTROM={fortran_re}", "mopac_gradient_norm"),
    ]

    for pattern, key in patterns:
        m = re.search(pattern, aux_content)
        if m:
            try:
                result[key] = parse_fortran_float(m.group(1))
            except (ValueError, IndexError):
                pass

    # Plain integer: NUMBER_SCF_CYCLES
    m = re.search(r"NUMBER_SCF_CYCLES=(\d+)", aux_content)
    if m:
        result["mopac_num_scf_cycles"] = int(m.group(1))

    # Plain string: POINT_GROUP
    m = re.search(r"POINT_GROUP=(\S+)", aux_content)
    if m:
        result["mopac_point_group"] = m.group(1)

    # Plain float: CPU_TIME
    m = re.search(r"CPU_TIME:SECONDS=\s*([\d.]+)", aux_content)
    if m:
        try:
            result["mopac_time_s"] = float(m.group(1))
        except ValueError:
            pass

    # HOMO/LUMO from NUM_ELECTRONS + EIGENVALUES
    m_electrons = re.search(r"NUM_ELECTRONS=(\d+)", aux_content)
    m_eigenvalues = re.search(
        r"EIGENVALUES\[\d+\]=\s*([\s\S]*?)(?:\n\s*[A-Z]|\Z)", aux_content
    )

    if m_electrons and m_eigenvalues:
        try:
            num_electrons = int(m_electrons.group(1))
            # Parse eigenvalue array (plain floats, whitespace-separated)
            eigenvalue_text = m_eigenvalues.group(1).strip()
            eigenvalues = [float(v) for v in eigenvalue_text.split()]

            homo_idx = (num_electrons // 2) - 1  # 0-indexed
            lumo_idx = num_electrons // 2

            if 0 <= homo_idx < len(eigenvalues):
                result["mopac_homo_ev"] = eigenvalues[homo_idx]
            if 0 <= lumo_idx < len(eigenvalues):
                result["mopac_lumo_ev"] = eigenvalues[lumo_idx]

            if not np.isnan(result["mopac_homo_ev"]) and not np.isnan(
                result["mopac_lumo_ev"]
            ):
                result["mopac_gap_ev"] = (
                    result["mopac_lumo_ev"] - result["mopac_homo_ev"]
                )
        except (ValueError, IndexError) as exc:
            LOG.debug(f"Failed to parse HOMO/LUMO from eigenvalues: {exc}")

    return result


def extract_mopac_descriptors(mopac_base_path: Path) -> dict[str, Any]:
    """Extract electronic descriptors from MOPAC .aux file.

    Args:
        mopac_base_path: Path to MOPAC output file (.out); .aux is derived from it

    Returns:
        Dictionary with 11 descriptor values (NaN/None if unavailable)
    """
    aux_path = mopac_base_path.with_suffix(".aux")

    if not aux_path.exists():
        LOG.warning(f"MOPAC .aux file not found: {aux_path}")
        return _empty_descriptors()

    try:
        content = aux_path.read_text(encoding="utf-8", errors="replace")
        descriptors = _parse_aux_file(content)
        LOG.info(f"Extracted MOPAC descriptors from {aux_path.name}")
        return descriptors
    except Exception as exc:
        LOG.warning(f"Failed to parse MOPAC .aux file {aux_path}: {exc}")
        return _empty_descriptors()
