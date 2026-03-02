"""MOPAC .out file descriptor parser.

Extracts electronic descriptors from MOPAC2016 .out (human-readable) output files.
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


def _parse_out_file(out_content: str) -> dict[str, Any]:
    """Parse MOPAC .out file content and extract descriptors.

    Args:
        out_content: Full text content of the .out file

    Returns:
        Dictionary with extracted descriptor values
    """
    result = _empty_descriptors()

    # HOF: "FINAL HEAT OF FORMATION = <value> KCAL/MOL"
    # (not extracted here — already handled by energy_extractor.py)

    # Dipole: "DIPOLE           =   X.XXXX DEBYE" (summary line, not component table)
    m = re.search(
        r"^\s*DIPOLE\s*=\s*([\d.]+)\s*DEBYE",
        out_content,
        re.MULTILINE,
    )
    if m:
        try:
            result["mopac_dipole_debye"] = float(m.group(1))
        except ValueError:
            pass

    # IP: "IONIZATION POTENTIAL    =  X.XXXXX EV"
    m = re.search(
        r"IONIZATION\s+POTENTIAL\s*=\s*([\d.]+)\s*EV",
        out_content,
        re.IGNORECASE,
    )
    if m:
        try:
            result["mopac_ionization_potential_ev"] = float(m.group(1))
        except ValueError:
            pass

    # HOMO/LUMO: three known format variants
    #   1) "HOMO LUMO ENERGIES (EV) =  -12.214   4.258"
    #   2) "HOMO LUMO ENERGIES (EV) = -12.214  4.258"
    #   3) "HOMO LUMO ENERGIES(EV) =  -12.214   4.258"
    m = re.search(
        r"HOMO\s+LUMO\s+ENERGIES\s*\(EV\)\s*=\s*"
        r"([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)",
        out_content,
        re.IGNORECASE,
    )
    if m:
        try:
            result["mopac_homo_ev"] = float(m.group(1))
            result["mopac_lumo_ev"] = float(m.group(2))
            result["mopac_gap_ev"] = result["mopac_lumo_ev"] - result["mopac_homo_ev"]
        except ValueError:
            pass

    # COSMO area: "COSMO AREA              =  X.XX SQUARE ANGSTROMS"
    m = re.search(
        r"COSMO\s+AREA\s*=\s*([\d.]+)\s*SQUARE\s+ANGSTROMS",
        out_content,
        re.IGNORECASE,
    )
    if m:
        try:
            result["mopac_cosmo_area_a2"] = float(m.group(1))
        except ValueError:
            pass

    # COSMO volume: "COSMO VOLUME            =  X.XX CUBIC ANGSTROMS"
    m = re.search(
        r"COSMO\s+VOLUME\s*=\s*([\d.]+)\s*CUBIC\s+ANGSTROMS",
        out_content,
        re.IGNORECASE,
    )
    if m:
        try:
            result["mopac_cosmo_volume_a3"] = float(m.group(1))
        except ValueError:
            pass

    # Gradient norm: "GRADIENT NORM =  X.XXXXX"
    m = re.search(
        r"GRADIENT\s+NORM\s*=\s*([\d.]+)",
        out_content,
        re.IGNORECASE,
    )
    if m:
        try:
            result["mopac_gradient_norm"] = float(m.group(1))
        except ValueError:
            pass

    # SCF cycles: "SCF CALCULATIONS        =         X"
    m = re.search(
        r"SCF\s+CALCULATIONS\s*=\s*(\d+)",
        out_content,
        re.IGNORECASE,
    )
    if m:
        result["mopac_num_scf_cycles"] = int(m.group(1))

    # Point group: "POINT GROUP:    C2v"
    m = re.search(
        r"POINT\s+GROUP:\s*(\S+)",
        out_content,
        re.IGNORECASE,
    )
    if m:
        result["mopac_point_group"] = m.group(1)

    # Wall-clock time: "WALL-CLOCK TIME         =  X.XXX SECONDS"
    m = re.search(
        r"WALL-CLOCK\s+TIME\s*=\s*([\d.]+)\s*SECONDS",
        out_content,
        re.IGNORECASE,
    )
    if m:
        try:
            result["mopac_time_s"] = float(m.group(1))
        except ValueError:
            pass

    return result


def extract_mopac_descriptors(mopac_out_path: Path) -> dict[str, Any]:
    """Extract electronic descriptors from MOPAC .out file.

    Args:
        mopac_out_path: Path to MOPAC output file (.out)

    Returns:
        Dictionary with 11 descriptor values (NaN/None if unavailable)
    """
    if not mopac_out_path.exists():
        LOG.warning(f"MOPAC .out file not found: {mopac_out_path}")
        return _empty_descriptors()

    try:
        content = mopac_out_path.read_text(encoding="utf-8", errors="replace")
        descriptors = _parse_out_file(content)
        LOG.info(f"Extracted MOPAC descriptors from {mopac_out_path.name}")
        return descriptors
    except Exception as exc:
        LOG.warning(f"Failed to parse MOPAC .out file {mopac_out_path}: {exc}")
        return _empty_descriptors()
