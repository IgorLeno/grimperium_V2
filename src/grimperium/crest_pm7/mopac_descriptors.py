"""MOPAC .out file descriptor parser.

Extracts electronic descriptors from MOPAC2016 .out (human-readable) output files.

Historical Note:
    Prior to v2.1, this module parsed `.aux` (auxiliary) files using Fortran D-notation
    conversion and eigenvalue array indexing to derive HOMO/LUMO energies.

    Migration to .out format (2026-03-02) for improved robustness:
    - HOMO/LUMO now read directly from MOPAC summary line (no index guessing)
    - Eliminates eigenvalue indexing errors when EIGENVALUES array is truncated
      (observed: only 20 eigenvalues for molecules with 50+ electrons)
    - Supports 3 output format variants (standard, SOMO, multiline)
    - ~50% code reduction vs previous .aux parser
    - Smaller file reads (.out ~5-8KB vs .aux ~15-20KB)
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
    #   1) Standard: "HOMO LUMO ENERGIES (EV) =  -10.096    -0.351"
    #   2) Open-shell/SOMO: "HOMO (SOMO) / LUMO ENERGIES (EV) = -10.096 / -0.351"
    #   3) Multiline: "HOMO ENERGY (EV) = -10.096" / "LUMO ENERGY (EV) = -0.351"
    _homo: float | None = None
    _lumo: float | None = None

    # Variant 1: standard joint line
    m = re.search(
        r"HOMO\s+LUMO\s+ENERGIES\s*\(EV\)\s*=\s*"
        r"([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)",
        out_content,
        re.IGNORECASE,
    )
    if m:
        try:
            _homo = float(m.group(1))
            _lumo = float(m.group(2))
        except ValueError:
            LOG.warning("Failed to parse HOMO/LUMO (variant 1): %s", m.group(0))

    # Variant 2: open-shell SOMO joint line
    if _homo is None or _lumo is None:
        m = re.search(
            r"HOMO.*?\(SOMO\).*?/.*?LUMO\s+ENERGIES?\s*\(EV\)\s*=\s*"
            r"([-+]?\d+\.?\d*)\s*/\s*([-+]?\d+\.?\d*)",
            out_content,
            re.IGNORECASE,
        )
        if m:
        if m:
            try:
                if _homo is None:
                    _homo = float(m.group(1))
                if _lumo is None:
                    _lumo = float(m.group(2))
            except ValueError:
                LOG.warning("Failed to parse HOMO/LUMO (variant 2): %s", m.group(0))
                LOG.warning("Failed to parse HOMO/LUMO (variant 2): %s", m.group(0))

    # Variant 3: separate HOMO ENERGY / LUMO ENERGY lines
    if _homo is None or _lumo is None:
        m_homo = re.search(
            r"HOMO\s+ENERG(?:Y|IES)\s*\(EV\)\s*=\s*([-+]?\d+\.?\d*)",
            out_content,
            re.IGNORECASE,
        )
        m_lumo = re.search(
            r"LUMO\s+ENERG(?:Y|IES)\s*\(EV\)\s*=\s*([-+]?\d+\.?\d*)",
            out_content,
            re.IGNORECASE,
        )
        try:
            if m_homo and _homo is None:
                _homo = float(m_homo.group(1))
        except ValueError:
            LOG.warning("Failed to parse HOMO (variant 3)")
        try:
            if m_lumo and _lumo is None:
                _lumo = float(m_lumo.group(1))
        except ValueError:
            LOG.warning("Failed to parse LUMO (variant 3)")

    # Store results and compute GAP with validation
    if _homo is not None:
        result["mopac_homo_ev"] = _homo
    if _lumo is not None:
        result["mopac_lumo_ev"] = _lumo
    if _homo is not None and _lumo is not None:
        if _lumo <= _homo:
            LOG.warning(
                "LUMO (%.3f eV) <= HOMO (%.3f eV); GAP set to NaN",
                _lumo,
                _homo,
            )
        else:
            result["mopac_gap_ev"] = _lumo - _homo

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
        try:
            result["mopac_num_scf_cycles"] = int(m.group(1))
        except ValueError:
            pass

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

    Parses the summary block of MOPAC2016 .out format to extract 11 descriptors:
    - mopac_dipole_debye: Dipole moment (Debye)
    - mopac_ionization_potential_ev: Ionization potential (eV)
    - mopac_homo_ev: HOMO energy (eV)
    - mopac_lumo_ev: LUMO energy (eV)
    - mopac_gap_ev: HOMO-LUMO gap (eV, computed as LUMO - HOMO)
    - mopac_cosmo_area_a2: COSMO surface area (Ų)
    - mopac_cosmo_volume_a3: COSMO volume (ų)
    - mopac_gradient_norm: Final gradient norm (kcal/mol/Å)
    - mopac_num_scf_cycles: Number of SCF cycles (int)
    - mopac_point_group: Molecular point group (str)
    - mopac_time_s: Wall-clock time (seconds)

    HOMO/LUMO parsing supports 3 format variants:
    1. Standard: ``HOMO LUMO ENERGIES (EV) = -10.096 -0.351``
    2. SOMO: ``HOMO (SOMO) / LUMO ENERGIES (EV) = -10.096 / -0.351``
    3. Multiline: Separate ``HOMO ENERGY (EV)`` and ``LUMO ENERGY (EV)`` lines

    Args:
        mopac_out_path: Path to MOPAC output file (.out)

    Returns:
        Dictionary with 11 descriptor values. Missing values are set to:
        - ``np.nan`` for numeric fields
        - ``None`` for mopac_point_group (string field)

    Notes:
        - This function reads .out files only. Prior to v2.1, it parsed .aux files.
        - If mopac_out_path does not exist, returns dict with all fields set to defaults.
        - Logs WARNING if file not found or if parsing raises an exception.

    Example:
        >>> desc = extract_mopac_descriptors(Path("mol_001_conf_0.out"))
        >>> desc["mopac_homo_ev"]
        -10.096
        >>> desc["mopac_gap_ev"]
        9.745
    """
    if not mopac_out_path.exists():
        LOG.warning("MOPAC .out file not found: %s", mopac_out_path)
        return _empty_descriptors()

    try:
        content = mopac_out_path.read_text(encoding="utf-8", errors="replace")
        descriptors = _parse_out_file(content)
        LOG.info("Extracted MOPAC descriptors from %s", mopac_out_path.name)
        return descriptors
    except Exception as exc:
        LOG.warning("Failed to parse MOPAC .out file %s: %s", mopac_out_path, exc)
        return _empty_descriptors()
