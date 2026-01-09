"""Energy extraction from MOPAC output files.

Implements 5 regex patterns for robust Heat of Formation (HOF) extraction
with confidence levels.
"""

import logging
import re
from typing import Optional

from .config import HOFConfidence

LOG = logging.getLogger("grimperium.crest_pm7.energy_extractor")

# HOF bounds (kcal/mol)
HOF_MIN_BOUND = -500.0
HOF_MAX_BOUND = +500.0

# 5 Regex patterns in order of confidence
PATTERNS = [
    # Pattern 1: FINAL HEAT OF FORMATION (HIGH)
    # Line: "FINAL HEAT OF FORMATION =      -17.93603 KCAL/MOL =      -75.04 KJ/MOL"
    {
        "name": "FINAL_HEAT_OF_FORMATION",
        "regex": re.compile(
            r"FINAL\s+HEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+\.?\d*)\s*KCAL/MOL",
            re.IGNORECASE,
        ),
        "confidence": HOFConfidence.HIGH,
        "use_finditer_last": False,
    },
    # Pattern 2: HEAT OF FORMATION without FINAL (HIGH)
    # Line: "HEAT OF FORMATION       =        -17.93603 KCAL/MOL"
    {
        "name": "HEAT_OF_FORMATION",
        "regex": re.compile(
            r"HEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+\.?\d*)\s*KCAL/MOL",
            re.IGNORECASE,
        ),
        "confidence": HOFConfidence.HIGH,
        "use_finditer_last": False,
    },
    # Pattern 3: Compact format (MEDIUM)
    # Line: "HEAT OF FORMATION=-17.93603"
    {
        "name": "HEAT_OF_FORMATION_COMPACT",
        "regex": re.compile(
            r"HEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+\.?\d+)",
            re.IGNORECASE,
        ),
        "confidence": HOFConfidence.MEDIUM,
        "use_finditer_last": False,
    },
    # Pattern 4: Thermal energy 298K (MEDIUM)
    # Line: "THERMAL ENERGY AT 298K = -12.345"
    {
        "name": "THERMAL_ENERGY_298K",
        "regex": re.compile(
            r"THERMAL\s+ENERGY\s+AT\s+298K?\s*=\s*([-+]?\d+\.?\d*)",
            re.IGNORECASE,
        ),
        "confidence": HOFConfidence.MEDIUM,
        "use_finditer_last": False,
    },
    # Pattern 5: Fallback - last occurrence (LOW)
    # Any "... -17.93603 KCAL/MOL ..."
    {
        "name": "FALLBACK_LAST_KCAL",
        "regex": re.compile(
            r"([-+]?\d+\.?\d+)\s+KCAL/MOL",
            re.IGNORECASE,
        ),
        "confidence": HOFConfidence.LOW,
        "use_finditer_last": True,
    },
]


def validate_hof(hof: float, nheavy: int) -> tuple[bool, str]:
    """Validate extracted HOF value.

    Args:
        hof: Heat of formation in kcal/mol
        nheavy: Number of heavy atoms in molecule

    Returns:
        Tuple of (valid, message)
    """
    # Hard fail if outside bounds
    if hof < HOF_MIN_BOUND or hof > HOF_MAX_BOUND:
        return False, f"HOF {hof:.2f} outside bounds [{HOF_MIN_BOUND}, {HOF_MAX_BOUND}]"

    # Soft warning if far from expected (~-10*nheavy)
    expected = -10.0 * nheavy
    if abs(hof - expected) > 300.0:
        LOG.warning(
            f"HOF {hof:.2f} deviates significantly from expected ~{expected:.0f}"
        )

    return True, "OK"


def extract_hof(
    content: str,
    nheavy: Optional[int] = None,
) -> tuple[Optional[float], Optional[str], Optional[HOFConfidence]]:
    """Extract Heat of Formation from MOPAC output.

    Tries 5 patterns in order of confidence, returning the first successful match.

    Args:
        content: MOPAC output file content
        nheavy: Number of heavy atoms (for validation)

    Returns:
        Tuple of (hof, method_name, confidence)
        - Success: (hof=-17.93, method="FINAL_HEAT_OF_FORMATION", confidence=HIGH)
        - Parse failure: (None, None, None)
        - Bounds failure: (None, method_name, confidence) - preserves diagnostic
    """
    for pattern in PATTERNS:
        name: str = pattern["name"]  # type: ignore[assignment]
        regex: re.Pattern[str] = pattern["regex"]  # type: ignore[assignment]
        confidence: HOFConfidence = pattern["confidence"]  # type: ignore[assignment]
        use_last: bool = pattern["use_finditer_last"]  # type: ignore[assignment]

        LOG.debug(f"Trying pattern: {name}")

        if use_last:
            # Use finditer and take last match
            matches = list(regex.finditer(content))
            if matches:
                match = matches[-1]
            else:
                match = None
        else:
            match = regex.search(content)

        if match:
            try:
                hof = float(match.group(1))
                LOG.debug(f"Pattern {name} matched: HOF={hof:.4f}")

                # Validate if nheavy provided
                if nheavy is not None:
                    valid, msg = validate_hof(hof, nheavy)
                    if not valid:
                        LOG.warning(f"HOF validation failed ({name}): {msg}")
                        return None, name, confidence

                return hof, name, confidence

            except (ValueError, IndexError) as e:
                LOG.debug(f"Pattern {name} parse error: {e}")
                continue

    LOG.warning("No HOF pattern matched")
    return None, None, None


def extract_hof_from_file(
    filepath: str,
    nheavy: Optional[int] = None,
) -> tuple[Optional[float], Optional[str], Optional[HOFConfidence]]:
    """Extract HOF from a MOPAC output file.

    Args:
        filepath: Path to MOPAC .out file
        nheavy: Number of heavy atoms (for validation)

    Returns:
        Same as extract_hof
    """
    try:
        with open(filepath, encoding="utf-8", errors="replace") as f:
            content = f.read()
        return extract_hof(content, nheavy)
    except OSError as e:
        LOG.warning(f"Cannot read file {filepath}: {e}")
        return None, None, None
