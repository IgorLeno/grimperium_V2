"""Classify MOPAC ``FORCE`` output without a naive sign-only rule.

MOPAC's manual documents two details that drive this parser:

* a molecule has five (linear) or six (non-linear) translational/rotational
  modes, commonly below about 30 cm-1; and
* the ``DESCRIPTION OF VIBRATIONS`` section reports the projected, non-trivial
  modes directly.

The latter is preferred.  When an adapter supplies a raw all-3N frequency
table, the expected trivial modes are removed by smallest absolute magnitude
only after the output/context identifies whether five or six are expected.
Small residual negative values above the configured imaginary threshold remain
auditable numerical low modes; they are not silently called a saddle.

Reference: https://openmopac.net/Manual/force.html and
https://openmopac.net/Manual/frame.html.
"""

from __future__ import annotations

import re
from collections.abc import Sequence

from semi_imperium.domain import VerificationSettings
from semi_imperium.mopac.models import (
    CandidateState,
    ForceClassification,
    FrequencyDiagnostics,
)

_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][-+]?\d+)?"
_VIBRATION_BLOCK = re.compile(
    r"\bVIBRATION\s+(?P<mode>\d+)\b" r"(?P<body>.*?)(?=\bVIBRATION\s+\d+\b|\Z)",
    re.IGNORECASE | re.DOTALL,
)
_VIBRATION_FREQUENCY = re.compile(
    rf"\bFREQ(?:UENCY)?\.?\s*(?:=\s*)?(?P<frequency>{_NUMBER})",
    re.IGNORECASE,
)
_FREQUENCY_LIST = re.compile(
    r"^\s*(?:ALL\s+3N\s+)?FREQUENC(?:Y|IES)(?:\s*\([^)]*\))?\s*[:=]"
    r"(?P<values>[^\n]+)$",
    re.IGNORECASE | re.MULTILINE,
)
_GRADIENT_NORM = re.compile(
    rf"\bGRADIENT\s+NORM\s*=\s*(?P<value>{_NUMBER})", re.IGNORECASE
)
_TRIVIAL_COUNT = re.compile(
    r"\b(?:EXPECTED\s+)?TRIVIAL\s+MODES?\s*[:=]\s*([356])\b",
    re.IGNORECASE,
)
_FORCE_FAILURE_PATTERNS = (
    re.compile(r"NOT\s+(?:AT|A)\s+(?:TRUE\s+)?STATIONARY\s+POINT", re.IGNORECASE),
    re.compile(r"GRADIENTS?\s+(?:ARE|IS|WERE)\s+(?:TOO\s+)?LARGE", re.IGNORECASE),
    re.compile(r"FORCE\s+CALCULATION\s+(?:FAILED|TERMINATED)", re.IGNORECASE),
    re.compile(r"ABNORMAL\s+TERMINATION", re.IGNORECASE),
    re.compile(r"SCF\s+(?:FAILED|NOT\s+CONVERGED)", re.IGNORECASE),
)
_STATIONARY_MARKERS = (
    re.compile(
        r"GRADIENTS?\s+(?:WERE|WAS|IS)\s+(?:INITIALLY\s+)?ACCEPTABLY\s+(?:LOW|SMALL)",
        re.I,
    ),
    re.compile(r"GRADIENT\s+NORM\s+IS\s+(?:ACCEPTABLY\s+)?(?:LOW|SMALL)", re.I),
    re.compile(r"STATIONARY\s+(?:POINT|GEOMETRY)\s+(?:IS\s+)?CONFIRMED", re.I),
)


def classify_force_output(
    output: str,
    settings: VerificationSettings,
    *,
    atom_count: int | None = None,
    is_linear: bool | None = None,
    frequencies_include_trivial_modes: bool | None = None,
    execution_error: str | None = None,
) -> ForceClassification:
    """Return a minimum, saddle or failure verdict with parsed evidence.

    ``frequencies_include_trivial_modes`` is an explicit adapter diagnostic.
    If omitted, a table headed ``ALL 3N FREQUENCIES`` or output containing
    ``TRIVIAL MODES INCLUDED`` is treated as raw; MOPAC's normal
    ``VIBRATION ... FREQ.`` descriptions are treated as already projected.
    """
    force_detected = bool(
        re.search(r"\b(?:START OF )?FORCE CALCULATION\b", output, re.IGNORECASE)
        or re.search(r"\bFORCE\s+-\s+FORCE CALCULATION SPECIFIED", output, re.I)
    )
    completed = bool(
        re.search(r"\b(?:MOPAC DONE|FORCE CALCULATION COMPLETED)\b", output, re.I)
        or "DESCRIPTION OF VIBRATIONS" in output.upper()
    )
    failure_marker = next(
        (
            pattern.pattern
            for pattern in _FORCE_FAILURE_PATTERNS
            if pattern.search(output)
        ),
        None,
    )
    stationary = any(marker.search(output) for marker in _STATIONARY_MARKERS)
    gradient = _parse_gradient_norm(output)
    notes: list[str] = []

    explicit_modes = _parse_vibration_descriptions(output)
    if explicit_modes:
        mode_numbers = tuple(mode for mode, _ in explicit_modes)
        frequencies = tuple(frequency for _, frequency in explicit_modes)
        source = "mopac_vibration_descriptions"
        includes_trivial = False
    else:
        frequencies = _parse_frequency_lists(output)
        mode_numbers = tuple(range(1, len(frequencies) + 1))
        source = "mopac_frequency_table" if frequencies else "none"
        detected_raw = bool(
            re.search(r"\bALL\s+3N\s+FREQUENC", output, re.IGNORECASE)
            or re.search(r"\bTRIVIAL\s+MODES?\s+INCLUDED\b", output, re.IGNORECASE)
        )
        includes_trivial = (
            detected_raw
            if frequencies_include_trivial_modes is None
            else frequencies_include_trivial_modes
        )

    expected_trivial = _expected_trivial_modes(output, atom_count, is_linear)
    nontrivial_pairs = tuple(zip(mode_numbers, frequencies))
    trivial: tuple[float, ...] = ()
    near_zero_count = sum(
        abs(frequency) <= settings.trivial_mode_cutoff_cm1 for frequency in frequencies
    )
    if includes_trivial and frequencies:
        if expected_trivial is None:
            notes.append("raw_3n_modes_without_linear_geometry_diagnostic")
        else:
            nontrivial_pairs, trivial = _remove_trivial_modes(
                nontrivial_pairs,
                expected_trivial,
                settings.trivial_mode_cutoff_cm1,
            )
            if len(trivial) != expected_trivial:
                notes.append(
                    "raw_3n_table_did_not_contain_the_expected_near_zero_modes"
                )
            elif near_zero_count != expected_trivial:
                notes.append("raw_3n_trivial_mode_assignment_is_ambiguous")

    nontrivial = tuple(frequency for _, frequency in nontrivial_pairs)
    numerical = tuple(
        value
        for value in nontrivial
        if settings.imaginary_frequency_threshold_cm1 <= value < 0.0
    )
    imaginary_pairs = tuple(
        (mode, value)
        for mode, value in nontrivial_pairs
        if value < settings.imaginary_frequency_threshold_cm1
    )
    imaginary = tuple(value for _, value in imaginary_pairs)
    imaginary_modes = tuple(mode for mode, _ in imaginary_pairs)

    if not stationary and force_detected and completed and frequencies:
        # FORCE refuses a genuinely non-stationary geometry unless overridden.
        # A completed normal-mode analysis is therefore itself a MOPAC
        # stationary-point diagnostic, even when a version omits the prose line.
        stationary = failure_marker is None
        if stationary:
            notes.append("stationary_point_inferred_from_completed_force_analysis")

    failure_reason = _failure_reason(
        execution_error=execution_error,
        force_detected=force_detected,
        completed=completed,
        stationary=stationary,
        frequencies=frequencies,
        failure_marker=failure_marker,
        includes_trivial=includes_trivial,
        expected_trivial=expected_trivial,
        trivial=trivial,
        near_zero_count=near_zero_count,
    )
    diagnostics = FrequencyDiagnostics(
        force_calculation_detected=force_detected,
        force_calculation_completed=completed,
        stationary_point_confirmed=stationary,
        frequencies_cm1=frequencies,
        nontrivial_frequencies_cm1=nontrivial,
        trivial_frequencies_cm1=trivial,
        numerical_low_frequencies_cm1=numerical,
        imaginary_frequencies_cm1=imaginary,
        imaginary_mode_numbers=imaginary_modes,
        expected_trivial_modes=expected_trivial,
        frequency_source=source,
        gradient_norm=gradient,
        failure_reason=failure_reason,
        notes=tuple(notes),
    )
    if failure_reason is not None:
        state = CandidateState.VERIFICATION_FAILED
    elif imaginary:
        state = CandidateState.SADDLE_DETECTED
    else:
        state = CandidateState.MINIMUM_VERIFIED
    return ForceClassification(state=state, diagnostics=diagnostics)


def _parse_vibration_descriptions(output: str) -> tuple[tuple[int, float], ...]:
    """Read projected MOPAC vibration descriptions in their printed order."""
    seen: set[int] = set()
    parsed: list[tuple[int, float]] = []
    for block in _VIBRATION_BLOCK.finditer(output):
        mode = int(block.group("mode"))
        if mode in seen:
            continue
        match = _VIBRATION_FREQUENCY.search(block.group("body"))
        if match is None:
            continue
        seen.add(mode)
        parsed.append((mode, _as_float(match.group("frequency"))))
    return tuple(parsed)


def _parse_frequency_lists(output: str) -> tuple[float, ...]:
    """Read compact adapter/raw frequency lists when descriptions are absent."""
    parsed: list[float] = []
    for match in _FREQUENCY_LIST.finditer(output):
        parsed.extend(
            _as_float(value) for value in re.findall(_NUMBER, match["values"])
        )
    return tuple(parsed)


def _parse_gradient_norm(output: str) -> float | None:
    matches = tuple(_GRADIENT_NORM.finditer(output))
    if not matches:
        return None
    return _as_float(matches[-1].group("value"))


def _expected_trivial_modes(
    output: str, atom_count: int | None, is_linear: bool | None
) -> int | None:
    explicit = _TRIVIAL_COUNT.search(output)
    if explicit:
        return int(explicit.group(1))
    if is_linear is not None:
        return 5 if is_linear else 6
    upper = output.upper()
    if "NON-LINEAR MOLECULE" in upper or "NONLINEAR MOLECULE" in upper:
        return 6
    if "LINEAR MOLECULE" in upper:
        return 5
    if atom_count == 1:
        return 3
    if atom_count == 2:
        return 5
    return None


def _remove_trivial_modes(
    pairs: Sequence[tuple[int, float]], expected: int, cutoff: float
) -> tuple[tuple[tuple[int, float], ...], tuple[float, ...]]:
    ordered = sorted(enumerate(pairs), key=lambda item: abs(item[1][1]))
    removable = [
        (position, pair)
        for position, pair in ordered[:expected]
        if abs(pair[1]) <= cutoff
    ]
    removed_positions = {position for position, _ in removable}
    kept = tuple(
        pair for position, pair in enumerate(pairs) if position not in removed_positions
    )
    trivial = tuple(pair[1] for _, pair in removable)
    return kept, trivial


def _failure_reason(
    *,
    execution_error: str | None,
    force_detected: bool,
    completed: bool,
    stationary: bool,
    frequencies: tuple[float, ...],
    failure_marker: str | None,
    includes_trivial: bool,
    expected_trivial: int | None,
    trivial: tuple[float, ...],
    near_zero_count: int,
) -> str | None:
    if execution_error:
        return execution_error
    if failure_marker:
        return f"MOPAC FORCE diagnostic matched {failure_marker!r}"
    if not force_detected:
        return "MOPAC FORCE calculation marker was not found"
    if not completed:
        return "MOPAC FORCE calculation did not reach a completed analysis"
    if not stationary:
        return "MOPAC did not confirm a stationary geometry"
    if not frequencies:
        return "MOPAC FORCE output contained no parseable frequencies"
    if includes_trivial and expected_trivial is None:
        return "Cannot classify an all-3N table without linearity diagnostics"
    if includes_trivial and len(trivial) != expected_trivial:
        return "MOPAC all-3N table lacks the expected near-zero trivial modes"
    if includes_trivial and near_zero_count != expected_trivial:
        return (
            "MOPAC all-3N table has an ambiguous number of near-zero modes; "
            "projected vibration diagnostics are required"
        )
    return None


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


__all__ = ["classify_force_output"]
