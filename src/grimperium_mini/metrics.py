"""Error metrics with explicit thermochemical units."""

from __future__ import annotations


def absolute_deviation(value: float | None, reference: float | None) -> float | None:
    """Return absolute deviation in the same unit as the inputs."""
    if value is None or reference is None:
        return None
    return abs(value - reference)


def relative_deviation_percent(
    value: float | None, reference: float | None
) -> float | None:
    """Return absolute relative deviation in percent."""
    if value is None or reference is None or reference == 0:
        return None
    return abs((value - reference) / reference) * 100.0
