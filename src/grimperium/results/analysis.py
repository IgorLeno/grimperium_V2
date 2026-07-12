"""Analysis helpers shared by results services and views."""

from __future__ import annotations

from typing import Final

import pandas as pd

from grimperium.results.models import DivergenceBucket

DIVERGENCE_THRESHOLDS: Final[tuple[tuple[str, float, float | None], ...]] = (
    ("LOW", 0.0, 10.0),
    ("MEDIUM", 10.0, 25.0),
    ("HIGH", 25.0, 50.0),
    ("CRITICAL", 50.0, None),
)


def divergence_distribution(scored_df: pd.DataFrame) -> list[DivergenceBucket]:
    """Build relative-error distribution buckets from an analyzed DataFrame."""
    if "relative_error_pct" not in scored_df.columns or scored_df.empty:
        return [
            DivergenceBucket(
                severity=severity,
                range_min=low,
                range_max=high,
                count=0,
                percentage=0.0,
            )
            for severity, low, high in DIVERGENCE_THRESHOLDS
        ]

    relative_error = scored_df["relative_error_pct"]
    total = int(len(relative_error))
    buckets: list[DivergenceBucket] = []
    for severity, low, high in DIVERGENCE_THRESHOLDS:
        if high is None:
            mask = relative_error >= low
        else:
            mask = (relative_error >= low) & (relative_error < high)
        count = int(mask.sum())
        percentage = count / total * 100.0 if total else 0.0
        buckets.append(
            DivergenceBucket(
                severity=severity,
                range_min=low,
                range_max=high,
                count=count,
                percentage=percentage,
            )
        )
    return buckets
