"""DTOs returned by the results analysis service."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from grimperium.ml.error_analysis import AnalysisResult


@dataclass(frozen=True)
class DivergenceBucket:
    """Relative divergence bucket for prediction-vs-reference analysis."""

    severity: str
    range_min: float
    range_max: float | None
    count: int
    percentage: float


@dataclass(frozen=True)
class ResultsAnalysisReport:
    """High-level analysis report for a dataset or run."""

    source_path: Path
    source_label: str
    analysis: AnalysisResult
    divergence_distribution: list[DivergenceBucket]
    model_label: str | None = None
    method_label: str | None = None
    run_id: str | None = None
    run_status: str | None = None

    @property
    def summary(self) -> dict[str, Any]:
        """Return the underlying analysis summary."""
        return self.analysis.summary

    @property
    def scored_df(self) -> pd.DataFrame:
        """Return the scored analysis rows."""
        return self.analysis.scored_df

    @property
    def outliers_df(self) -> pd.DataFrame:
        """Return outlier rows."""
        return self.analysis.outliers_df

    @property
    def top_n_df(self) -> pd.DataFrame:
        """Return top-error rows."""
        return self.analysis.top_n_df


@dataclass(frozen=True)
class RunComparisonRow:
    """One run summarized for comparison."""

    run_id: str
    method_label: str
    model_label: str | None
    status: str
    n_molecules: int
    mae: float
    rmse: float
    r2: float


@dataclass(frozen=True)
class CompareRunsReport:
    """Basic comparison across runs."""

    rows: list[RunComparisonRow]
    best_run_id: str | None
