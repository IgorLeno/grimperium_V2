"""DTOs returned by the results analysis service."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

import pandas as pd

from grimperium.ml.error_analysis import AnalysisResult


class ResultsAnalysisMode(str, Enum):
    """Capability mode describing what kind of analysis a run supports."""

    PREDICTION_WITH_REFERENCE = "prediction_with_reference"
    BASELINE_WITH_REFERENCE = "baseline_with_reference"
    SCIENTIFIC_SUMMARY_ONLY = "scientific_summary_only"


@dataclass(frozen=True)
class DivergenceBucket:
    """Relative divergence bucket for prediction-vs-reference analysis."""

    severity: str
    range_min: float
    range_max: float | None
    count: int
    percentage: float


@dataclass(frozen=True)
class ScientificRunSummary:
    """Scientific overview available for every analyzable run output."""

    n_molecules: int
    n_estimates: int
    roles: tuple[str, ...]
    hamiltonians: tuple[str, ...]
    value_min: float | None
    value_max: float | None
    value_mean: float | None
    value_median: float | None
    run_status: str | None = None
    method_id: str | None = None
    model_label: str | None = None
    warnings: tuple[str, ...] = ()
    comparison_label: str | None = None
    started_at: str | None = None
    completed_at: str | None = None
    duration_seconds: float | None = None


@dataclass(frozen=True)
class ResultsAnalysisReport:
    """High-level analysis report for a dataset or run."""

    source_path: Path
    source_label: str
    analysis_mode: ResultsAnalysisMode
    scientific_summary: ScientificRunSummary
    analysis: AnalysisResult | None = None
    divergence_distribution: list[DivergenceBucket] = field(default_factory=list)
    model_label: str | None = None
    method_label: str | None = None
    run_id: str | None = None
    run_status: str | None = None

    @property
    def summary(self) -> dict[str, Any]:
        """Return comparative metrics when available, else a non-metric summary."""
        if self.analysis is not None:
            return self.analysis.summary
        sci = self.scientific_summary
        return {
            "n_molecules": sci.n_molecules,
            "n_estimates": sci.n_estimates,
            "roles": list(sci.roles),
            "hamiltonians": list(sci.hamiltonians),
            "value_min": sci.value_min,
            "value_max": sci.value_max,
            "value_mean": sci.value_mean,
            "value_median": sci.value_median,
            "mae": None,
            "rmse": None,
            "r2": None,
            "bias": None,
            "max_error": None,
            "pearson_r": None,
            "pct_within_1": None,
            "pct_within_2": None,
            "pct_within_5": None,
            "comparison_label": sci.comparison_label,
            "warnings": list(sci.warnings),
        }

    @property
    def scored_df(self) -> pd.DataFrame:
        """Return the scored analysis rows (empty when no comparative analysis)."""
        if self.analysis is None:
            return pd.DataFrame()
        return self.analysis.scored_df

    @property
    def outliers_df(self) -> pd.DataFrame:
        """Return outlier rows (empty when no comparative analysis)."""
        if self.analysis is None:
            return pd.DataFrame()
        return self.analysis.outliers_df

    @property
    def top_n_df(self) -> pd.DataFrame:
        """Return top-error rows (empty when no comparative analysis)."""
        if self.analysis is None:
            return pd.DataFrame()
        return self.analysis.top_n_df

    @property
    def has_comparative_metrics(self) -> bool:
        """Whether MAE/RMSE-style metrics are available."""
        return self.analysis is not None


@dataclass(frozen=True)
class RunComparisonRow:
    """One run summarized for comparison."""

    run_id: str
    property_id: str
    reference_label: str | None
    analysis_mode: str
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
