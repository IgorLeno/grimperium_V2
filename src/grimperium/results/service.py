"""Service contracts for results analysis and run comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import pandas as pd

from grimperium.ml import error_analysis
from grimperium.ml.persistence import load_model_metadata
from grimperium.results.analysis import divergence_distribution
from grimperium.results.loaders import load_analysis_dataframe
from grimperium.results.models import (
    CompareRunsReport,
    ResultsAnalysisReport,
    RunComparisonRow,
)
from grimperium.runs.models import RunManifest
from grimperium.runs.service import RunService


class ResultsService:
    """Analyze datasets and persisted runs without depending on CLI views."""

    def __init__(self, run_service: RunService | None = None) -> None:
        self.run_service = run_service or RunService.from_environment()

    def analyze_dataset(
        self, dataset_path: Path | str | object
    ) -> ResultsAnalysisReport:
        """Analyze a CSV path or dataset-like object with a ``path`` attribute."""
        path = _source_path(dataset_path)
        df = load_analysis_dataframe(path)
        result = error_analysis.analyze(df)
        scored_df = getattr(result, "scored_df", None)
        distribution = (
            divergence_distribution(cast(pd.DataFrame, scored_df))
            if scored_df is not None
            else []
        )
        return ResultsAnalysisReport(
            source_path=path,
            source_label=_source_label(dataset_path, path),
            analysis=result,
            divergence_distribution=distribution,
        )

    def analyze_run(self, run: str | RunManifest) -> ResultsAnalysisReport:
        """Analyze the canonical result output attached to a run manifest."""
        manifest = self.run_service.get_run(run) if isinstance(run, str) else run
        output_path = _analysis_output_path(manifest)
        report = self.analyze_dataset(output_path)
        return ResultsAnalysisReport(
            source_path=report.source_path,
            source_label=manifest.run_id,
            analysis=report.analysis,
            divergence_distribution=report.divergence_distribution,
            model_label=_model_label(manifest.model_ref),
            method_label=_method_label(manifest),
            run_id=manifest.run_id,
            run_status=manifest.status.value,
        )

    def compare_runs(self, run_ids: list[str]) -> CompareRunsReport:
        """Return a basic MAE-sorted comparison for selected runs."""
        rows: list[RunComparisonRow] = []
        for run_id in run_ids:
            report = self.analyze_run(run_id)
            summary = report.summary
            rows.append(
                RunComparisonRow(
                    run_id=run_id,
                    method_label=report.method_label or "-",
                    model_label=report.model_label,
                    status=report.run_status or "-",
                    n_molecules=int(summary["n_molecules"]),
                    mae=float(summary["mae"]),
                    rmse=float(summary["rmse"]),
                    r2=float(summary["r2"]),
                )
            )
        best_run_id = min(rows, key=lambda row: row.mae).run_id if rows else None
        return CompareRunsReport(rows=rows, best_run_id=best_run_id)

    def model_metadata(
        self,
        model_path: Path | None,
        *,
        model_name: str | None = None,
    ) -> dict[str, object] | None:
        """Load display metadata for the selected model when available."""
        if model_path is None:
            return None
        try:
            metadata = load_model_metadata(model_path)
        except (FileNotFoundError, KeyError, OSError, ValueError):
            return {
                "model_label": model_name or model_path.stem,
                "algorithm": "metadata unavailable",
                "mae": None,
                "r2": None,
                "status": "Unavailable",
            }
        test_metrics = metadata.get("metrics", {}).get("test", {})
        algorithm = metadata.get("algorithm") or metadata.get("feature_schema_id")
        return {
            "model_label": model_name or model_path.stem,
            "algorithm": str(algorithm or "model bundle"),
            "mae": test_metrics.get("mae"),
            "r2": test_metrics.get("r2"),
            "status": "Ready",
        }


def _source_path(source: Path | str | object) -> Path:
    if isinstance(source, str | Path):
        return Path(source)
    path = getattr(source, "path", None)
    if path is None:
        raise TypeError("analysis source must be a path or expose a path attribute")
    return Path(path)


def _source_label(source: Path | str | object, path: Path) -> str:
    if isinstance(source, str | Path):
        return path.stem
    for attr in ("alias", "name"):
        value = getattr(source, attr, None)
        if value:
            return str(value)
    return path.stem


def _analysis_output_path(manifest: RunManifest) -> Path:
    preferred_keys = (
        "results_csv",
        "calculation_results_csv",
        "canonical_results_csv",
        "csv",
    )
    for key in preferred_keys:
        path = manifest.output_paths.get(key)
        if path is not None:
            return path
    for path in manifest.output_paths.values():
        if path.suffix == ".csv":
            return path
    raise ValueError(f"Run {manifest.run_id} has no CSV output path")


def _model_label(model_ref: dict[str, Any] | None) -> str | None:
    if model_ref is None:
        return None
    name = model_ref.get("name")
    if name:
        return str(name)
    path = model_ref.get("path")
    if path:
        return Path(str(path)).stem
    return None


def _method_label(manifest: RunManifest) -> str:
    display_name = manifest.method_snapshot.get("display_name")
    return str(display_name or manifest.method_id)
