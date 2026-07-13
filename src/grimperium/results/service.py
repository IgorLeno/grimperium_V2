"""Service contracts for results analysis and run comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import pandas as pd

from grimperium.calculation.contracts.quantity import KJ_PER_KCAL
from grimperium.ml import error_analysis
from grimperium.ml.persistence import load_model_metadata
from grimperium.results.analysis import divergence_distribution
from grimperium.results.loaders import (
    CANONICAL_MARKER_COLUMNS,
    join_reference_from_dataset,
    load_analysis_dataframe,
    load_canonical_long_form,
)
from grimperium.results.models import (
    CompareRunsReport,
    ResultsAnalysisMode,
    ResultsAnalysisReport,
    RunComparisonRow,
    ScientificRunSummary,
)
from grimperium.runs.models import RunManifest
from grimperium.runs.service import RunService


class ResultsService:
    """Analyze datasets and persisted runs without depending on CLI views."""

    def __init__(self, run_service: RunService | None = None) -> None:
        self.run_service = run_service or RunService.from_environment()

    def analyze_dataset(
        self,
        dataset_path: Path | str | object,
        *,
        dataset_ref: dict[str, Any] | None = None,
        run_meta: RunManifest | None = None,
    ) -> ResultsAnalysisReport:
        """Analyze a CSV path or dataset-like object with a ``path`` attribute."""
        path = _source_path(dataset_path)
        analysis_df = load_analysis_dataframe(path)
        warnings: list[str] = []

        if (
            "H298_cbs" not in analysis_df.columns
            or analysis_df["H298_cbs"].isna().all()
        ):
            ref_path = _dataset_ref_path(dataset_ref)
            if ref_path is not None:
                analysis_df, join_warnings = join_reference_from_dataset(
                    analysis_df, ref_path
                )
                warnings.extend(join_warnings)

        mode = _detect_analysis_mode(analysis_df)
        scientific = _build_scientific_summary(
            path,
            analysis_df,
            run_meta=run_meta,
            warnings=warnings,
            comparison_label=_comparison_label(mode),
        )

        if mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY:
            return ResultsAnalysisReport(
                source_path=path,
                source_label=_source_label(dataset_path, path),
                analysis_mode=mode,
                scientific_summary=scientific,
                analysis=None,
                divergence_distribution=[],
                model_label=_model_label(run_meta.model_ref) if run_meta else None,
                method_label=_method_label(run_meta) if run_meta else None,
                run_id=run_meta.run_id if run_meta else None,
                run_status=run_meta.status.value if run_meta else None,
            )

        stats_df = analysis_df
        if mode is ResultsAnalysisMode.BASELINE_WITH_REFERENCE:
            stats_df = _baseline_as_predicted_for_stats(analysis_df)

        result = error_analysis.analyze(stats_df)
        scored_df = getattr(result, "scored_df", None)
        distribution = (
            divergence_distribution(cast(pd.DataFrame, scored_df))
            if scored_df is not None
            else []
        )
        return ResultsAnalysisReport(
            source_path=path,
            source_label=_source_label(dataset_path, path),
            analysis_mode=mode,
            scientific_summary=scientific,
            analysis=result,
            divergence_distribution=distribution,
            model_label=_model_label(run_meta.model_ref) if run_meta else None,
            method_label=_method_label(run_meta) if run_meta else None,
            run_id=run_meta.run_id if run_meta else None,
            run_status=run_meta.status.value if run_meta else None,
        )

    def analyze_run(self, run: str | RunManifest) -> ResultsAnalysisReport:
        """Analyze the canonical result output attached to a run manifest."""
        manifest = self.run_service.get_run(run) if isinstance(run, str) else run
        try:
            output_path = _analysis_output_path(manifest)
        except ValueError as exc:
            raise ValueError(f"Run incompatible for Results: {exc}") from exc
        if not output_path.exists():
            raise FileNotFoundError(
                f"Output ausente for run {manifest.run_id}: {output_path}"
            )
        try:
            report = self.analyze_dataset(
                output_path,
                dataset_ref=manifest.dataset_ref,
                run_meta=manifest,
            )
        except ValueError as exc:
            message = str(exc)
            if "H298" in message or "column" in message.lower():
                raise ValueError(f"Arquivo inválido for Results: {message}") from exc
            raise
        return ResultsAnalysisReport(
            source_path=report.source_path,
            source_label=manifest.run_id,
            analysis_mode=report.analysis_mode,
            scientific_summary=report.scientific_summary,
            analysis=report.analysis,
            divergence_distribution=report.divergence_distribution,
            model_label=report.model_label or _model_label(manifest.model_ref),
            method_label=_method_label(manifest),
            run_id=manifest.run_id,
            run_status=manifest.status.value,
        )

    def compare_runs(self, run_ids: list[str]) -> CompareRunsReport:
        """Return a basic MAE-sorted comparison for selected runs."""
        rows: list[RunComparisonRow] = []
        expected_property_id: str | None = None
        expected_reference_label: str | None = None
        expected_analysis_mode: ResultsAnalysisMode | None = None
        for run_id in run_ids:
            manifest = self.run_service.get_run(run_id)
            report = self.analyze_run(manifest)
            if not report.has_comparative_metrics:
                raise ValueError(
                    f"Run {run_id} has no comparative metrics "
                    f"(mode={report.analysis_mode.value})"
                )
            reference_label = _reference_label(manifest.dataset_ref)
            if expected_property_id is None:
                expected_property_id = manifest.property_id
                expected_reference_label = reference_label
                expected_analysis_mode = report.analysis_mode
            elif manifest.property_id != expected_property_id:
                raise ValueError(
                    "Cannot compare runs with incompatible property_id values: "
                    f"{expected_property_id} != {manifest.property_id}"
                )
            elif reference_label != expected_reference_label:
                raise ValueError(
                    "Cannot compare runs with incompatible reference values: "
                    f"{expected_reference_label} != {reference_label}"
                )
            elif report.analysis_mode is not expected_analysis_mode:
                raise ValueError(
                    "Cannot compare runs with incompatible analysis_mode values: "
                    f"{expected_analysis_mode.value if expected_analysis_mode else None} "
                    f"!= {report.analysis_mode.value}"
                )
            summary = report.summary
            rows.append(
                RunComparisonRow(
                    run_id=run_id,
                    property_id=manifest.property_id,
                    reference_label=reference_label,
                    analysis_mode=report.analysis_mode.value,
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


def _reference_label(dataset_ref: dict[str, Any] | None) -> str | None:
    if not dataset_ref:
        return None
    for key in ("database_id", "alias", "name", "path"):
        value = dataset_ref.get(key)
        if value:
            return str(value)
    return None


def _method_label(manifest: RunManifest) -> str:
    display_name = manifest.method_snapshot.get("display_name")
    return str(display_name or manifest.method_id)


def _dataset_ref_path(dataset_ref: dict[str, Any] | None) -> Path | None:
    if not dataset_ref:
        return None
    path = dataset_ref.get("path")
    if path is None:
        return None
    resolved = Path(str(path))
    return resolved if resolved.exists() else resolved


def _detect_analysis_mode(df: pd.DataFrame) -> ResultsAnalysisMode:
    has_ref = _has_usable_column(df, "H298_cbs")
    has_pred = _has_usable_column(df, "H298_predicted")
    has_baseline = _has_usable_column(df, "H298_pm7")
    if has_ref and has_pred:
        return ResultsAnalysisMode.PREDICTION_WITH_REFERENCE
    if has_ref and has_baseline:
        return ResultsAnalysisMode.BASELINE_WITH_REFERENCE
    return ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY


def _has_usable_column(df: pd.DataFrame, column: str) -> bool:
    if column not in df.columns:
        return False
    series = df[column]
    return bool(series.notna().any())


def _baseline_as_predicted_for_stats(df: pd.DataFrame) -> pd.DataFrame:
    """In-memory adaptation only: map baseline → predicted for error_analysis."""
    adapted = df.copy()
    adapted["H298_cbs"] = pd.to_numeric(adapted["H298_cbs"], errors="coerce")
    adapted["H298_pm7"] = pd.to_numeric(adapted["H298_pm7"], errors="coerce")
    adapted["H298_predicted"] = adapted["H298_pm7"]
    return adapted


def _comparison_label(mode: ResultsAnalysisMode) -> str | None:
    if mode is ResultsAnalysisMode.BASELINE_WITH_REFERENCE:
        return "PM7 baseline vs reference"
    if mode is ResultsAnalysisMode.PREDICTION_WITH_REFERENCE:
        return "prediction vs reference"
    return None


def _build_scientific_summary(
    path: Path,
    analysis_df: pd.DataFrame,
    *,
    run_meta: RunManifest | None,
    warnings: list[str],
    comparison_label: str | None,
) -> ScientificRunSummary:
    canonical = load_canonical_long_form(path)
    if canonical is not None and CANONICAL_MARKER_COLUMNS.issubset(canonical.columns):
        roles = tuple(
            sorted({str(r).lower() for r in canonical["role"].dropna().unique()})
        )
        hamiltonians = tuple(
            sorted(
                {
                    str(h).upper()
                    for h in canonical.get("hamiltonian", pd.Series(dtype=object))
                    .dropna()
                    .unique()
                    if str(h).strip()
                }
            )
        )
        values = _canonical_values_kcal(canonical)
        n_estimates = len(canonical)
        n_molecules = (
            int(analysis_df["mol_id"].nunique())
            if "mol_id" in analysis_df.columns and not analysis_df.empty
            else int(canonical.get("molecule_name", canonical.index).nunique())
        )
    else:
        roles = tuple(
            sorted(
                col.replace("H298_", "").lower()
                for col in (
                    "H298_pm7",
                    "H298_predicted",
                    "H298_cbs",
                    "delta_correction",
                )
                if col in analysis_df.columns and analysis_df[col].notna().any()
            )
        )
        hamiltonians = ()
        value_cols = [
            c
            for c in ("H298_predicted", "H298_pm7", "H298_cbs")
            if c in analysis_df.columns
        ]
        values = (
            pd.concat([analysis_df[c].dropna() for c in value_cols], ignore_index=True)
            if value_cols
            else pd.Series(dtype=float)
        )
        n_estimates = int(sum(analysis_df[c].notna().sum() for c in value_cols))
        n_molecules = len(analysis_df)

    return ScientificRunSummary(
        n_molecules=n_molecules,
        n_estimates=n_estimates,
        roles=roles,
        hamiltonians=hamiltonians,
        value_min=float(values.min()) if len(values) else None,
        value_max=float(values.max()) if len(values) else None,
        value_mean=float(values.mean()) if len(values) else None,
        value_median=float(values.median()) if len(values) else None,
        run_status=run_meta.status.value if run_meta else None,
        method_id=run_meta.method_id if run_meta else None,
        model_label=_model_label(run_meta.model_ref) if run_meta else None,
        warnings=tuple(warnings),
        comparison_label=comparison_label,
        started_at=(
            run_meta.started_at.isoformat()
            if run_meta is not None and run_meta.started_at is not None
            else None
        ),
        completed_at=(
            run_meta.completed_at.isoformat()
            if run_meta is not None and run_meta.completed_at is not None
            else None
        ),
        duration_seconds=_run_duration_seconds(run_meta),
    )


def _run_duration_seconds(run_meta: RunManifest | None) -> float | None:
    if run_meta is None or run_meta.started_at is None or run_meta.completed_at is None:
        return None
    return (run_meta.completed_at - run_meta.started_at).total_seconds()


def _canonical_values_kcal(canonical: pd.DataFrame) -> pd.Series:
    values: list[float] = []
    for _, row in canonical.iterrows():
        value_kcal = row.get("value_kcal_mol")
        if pd.notna(value_kcal) and str(value_kcal) != "":
            values.append(float(value_kcal))
            continue
        if pd.isna(row.get("canonical_value")):
            continue
        value = float(row["canonical_value"])
        unit = str(row.get("canonical_unit", "kcal/mol"))
        if unit == "kJ/mol":
            value = value / KJ_PER_KCAL
        values.append(value)
    return pd.Series(values, dtype=float)
