"""Single source of truth for post-prediction error analysis and outlier detection.

All consumers (ResultsView, html_report, CLI) call analyze() exactly once
and read from AnalysisResult — none recalculate statistics independently.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import pandas as pd

from grimperium.core.metrics import compute_all_metrics

_REQUIRED_COLS: tuple[str, ...] = ("H298_cbs", "H298_predicted")


@dataclass
class ErrorAnalysisConfig:
    threshold_kcalmol: float = 10.0
    zscore_threshold: float = 3.0
    percentile_outlier: float = 95.0
    extreme_iqr_factor: float = 3.0
    top_n: int = 50


@dataclass
class AnalysisResult:
    scored_df: pd.DataFrame
    summary: dict[str, Any]
    outliers_df: pd.DataFrame
    top_n_df: pd.DataFrame
    warnings: list[str]


def analyze(
    df: pd.DataFrame,
    config: ErrorAnalysisConfig | None = None,
) -> AnalysisResult:
    """Analyze prediction errors and flag outliers.

    Args:
        df: DataFrame with at least H298_cbs and H298_predicted columns.
        config: Optional configuration; uses defaults if None.

    Returns:
        AnalysisResult containing scored DataFrame, summary statistics,
        outlier subset, top-N ranking, and data-quality warnings.

    Raises:
        ValueError: If required columns are missing or no valid rows remain
                    after dropping NaN.
    """
    if config is None:
        config = ErrorAnalysisConfig()

    _validate_columns(df)

    valid = df.dropna(subset=list(_REQUIRED_COLS)).copy()
    warnings = _build_warnings(df, valid)

    if len(valid) == 0:
        raise ValueError(
            "No valid rows: all rows have NaN in H298_cbs or H298_predicted."
        )

    scored = _compute_scored(valid)
    scored = _detect_outliers(scored, config)

    y_cbs = scored["H298_cbs"].to_numpy(dtype=np.float64)
    y_pred = scored["H298_predicted"].to_numpy(dtype=np.float64)
    errors = scored["signed_error"].to_numpy(dtype=np.float64)
    abs_errs = scored["abs_error"].to_numpy(dtype=np.float64)
    n = len(scored)

    base = compute_all_metrics(y_cbs, y_pred)

    bias = float(np.mean(errors))
    median_error = float(np.median(errors))
    std_error = float(np.std(errors))

    if np.std(y_cbs) == 0.0 or np.std(y_pred) == 0.0:
        pearson_r: float = float("nan")
    else:
        pearson_r = float(np.corrcoef(y_cbs, y_pred)[0, 1])

    pct_within_1 = float(np.sum(abs_errs <= 1.0) / n * 100)
    pct_within_2 = float(np.sum(abs_errs <= 2.0) / n * 100)
    pct_within_5 = float(np.sum(abs_errs <= 5.0) / n * 100)

    n_outliers = int((scored["outlier_score"] > 0).sum())

    severity_counts: dict[str, int] = {
        f"n_{s}": int((scored["severity"] == s).sum())
        for s in ("normal", "elevated", "high", "critical", "extreme")
    }

    summary: dict[str, Any] = {
        **base,
        "bias": bias,
        "median_error": median_error,
        "std_error": std_error,
        "pearson_r": pearson_r,
        "pct_within_1": pct_within_1,
        "pct_within_2": pct_within_2,
        "pct_within_5": pct_within_5,
        "n_molecules": n,
        "n_outliers": n_outliers,
        "outlier_ratio": n_outliers / n,
        "p50_abs_error": float(np.percentile(abs_errs, 50)),
        "p75_abs_error": float(np.percentile(abs_errs, 75)),
        "p90_abs_error": float(np.percentile(abs_errs, 90)),
        "p95_abs_error": float(np.percentile(abs_errs, 95)),
        "p99_abs_error": float(np.percentile(abs_errs, 99)),
        **severity_counts,
    }

    outliers_df = scored[scored["outlier_score"] > 0].copy()

    top_n_df = scored.nlargest(config.top_n, "abs_error").reset_index(drop=True).copy()
    top_n_df.insert(0, "rank", range(1, len(top_n_df) + 1))

    return AnalysisResult(
        scored_df=scored,
        summary=summary,
        outliers_df=outliers_df,
        top_n_df=top_n_df,
        warnings=warnings,
    )


# ── Internal helpers ──────────────────────────────────────────────────────────


def _validate_columns(df: pd.DataFrame) -> None:
    missing = [c for c in _REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"DataFrame missing required columns: {missing}")


def _severity_label(abs_error: float) -> str:
    if abs_error < 5.0:
        return "normal"
    if abs_error < 10.0:
        return "elevated"
    if abs_error < 25.0:
        return "high"
    if abs_error < 50.0:
        return "critical"
    return "extreme"


def _compute_scored(df: pd.DataFrame) -> pd.DataFrame:
    scored: pd.DataFrame = cast(pd.DataFrame, df.copy())

    scored["signed_error"] = scored["H298_predicted"] - scored["H298_cbs"]
    scored["abs_error"] = scored["signed_error"].abs()
    scored["squared_error"] = scored["signed_error"] ** 2

    cbs_abs = scored["H298_cbs"].abs()
    safe_cbs = cbs_abs.where(cbs_abs > 1e-10, other=1e-10)
    scored["relative_error_pct"] = (scored["abs_error"] / safe_cbs) * 100.0

    abs_series = scored["abs_error"]
    std = float(abs_series.std())
    mean_val = float(abs_series.mean())
    if std > 0.0:
        scored["abs_error_zscore"] = (abs_series - mean_val) / std
    else:
        scored["abs_error_zscore"] = pd.Series(0.0, index=scored.index)

    scored["severity"] = scored["abs_error"].apply(_severity_label)

    return scored


def _detect_outliers(df: pd.DataFrame, config: ErrorAnalysisConfig) -> pd.DataFrame:
    scored: pd.DataFrame = cast(pd.DataFrame, df.copy())
    abs_series = scored["abs_error"]

    q1 = float(abs_series.quantile(0.25))
    q3 = float(abs_series.quantile(0.75))
    iqr = q3 - q1
    iqr_upper = q3 + 1.5 * iqr
    extreme_upper = q3 + config.extreme_iqr_factor * iqr
    p_thresh = float(abs_series.quantile(config.percentile_outlier / 100.0))

    scored["iqr_outlier"] = abs_series > iqr_upper
    scored["extreme_iqr_outlier"] = abs_series > extreme_upper
    scored["percentile_outlier"] = abs_series > p_thresh
    scored["threshold_outlier"] = abs_series > config.threshold_kcalmol

    scored["outlier_score"] = (
        scored["threshold_outlier"].astype(int)
        + scored["iqr_outlier"].astype(int)
        + scored["extreme_iqr_outlier"].astype(int) * 2
        + (scored["abs_error_zscore"] >= config.zscore_threshold).astype(int)
        + scored["percentile_outlier"].astype(int)
        + (scored["severity"] == "critical").astype(int)
        + (scored["severity"] == "extreme").astype(int) * 2
    )

    return scored


def _build_warnings(original: pd.DataFrame, valid: pd.DataFrame) -> list[str]:
    warnings: list[str] = []

    n_dropped = len(original) - len(valid)
    if n_dropped > 0:
        warnings.append(
            f"{n_dropped} row(s) dropped: NaN in H298_cbs or H298_predicted."
        )

    if "H298_pm7" in original.columns:
        n_pm7_nan = int(original["H298_pm7"].isna().sum())
        if n_pm7_nan > 0:
            warnings.append(f"{n_pm7_nan} row(s) with NaN in H298_pm7.")

    if "status" in valid.columns:
        n_bad = int((valid["status"] != "OK").sum())
        if n_bad > 0:
            warnings.append(f"{n_bad} molecule(s) with status != 'OK'.")

    if "cbs_quality_flag" in valid.columns:
        n_bad_flag = int((valid["cbs_quality_flag"] != "OK").sum())
        if n_bad_flag > 0:
            warnings.append(f"{n_bad_flag} molecule(s) with cbs_quality_flag != 'OK'.")

    if len(valid) > 1:
        abs_errs = (valid["H298_predicted"] - valid["H298_cbs"]).abs()
        max_err = float(abs_errs.max())
        rest = abs_errs[abs_errs < max_err]
        if len(rest) > 0:
            rest_max = float(rest.max())
            if rest_max > 0.0 and max_err > 3.0 * rest_max:
                warnings.append(
                    f"Max abs_error ({max_err:.2f}) is >3× the next highest "
                    f"({rest_max:.2f}). Possible anomaly."
                )

    if len(valid) > 0:
        cbs_abs = valid["H298_cbs"].abs()
        safe_cbs = cbs_abs.where(cbs_abs > 1e-10, other=1e-10)
        rel_err = (valid["H298_predicted"] - valid["H298_cbs"]).abs() / safe_cbs * 100.0
        n_absurd = int((rel_err > 500.0).sum())
        if n_absurd > 0:
            warnings.append(f"{n_absurd} molecule(s) with relative_error_pct > 500%.")

    return warnings
