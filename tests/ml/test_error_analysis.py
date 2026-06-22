"""Tests for src/grimperium/ml/error_analysis.py.

All tests use synthetic DataFrames — never read data/thermo_pm7.csv.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from grimperium.ml.error_analysis import (
    AnalysisResult,
    ErrorAnalysisConfig,
    analyze,
)

# ── Fixtures ──────────────────────────────────────────────────────────────────


def _make_df(
    cbs: list[float],
    pred: list[float],
    **extra: object,
) -> pd.DataFrame:
    data: dict[str, object] = {"H298_cbs": cbs, "H298_predicted": pred}
    data.update(extra)
    return pd.DataFrame(data)


@pytest.fixture()
def simple_df() -> pd.DataFrame:
    """5 molecules with known errors: [1, 2, 3, 1, 5] kcal/mol."""
    return _make_df(
        [-100.0, -200.0, -300.0, -400.0, -500.0],
        [-101.0, -198.0, -303.0, -399.0, -505.0],
    )


@pytest.fixture()
def outlier_df() -> pd.DataFrame:
    """9 normal molecules + 1 extreme outlier (abs_error=100)."""
    cbs = [-100.0] * 10
    pred = [-99.0] * 9 + [0.0]
    return pd.DataFrame({"H298_cbs": cbs, "H298_predicted": pred})


# ── Column derivation ─────────────────────────────────────────────────────────


def test_signed_error_computed_correctly(simple_df: pd.DataFrame) -> None:
    result = analyze(simple_df)
    scored = result.scored_df
    expected = simple_df["H298_predicted"] - simple_df["H298_cbs"]
    pd.testing.assert_series_equal(
        scored["signed_error"].reset_index(drop=True),
        expected.reset_index(drop=True),
        check_names=False,
    )


def test_abs_error_is_positive(simple_df: pd.DataFrame) -> None:
    result = analyze(simple_df)
    assert (result.scored_df["abs_error"] >= 0).all()


def test_squared_error_computed_correctly(simple_df: pd.DataFrame) -> None:
    result = analyze(simple_df)
    scored = result.scored_df
    np.testing.assert_allclose(
        scored["squared_error"].to_numpy(),
        scored["signed_error"].to_numpy() ** 2,
    )


def test_relative_error_pct_near_zero_ref_is_finite() -> None:
    """H298_cbs near zero should not produce inf or NaN."""
    df = _make_df([0.0, 1e-15, -1e-15], [1.0, 1.0, -1.0])
    result = analyze(df)
    rel = result.scored_df["relative_error_pct"]
    assert rel.notna().all(), "relative_error_pct must not contain NaN"
    assert np.isfinite(rel.to_numpy()).all(), "relative_error_pct must be finite"


# ── Severity boundaries ───────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "abs_error,expected_severity",
    [
        (4.9, "normal"),
        (5.0, "elevated"),
        (9.9, "elevated"),
        (10.0, "high"),
        (24.9, "high"),
        (25.0, "critical"),
        (49.9, "critical"),
        (50.0, "extreme"),
    ],
)
def test_severity_exact_boundaries(abs_error: float, expected_severity: str) -> None:
    """Exact boundary values map to the correct severity label."""
    df = _make_df([0.0], [abs_error])
    result = analyze(df)
    assert result.scored_df.iloc[0]["severity"] == expected_severity


# ── Outlier detection flags ───────────────────────────────────────────────────


def test_iqr_outlier_flag(outlier_df: pd.DataFrame) -> None:
    """The extreme molecule is flagged as IQR outlier."""
    result = analyze(outlier_df)
    scored = result.scored_df
    extreme_idx = scored["abs_error"].idxmax()
    assert scored.loc[extreme_idx, "iqr_outlier"] is True or bool(
        scored.loc[extreme_idx, "iqr_outlier"]
    )


def test_extreme_iqr_outlier_flag(outlier_df: pd.DataFrame) -> None:
    """The extreme molecule is flagged as extreme-IQR outlier."""
    result = analyze(outlier_df)
    scored = result.scored_df
    extreme_idx = scored["abs_error"].idxmax()
    assert bool(scored.loc[extreme_idx, "extreme_iqr_outlier"])


def test_threshold_outlier_flag() -> None:
    """Molecules above the threshold are flagged."""
    df = _make_df([-100.0, -200.0], [-90.0, -199.0])  # errors: 10, 1
    config = ErrorAnalysisConfig(threshold_kcalmol=5.0)
    result = analyze(df, config)
    scored = result.scored_df
    above = scored[scored["abs_error"] > 5.0]
    assert bool(above["threshold_outlier"].all())
    below = scored[scored["abs_error"] <= 5.0]
    assert not bool(below["threshold_outlier"].any())


def test_percentile_outlier_flag(outlier_df: pd.DataFrame) -> None:
    """The extreme molecule is above the 95th percentile threshold."""
    config = ErrorAnalysisConfig(percentile_outlier=95.0)
    result = analyze(outlier_df, config)
    scored = result.scored_df
    extreme_idx = scored["abs_error"].idxmax()
    assert bool(scored.loc[extreme_idx, "percentile_outlier"])


def test_zscore_flag_on_extreme() -> None:
    """With 19 near-zero-error molecules + 1 extreme (abs_error=100), zscore>=3.

    10 molecules are not enough: sample zscore is bounded by (n-1)/sqrt(n),
    which for n=10 gives ~2.85. 20 molecules gives ~4.25, safely above 3.
    """
    cbs = [-100.0 * (i + 1) for i in range(20)]
    pred = list(cbs[:-1]) + [cbs[-1] + 100.0]
    df = pd.DataFrame({"H298_cbs": cbs, "H298_predicted": pred})
    result = analyze(df)
    scored = result.scored_df
    extreme_idx = scored["abs_error"].idxmax()
    assert float(scored.loc[extreme_idx, "abs_error_zscore"]) >= 3.0


# ── Outlier score ─────────────────────────────────────────────────────────────


def test_outlier_score_equals_weighted_sum(simple_df: pd.DataFrame) -> None:
    """outlier_score equals the weighted sum of individual criterion flags."""
    config = ErrorAnalysisConfig(zscore_threshold=3.0)
    result = analyze(simple_df, config)
    scored = result.scored_df

    expected = (
        scored["threshold_outlier"].astype(int)
        + scored["iqr_outlier"].astype(int)
        + scored["extreme_iqr_outlier"].astype(int) * 2
        + (scored["abs_error_zscore"] >= config.zscore_threshold).astype(int)
        + scored["percentile_outlier"].astype(int)
        + (scored["severity"] == "critical").astype(int)
        + (scored["severity"] == "extreme").astype(int) * 2
    )

    pd.testing.assert_series_equal(
        scored["outlier_score"].astype(int).reset_index(drop=True),
        expected.astype(int).reset_index(drop=True),
        check_names=False,
    )


# ── Empty outlier path ────────────────────────────────────────────────────────


def test_no_outliers_returns_empty_df() -> None:
    """Perfect predictions give abs_error=0 everywhere — outliers_df is empty."""
    df = _make_df(
        [-100.0, -200.0, -300.0],
        [-100.0, -200.0, -300.0],  # identical → abs_error = 0 for all
    )
    result = analyze(df)
    assert result.outliers_df.empty


# ── Optional columns ──────────────────────────────────────────────────────────


def test_optional_columns_absent_no_keyerror(simple_df: pd.DataFrame) -> None:
    """analyze() works without optional columns — no KeyError."""
    result = analyze(simple_df)
    assert isinstance(result, AnalysisResult)


def test_optional_columns_present_carried_through() -> None:
    """Optional columns survive into scored_df when present in input."""
    df = _make_df(
        [-100.0],
        [-101.0],
        batch_id=["b1"],
        rdkit_mol_weight=[100.0],
    )
    result = analyze(df)
    assert "batch_id" in result.scored_df.columns
    assert "rdkit_mol_weight" in result.scored_df.columns


# ── NaN handling ──────────────────────────────────────────────────────────────


def test_nan_rows_ignored_no_crash() -> None:
    """Rows with NaN in required columns are silently dropped."""
    df = pd.DataFrame(
        {
            "H298_cbs": [-100.0, -200.0, np.nan, -400.0],
            "H298_predicted": [-100.0, np.nan, -300.0, -400.0],
        }
    )
    result = analyze(df)
    # Only rows 0 and 3 are valid — both have perfect predictions
    assert int(result.summary["n_molecules"]) == 2
    assert abs(result.summary["mae"]) < 1e-10


def test_nan_rows_appear_in_warnings() -> None:
    """Dropped NaN rows are reported in warnings."""
    df = pd.DataFrame(
        {
            "H298_cbs": [-100.0, np.nan],
            "H298_predicted": [-101.0, -200.0],
        }
    )
    result = analyze(df)
    assert any("dropped" in w.lower() or "nan" in w.lower() for w in result.warnings)


def test_all_nan_raises_value_error() -> None:
    """All-NaN DataFrame raises ValueError after filtering."""
    df = pd.DataFrame(
        {
            "H298_cbs": [np.nan, np.nan],
            "H298_predicted": [np.nan, np.nan],
        }
    )
    with pytest.raises(ValueError, match="No valid rows"):
        analyze(df)


# ── Warnings ─────────────────────────────────────────────────────────────────


def test_warnings_bad_status() -> None:
    """status != 'OK' triggers a warning."""
    df = _make_df(
        [-100.0, -200.0],
        [-101.0, -199.0],
        status=["OK", "FAILED"],
    )
    result = analyze(df)
    assert any("status" in w for w in result.warnings)


def test_warnings_bad_cbs_quality_flag() -> None:
    """cbs_quality_flag != 'OK' triggers a warning."""
    df = _make_df(
        [-100.0, -200.0],
        [-101.0, -199.0],
        cbs_quality_flag=["OK", "SUSPECT"],
    )
    result = analyze(df)
    assert any("cbs_quality_flag" in w for w in result.warnings)


def test_warnings_3x_max_anomaly() -> None:
    """When one molecule's error is >3× the next highest, a warning is emitted."""
    # 4 molecules at ~1 kcal/mol error, 1 at 100 kcal/mol
    df = _make_df(
        [-100.0, -200.0, -300.0, -400.0, -500.0],
        [-101.0, -201.0, -301.0, -401.0, -400.0],  # last has abs_error=100
    )
    result = analyze(df)
    assert any("3" in w or "anomaly" in w.lower() for w in result.warnings)


def test_warnings_absurd_relative_error() -> None:
    """relative_error_pct > 500% triggers a warning."""
    # H298_cbs = -1.0, H298_predicted = -1000.0 → abs_err=999, rel=99900%
    df = _make_df([-1.0, -100.0], [-1000.0, -101.0])
    result = analyze(df)
    assert any("500%" in w or "relative_error" in w for w in result.warnings)


# ── Missing required columns ──────────────────────────────────────────────────


def test_missing_required_column_raises() -> None:
    """ValueError is raised if H298_cbs or H298_predicted is absent."""
    df = pd.DataFrame({"H298_cbs": [-100.0], "other": [1.0]})
    with pytest.raises(ValueError, match="missing required columns"):
        analyze(df)


# ── Summary values match compute_all_metrics ─────────────────────────────────


def test_summary_uses_compute_all_metrics(simple_df: pd.DataFrame) -> None:
    """summary metrics match compute_all_metrics for the same data."""
    from grimperium.core.metrics import compute_all_metrics

    result = analyze(simple_df)
    s = result.summary

    y_cbs = simple_df["H298_cbs"].to_numpy(dtype=np.float64)
    y_pred = simple_df["H298_predicted"].to_numpy(dtype=np.float64)
    expected = compute_all_metrics(y_cbs, y_pred)

    assert abs(s["mae"] - expected["mae"]) < 1e-10
    assert abs(s["rmse"] - expected["rmse"]) < 1e-10
    assert abs(s["r2"] - expected["r2"]) < 1e-10
    assert abs(s["mape"] - expected["mape"]) < 1e-10
    assert abs(s["max_error"] - expected["max_error"]) < 1e-10


def test_summary_contains_additional_keys(simple_df: pd.DataFrame) -> None:
    """Summary contains bias, pearson_r, percentiles, severity counts."""
    result = analyze(simple_df)
    s = result.summary
    for key in (
        "bias",
        "median_error",
        "std_error",
        "pearson_r",
        "pct_within_1",
        "pct_within_2",
        "pct_within_5",
        "n_molecules",
        "n_outliers",
        "outlier_ratio",
        "p50_abs_error",
        "p95_abs_error",
        "n_normal",
        "n_elevated",
        "n_high",
        "n_critical",
        "n_extreme",
    ):
        assert key in s, f"Missing summary key: {key}"


# ── top_n_df ──────────────────────────────────────────────────────────────────


def test_top_n_df_has_rank_column(simple_df: pd.DataFrame) -> None:
    result = analyze(simple_df)
    assert "rank" in result.top_n_df.columns


def test_top_n_df_sorted_descending(simple_df: pd.DataFrame) -> None:
    result = analyze(simple_df)
    abs_errors = result.top_n_df["abs_error"].to_numpy()
    assert (abs_errors[:-1] >= abs_errors[1:]).all()


def test_top_n_df_max_rows() -> None:
    """top_n_df contains at most config.top_n rows."""
    df = _make_df(list(range(-20, 0)), list(range(-21, -1)))
    config = ErrorAnalysisConfig(top_n=5)
    result = analyze(df, config)
    assert len(result.top_n_df) <= 5
