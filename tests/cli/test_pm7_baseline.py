# tests/cli/test_pm7_baseline.py
"""Tests for _compute_pm7_stats — PM7 baseline absolute-error distribution."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from grimperium.cli.views.databases_view import _compute_pm7_stats


def _make_df(cbs: list[float], pm7: list[float]) -> pd.DataFrame:
    return pd.DataFrame({"H298_cbs": cbs, "H298_pm7": pm7})


class TestComputePm7StatsAbsoluteMetrics:
    def test_mare_is_mean_absolute_error(self) -> None:
        # Errors: 2, 3, 4 → MARE = 3.0
        df = _make_df([5.0, 5.0, 5.0], [3.0, 2.0, 1.0])
        result = _compute_pm7_stats(df)
        assert result["mare"] == pytest.approx(3.0)

    def test_bias_is_mean_pm7_minus_cbs(self) -> None:
        # pm7 - cbs: -2, -3, -4 → bias = -3.0
        df = _make_df([5.0, 5.0, 5.0], [3.0, 2.0, 1.0])
        result = _compute_pm7_stats(df)
        assert result["bias"] == pytest.approx(-3.0)

    def test_r2_is_coefficient_of_determination(self) -> None:
        # Perfect predictor: pm7 == cbs → r2 = 1.0
        df = _make_df([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
        result = _compute_pm7_stats(df)
        assert result["r2"] == pytest.approx(1.0)

    def test_n_is_row_count(self) -> None:
        df = _make_df([1.0, 2.0, 3.0], [0.5, 1.5, 2.5])
        result = _compute_pm7_stats(df)
        assert result["n"] == 3


class TestComputePm7StatsNoRelativeError:
    def test_mre_pct_key_absent(self) -> None:
        df = _make_df([0.1, 0.2], [-3.9, -4.8])  # H298_cbs near zero → RE% would be huge
        result = _compute_pm7_stats(df)
        assert "mre_pct" not in result

    def test_mdre_pct_key_absent(self) -> None:
        df = _make_df([0.1, 0.2], [-3.9, -4.8])
        assert "mdre_pct" not in _compute_pm7_stats(df)

    def test_max_re_pct_key_absent(self) -> None:
        df = _make_df([0.1, 0.2], [-3.9, -4.8])
        assert "max_re_pct" not in _compute_pm7_stats(df)

    def test_std_re_pct_key_absent(self) -> None:
        df = _make_df([0.1, 0.2], [-3.9, -4.8])
        assert "std_re_pct" not in _compute_pm7_stats(df)


class TestComputePm7StatsPercentiles:
    def test_p50_is_median_absolute_error(self) -> None:
        # Errors: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        df = _make_df(list(range(10)), [0.0] * 10)
        result = _compute_pm7_stats(df)
        assert result["p50"] == pytest.approx(np.median(range(10)))

    def test_p90_is_90th_percentile_absolute_error(self) -> None:
        df = _make_df(list(range(10)), [0.0] * 10)
        result = _compute_pm7_stats(df)
        assert result["p90"] == pytest.approx(np.percentile(range(10), 90))

    def test_p95_is_95th_percentile_absolute_error(self) -> None:
        df = _make_df(list(range(10)), [0.0] * 10)
        result = _compute_pm7_stats(df)
        assert result["p95"] == pytest.approx(np.percentile(range(10), 95))


class TestComputePm7StatsAbsoluteThresholds:
    def test_pct_lt_1_counts_errors_below_1_kcalmol(self) -> None:
        # Errors: 0.5, 1.5, 2.5, 4.5, 6.5 → only 0.5 < 1 → 20%
        df = _make_df([0.0, 0.0, 0.0, 0.0, 0.0], [-0.5, -1.5, -2.5, -4.5, -6.5])
        result = _compute_pm7_stats(df)
        assert result["pct_lt_1"] == pytest.approx(20.0)

    def test_pct_lt_2_counts_errors_below_2_kcalmol(self) -> None:
        # Errors: 0.5, 1.5, 2.5, 4.5, 6.5 → 0.5 and 1.5 < 2 → 40%
        df = _make_df([0.0, 0.0, 0.0, 0.0, 0.0], [-0.5, -1.5, -2.5, -4.5, -6.5])
        result = _compute_pm7_stats(df)
        assert result["pct_lt_2"] == pytest.approx(40.0)

    def test_pct_lt_5_counts_errors_below_5_kcalmol(self) -> None:
        # Errors: 0.5, 1.5, 2.5, 4.5, 6.5 → first four < 5 → 80%
        df = _make_df([0.0, 0.0, 0.0, 0.0, 0.0], [-0.5, -1.5, -2.5, -4.5, -6.5])
        result = _compute_pm7_stats(df)
        assert result["pct_lt_5"] == pytest.approx(80.0)

    def test_near_zero_h298_cbs_does_not_distort_absolute_metrics(self) -> None:
        """Near-zero H298_cbs that caused RE% blowup is harmless for absolute stats."""
        # mol_00001-style: H298_cbs=0.152, H298_pm7=-3.83 → |error|=3.98
        df = _make_df([0.152, -12.641], [-3.830, -14.660])
        result = _compute_pm7_stats(df)
        assert result["mare"] == pytest.approx((3.982 + 2.019) / 2, abs=0.01)
        assert "mre_pct" not in result
