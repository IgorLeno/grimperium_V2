# tests/ml/test_gate.py
"""Tests for grimperium.ml.gate — quality-gate evaluation."""

from __future__ import annotations

import pytest

from grimperium.ml.gate import GATE_MAE_MAX, GATE_R2_MIN, GATE_RMSE_MAX, evaluate_gate

_GOOD = {"mae": 2.5, "r2": 0.98, "rmse": 3.5, "mape": 8.0, "max_error": 15.0}


class TestGateThresholdConstants:
    def test_mae_max(self) -> None:
        assert GATE_MAE_MAX == 3.5

    def test_r2_min(self) -> None:
        assert GATE_R2_MIN == 0.97

    def test_rmse_max(self) -> None:
        assert GATE_RMSE_MAX == 5.0


class TestEvaluateGate:
    def test_passes_when_all_criteria_met(self) -> None:
        assert evaluate_gate(_GOOD) is True

    def test_fails_when_mae_exceeds_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "mae": 3.51}) is False

    def test_passes_when_mae_at_exact_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "mae": 3.5}) is True  # ≤, not <

    def test_fails_when_r2_below_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "r2": 0.969}) is False

    def test_passes_when_r2_at_exact_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "r2": 0.97}) is True  # ≥, not >

    def test_fails_when_rmse_exceeds_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "rmse": 5.01}) is False

    def test_passes_when_rmse_at_exact_threshold(self) -> None:
        assert evaluate_gate({**_GOOD, "rmse": 5.0}) is True  # ≤, not <

    def test_high_max_error_does_not_fail_gate(self) -> None:
        """max_error excluded — organic radicals drive it, not model quality."""
        assert evaluate_gate({**_GOOD, "max_error": 999.0}) is True

    def test_high_mape_does_not_fail_gate(self) -> None:
        """MAPE excluded — undefined for H298 near zero."""
        assert evaluate_gate({**_GOOD, "mape": 9999.0}) is True

    def test_all_criteria_fail_returns_false(self) -> None:
        assert evaluate_gate({"mae": 5.0, "r2": 0.90, "rmse": 7.0}) is False

    def test_missing_keys_treated_as_worst_case(self) -> None:
        """Missing mae/r2/rmse default to inf/0.0/inf → gate fails."""
        assert evaluate_gate({}) is False

    def test_current_model_passes(self) -> None:
        """Regression guard: model with MAE=2.62, R²=0.9968, RMSE=3.67 passes."""
        assert evaluate_gate({"mae": 2.62, "r2": 0.9968, "rmse": 3.67}) is True
