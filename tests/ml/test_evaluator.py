"""Tests for grimperium.ml.evaluator — Model evaluation utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from grimperium.ml.evaluator import evaluate
from grimperium.ml.trainer import train


class TestEvaluate:
    """Tests for the evaluate() function."""

    def test_returns_metrics_dict(self, synthetic_csv_path: Path) -> None:
        """evaluate returns metric dict with expected keys.

        Args:
            synthetic_csv_path: Path to synthetic thermo_pm7-like CSV fixture.

        Returns:
            None: Asserts key set and finite metric values.
        """
        learner, _, _, pipeline = train(
            synthetic_csv_path,
            test_size=0.2,
            random_state=42,
            return_pipeline=True,
        )
        metrics = evaluate(learner, synthetic_csv_path, pipeline=pipeline)

        expected_keys = {
            "rmse",
            "mae",
            "r2",
            "mape",
            "max_error",
            "gate_pass",
            "rmse_gate",
            "r2_gate",
        }
        assert set(metrics.keys()) == expected_keys
        for key in ("rmse", "mae", "r2", "mape", "max_error", "rmse_gate", "r2_gate"):
            assert np.isfinite(metrics[key]), f"{key} is not finite"
        assert isinstance(metrics["gate_pass"], bool)

    def test_full_dataset_evaluation(self, synthetic_csv_path: Path) -> None:
        """evaluate() on full dataset produces lower/equal RMSE than test-only."""
        learner, _, test_metrics, pipeline = train(
            synthetic_csv_path,
            test_size=0.2,
            random_state=42,
            return_pipeline=True,
        )
        full_metrics = evaluate(learner, synthetic_csv_path, pipeline=pipeline)

        # Full dataset includes training data, so RMSE should be lower/equal.
        assert full_metrics["rmse"] <= test_metrics["rmse"]
        assert np.isfinite(full_metrics["r2"])

    def test_raises_without_prefitted_pipeline(self, synthetic_csv_path: Path) -> None:
        """evaluate requires a pre-fitted FeaturePipeline to avoid leakage.

        Args:
            synthetic_csv_path: Path to synthetic thermo_pm7-like CSV fixture.

        Returns:
            None: Asserts leakage-prevention guard raises ValueError.
        """
        learner, _, _ = train(synthetic_csv_path, test_size=0.2, random_state=42)

        with pytest.raises(ValueError, match="pre-fitted FeaturePipeline"):
            evaluate(learner, synthetic_csv_path)
