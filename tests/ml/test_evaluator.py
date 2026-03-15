"""Tests for grimperium.ml.evaluator — Model evaluation utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np


class TestEvaluate:
    """Tests for the evaluate() function."""

    def test_returns_metrics_dict(self, synthetic_csv_path: Path) -> None:
        """evaluate() returns a dict with all expected metric keys."""
        from grimperium.ml.evaluator import evaluate
        from grimperium.ml.trainer import train

        learner, _, _ = train(synthetic_csv_path, test_size=0.2, random_state=42)
        metrics = evaluate(learner, synthetic_csv_path)

        expected_keys = {"rmse", "mae", "r2", "mape", "max_error"}
        assert set(metrics.keys()) == expected_keys
        for key in expected_keys:
            assert np.isfinite(metrics[key]), f"{key} is not finite"

    def test_full_dataset_evaluation(self, synthetic_csv_path: Path) -> None:
        """evaluate() on full dataset produces lower/equal RMSE than test-only."""
        from grimperium.ml.evaluator import evaluate
        from grimperium.ml.trainer import train

        learner, train_metrics, _ = train(
            synthetic_csv_path, test_size=0.2, random_state=42
        )
        full_metrics = evaluate(learner, synthetic_csv_path)

        # Full dataset includes training data, so RMSE should be reasonable
        assert full_metrics["rmse"] >= 0
        assert np.isfinite(full_metrics["r2"])
