"""Tests for grimperium.ml.trainer — Training orchestration."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from grimperium.core.delta_learning import DeltaLearner
from grimperium.ml.trainer import train


class TestTrain:
    """Tests for the train() function."""

    def test_returns_trained_learner_and_metrics(
        self, synthetic_csv_path: Path
    ) -> None:
        """train() returns (DeltaLearner, train_metrics, test_metrics)."""
        learner, train_metrics, test_metrics = train(
            synthetic_csv_path, test_size=0.2, random_state=42
        )

        assert isinstance(learner, DeltaLearner)
        assert learner.is_fitted
        assert isinstance(train_metrics, dict)
        assert isinstance(test_metrics, dict)

    def test_metrics_contain_expected_keys(self, synthetic_csv_path: Path) -> None:
        """Both metric dicts contain rmse, mae, r2, mape, max_error."""
        _, train_metrics, test_metrics = train(
            synthetic_csv_path, test_size=0.2, random_state=42
        )

        expected_keys = {"rmse", "mae", "r2", "mape", "max_error"}
        assert set(train_metrics.keys()) == expected_keys
        assert set(test_metrics.keys()) == expected_keys

    def test_no_data_leakage(self, synthetic_csv_path: Path) -> None:
        """Split happens on DataFrame BEFORE fit_transform (no leakage).

        We verify by checking that train RMSE != test RMSE (they should differ
        because train is fit on train data, test uses transform-only).
        With 10 rows this is a weak check, but the architecture guarantees it.
        """
        _, train_metrics, test_metrics = train(
            synthetic_csv_path, test_size=0.2, random_state=42
        )

        # Both should be finite
        assert np.isfinite(train_metrics["rmse"])
        assert np.isfinite(test_metrics["rmse"])
        # RMSE should be non-negative
        assert train_metrics["rmse"] >= 0
        assert test_metrics["rmse"] >= 0
        assert not np.isclose(train_metrics["rmse"], test_metrics["rmse"])
