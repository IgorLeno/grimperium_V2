"""Tests for grimperium.ml.trainer — Training orchestration."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from sklearn.model_selection import train_test_split

from grimperium.core.delta_learning import DeltaLearner
from grimperium.ml.data_loader import load_ml_data
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
        because train is fit on train data, test uses transform-only), and by
        asserting NaN values on held-out rows are imputed from TRAIN medians.
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

        # Controlled check: ensure held-out NaN rows are transformed with
        # medians from train split only (no fit on full dataset).
        _, _, _, pipeline = train(
            synthetic_csv_path,
            test_size=0.2,
            random_state=0,
            return_pipeline=True,
        )
        df, y_cbs, y_pm7 = load_ml_data(synthetic_csv_path)
        df_train, df_test, *_ = train_test_split(
            df,
            y_cbs,
            y_pm7,
            test_size=0.2,
            random_state=0,
        )

        mopac_cols = ["mopac_homo_ev", "mopac_lumo_ev", "mopac_gap_ev"]
        test_nan_mask = df_test[mopac_cols].isna()
        assert test_nan_mask.any().any(), "Fixture must place NaN in held-out rows."

        train_medians = df_train[mopac_cols].median()
        transformed_test = pipeline.transform(df_test)
        feature_index = {col: idx for idx, col in enumerate(pipeline.feature_cols)}

        for col in mopac_cols:
            col_nan_rows = np.where(test_nan_mask[col].to_numpy())[0]
            if len(col_nan_rows) == 0:
                continue
            col_idx = feature_index[col]
            for row_idx in col_nan_rows:
                assert np.isclose(transformed_test[row_idx, col_idx], train_medians[col])

    def test_raises_on_invalid_test_size(self, synthetic_csv_path: Path) -> None:
        """train validates test_size bounds before train_test_split call."""
        with pytest.raises(ValueError, match="test_size must be between 0 and 1"):
            train(synthetic_csv_path, test_size=1.0, random_state=42)
