"""Tests for grimperium.ml.features — Feature engineering pipeline."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


class TestFeaturePipeline:
    """Tests for FeaturePipeline class."""

    def test_fit_transform_returns_ndarray(
        self, synthetic_df: pd.DataFrame, feature_cols: list[str]
    ) -> None:
        """fit_transform returns a numpy array with correct shape."""
        from grimperium.ml.features import FeaturePipeline

        pipeline = FeaturePipeline(feature_cols)
        X = pipeline.fit_transform(synthetic_df)

        assert isinstance(X, np.ndarray)
        assert X.shape == (10, len(feature_cols))

    def test_imputes_nan_with_median(
        self, synthetic_df: pd.DataFrame, feature_cols: list[str]
    ) -> None:
        """NaN values in mopac columns are imputed (no NaN in output)."""
        from grimperium.ml.features import FeaturePipeline

        # Verify NaNs exist in input
        assert (
            synthetic_df[["mopac_homo_ev", "mopac_lumo_ev", "mopac_gap_ev"]]
            .isna()
            .any()
            .any()
        )

        pipeline = FeaturePipeline(feature_cols)
        X = pipeline.fit_transform(synthetic_df)

        assert not np.isnan(X).any(), "Output should have no NaN after imputation"

    def test_transform_uses_train_statistics(
        self, synthetic_df: pd.DataFrame, feature_cols: list[str]
    ) -> None:
        """transform() applies statistics learned during fit_transform (no leakage)."""
        from grimperium.ml.features import FeaturePipeline

        pipeline = FeaturePipeline(feature_cols)

        # Fit on first 8 rows (train)
        train_df = synthetic_df.iloc[:8].copy()
        pipeline.fit_transform(train_df)

        # Transform test rows (last 2 — these have NaN in mopac columns)
        test_df = synthetic_df.iloc[8:].copy()
        X_test = pipeline.transform(test_df)

        assert X_test.shape == (2, len(feature_cols))
        assert not np.isnan(X_test).any(), "Test imputation should use train medians"

    def test_raises_on_missing_column(self, synthetic_df: pd.DataFrame) -> None:
        """Raises KeyError if a required feature column is missing."""
        from grimperium.ml.features import FeaturePipeline

        pipeline = FeaturePipeline(["nheavy", "nonexistent_column"])
        with pytest.raises(KeyError):
            pipeline.fit_transform(synthetic_df)
