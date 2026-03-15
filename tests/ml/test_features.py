"""Tests for grimperium.ml.features — Feature engineering pipeline."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from grimperium.ml.features import FeaturePipeline


class TestFeaturePipeline:
    """Tests for FeaturePipeline class."""

    def test_fit_transform_returns_ndarray(
        self, synthetic_df: pd.DataFrame, feature_cols: list[str]
    ) -> None:
        """fit_transform returns ndarray with expected shape.

        Args:
            synthetic_df: Synthetic tabular dataset used as pipeline input.
            feature_cols: Feature names extracted from the DataFrame.

        Returns:
            None: Asserts output type and matrix dimensions.
        """
        pipeline = FeaturePipeline(feature_cols)
        X = pipeline.fit_transform(synthetic_df)

        assert isinstance(X, np.ndarray)
        assert X.shape == (len(synthetic_df), len(feature_cols))

    def test_imputes_nan_with_median(
        self, synthetic_df: pd.DataFrame, feature_cols: list[str]
    ) -> None:
        """Median imputation removes NaN values from transformed matrix.

        Args:
            synthetic_df: Synthetic dataset containing NaN in mopac columns.
            feature_cols: Feature names extracted from the DataFrame.

        Returns:
            None: Asserts transformed matrix has no NaN values.
        """
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
        """transform uses statistics learned from the train split only.

        Args:
            synthetic_df: Synthetic dataset split into train/test slices.
            feature_cols: Feature names extracted from the DataFrame.

        Returns:
            None: Asserts test transformation uses train-fitted imputation stats.
        """
        pipeline = FeaturePipeline(feature_cols)

        mopac_columns = ["mopac_homo_ev", "mopac_lumo_ev", "mopac_gap_ev"]
        test_slice = synthetic_df.iloc[8:][mopac_columns]
        assert test_slice.isna().any().any(), (
            "Fixture assumption violated: rows used as test split must contain "
            "NaN in at least one mopac column."
        )

        # Fit on first 8 rows (train)
        train_df = synthetic_df.iloc[:8].copy()
        pipeline.fit_transform(train_df)

        # Transform test rows (last 2 — these have NaN in mopac columns)
        test_df = synthetic_df.iloc[8:].copy()
        X_test = pipeline.transform(test_df)

        assert X_test.shape == (len(test_df), len(feature_cols))
        assert not np.isnan(X_test).any(), "Test imputation should use train medians"

    def test_raises_on_missing_column(self, synthetic_df: pd.DataFrame) -> None:
        """Raises KeyError when a required feature column is missing.

        Args:
            synthetic_df: Synthetic dataset missing no required columns initially.

        Returns:
            None: Asserts fit_transform fails with unknown column names.
        """
        pipeline = FeaturePipeline(["nheavy", "nonexistent_column"])
        with pytest.raises(KeyError):
            pipeline.fit_transform(synthetic_df)

    def test_train_persists_metadata(self, synthetic_df: pd.DataFrame) -> None:
        """train fits imputer state and stores metadata for reuse.

        Args:
            synthetic_df: Synthetic dataset used to fit the pipeline.

        Returns:
            None: Asserts train output shape and metadata persistence.
        """
        pipeline = FeaturePipeline(["nheavy", "mopac_homo_ev"])
        X_train = pipeline.train(synthetic_df)

        assert X_train.shape == (len(synthetic_df), 2)
        assert pipeline.n_features_in_ == 2
        assert pipeline.n_samples_seen_ == len(synthetic_df)
