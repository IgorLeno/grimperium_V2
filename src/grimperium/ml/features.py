"""Feature engineering pipeline with NaN imputation."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer

from grimperium import MatrixFloat


class FeaturePipeline:
    """Extract and impute feature columns from a DataFrame.

    Uses median imputation for NaN values (mopac_homo/lumo/gap_ev have ~12.4% NaN).
    The imputer is fitted on training data only to prevent data leakage.

    Parameters:
        feature_cols: Column names to extract as features.
    """

    def __init__(self, feature_cols: list[str]) -> None:
        self.feature_cols = feature_cols
        self._imputer = SimpleImputer(strategy="median")
        self._is_fitted = False

    def fit_transform(self, df: pd.DataFrame) -> MatrixFloat:
        """Fit imputer on df and return transformed feature matrix.

        Parameters:
            df: Training DataFrame with feature columns.

        Returns:
            Feature matrix (n_samples, n_features) with no NaN.

        Raises:
            KeyError: If a required feature column is missing.
        """
        raw = df[self.feature_cols].to_numpy(dtype=np.float64)
        result: MatrixFloat = self._imputer.fit_transform(raw)
        self._is_fitted = True
        return result

    def transform(self, df: pd.DataFrame) -> MatrixFloat:
        """Transform df using statistics learned during fit_transform.

        Parameters:
            df: Test DataFrame with feature columns.

        Returns:
            Feature matrix (n_samples, n_features) with no NaN.

        Raises:
            RuntimeError: If called before fit_transform.
            KeyError: If a required feature column is missing.
        """
        if not self._is_fitted:
            msg = "FeaturePipeline not fitted. Call fit_transform() first."
            raise RuntimeError(msg)
        raw = df[self.feature_cols].to_numpy(dtype=np.float64)
        result: MatrixFloat = self._imputer.transform(raw)
        return result
