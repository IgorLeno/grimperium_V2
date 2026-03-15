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

    Args:
        feature_cols: Column names to extract as features.

    Example:
        >>> import pandas as pd
        >>> df_train = pd.DataFrame(
        ...     {"nheavy": [1, 2], "mopac_homo_ev": [-10.5, None]}
        ... )
        >>> pipeline = FeaturePipeline(feature_cols=["nheavy", "mopac_homo_ev"])
        >>> pipeline.fit(df_train)  # Fit somente no treino para evitar leakage
        FeaturePipeline(...)
        >>> X_train = pipeline.transform(df_train)
        >>> X_train.shape
        (2, 2)
    """

    def __init__(self, feature_cols: list[str]) -> None:
        self.feature_cols: list[str] = list(feature_cols)
        self._imputer: SimpleImputer = SimpleImputer(strategy="median")
        self._is_fitted: bool = False
        self.n_features_in_: int = 0
        self.n_samples_seen_: int = 0
        self.feature_names_in_: tuple[str, ...] = ()

    def _extract_raw(self, df: pd.DataFrame) -> MatrixFloat:
        """Extract configured feature columns as float64 matrix."""
        raw: MatrixFloat = df[self.feature_cols].to_numpy(dtype=np.float64)
        return raw

    def fit(self, df: pd.DataFrame) -> FeaturePipeline:
        """Fit internal transforms using training data only.

        Args:
            df: Training DataFrame with all configured feature columns.

        Returns:
            FeaturePipeline: Fitted pipeline instance.

        Raises:
            KeyError: If a required feature column is missing.
        """
        raw = self._extract_raw(df)
        self._imputer.fit(raw)
        self._is_fitted = True
        self.n_features_in_ = int(raw.shape[1])
        self.n_samples_seen_ = int(raw.shape[0])
        self.feature_names_in_ = tuple(self.feature_cols)
        return self

    def train(self, df: pd.DataFrame) -> MatrixFloat:
        """Fit pipeline state and return transformed training matrix.

        Args:
            df: Training DataFrame with feature columns.

        Returns:
            Feature matrix (n_samples, n_features) with imputed values.
        """
        # TODO(spec-gap): `FeaturePipeline.train` foi truncado no prompt original
        # da Fase D; esta reconstrução preserva fit + metadados + transform.
        self.fit(df)
        return self.transform(df)

    def fit_transform(self, df: pd.DataFrame) -> MatrixFloat:
        """Fit imputer on df and return transformed feature matrix.

        Args:
            df: Training DataFrame with feature columns.

        Returns:
            Feature matrix (n_samples, n_features) with no NaN.

        Raises:
            KeyError: If a required feature column is missing.
        """
        return self.train(df)

    def transform(self, df: pd.DataFrame) -> MatrixFloat:
        """Transform df using statistics learned during fit_transform.

        Args:
            df: Test DataFrame with feature columns.

        Returns:
            Feature matrix (n_samples, n_features) with no NaN.

        Raises:
            RuntimeError: If called before fit_transform.
            KeyError: If a required feature column is missing.
        """
        if not self._is_fitted:
            msg = "FeaturePipeline not fitted. Call fit(), train(), or fit_transform() first."
            raise RuntimeError(msg)
        raw = self._extract_raw(df)
        result: MatrixFloat = self._imputer.transform(raw)
        return result
