"""
High-level API for Grimperium.

This module provides a simple, user-friendly interface for:
    - Training delta-learning models
    - Making predictions on new molecules
    - Evaluating model performance

Example:
    >>> from grimperium import GrimperiumAPI
    >>> api = GrimperiumAPI()
    >>> api.load_data("chemperium.csv")
    >>> api.train()
    >>> predictions = api.predict(["CCO", "CC(=O)O"])

"""

from pathlib import Path

import pandas as pd

from grimperium import MatrixFloat
from grimperium.config import GrimperiumConfig


class GrimperiumAPI:
    """
    High-level API for Grimperium delta-learning framework.

    This class orchestrates the entire workflow:
        1. Data loading (Chemperium + PM7)
        2. Feature engineering (Morgan FP + RDKit + tabular)
        3. Model training (KRR + XGBoost ensemble)
        4. Prediction and evaluation

    Attributes:
        config: Configuration object
        is_fitted: Whether the model has been trained

    Example:
        >>> api = GrimperiumAPI()
        >>> api.load_data("chemperium.csv")
        >>> api.train()
        >>> predictions = api.predict(["CCO", "CC(=O)O"])
        >>> metrics = api.evaluate()

    """

    def __init__(self, config: GrimperiumConfig | None = None) -> None:
        """
        Initialize GrimperiumAPI.

        Args:
            config: Configuration object. If None, uses defaults.

        """
        self.config = config or GrimperiumConfig()
        self.is_fitted = False
        self._data: pd.DataFrame | None = None
        self._features: MatrixFloat | None = None
        self._targets: MatrixFloat | None = None
        self._model = None

    def load_data(
        self,
        chemperium_path: str | Path,
        pm7_path: str | Path | None = None,
    ) -> "GrimperiumAPI":
        """
        Load and fuse Chemperium and PM7 datasets.

        Args:
            chemperium_path: Path to Chemperium CSV/parquet
            pm7_path: Optional path to precomputed PM7 data

        Returns:
            self for method chaining

        Raises:
            FileNotFoundError: If data files not found
            ValueError: If data validation fails

        """
        raise NotImplementedError("Will be implemented in Batch 2")

    def compute_features(self) -> "GrimperiumAPI":
        """
        Compute features from loaded data.

        Generates:
            - Morgan Fingerprints (256 bits by default)
            - RDKit descriptors (MolWt, TPSA, LogP, etc.)
            - Tabular features (nheavy, charge, multiplicity)

        Returns:
            self for method chaining

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def train(
        self,
        X: MatrixFloat | None = None,
        y: MatrixFloat | None = None,
    ) -> "GrimperiumAPI":
        """
        Train the delta-learning ensemble model.

        Args:
            X: Optional feature matrix (uses computed features if None)
            y: Optional target deltas (uses computed deltas if None)

        Returns:
            self for method chaining

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def predict(
        self,
        smiles: list[str] | None = None,
        X: MatrixFloat | None = None,
        h298_pm7: MatrixFloat | None = None,
    ) -> MatrixFloat:
        """
        Predict H298_CBS for new molecules.

        Args:
            smiles: List of SMILES strings
            X: Optional precomputed features
            h298_pm7: Optional precomputed PM7 values

        Returns:
            Predicted H298_CBS values

        Raises:
            ValueError: If model not fitted
            ValueError: If neither smiles nor X provided

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def evaluate(
        self,
        X_test: MatrixFloat | None = None,
        y_test: MatrixFloat | None = None,
    ) -> dict[str, float]:
        """
        Evaluate model performance.

        Args:
            X_test: Test features (uses hold-out set if None)
            y_test: Test targets (uses hold-out set if None)

        Returns:
            Dictionary with metrics: RMSE, MAE, R²

        """
        raise NotImplementedError("Will be implemented in Batch 4")

    def cross_validate(self, n_folds: int = 5) -> dict[str, list[float]]:
        """
        Perform k-fold cross-validation.

        Args:
            n_folds: Number of CV folds

        Returns:
            Dictionary with per-fold metrics

        """
        raise NotImplementedError("Will be implemented in Batch 4")

    def save(self, path: str | Path) -> None:
        """Save trained model to disk."""
        raise NotImplementedError("Will be implemented in Batch 5")

    def load(self, path: str | Path) -> "GrimperiumAPI":
        """Load trained model from disk."""
        raise NotImplementedError("Will be implemented in Batch 5")
