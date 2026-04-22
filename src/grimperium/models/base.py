"""
Base model abstract class for Grimperium.

This module defines the BaseModel abstract class that all ML models
must inherit from, ensuring consistent interface across:
    - KernelRidgeRegressor
    - XGBoostRegressor
    - DeltaLearningEnsemble

The interface follows scikit-learn conventions with fit/predict methods.

Example:
    >>> class MyModel(BaseModel):
    ...     def fit(self, X, y):
    ...         # implementation
    ...     def predict(self, X):
    ...         # implementation

"""

from abc import ABC, abstractmethod
from typing import Any, Optional

import numpy as np

from grimperium import MatrixFloat


class BaseModel(ABC):
    """
    Abstract base class for all Grimperium ML models.

    Defines the standard interface for delta-learning models
    following scikit-learn conventions.

    Attributes:
        is_fitted: Whether model has been trained
        n_features_in_: Number of features seen during fit
        feature_names_in_: Feature names (if available)

    Methods:
        fit: Train the model
        predict: Make predictions
        score: Evaluate model performance

    """

    def __init__(self) -> None:
        """Initialize BaseModel."""
        self.is_fitted: bool = False
        self.n_features_in_: Optional[int] = None
        self.feature_names_in_: Optional[MatrixFloat] = None

    @abstractmethod
    def fit(
        self,
        X: MatrixFloat,
        y: MatrixFloat,
        sample_weight: Optional[MatrixFloat] = None,
    ) -> "BaseModel":
        """
        Fit the model to training data.

        Args:
            X: Training features of shape (n_samples, n_features)
            y: Target values of shape (n_samples,)
            sample_weight: Optional sample weights

        Returns:
            self for method chaining

        Raises:
            ValueError: If X and y have incompatible shapes

        """
        pass

    @abstractmethod
    def predict(self, X: MatrixFloat) -> MatrixFloat:
        """
        Predict target values for X.

        Args:
            X: Features of shape (n_samples, n_features)

        Returns:
            Predicted values of shape (n_samples,)

        Raises:
            ValueError: If model not fitted
            ValueError: If X has wrong number of features

        """
        pass

    def score(
        self,
        X: MatrixFloat,
        y: MatrixFloat,
        sample_weight: Optional[MatrixFloat] = None,
    ) -> float:
        """
        Compute R² score for predictions.

        Args:
            X: Test features
            y: True target values
            sample_weight: Optional sample weights

        Returns:
            R² score (higher is better, max 1.0)

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def get_params(self, deep: bool = True) -> dict[str, Any]:
        """
        Get model parameters (sklearn-compatible).

        Args:
            deep: Whether to return nested parameters

        Returns:
            Dictionary of parameters

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def set_params(self, **params: Any) -> "BaseModel":
        """
        Set model parameters (sklearn-compatible).

        Args:
            **params: Parameters to set

        Returns:
            self for method chaining

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def _check_is_fitted(self) -> None:
        """Check if model is fitted, raise if not."""
        if not self.is_fitted:
            raise ValueError(
                f"{self.__class__.__name__} is not fitted. "
                "Call fit() before predict()."
            )

    def _validate_X(self, X: MatrixFloat, reset: bool = False) -> MatrixFloat:
        """Validate input features."""
        raise NotImplementedError("Will be implemented in Batch 3")

    def __repr__(self) -> str:
        """String representation."""
        status = "fitted" if self.is_fitted else "not fitted"
        return f"{self.__class__.__name__}({status})"
