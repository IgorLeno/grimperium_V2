"""
Kernel Ridge Regression model for Grimperium.

This module provides KernelRidgeRegressor, a wrapper around
scikit-learn's KernelRidge with optimizations for delta-learning:
    - RBF kernel (default) for smooth corrections
    - Automatic hyperparameter tuning
    - Integration with Grimperium pipeline

KRR is particularly good for learning smooth, continuous corrections
to semiempirical energies.

Example:
    >>> from grimperium.models import KernelRidgeRegressor
    >>> model = KernelRidgeRegressor(alpha=1.0, kernel="rbf")
    >>> model.fit(X_train, y_train)
    >>> predictions = model.predict(X_test)

"""

from typing import Any, Optional, Union

import numpy as np
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler

from grimperium import DictStrAny, MatrixFloat
from grimperium.models.base import BaseModel


class KernelRidgeRegressor(BaseModel):
    """
    Kernel Ridge Regression for delta-learning.

    Wraps scikit-learn's KernelRidge with optimizations for
    molecular property prediction.

    Attributes:
        alpha: Regularization strength
        kernel: Kernel type ('rbf', 'laplacian', 'polynomial')
        gamma: Kernel coefficient (None for auto)
        degree: Polynomial degree (only for polynomial kernel)

    Example:
        >>> krr = KernelRidgeRegressor(alpha=0.1)
        >>> krr.fit(X_train, y_delta)
        >>> delta_pred = krr.predict(X_test)
        >>> h298_pred = h298_pm7_test + delta_pred

    Notes:
        - RBF kernel works best for energy corrections
        - gamma='scale' uses 1/(n_features * X.var())
        - alpha should be tuned via cross-validation

    """

    def __init__(
        self,
        alpha: float = 1.0,
        kernel: str = "rbf",
        gamma: Optional[float] = None,
        degree: int = 3,
        coef0: float = 1.0,
    ) -> None:
        """
        Initialize KernelRidgeRegressor.

        Args:
            alpha: Regularization strength (higher = more regularization)
            kernel: Kernel type ('rbf', 'laplacian', 'polynomial', 'linear')
            gamma: Kernel coefficient (None for 1/n_features)
            degree: Polynomial degree (only for polynomial kernel)
            coef0: Independent term in polynomial kernel

        """
        super().__init__()
        self.alpha = alpha
        self.kernel = kernel
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.scaler = StandardScaler()
        self._model: Optional[KernelRidge] = None

    def fit(
        self,
        X: MatrixFloat,
        y: MatrixFloat,
        sample_weight: Optional[MatrixFloat] = None,
    ) -> "KernelRidgeRegressor":
        """
        Fit the KRR model.

        Args:
            X: Training features of shape (n_samples, n_features)
            y: Target deltas of shape (n_samples,)
            sample_weight: Optional sample weights (not used in KRR)

        Returns:
            self for method chaining

        """
        # Scale features
        X_scaled = self.scaler.fit_transform(X)

        # Store number of features
        self.n_features_in_ = X.shape[1]

        # Initialize and fit KernelRidge
        model: KernelRidge = KernelRidge(
            alpha=self.alpha,
            kernel=self.kernel,
            gamma=self.gamma,
            degree=self.degree,
            coef0=self.coef0,
        )
        model.fit(X_scaled, y)
        self._model = model

        self.is_fitted = True
        return self

    def predict(self, X: MatrixFloat) -> MatrixFloat:
        """
        Predict delta values.

        Args:
            X: Features of shape (n_samples, n_features)

        Returns:
            Predicted deltas of shape (n_samples,)

        """
        self._check_is_fitted()

        # Scale features using fitted scaler
        X_scaled = self.scaler.transform(X)

        assert self._model is not None, "Model should be fitted"
        predictions: MatrixFloat = self._model.predict(X_scaled)
        return predictions

    def tune_hyperparameters(
        self,
        X: MatrixFloat,
        y: MatrixFloat,
        param_grid: Optional[DictStrAny] = None,
        cv: int = 5,
        scoring: str = "neg_root_mean_squared_error",
    ) -> DictStrAny:
        """
        Tune hyperparameters via grid search.

        Args:
            X: Training features
            y: Target values
            param_grid: Parameter grid (uses default if None)
            cv: Number of CV folds
            scoring: Scoring metric

        Returns:
            Best parameters found

        """
        raise NotImplementedError("Will be implemented in Batch 3")

    def get_params(self, deep: bool = True) -> DictStrAny:
        """Get model parameters."""
        return {
            "alpha": self.alpha,
            "kernel": self.kernel,
            "gamma": self.gamma,
            "degree": self.degree,
            "coef0": self.coef0,
        }

    def set_params(self, **params: Any) -> "KernelRidgeRegressor":
        """Set model parameters."""
        for key, value in params.items():
            setattr(self, key, value)
        return self

    def __repr__(self) -> str:
        """String representation."""
        status = "fitted" if self.is_fitted else "not fitted"
        return (
            f"KernelRidgeRegressor(alpha={self.alpha}, "
            f"kernel='{self.kernel}', {status})"
        )
