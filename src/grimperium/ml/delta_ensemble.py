from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from grimperium.models.kernel_ridge import KernelRidgeRegressor
from grimperium.models.xgboost_model import XGBoostRegressor

DictStrAny = dict[str, Any]
MatrixFloat = NDArray[np.floating[Any]]


class DeltaLearningEnsemble:
    """
    Simple ensemble combining KRR and XGBoost.

    Week 1 Strategy: Weighted average (w_krr=0.5, w_xgb=0.5)
    - Simple, fast validation
    - If gate passes, keep simple
    - If fails, optimize weights in Week 2
    """

    def __init__(
        self,
        weights: tuple[float, float] = (0.5, 0.5),
        krr_params: DictStrAny | None = None,
        xgb_params: DictStrAny | None = None,
        *,
        w_krr: float | None = None,
        w_xgb: float | None = None,
    ):
        """
        Initialize ensemble with two models.

        Parameters:
            w_krr: Weight for KRR predictions (0 to 1)
            w_xgb: Weight for XGB predictions (0 to 1)
                   Note: w_krr + w_xgb should equal 1.0
            krr_params: Dict of KRR hyperparameters (alpha, gamma)
            xgb_params: Dict of XGB hyperparameters (n_estimators, max_depth, etc)
        """
        if (w_krr is not None or w_xgb is not None) and weights != (0.5, 0.5):
            raise ValueError("Use apenas `weights` ou `w_krr`/`w_xgb`, não ambos.")

        if w_krr is not None or w_xgb is not None:
            if w_krr is None or w_xgb is None:
                raise ValueError("Informe ambos `w_krr` e `w_xgb`.")
            weights = (float(w_krr), float(w_xgb))

        if len(weights) != 2:
            raise ValueError("`weights` deve ter exatamente 2 valores (krr, xgb).")
        if not np.isclose(weights[0] + weights[1], 1.0):
            raise ValueError("`weights` deve sum to 1.")

        self.weights: tuple[float, float] = (float(weights[0]), float(weights[1]))
        self.w_krr, self.w_xgb = self.weights

        # Week 1 defaults
        krr_defaults = {"alpha": 1.0, "gamma": 0.1}
        xgb_defaults = {
            "n_estimators": 100,
            "max_depth": 6,
            "learning_rate": 0.1,
            "early_stopping_rounds": 10,
        }

        # Merge user params with defaults
        self.krr_params: dict[str, Any] = {**krr_defaults, **(krr_params or {})}
        self.xgb_params: dict[str, Any] = {**xgb_defaults, **(xgb_params or {})}

        # Create models
        self.krr = KernelRidgeRegressor(**self.krr_params)
        self.xgb = XGBoostRegressor(**self.xgb_params)
        self.is_fitted = False

    def fit(self, X: MatrixFloat, y: MatrixFloat) -> DeltaLearningEnsemble:
        """
        Train both models on the same target (delta).

        Parameters:
            X: Features (n_samples, n_features)
            y: Target (y_delta) (n_samples,)

        Returns:
            self
        """
        self.krr.fit(X, y)
        self.xgb.fit(X, y)
        self.is_fitted = True
        return self

    def predict(self, X: MatrixFloat) -> MatrixFloat:
        """
        Predict using weighted average ensemble.

        Formula: y_pred = w_krr * krr_pred + w_xgb * xgb_pred

        Parameters:
            X: Features (n_samples, n_features)

        Returns:
            Ensemble predictions (n_samples,)
        """
        if not self.is_fitted:
            raise ValueError("Ensemble not fitted. Call fit() first.")

        krr_pred: MatrixFloat = self.krr.predict(X)
        xgb_pred: MatrixFloat = self.xgb.predict(X)

        ensemble_pred: MatrixFloat = self.w_krr * krr_pred + self.w_xgb * xgb_pred
        return ensemble_pred

    def predict_individual(self, X: MatrixFloat) -> tuple[MatrixFloat, MatrixFloat]:
        """
        Return individual model predictions (for analysis/debugging).

        Parameters:
            X: Features (n_samples, n_features)

        Returns:
            (krr_predictions, xgb_predictions)
        """
        if not self.is_fitted:
            raise ValueError("Ensemble not fitted. Call fit() first.")

        return self.krr.predict(X), self.xgb.predict(X)

    def get_params(self, deep: bool = True) -> dict[str, Any]:
        """Get ensemble parameters (sklearn-compatible)."""
        params: dict[str, Any] = {
            "weights": self.weights,
            "krr_params": dict(self.krr_params),
            "xgb_params": dict(self.xgb_params),
        }

        if deep:
            params.update({f"krr__{k}": v for k, v in self.krr.get_params().items()})
            params.update({f"xgb__{k}": v for k, v in self.xgb.get_params().items()})

        return params

    def __repr__(self) -> str:
        status = "fitted" if self.is_fitted else "not fitted"
        return f"DeltaLearningEnsemble(weights={self.weights}, {status})"
