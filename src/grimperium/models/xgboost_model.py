from typing import Any

import xgboost as xgb
from sklearn.preprocessing import StandardScaler

from grimperium import MatrixFloat
from grimperium.models.base import BaseModel


class XGBoostRegressor(BaseModel):
    """
    XGBoost regressor with encapsulated StandardScaler.

    Same interface as KernelRidgeRegressor for consistency.
    Week 1 defaults: max_depth=5, learning_rate=0.1, early_stopping=10
    """

    def __init__(
        self,
        n_estimators: int = 100,
        max_depth: int = 6,
        learning_rate: float = 0.1,
        early_stopping_rounds: int = 10,
        random_state: int = 42,
    ):
        """
        Initialize XGBoost regressor.

        Parameters:
            n_estimators: Number of boosting rounds
            max_depth: Maximum tree depth (5 is shallow, prevents overfitting)
            learning_rate: Learning rate for boosting (0.1 is standard)
            early_stopping_rounds: Stop if no improvement for N rounds
            random_state: For reproducibility
        """
        super().__init__()
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.early_stopping_rounds = early_stopping_rounds
        self.random_state = random_state

        self.scaler = StandardScaler()  # ← Encapsulated!
        self._model: xgb.XGBRegressor | None = None
        self.is_fitted = False

    def fit(
        self,
        X: MatrixFloat,
        y: MatrixFloat,
        sample_weight: MatrixFloat | None = None,
        eval_set: tuple[MatrixFloat, MatrixFloat] | None = None,
        **kwargs: Any,
    ) -> "XGBoostRegressor":
        """
        Train XGBoost on scaled features.

        Parameters:
            X: Training features (n_samples, n_features)
            y: Training target (n_samples,)
            sample_weight: Optional sample weights (for BaseModel compatibility)
            eval_set: Optional validation set for early stopping
                     Format: (X_val, y_val) or [(X_val, y_val)]

        Returns:
            self (for method chaining)
        """
        _ = sample_weight
        _ = kwargs

        # Scale features
        X_scaled = self.scaler.fit_transform(X)

        # Create model
        model: xgb.XGBRegressor = xgb.XGBRegressor(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=self.learning_rate,
            random_state=self.random_state,
            verbosity=0,
        )

        # Prepare fit parameters
        fit_params: dict[str, Any] = {}
        if eval_set is not None:
            # Handle both (X, y) and [(X, y)] formats
            if isinstance(eval_set, tuple) and len(eval_set) == 2:
                X_eval_scaled = self.scaler.transform(eval_set[0])
                fit_params["eval_set"] = [(X_eval_scaled, eval_set[1])]
            else:
                X_eval_scaled = self.scaler.transform(eval_set[0][0])
                fit_params["eval_set"] = [(X_eval_scaled, eval_set[0][1])]

        # Train
        model.fit(X_scaled, y, **fit_params)
        self._model = model
        self.is_fitted = True
        return self

    def predict(self, X: MatrixFloat) -> MatrixFloat:
        """
        Predict using trained model.

        Parameters:
            X: Features to predict on (n_samples, n_features)

        Returns:
            Predictions (n_samples,)
        """
        if not self.is_fitted:
            raise ValueError("Model not fitted. Call fit() first.")

        if self._model is None:
            raise ValueError("Model not fitted. Call fit() first.")
        X_scaled = self.scaler.transform(X)
        predictions: MatrixFloat = self._model.predict(X_scaled)
        return predictions

    def feature_importances(self) -> MatrixFloat:
        """Return feature importances from XGBoost trees."""
        if not self.is_fitted:
            raise ValueError("Model not fitted. Call fit() first.")
        if self._model is None:
            raise ValueError("Model not fitted. Call fit() first.")
        importances: MatrixFloat = self._model.feature_importances_
        return importances

    def get_params(self, deep: bool = True) -> dict[str, Any]:
        """Get model parameters (sklearn-compatible)."""
        _ = deep
        return {
            "n_estimators": self.n_estimators,
            "max_depth": self.max_depth,
            "learning_rate": self.learning_rate,
            "early_stopping_rounds": self.early_stopping_rounds,
            "random_state": self.random_state,
        }

    def set_params(self, **params: Any) -> "XGBoostRegressor":
        """Set model parameters (sklearn-compatible)."""
        for key, value in params.items():
            setattr(self, key, value)
        return self

    def __repr__(self) -> str:
        status = "fitted" if self.is_fitted else "not fitted"
        return (
            f"XGBoostRegressor(n_estimators={self.n_estimators}, "
            f"max_depth={self.max_depth}, learning_rate={self.learning_rate}, {status})"
        )
