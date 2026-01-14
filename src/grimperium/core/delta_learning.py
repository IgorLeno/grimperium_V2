import numpy as np

from grimperium.config import GrimperiumConfig
from grimperium.core.metrics import compute_all_metrics
from grimperium.models.delta_ensemble import DeltaLearningEnsemble


class DeltaLearner:
    """
    Delta-Learning orchestrator.

    Core hypothesis: Learning delta (correction) is easier than learning y_cbs directly.

    Architecture:
    - Input: X (features), y_cbs (CBS results), y_pm7 (PM7 approximation)
    - Process:
        y_delta = y_cbs - y_pm7  ← EXPLICIT!
        ensemble.fit(X, y_delta)  ← Train on DELTA!
    - Prediction:
        delta_pred = ensemble.predict(X)
        y_pred = y_pm7 + delta_pred  ← Composition! EXPLICIT!
    """

    def __init__(
        self,
        w_krr: float = 0.5,
        w_xgb: float = 0.5,
        krr_params: dict | None = None,
        xgb_params: dict | None = None,
        *,
        config: GrimperiumConfig | None = None,
    ):
        """
        Initialize DeltaLearner with ensemble.

        Parameters:
            w_krr: Weight for KRR in ensemble
            w_xgb: Weight for XGB in ensemble
            krr_params: KRR hyperparameters
            xgb_params: XGB hyperparameters
            config: Configuração global do Grimperium (opcional).
        """
        self.config: GrimperiumConfig = config or GrimperiumConfig()

        self.w_krr = w_krr
        self.w_xgb = w_xgb
        self.ensemble = DeltaLearningEnsemble(
            w_krr=w_krr, w_xgb=w_xgb, krr_params=krr_params, xgb_params=xgb_params
        )

        # Estado interno esperado pelos testes/unit (mantido por compatibilidade).
        self._data: object | None = None
        self._is_trained: bool = False
        self._n_samples: int = 0

        # Alias legado usado no restante do código.
        self.is_fitted: bool = False

    def __repr__(self) -> str:
        status = "trained" if self._is_trained else "not trained"
        return f"DeltaLearner(n_samples={self._n_samples}, {status})"

    def fit(
        self, X: np.ndarray, y_cbs: np.ndarray, y_pm7: np.ndarray
    ) -> "DeltaLearner":
        """
        Fit DeltaLearner.

        EXPLICIT PROCESS:
        1. Calculate delta: y_delta = y_cbs - y_pm7
        2. Train ensemble on y_delta (not on y_cbs!)

        Parameters:
            X: Features (n_samples, n_features)
            y_cbs: CBS results (ground truth) (n_samples,)
            y_pm7: PM7 approximation (n_samples,)

        Returns:
            self

        Example:
            learner = DeltaLearner()
            learner.fit(X_train, y_cbs_train, y_pm7_train)
        """
        # STEP 1: Calculate delta EXPLICITLY
        y_delta = y_cbs - y_pm7

        # Sanity checks
        assert X.shape[0] == len(y_cbs), "X and y_cbs must have same n_samples"
        assert X.shape[0] == len(y_pm7), "X and y_pm7 must have same n_samples"

        # STEP 2: Train ensemble on y_delta
        self.ensemble.fit(X, y_delta)

        self._n_samples = int(X.shape[0])
        self._is_trained = True
        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray, y_pm7: np.ndarray) -> np.ndarray:
        """
        Predict H298_CBS using delta-learning.

        EXPLICIT PROCESS:
        1. Predict delta correction: delta_pred = ensemble.predict(X)
        2. Compose with PM7: y_pred = y_pm7 + delta_pred

        Parameters:
            X: Features (n_samples, n_features)
            y_pm7: PM7 approximation (n_samples,)

        Returns:
            Predicted H298_CBS (n_samples,)

        Example:
            y_pred = learner.predict(X_test, y_pm7_test)
        """
        if not self.is_fitted:
            raise ValueError("DeltaLearner not fitted. Call fit() first.")

        # STEP 1: Predict delta
        delta_pred: np.ndarray = self.ensemble.predict(X)

        # STEP 2: Compose with PM7 EXPLICITLY
        y_pred: np.ndarray = y_pm7 + delta_pred

        return y_pred

    def evaluate(self, X: np.ndarray, y_cbs: np.ndarray, y_pm7: np.ndarray) -> dict:
        """
        Evaluate model and return all metrics.

        Parameters:
            X: Features
            y_cbs: Ground truth CBS
            y_pm7: PM7 approximation

        Returns:
            Dict with RMSE, MAE, R², MAPE, max_error

        Example:
            metrics = learner.evaluate(X_test, y_cbs_test, y_pm7_test)
            print(f"RMSE: {metrics['rmse']:.2f}")
        """
        y_pred = self.predict(X, y_pm7)
        return compute_all_metrics(y_cbs, y_pred)
