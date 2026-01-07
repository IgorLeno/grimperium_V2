"""
Unit tests for ML models.

Tests cover:
    - BaseModel interface
    - KernelRidgeRegressor
    - XGBoostRegressor
    - DeltaLearningEnsemble

"""

import numpy as np
import pytest

from grimperium.models.base import BaseModel
from grimperium.models.delta_ensemble import DeltaLearningEnsemble
from grimperium.models.kernel_ridge import KernelRidgeRegressor
from grimperium.models.xgboost_model import XGBoostRegressor


class TestBaseModel:
    """Tests for BaseModel abstract class."""

    def test_cannot_instantiate(self):
        """Test that BaseModel cannot be instantiated directly."""
        with pytest.raises(TypeError):
            BaseModel()

    def test_check_is_fitted_raises(self):
        """Test that _check_is_fitted raises when not fitted."""

        class ConcreteModel(BaseModel):
            def fit(self, X, y, sample_weight=None):
                pass

            def predict(self, X):
                pass

        model = ConcreteModel()
        with pytest.raises(ValueError, match="not fitted"):
            model._check_is_fitted()


class TestKernelRidgeRegressor:
    """Tests for KernelRidgeRegressor."""

    def test_init_default(self):
        """Test default initialization."""
        krr = KernelRidgeRegressor()

        assert krr.alpha == 1.0
        assert krr.kernel == "rbf"
        assert krr.gamma is None
        assert krr.is_fitted is False

    def test_init_custom_params(self):
        """Test initialization with custom parameters."""
        krr = KernelRidgeRegressor(alpha=0.5, kernel="laplacian", gamma=0.1)

        assert krr.alpha == 0.5
        assert krr.kernel == "laplacian"
        assert krr.gamma == 0.1

    def test_get_params(self):
        """Test get_params method."""
        krr = KernelRidgeRegressor(alpha=0.5, kernel="rbf")
        params = krr.get_params()

        assert params["alpha"] == 0.5
        assert params["kernel"] == "rbf"

    def test_set_params(self):
        """Test set_params method."""
        krr = KernelRidgeRegressor()
        krr.set_params(alpha=0.1, kernel="polynomial")

        assert krr.alpha == 0.1
        assert krr.kernel == "polynomial"

    def test_repr(self):
        """Test string representation."""
        krr = KernelRidgeRegressor(alpha=0.5)
        repr_str = repr(krr)

        assert "KernelRidgeRegressor" in repr_str
        assert "alpha=0.5" in repr_str
        assert "rbf" in repr_str
        assert "not fitted" in repr_str

    @pytest.mark.skip(reason="Fit not implemented yet")
    def test_fit(self, train_test_split):
        """Test fit method."""
        X_train, _, y_train, _ = train_test_split
        krr = KernelRidgeRegressor()
        krr.fit(X_train, y_train)

        assert krr.is_fitted is True

    @pytest.mark.skip(reason="Predict not implemented yet")
    def test_predict_before_fit_raises(self):
        """Test that predict before fit raises error."""
        krr = KernelRidgeRegressor()
        X = np.random.randn(5, 10)

        with pytest.raises(ValueError, match="not fitted"):
            krr.predict(X)


class TestXGBoostRegressor:
    """Tests for XGBoostRegressor."""

    def test_init_default(self):
        """Test default initialization."""
        xgb = XGBoostRegressor()

        assert xgb.n_estimators == 100
        assert xgb.max_depth == 6
        assert xgb.learning_rate == 0.1
        assert xgb.is_fitted is False

    def test_init_custom_params(self):
        """Test initialization with custom parameters."""
        xgb = XGBoostRegressor(
            n_estimators=200,
            max_depth=8,
            learning_rate=0.05,
        )

        assert xgb.n_estimators == 200
        assert xgb.max_depth == 8
        assert xgb.learning_rate == 0.05

    def test_get_params(self):
        """Test get_params method."""
        xgb = XGBoostRegressor(n_estimators=50)
        params = xgb.get_params()

        assert params["n_estimators"] == 50
        assert "max_depth" in params
        assert "learning_rate" in params

    def test_repr(self):
        """Test string representation."""
        xgb = XGBoostRegressor(n_estimators=50, max_depth=4)
        repr_str = repr(xgb)

        assert "XGBoostRegressor" in repr_str
        assert "n_estimators=50" in repr_str
        assert "max_depth=4" in repr_str


class TestDeltaLearningEnsemble:
    """Tests for DeltaLearningEnsemble."""

    def test_init_default(self):
        """Test default initialization."""
        ensemble = DeltaLearningEnsemble()

        assert ensemble.weights == (0.5, 0.5)
        assert isinstance(ensemble.krr, KernelRidgeRegressor)
        assert isinstance(ensemble.xgb, XGBoostRegressor)
        assert ensemble.is_fitted is False

    def test_init_custom_weights(self):
        """Test initialization with custom weights."""
        ensemble = DeltaLearningEnsemble(weights=(0.7, 0.3))
        assert ensemble.weights == (0.7, 0.3)

    def test_init_invalid_weights_raises(self):
        """Test that invalid weights raise error."""
        with pytest.raises(ValueError, match="sum to 1"):
            DeltaLearningEnsemble(weights=(0.5, 0.6))

    def test_init_with_model_params(self):
        """Test initialization with model parameters."""
        ensemble = DeltaLearningEnsemble(
            krr_params={"alpha": 0.1},
            xgb_params={"n_estimators": 50},
        )

        assert ensemble.krr.alpha == 0.1
        assert ensemble.xgb.n_estimators == 50

    def test_get_params(self):
        """Test get_params method."""
        ensemble = DeltaLearningEnsemble(weights=(0.6, 0.4))
        params = ensemble.get_params(deep=True)

        assert params["weights"] == (0.6, 0.4)
        assert "krr_params" in params
        assert "xgb_params" in params

    def test_repr(self):
        """Test string representation."""
        ensemble = DeltaLearningEnsemble(weights=(0.6, 0.4))
        repr_str = repr(ensemble)

        assert "DeltaLearningEnsemble" in repr_str
        assert "(0.6, 0.4)" in repr_str
        assert "not fitted" in repr_str

    @pytest.mark.skip(reason="Fit not implemented yet")
    def test_fit(self, train_test_split):
        """Test fit method."""
        X_train, _, y_train, _ = train_test_split
        ensemble = DeltaLearningEnsemble()
        ensemble.fit(X_train, y_train)

        assert ensemble.is_fitted is True
        assert ensemble.krr.is_fitted is True
        assert ensemble.xgb.is_fitted is True
