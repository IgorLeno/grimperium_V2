"""
Unit tests for DeltaLearner.

Tests cover:
    - Initialization
    - Data loading interface
    - Delta computation
    - Training interface
    - Evaluation interface

"""

import numpy as np
import pytest

from grimperium.config import GrimperiumConfig
from grimperium.core.delta_learning import DeltaLearner


class TestDeltaLearner:
    """Tests for DeltaLearner class."""

    def test_init_default(self):
        """Test default initialization."""
        learner = DeltaLearner()

        assert isinstance(learner.config, GrimperiumConfig)
        assert learner._data is None
        assert learner._is_trained is False

    def test_init_with_config(self, default_config):
        """Test initialization with custom config."""
        learner = DeltaLearner(config=default_config)
        assert learner.config is default_config

    def test_repr(self):
        """Test string representation."""
        learner = DeltaLearner()
        repr_str = repr(learner)

        assert "DeltaLearner" in repr_str
        assert "n_samples=0" in repr_str
        assert "not trained" in repr_str

    @pytest.mark.skip(reason="Load data not implemented yet")
    def test_load_data(self, tmp_path, mock_chemperium_df, mock_pm7_df):
        """Test data loading."""
        # Save mock data
        chem_path = tmp_path / "chemperium.csv"
        pm7_path = tmp_path / "pm7.csv"
        mock_chemperium_df.to_csv(chem_path, index=False)
        mock_pm7_df.to_csv(pm7_path, index=False)

        learner = DeltaLearner()
        learner.load_data(chem_path, pm7_path)

        assert learner._data is not None

    @pytest.mark.skip(reason="Compute features not implemented yet")
    def test_compute_features(self, mock_merged_df):
        """Test feature computation."""
        learner = DeltaLearner()
        learner._data = mock_merged_df
        learner.compute_features()

        assert learner._features is not None
        assert learner._features.shape[0] == len(mock_merged_df)

    @pytest.mark.skip(reason="Train not implemented yet")
    def test_train(self, mock_features, mock_targets):
        """Test training."""
        learner = DeltaLearner()
        learner._features = mock_features
        learner._targets = mock_targets
        learner.train()

        assert learner._is_trained is True
        assert learner._model is not None


class TestDeltaLearnerMetrics:
    """Tests for DeltaLearner evaluation."""

    @pytest.mark.skip(reason="Evaluate not implemented yet")
    def test_evaluate_returns_metrics(self, mock_features, mock_targets):
        """Test that evaluate returns expected metrics."""
        learner = DeltaLearner()
        learner._features = mock_features
        learner._targets = mock_targets
        learner._is_trained = True

        metrics = learner.evaluate()

        assert "rmse" in metrics
        assert "mae" in metrics
        assert "r2" in metrics

    @pytest.mark.skip(reason="Cross validate not implemented yet")
    def test_cross_validate(self, mock_features, mock_targets):
        """Test cross-validation."""
        learner = DeltaLearner()
        learner._features = mock_features
        learner._targets = mock_targets

        cv_results = learner.cross_validate(n_folds=3)

        assert "rmse" in cv_results
        assert len(cv_results["rmse"]) == 3


class TestDeltaLearningConcept:
    """Tests documenting the delta-learning concept."""

    def test_delta_formula(self, mock_merged_df):
        """
        Document the delta-learning formula:
        delta = H298_CBS - H298_PM7

        The ML model learns to predict this delta.
        Final prediction: H298_CBS ≈ H298_PM7 + delta_ML
        """
        delta = mock_merged_df["H298_cbs"] - mock_merged_df["H298_pm7"]

        # Delta represents the error in PM7 relative to CBS
        assert delta is not None
        assert len(delta) == len(mock_merged_df)

    def test_prediction_reconstruction(self, mock_merged_df):
        """
        Test that original CBS can be reconstructed from PM7 + delta.

        H298_CBS = H298_PM7 + delta
        """
        h298_cbs = mock_merged_df["H298_cbs"]
        h298_pm7 = mock_merged_df["H298_pm7"]
        delta = h298_cbs - h298_pm7

        reconstructed = h298_pm7 + delta
        np.testing.assert_array_almost_equal(h298_cbs, reconstructed)

    def test_delta_learning_reduces_error(self):
        """
        Document expected behavior: delta-learning should reduce error.

        - PM7 raw error vs CBS: ~3-5 kcal/mol typical
        - PM7 + delta_ML error vs CBS: ~0.5-2 kcal/mol target
        """
        # This test documents the expected improvement
        pm7_typical_mae = 4.0  # kcal/mol
        target_mae = 1.5  # kcal/mol with ML correction

        expected_improvement = pm7_typical_mae - target_mae
        assert expected_improvement > 2  # At least 2 kcal/mol improvement
