"""Tests for model serialization and deserialization (Phase E)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest


class TestSaveModel:
    """Tests for save_model function."""

    def test_save_creates_joblib_file(
        self, trained_model_fixture: dict[str, Any], tmp_path: Path
    ) -> None:
        """Verify that save_model creates a .joblib file at the given path.

        Args:
            trained_model_fixture: Bundle dict with learner, pipeline, metrics.
            tmp_path: Temporary directory provided by pytest.
        """
        from grimperium.ml.persistence import save_model

        path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, path)
        assert path.exists()
        assert path.suffix == ".joblib"

    def test_save_creates_parent_dirs(
        self, trained_model_fixture: dict[str, Any], tmp_path: Path
    ) -> None:
        """Verify that save_model creates parent directories when missing.

        Args:
            trained_model_fixture: Bundle dict with learner, pipeline, metrics.
            tmp_path: Temporary directory provided by pytest.
        """
        from grimperium.ml.persistence import save_model

        path = tmp_path / "nested" / "deep" / "model.joblib"
        save_model(trained_model_fixture, path)
        assert path.exists()

    def test_save_includes_metadata(
        self, trained_model_fixture: dict[str, Any], tmp_path: Path
    ) -> None:
        """Verify saved bundle contains all expected keys.

        Args:
            trained_model_fixture: Bundle dict with learner, pipeline, metrics.
            tmp_path: Temporary directory provided by pytest.
        """
        import joblib

        from grimperium.ml.persistence import save_model

        path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, path)
        bundle = joblib.load(path)
        assert "learner" in bundle
        assert "pipeline" in bundle
        assert "metrics" in bundle
        assert "version" in bundle
        assert "trained_at" in bundle


class TestLoadModel:
    """Tests for load_model function."""

    def test_load_returns_learner_and_pipeline(
        self, trained_model_fixture: dict[str, Any], tmp_path: Path
    ) -> None:
        """Verify load_model returns a (DeltaLearner, FeaturePipeline) tuple.

        Args:
            trained_model_fixture: Bundle dict with learner, pipeline, metrics.
            tmp_path: Temporary directory provided by pytest.
        """
        from grimperium.core.delta_learning import DeltaLearner
        from grimperium.ml.features import FeaturePipeline
        from grimperium.ml.persistence import load_model, save_model

        path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, path)
        learner, pipeline = load_model(path)
        assert isinstance(learner, DeltaLearner)
        assert isinstance(pipeline, FeaturePipeline)

    def test_load_model_is_fitted(
        self, trained_model_fixture: dict[str, Any], tmp_path: Path
    ) -> None:
        """Verify loaded model is fitted and ready for prediction.

        Args:
            trained_model_fixture: Bundle dict with learner, pipeline, metrics.
            tmp_path: Temporary directory provided by pytest.
        """
        from grimperium.ml.persistence import load_model, save_model

        path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, path)
        learner, _ = load_model(path)
        assert learner.is_fitted

    def test_load_raises_on_missing_file(self, tmp_path: Path) -> None:
        """Verify load_model raises FileNotFoundError for non-existent file.

        Args:
            tmp_path: Temporary directory provided by pytest.
        """
        from grimperium.ml.persistence import load_model

        with pytest.raises(FileNotFoundError, match="Model file not found"):
            load_model(tmp_path / "nonexistent.joblib")

    def test_roundtrip_predictions_identical(
        self,
        trained_model_fixture: dict[str, Any],
        synthetic_csv_path: Path,
        tmp_path: Path,
    ) -> None:
        """Verify predictions before and after save/load are identical.

        Args:
            trained_model_fixture: Bundle dict with learner, pipeline, metrics.
            synthetic_csv_path: Path to synthetic CSV fixture.
            tmp_path: Temporary directory provided by pytest.
        """
        import numpy as np

        from grimperium.ml.data_loader import load_ml_data
        from grimperium.ml.persistence import load_model, save_model

        path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, path)
        learner_orig = trained_model_fixture["learner"]
        pipeline_orig = trained_model_fixture["pipeline"]
        learner_loaded, pipeline_loaded = load_model(path)

        df, _, y_pm7 = load_ml_data(synthetic_csv_path)
        X_orig = pipeline_orig.transform(df)
        X_loaded = pipeline_loaded.transform(df)
        preds_orig = learner_orig.predict(X_orig, y_pm7)
        preds_loaded = learner_loaded.predict(X_loaded, y_pm7)
        np.testing.assert_array_almost_equal(preds_orig, preds_loaded)
