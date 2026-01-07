"""
Integration tests for Grimperium pipeline.

Tests cover end-to-end workflows:
    - Data loading → Feature engineering → Training → Evaluation
    - API usage patterns
    - Full delta-learning pipeline

These tests use mock data to avoid external dependencies.

"""

import numpy as np
import pandas as pd
import pytest

from tests.fixtures.mock_data import (
    generate_training_data,
    make_chemperium_with_pm7_df,
    make_small_chemperium_df,
)


class TestFullPipeline:
    """End-to-end pipeline tests."""

    @pytest.fixture
    def training_data(self):
        """Generate mock training data."""
        return generate_training_data(n_samples=100, test_size=0.2)

    @pytest.mark.skip(reason="Full pipeline not implemented yet")
    def test_full_pipeline(self, training_data):
        """Test complete delta-learning pipeline."""
        from grimperium import GrimperiumAPI

        api = GrimperiumAPI()

        # Load data
        # api.load_data(...)

        # Compute features
        # api.compute_features()

        # Train
        # api.train()

        # Evaluate
        # metrics = api.evaluate()

        # Assert improvement over baseline
        # assert metrics["rmse"] < baseline_rmse

    @pytest.mark.skip(reason="API not fully implemented yet")
    def test_api_predict(self, training_data):
        """Test API prediction."""
        from grimperium import GrimperiumAPI

        api = GrimperiumAPI()
        # Predictions should match expected shape
        pass


class TestDataPipeline:
    """Tests for data loading and fusion pipeline."""

    def test_load_and_fuse(self, tmp_path):
        """Test data loading and fusion end-to-end."""
        from grimperium.data import ChemperiumLoader, DataFusion

        # Create mock data files
        chem_df = make_small_chemperium_df(n=50)
        pm7_df = pd.DataFrame(
            {
                "smiles": chem_df["smiles"].tolist(),
                "H298_pm7": chem_df["H298_cbs"]
                + np.random.normal(3.0, 1.5, len(chem_df)),
            }
        )

        # Save to files
        chem_path = tmp_path / "chemperium.csv"
        pm7_path = tmp_path / "pm7.csv"
        chem_df.to_csv(chem_path, index=False)
        pm7_df.to_csv(pm7_path, index=False)

        # Load Chemperium
        loader = ChemperiumLoader()
        loaded_df = loader.load(chem_path)
        assert len(loaded_df) == 50

        # Load PM7 manually (semiempirical handler not implemented yet)
        pm7_loaded = pd.read_csv(pm7_path)

        # Fuse datasets
        fusion = DataFusion()
        merged = fusion.merge(loaded_df, pm7_loaded, on="smiles")
        assert "H298_cbs" in merged.columns
        assert "H298_pm7" in merged.columns

        # Compute deltas
        result = fusion.compute_deltas(merged)
        assert "delta_pm7" in result.columns

        # Analyze deltas
        stats = fusion.analyze_deltas(result)
        assert "mean" in stats
        assert abs(stats["mean"]) < 10  # Realistic delta range

    def test_loader_to_task_view(self, tmp_path):
        """Test loading data and creating task views."""
        from grimperium.data import ChemperiumLoader, DataFusion

        # Create and save mock data
        df = make_small_chemperium_df(n=30)
        csv_path = tmp_path / "test_data.csv"
        df.to_csv(csv_path, index=False)

        # Load data
        loader = ChemperiumLoader()
        loaded_df = loader.load(csv_path)

        # Create task views
        fusion = DataFusion()

        # Enthalpy task
        X_enth, y_enth = fusion.select_task_view(loaded_df, task="enthalpy")
        assert len(X_enth) == 30
        assert len(y_enth) == 30
        assert y_enth.name == "H298_cbs"

        # Entropy task
        X_ent, y_ent = fusion.select_task_view(loaded_df, task="entropy")
        assert len(X_ent) == 30
        assert y_ent.name == "S298"

        # Heat capacity task (multioutput)
        X_cp, Y_cp = fusion.select_task_view(loaded_df, task="heat_capacity")
        assert len(X_cp) == 30
        assert Y_cp.shape[1] == 45  # 45 Cp columns

    def test_full_data_pipeline(self, tmp_path):
        """Test full data pipeline: load -> filter -> split -> task view."""
        from grimperium.data import ChemperiumLoader, DataFusion

        # Create mock data with varied nheavy values
        df = make_small_chemperium_df(n=100)
        csv_path = tmp_path / "full_pipeline.csv"
        df.to_csv(csv_path, index=False)

        # Load with filter
        loader = ChemperiumLoader()
        filtered_df = loader.load(csv_path, max_nheavy=10)
        assert all(filtered_df["nheavy"] <= 10)

        # 3-way split
        train_df, val_df, test_df = loader.train_val_test_split(
            filtered_df, test_size=0.2, val_size=0.1, random_state=42
        )
        assert len(train_df) + len(val_df) + len(test_df) == len(filtered_df)

        # Create task view for training
        fusion = DataFusion()
        X_train, y_train = fusion.select_task_view(train_df, task="enthalpy")
        X_val, y_val = fusion.select_task_view(val_df, task="enthalpy")
        X_test, y_test = fusion.select_task_view(test_df, task="enthalpy")

        # Verify shapes are consistent
        assert X_train.shape[1] == X_val.shape[1] == X_test.shape[1]

    def test_delta_learning_data_preparation(self):
        """Test preparing data for delta learning."""
        from grimperium.data import DataFusion

        # Create mock data with PM7
        df = make_chemperium_with_pm7_df(n=50)

        # Create fusion and compute deltas
        fusion = DataFusion()
        fusion.merged_data = df
        result = fusion.compute_deltas()

        # Verify delta values
        expected_delta = df["H298_cbs"] - df["H298_pm7"]
        np.testing.assert_array_almost_equal(
            result["delta_pm7"].values, expected_delta.values
        )

        # Get training data
        features, deltas = fusion.get_training_data()
        assert len(features) == 50
        assert len(deltas) == 50


class TestModelPipeline:
    """Tests for model training pipeline."""

    @pytest.fixture
    def training_data(self):
        """Generate mock training data."""
        return generate_training_data(n_samples=100, test_size=0.2)

    @pytest.mark.skip(reason="Model pipeline not implemented yet")
    def test_ensemble_training(self, training_data):
        """Test ensemble model training."""
        from grimperium.models import DeltaLearningEnsemble

        X_train = training_data["X_train"]
        y_train = training_data["y_train"]

        ensemble = DeltaLearningEnsemble()
        # ensemble.fit(X_train, y_train)

        # assert ensemble.is_fitted

    @pytest.mark.skip(reason="Model pipeline not implemented yet")
    def test_model_comparison(self, training_data):
        """Test comparing KRR vs XGB vs Ensemble."""

        # Train each model and compare metrics
        pass


class TestFeaturePipeline:
    """Tests for feature engineering pipeline."""

    @pytest.mark.skip(reason="Feature pipeline not implemented yet")
    def test_feature_engineering(self, sample_smiles):
        """Test feature engineering from SMILES."""
        from grimperium.utils import FeatureEngineer

        fe = FeatureEngineer(morgan_bits=256)
        # features = fe.fit_transform(sample_smiles)

        # assert features.shape[0] == len(sample_smiles)
        # assert features.shape[1] > 0

    @pytest.mark.skip(reason="Feature pipeline not implemented yet")
    def test_morgan_fingerprints(self, sample_smiles):
        """Test Morgan fingerprint computation."""

        # fps = compute_morgan_fingerprints(sample_smiles)
        # assert fps.shape == (len(sample_smiles), 256)


class TestEvaluationPipeline:
    """Tests for evaluation and metrics."""

    @pytest.mark.skip(reason="Evaluation pipeline not implemented yet")
    def test_metrics_computation(self):
        """Test metrics computation."""

        y_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y_pred = np.array([1.1, 2.0, 2.9, 4.1, 5.0])

        # rmse_val = rmse(y_true, y_pred)
        # mae_val = mae(y_true, y_pred)
        # r2_val = r2_score(y_true, y_pred)

        # assert rmse_val < 0.1
        # assert mae_val < 0.1
        # assert r2_val > 0.99

    @pytest.mark.skip(reason="Evaluation pipeline not implemented yet")
    def test_delta_improvement(self, training_data):
        """Test that delta-learning improves over PM7 baseline."""
        training_data = generate_training_data(n_samples=100, test_size=0.2)
        h298_cbs = training_data["h298_cbs_test"]

        # Baseline: PM7 raw error
        # pm7_error = np.abs(h298_pm7 - h298_cbs).mean()

        # Delta-corrected error should be lower
        # delta_pred = model.predict(X_test)
        # corrected = h298_pm7 + delta_pred
        # corrected_error = np.abs(corrected - h298_cbs).mean()

        # assert corrected_error < pm7_error
