"""
Unit tests for DataFusion.

Tests cover:
    - Merging Chemperium and PM7 datasets
    - Delta computation
    - Delta statistics
    - Training data extraction
    - Task-specific views

"""

import numpy as np
import pandas as pd
import pytest

from grimperium.data.fusion import DataFusion


class TestDataFusion:
    """Tests for DataFusion class."""

    def test_init_default(self):
        """Test default initialization."""
        fusion = DataFusion()

        assert fusion.target_column == "H298_cbs"
        assert fusion.semiempirical_column == "H298_pm7"
        assert fusion.delta_column == "delta_pm7"
        assert fusion.merged_data is None

    def test_init_custom_columns(self):
        """Test initialization with custom column names."""
        fusion = DataFusion(
            target_column="custom_target",
            semiempirical_column="custom_semiemp",
            delta_column="custom_delta",
        )

        assert fusion.target_column == "custom_target"
        assert fusion.semiempirical_column == "custom_semiemp"
        assert fusion.delta_column == "custom_delta"

    def test_feature_columns_defined(self):
        """Test that feature columns are defined."""
        assert "nheavy" in DataFusion.FEATURE_COLUMNS
        assert "charge" in DataFusion.FEATURE_COLUMNS
        assert "multiplicity" in DataFusion.FEATURE_COLUMNS

    def test_repr(self):
        """Test string representation."""
        fusion = DataFusion()
        repr_str = repr(fusion)

        assert "DataFusion" in repr_str
        assert "n_merged=0" in repr_str
        assert "delta_pm7" in repr_str

    def test_merge_inner(self, mock_chemperium_df, mock_pm7_df):
        """Test inner merge of datasets."""
        fusion = DataFusion()
        merged = fusion.merge(mock_chemperium_df, mock_pm7_df, on="smiles")

        assert len(merged) <= len(mock_chemperium_df)
        assert "H298_cbs" in merged.columns
        assert "H298_pm7" in merged.columns
        assert fusion.merged_data is not None

    def test_merge_stores_data(self, mock_chemperium_df, mock_pm7_df):
        """Test that merge stores data in merged_data attribute."""
        fusion = DataFusion()
        merged = fusion.merge(mock_chemperium_df, mock_pm7_df, on="smiles")

        assert fusion.merged_data is not None
        pd.testing.assert_frame_equal(fusion.merged_data, merged)

    def test_merge_empty_raises(self):
        """Test that merging with no common keys raises error."""
        fusion = DataFusion()
        df1 = pd.DataFrame({"smiles": ["A", "B"], "H298_cbs": [1.0, 2.0]})
        df2 = pd.DataFrame({"smiles": ["C", "D"], "H298_pm7": [1.0, 2.0]})

        with pytest.raises(ValueError, match="empty result"):
            fusion.merge(df1, df2, on="smiles")

    def test_compute_deltas(self, mock_merged_df):
        """Test delta computation."""
        fusion = DataFusion()
        fusion.merged_data = mock_merged_df.copy()

        result = fusion.compute_deltas()

        assert fusion.delta_column in result.columns
        expected_delta = mock_merged_df["H298_cbs"] - mock_merged_df["H298_pm7"]
        np.testing.assert_array_almost_equal(
            result[fusion.delta_column].values,
            expected_delta.values,
        )

    def test_compute_deltas_updates_merged_data(self, mock_merged_df):
        """Test that compute_deltas updates merged_data attribute."""
        fusion = DataFusion()
        fusion.merged_data = mock_merged_df.copy()

        fusion.compute_deltas()

        assert fusion.delta_column in fusion.merged_data.columns

    def test_compute_deltas_no_data_raises(self):
        """Test compute_deltas raises error when no data."""
        fusion = DataFusion()

        with pytest.raises(ValueError, match="No data available"):
            fusion.compute_deltas()

    def test_compute_deltas_missing_column_raises(self, mock_chemperium_df):
        """Test compute_deltas raises error when semiempirical column missing."""
        fusion = DataFusion()
        fusion.merged_data = mock_chemperium_df  # Has H298_cbs but not H298_pm7

        with pytest.raises(ValueError, match="not found"):
            fusion.compute_deltas()

    def test_analyze_deltas(self, mock_merged_df):
        """Test delta statistics."""
        fusion = DataFusion()
        # Add delta column
        df = mock_merged_df.copy()
        df["delta_pm7"] = df["H298_cbs"] - df["H298_pm7"]
        fusion.merged_data = df

        stats = fusion.analyze_deltas()

        assert "mean" in stats
        assert "std" in stats
        assert "min" in stats
        assert "max" in stats
        assert "median" in stats
        # Check values are floats
        assert isinstance(stats["mean"], float)
        assert isinstance(stats["std"], float)

    def test_analyze_deltas_no_delta_column_raises(self, mock_chemperium_df, mock_pm7_df):
        """Test analyze_deltas raises error when delta column missing."""
        fusion = DataFusion()
        # Create merged df WITHOUT delta column
        merged_no_delta = mock_chemperium_df.merge(mock_pm7_df, on="smiles")
        fusion.merged_data = merged_no_delta

        with pytest.raises(ValueError, match="not found"):
            fusion.analyze_deltas()

    def test_get_training_data(self, mock_merged_df):
        """Test getting training data."""
        fusion = DataFusion()
        # Add delta column
        df = mock_merged_df.copy()
        df["delta_pm7"] = df["H298_cbs"] - df["H298_pm7"]
        fusion.merged_data = df

        features, deltas = fusion.get_training_data()

        assert isinstance(features, pd.DataFrame)
        assert isinstance(deltas, np.ndarray)
        assert len(features) == len(deltas)


class TestDataFusionDeltaValues:
    """Tests for delta value correctness."""

    def test_delta_sign_convention(self, mock_merged_df):
        """Test that delta = CBS - PM7 (positive when PM7 underestimates)."""
        # Manually compute expected deltas
        expected = mock_merged_df["H298_cbs"] - mock_merged_df["H298_pm7"]

        # In our mock data, PM7 tends to underestimate binding
        # so delta should be mostly negative (CBS more negative than PM7)
        # This matches real PM7 behavior
        assert expected.mean() < 0  # PM7 gives less negative values than CBS

    def test_delta_range_realistic(self, mock_merged_df):
        """Test that delta values are in realistic range."""
        delta = mock_merged_df["H298_cbs"] - mock_merged_df["H298_pm7"]

        # Typical PM7 errors are in the range of ±10 kcal/mol
        assert delta.min() > -20
        assert delta.max() < 20
        assert abs(delta.mean()) < 10


class TestDataFusionTaskViews:
    """Tests for task-specific data views."""

    def test_select_task_view_enthalpy(self, mock_chemperium_df):
        """Test enthalpy task view."""
        fusion = DataFusion()
        X, y = fusion.select_task_view(mock_chemperium_df, task="enthalpy")

        assert isinstance(X, pd.DataFrame)
        assert isinstance(y, pd.Series)
        assert len(X) == len(y)
        assert "H298_cbs" not in X.columns  # Target should not be in features
        assert y.name == "H298_cbs"

    def test_select_task_view_entropy(self, mock_chemperium_df):
        """Test entropy task view."""
        fusion = DataFusion()
        X, y = fusion.select_task_view(mock_chemperium_df, task="entropy")

        assert isinstance(X, pd.DataFrame)
        assert isinstance(y, pd.Series)
        assert len(X) == len(y)
        assert "S298" not in X.columns
        assert y.name == "S298"

    def test_select_task_view_heat_capacity(self, mock_chemperium_df):
        """Test heat capacity task view (multioutput)."""
        fusion = DataFusion()
        X, Y = fusion.select_task_view(mock_chemperium_df, task="heat_capacity")

        assert isinstance(X, pd.DataFrame)
        assert isinstance(Y, pd.DataFrame)  # Multi-output returns DataFrame
        assert len(X) == len(Y)
        # Y should have cp_* columns
        cp_cols = [c for c in Y.columns if c.startswith("cp_")]
        assert len(cp_cols) > 0

    def test_select_task_view_invalid_task(self, mock_chemperium_df):
        """Test that invalid task raises error."""
        fusion = DataFusion()

        with pytest.raises(ValueError, match="Invalid task"):
            fusion.select_task_view(mock_chemperium_df, task="invalid")

    def test_select_task_view_missing_column_raises(self):
        """Test that missing target column raises error."""
        fusion = DataFusion()
        df = pd.DataFrame({
            "smiles": ["C", "CC"],
            "nheavy": [1, 2],
            "charge": [0, 0],
            "multiplicity": [1, 1],
        })  # Missing H298_cbs

        with pytest.raises(ValueError, match="required for enthalpy"):
            fusion.select_task_view(df, task="enthalpy")

    def test_select_task_view_no_cp_columns_raises(self):
        """Test heat_capacity task fails without cp columns."""
        fusion = DataFusion()
        df = pd.DataFrame({
            "smiles": ["C", "CC"],
            "nheavy": [1, 2],
            "charge": [0, 0],
            "multiplicity": [1, 1],
            "H298_cbs": [1.0, 2.0],
        })  # No cp_* columns

        with pytest.raises(ValueError, match="No cp_\\* columns"):
            fusion.select_task_view(df, task="heat_capacity")

    def test_select_task_view_features_correct(self, mock_chemperium_df):
        """Test that features include expected columns."""
        fusion = DataFusion()
        X, _ = fusion.select_task_view(mock_chemperium_df, task="enthalpy")

        # Should include tabular features
        expected_features = ["nheavy", "charge", "multiplicity"]
        for col in expected_features:
            if col in mock_chemperium_df.columns:
                assert col in X.columns
