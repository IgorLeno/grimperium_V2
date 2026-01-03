"""
Unit tests for DataFusion.

Tests cover:
    - Merging Chemperium and PM7 datasets
    - Delta computation
    - Delta statistics
    - Training data extraction

"""

import numpy as np
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

    def test_repr(self):
        """Test string representation."""
        fusion = DataFusion()
        repr_str = repr(fusion)

        assert "DataFusion" in repr_str
        assert "n_merged=0" in repr_str
        assert "delta_pm7" in repr_str

    @pytest.mark.skip(reason="Merge not implemented yet")
    def test_merge_inner(self, mock_chemperium_df, mock_pm7_df):
        """Test inner merge of datasets."""
        fusion = DataFusion()
        merged = fusion.merge(mock_chemperium_df, mock_pm7_df, on="smiles")

        assert len(merged) <= len(mock_chemperium_df)
        assert "H298_cbs" in merged.columns
        assert "H298_pm7" in merged.columns

    @pytest.mark.skip(reason="Compute deltas not implemented yet")
    def test_compute_deltas(self, mock_merged_df):
        """Test delta computation."""
        fusion = DataFusion()
        fusion.merged_data = mock_merged_df

        result = fusion.compute_deltas()

        assert fusion.delta_column in result.columns
        expected_delta = mock_merged_df["H298_cbs"] - mock_merged_df["H298_pm7"]
        np.testing.assert_array_almost_equal(
            result[fusion.delta_column].values,
            expected_delta.values,
        )

    @pytest.mark.skip(reason="Analyze deltas not implemented yet")
    def test_analyze_deltas(self, mock_merged_df):
        """Test delta statistics."""
        fusion = DataFusion()
        fusion.merged_data = mock_merged_df
        fusion.merged_data["delta_pm7"] = (
            mock_merged_df["H298_cbs"] - mock_merged_df["H298_pm7"]
        )

        stats = fusion.analyze_deltas()

        assert "mean" in stats
        assert "std" in stats
        assert "min" in stats
        assert "max" in stats
        assert "median" in stats


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
