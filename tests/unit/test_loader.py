"""
Unit tests for ChemperiumLoader.

Tests cover:
    - Loading from CSV and Parquet
    - Column validation
    - Train/test splitting
    - Feature extraction
    - Filtering

"""

import pandas as pd
import pytest

from grimperium.data.loader import ChemperiumLoader


class TestChemperiumLoader:
    """Tests for ChemperiumLoader class."""

    def test_init_default(self):
        """Test default initialization."""
        loader = ChemperiumLoader()
        assert loader.validate is True
        assert loader.data is None

    def test_init_no_validate(self):
        """Test initialization with validation disabled."""
        loader = ChemperiumLoader(validate=False)
        assert loader.validate is False

    def test_required_columns_defined(self):
        """Test that required columns are defined."""
        assert "smiles" in ChemperiumLoader.REQUIRED_COLUMNS
        assert "H298_cbs" in ChemperiumLoader.REQUIRED_COLUMNS
        assert "charge" in ChemperiumLoader.REQUIRED_COLUMNS

    def test_cp_columns_defined(self):
        """Test that Cp columns are properly defined."""
        assert len(ChemperiumLoader.CP_COLUMNS) == 45
        assert ChemperiumLoader.CP_COLUMNS[0] == "cp_1"
        assert ChemperiumLoader.CP_COLUMNS[-1] == "cp_45"

    def test_feature_columns_defined(self):
        """Test that feature columns are defined."""
        assert "nheavy" in ChemperiumLoader.FEATURE_COLUMNS
        assert "charge" in ChemperiumLoader.FEATURE_COLUMNS
        assert "multiplicity" in ChemperiumLoader.FEATURE_COLUMNS

    def test_repr(self):
        """Test string representation."""
        loader = ChemperiumLoader()
        assert "ChemperiumLoader" in repr(loader)
        assert "n_molecules=0" in repr(loader)

    def test_load_csv(self, tmp_path, mock_chemperium_df):
        """Test loading from CSV file."""
        # Save mock data
        csv_path = tmp_path / "test.csv"
        mock_chemperium_df.to_csv(csv_path, index=False)

        # Load
        loader = ChemperiumLoader()
        df = loader.load(csv_path)

        assert len(df) == len(mock_chemperium_df)
        assert "smiles" in df.columns
        assert "H298_cbs" in df.columns
        assert loader.data is not None

    def test_load_parquet(self, tmp_path, mock_chemperium_df):
        """Test loading from Parquet file."""
        pytest.importorskip("pyarrow", reason="pyarrow required for parquet tests")

        # Save mock data
        parquet_path = tmp_path / "test.parquet"
        mock_chemperium_df.to_parquet(parquet_path, index=False)

        # Load
        loader = ChemperiumLoader()
        df = loader.load(parquet_path)

        assert len(df) == len(mock_chemperium_df)
        assert "smiles" in df.columns
        assert "H298_cbs" in df.columns

    def test_load_missing_file(self):
        """Test loading non-existent file raises error."""
        loader = ChemperiumLoader()
        with pytest.raises(FileNotFoundError):
            loader.load("nonexistent.csv")

    def test_load_unsupported_format(self, tmp_path):
        """Test loading unsupported file format raises error."""
        # Create a file with unsupported extension
        txt_path = tmp_path / "test.txt"
        txt_path.write_text("some data")

        loader = ChemperiumLoader()
        with pytest.raises(ValueError, match="Unsupported file type"):
            loader.load(txt_path)

    def test_load_with_max_nheavy_filter(self, tmp_path, mock_chemperium_df):
        """Test loading with max_nheavy filter."""
        # Save mock data
        csv_path = tmp_path / "test.csv"
        mock_chemperium_df.to_csv(csv_path, index=False)

        # Load with filter
        loader = ChemperiumLoader()
        df = loader.load(csv_path, max_nheavy=5)

        # All loaded molecules should have nheavy <= 5
        assert all(df["nheavy"] <= 5)
        # Should have fewer molecules than original (if any had nheavy > 5)
        assert len(df) <= len(mock_chemperium_df)

    def test_split(self, mock_chemperium_df):
        """Test train/test split."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        train, test = loader.split(test_size=0.2, random_state=42)

        assert len(train) + len(test) == len(mock_chemperium_df)
        # With 10 samples, 20% test = 2 samples
        assert len(test) == 2

    def test_split_reproducibility(self, mock_chemperium_df):
        """Test that split with same seed produces same results."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        train1, test1 = loader.split(test_size=0.2, random_state=42)
        train2, test2 = loader.split(test_size=0.2, random_state=42)

        pd.testing.assert_frame_equal(train1.reset_index(drop=True), train2.reset_index(drop=True))
        pd.testing.assert_frame_equal(test1.reset_index(drop=True), test2.reset_index(drop=True))

    def test_train_val_test_split(self, mock_chemperium_df):
        """Test 3-way train/val/test split."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        train, val, test = loader.train_val_test_split(
            test_size=0.2, val_size=0.1, random_state=42
        )

        # Total should equal original
        assert len(train) + len(val) + len(test) == len(mock_chemperium_df)
        # Check proportions (approximate due to small sample size)
        assert len(test) == 2  # 20% of 10
        assert len(val) == 1   # 10% of 10 (approximately)

    def test_train_val_test_split_invalid_sizes(self, mock_chemperium_df):
        """Test that invalid split sizes raise error."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        with pytest.raises(ValueError, match="must be less than 1.0"):
            loader.train_val_test_split(test_size=0.6, val_size=0.5)

    def test_get_features(self, mock_chemperium_df):
        """Test feature extraction."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        features = loader.get_features()

        assert "nheavy" in features.columns
        assert "charge" in features.columns
        assert "H298_cbs" not in features.columns

    def test_get_features_with_cp(self, mock_chemperium_df):
        """Test feature extraction including Cp columns."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        features = loader.get_features(include_cp=True)

        assert "nheavy" in features.columns
        # Should include cp columns that exist in the data
        cp_cols = [c for c in features.columns if c.startswith("cp_")]
        assert len(cp_cols) > 0

    def test_get_targets(self, mock_chemperium_df):
        """Test target extraction."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        targets = loader.get_targets()

        assert isinstance(targets, pd.Series)
        assert len(targets) == len(mock_chemperium_df)

    def test_get_targets_custom_column(self, mock_chemperium_df):
        """Test target extraction with custom column."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        targets = loader.get_targets(target="S298")

        assert isinstance(targets, pd.Series)
        assert len(targets) == len(mock_chemperium_df)

    def test_get_targets_missing_column(self, mock_chemperium_df):
        """Test target extraction fails for missing column."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        with pytest.raises(ValueError, match="not found"):
            loader.get_targets(target="nonexistent")

    def test_split_no_data_raises(self):
        """Test split raises error when no data loaded."""
        loader = ChemperiumLoader()

        with pytest.raises(ValueError, match="No data available"):
            loader.split()

    def test_get_features_no_data_raises(self):
        """Test get_features raises error when no data loaded."""
        loader = ChemperiumLoader()

        with pytest.raises(ValueError, match="No data available"):
            loader.get_features()


class TestChemperiumLoaderValidation:
    """Tests for ChemperiumLoader validation."""

    def test_validate_missing_columns(self):
        """Test validation fails with missing columns."""
        loader = ChemperiumLoader()
        df = pd.DataFrame({"smiles": ["CCO"]})  # Missing required columns

        with pytest.raises(ValueError, match="Missing required columns"):
            loader._validate_dataframe(df)

    def test_validate_valid_dataframe(self, mock_chemperium_df):
        """Test validation passes with valid DataFrame."""
        loader = ChemperiumLoader()
        # Should not raise
        loader._validate_dataframe(mock_chemperium_df)

    def test_load_with_validation_error(self, tmp_path):
        """Test loading file with missing columns raises validation error."""
        # Create CSV with incomplete data
        incomplete_df = pd.DataFrame({"smiles": ["CCO", "CC"], "H298_cbs": [1.0, 2.0]})
        csv_path = tmp_path / "incomplete.csv"
        incomplete_df.to_csv(csv_path, index=False)

        loader = ChemperiumLoader(validate=True)
        with pytest.raises(ValueError, match="Missing required columns"):
            loader.load(csv_path)

    def test_load_without_validation(self, tmp_path):
        """Test loading file without validation skips column check."""
        # Create CSV with incomplete data
        incomplete_df = pd.DataFrame({"smiles": ["CCO", "CC"], "H298_cbs": [1.0, 2.0]})
        csv_path = tmp_path / "incomplete.csv"
        incomplete_df.to_csv(csv_path, index=False)

        loader = ChemperiumLoader(validate=False)
        df = loader.load(csv_path)

        # Should load without error
        assert len(df) == 2
