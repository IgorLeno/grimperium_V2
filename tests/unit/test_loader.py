"""
Unit tests for ChemperiumLoader.

Tests cover:
    - Loading from CSV and Parquet
    - Column validation
    - Train/test splitting
    - Feature extraction

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

    def test_repr(self):
        """Test string representation."""
        loader = ChemperiumLoader()
        assert "ChemperiumLoader" in repr(loader)
        assert "n_molecules=0" in repr(loader)

    @pytest.mark.skip(reason="Load not implemented yet")
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

    @pytest.mark.skip(reason="Load not implemented yet")
    def test_load_missing_file(self):
        """Test loading non-existent file raises error."""
        loader = ChemperiumLoader()
        with pytest.raises(FileNotFoundError):
            loader.load("nonexistent.csv")

    @pytest.mark.skip(reason="Split not implemented yet")
    def test_split(self, mock_chemperium_df):
        """Test train/test split."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        train, test = loader.split(test_size=0.2, random_state=42)

        assert len(train) + len(test) == len(mock_chemperium_df)
        assert len(test) == int(len(mock_chemperium_df) * 0.2)

    @pytest.mark.skip(reason="Get features not implemented yet")
    def test_get_features(self, mock_chemperium_df):
        """Test feature extraction."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        features = loader.get_features()

        assert "nheavy" in features.columns
        assert "charge" in features.columns
        assert "H298_cbs" not in features.columns

    @pytest.mark.skip(reason="Get targets not implemented yet")
    def test_get_targets(self, mock_chemperium_df):
        """Test target extraction."""
        loader = ChemperiumLoader()
        loader.data = mock_chemperium_df

        targets = loader.get_targets()

        assert isinstance(targets, pd.Series)
        assert len(targets) == len(mock_chemperium_df)


class TestChemperiumLoaderValidation:
    """Tests for ChemperiumLoader validation."""

    @pytest.mark.skip(reason="Validation not implemented yet")
    def test_validate_missing_columns(self):
        """Test validation fails with missing columns."""
        loader = ChemperiumLoader()
        df = pd.DataFrame({"smiles": ["CCO"]})  # Missing required columns

        with pytest.raises(ValueError, match="missing"):
            loader._validate_dataframe(df)

    @pytest.mark.skip(reason="Validation not implemented yet")
    def test_validate_valid_dataframe(self, mock_chemperium_df):
        """Test validation passes with valid DataFrame."""
        loader = ChemperiumLoader()
        # Should not raise
        loader._validate_dataframe(mock_chemperium_df)
