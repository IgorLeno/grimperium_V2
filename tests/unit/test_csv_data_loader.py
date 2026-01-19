"""
Tests for CSVDataLoader and BatchDataManager.

Includes tests for Ajuste #3 (duplicate check integration).
"""

from pathlib import Path

import pytest

from grimperium.core.csv_data_loader import (
    BatchDataManager,
    CSVDataLoader,
    CSVDataLoaderError,
)


@pytest.fixture
def valid_csv_content():
    """Minimal valid CSV content."""
    return (
        "mol_id,smiles,nheavy,status\n"
        "mol_001,C,1,Pending\n"
        "mol_002,CC,2,Pending\n"
        "mol_003,CCC,3,OK\n"
    )


@pytest.fixture
def csv_with_duplicates():
    """CSV with duplicate mol_ids."""
    return (
        "mol_id,smiles,nheavy,status\n"
        "mol_001,C,1,Pending\n"
        "mol_002,CC,2,Pending\n"
        "mol_001,C,1,Pending\n"  # Duplicate!
    )


@pytest.fixture
def csv_with_invalid_rows():
    """CSV with some invalid rows."""
    return (
        "mol_id,smiles,nheavy,status\n"
        "mol_001,C,1,Pending\n"
        ",CC,2,Pending\n"  # Missing mol_id
        "mol_003,CCC,3,OK\n"
        "mol_004,,4,Pending\n"  # Missing smiles
        "mol_005,CCCCC,5,Pending\n"
    )


class TestCSVDataLoaderBasic:
    """Basic tests for CSVDataLoader."""

    def test_load_valid_csv(self, tmp_path, valid_csv_content):
        """Load valid CSV successfully."""
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text(valid_csv_content)

        loader = CSVDataLoader(csv_file, strict=True)
        df = loader.load_dataframe()

        assert len(df) == 3
        assert "mol_001" in df["mol_id"].values

    def test_file_not_found(self, tmp_path):
        """Raise error for missing file."""
        csv_file = tmp_path / "nonexistent.csv"

        loader = CSVDataLoader(csv_file, strict=True)
        with pytest.raises(CSVDataLoaderError) as exc:
            loader.load_dataframe()

        assert "not found" in str(exc.value)

    def test_missing_required_columns(self, tmp_path):
        """Raise error for missing required columns."""
        csv_file = tmp_path / "bad.csv"
        csv_file.write_text("col1,col2\n1,2\n")

        loader = CSVDataLoader(csv_file, strict=True)
        with pytest.raises(CSVDataLoaderError) as exc:
            loader.load_dataframe()

        assert "Missing required columns" in str(exc.value)

    def test_adds_optional_columns(self, tmp_path, valid_csv_content):
        """Optional columns are added with defaults."""
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text(valid_csv_content)

        loader = CSVDataLoader(csv_file, strict=True)
        df = loader.load_dataframe()

        # Check optional columns were added
        assert "charge" in df.columns
        assert "multiplicity" in df.columns
        assert "H298_cbs" in df.columns


class TestCSVDataLoaderAjuste3:
    """Tests for Ajuste #3 - duplicate check integration."""

    def test_strict_mode_checks_duplicates_first(
        self, tmp_path, csv_with_duplicates
    ):
        """AJUSTE #3: Strict mode checks duplicates first."""
        csv_file = tmp_path / "duplicates.csv"
        csv_file.write_text(csv_with_duplicates)

        loader = CSVDataLoader(csv_file, strict=True)
        with pytest.raises(CSVDataLoaderError) as exc:
            loader.load_dataframe()

        assert "Duplicate mol_id" in str(exc.value)

    def test_permissive_mode_also_fails_on_duplicates(
        self, tmp_path, csv_with_duplicates
    ):
        """AJUSTE #3: Permissive mode also fails on duplicates (critical)."""
        csv_file = tmp_path / "duplicates.csv"
        csv_file.write_text(csv_with_duplicates)

        loader = CSVDataLoader(csv_file, strict=False)
        with pytest.raises(CSVDataLoaderError) as exc:
            loader.load_dataframe()

        # Permissive mode still fails on duplicates
        assert "Duplicate mol_id" in str(exc.value)

    def test_duplicate_error_message_helpful(
        self, tmp_path, csv_with_duplicates
    ):
        """AJUSTE #3: Error message includes recovery steps."""
        csv_file = tmp_path / "duplicates.csv"
        csv_file.write_text(csv_with_duplicates)

        loader = CSVDataLoader(csv_file, strict=True)
        with pytest.raises(CSVDataLoaderError) as exc:
            loader.load_dataframe()

        error_msg = str(exc.value)
        assert "Backup CSV" in error_msg or "data corruption" in error_msg


class TestCSVDataLoaderPermissive:
    """Tests for permissive validation mode."""

    def test_permissive_skips_invalid_rows(
        self, tmp_path, csv_with_invalid_rows
    ):
        """Permissive mode skips invalid rows."""
        csv_file = tmp_path / "invalid.csv"
        csv_file.write_text(csv_with_invalid_rows)

        loader = CSVDataLoader(csv_file, strict=False)
        df = loader.load_dataframe()

        # Only valid rows are returned
        assert len(df) == 3  # mol_001, mol_003, mol_005

    def test_permissive_reports_skipped_rows(
        self, tmp_path, csv_with_invalid_rows
    ):
        """Permissive mode reports skipped rows."""
        csv_file = tmp_path / "invalid.csv"
        csv_file.write_text(csv_with_invalid_rows)

        loader = CSVDataLoader(csv_file, strict=False)
        loader.load_dataframe()

        report = loader.get_validation_report()
        assert report.total_errors == 2  # Two invalid rows

    def test_strict_fails_on_first_invalid_row(
        self, tmp_path, csv_with_invalid_rows
    ):
        """Strict mode fails on first invalid row."""
        csv_file = tmp_path / "invalid.csv"
        csv_file.write_text(csv_with_invalid_rows)

        loader = CSVDataLoader(csv_file, strict=True)
        with pytest.raises(CSVDataLoaderError) as exc:
            loader.load_dataframe()

        assert "mol_id" in str(exc.value) or "smiles" in str(exc.value)


class TestBatchDataManager:
    """Tests for BatchDataManager."""

    def test_load_batch_creates_molecules(self, tmp_path, valid_csv_content):
        """load_batch creates Molecule objects."""
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text(valid_csv_content)

        manager = BatchDataManager(csv_file, strict=True)
        molecules = manager.load_batch()

        assert len(molecules) == 3
        assert molecules[0].mol_id == "mol_001"
        assert molecules[0].smiles == "C"

    def test_count_by_status(self, tmp_path, valid_csv_content):
        """count_by_status returns status counts."""
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text(valid_csv_content)

        manager = BatchDataManager(csv_file, strict=True)
        manager.load_batch()

        counts = manager.count_by_status()
        assert counts.get("pending", 0) == 2
        assert counts.get("ok", 0) == 1

    def test_get_dataframe(self, tmp_path, valid_csv_content):
        """get_dataframe returns the loaded DataFrame."""
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text(valid_csv_content)

        manager = BatchDataManager(csv_file, strict=True)
        manager.load_batch()

        df = manager.get_dataframe()
        assert df is not None
        assert len(df) == 3
