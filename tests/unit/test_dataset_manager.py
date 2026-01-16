"""Unit tests for DatasetManager.

Tests the source-of-truth CSV synchronization functionality.
"""

from pathlib import Path

import pandas as pd
import pytest

from grimperium.cli.dataset_manager import DatasetManager


class TestDatasetManager:
    """Tests for DatasetManager class."""

    @pytest.fixture
    def source_csv(self, tmp_path: Path) -> Path:
        """Create a mock source CSV for testing."""
        csv_path = tmp_path / "source.csv"
        pd.DataFrame(
            {
                "smiles": ["CCO", "CC(C)O", "CC(C)C(O)C", "C1CCCCC1"],
                "nheavy": [3, 4, 6, 6],
                "multiplicity": [1, 1, 1, 1],
                "charge": [0, 0, 0, 0],
                "H298_cbs": [0.1, 0.2, 0.3, 0.4],
                "H298_b3": [0.05, 0.15, 0.25, 0.35],
            }
        ).to_csv(csv_path, index=False)
        return csv_path

    @pytest.fixture
    def manager(self, source_csv: Path, tmp_path: Path) -> DatasetManager:
        """Create a DatasetManager instance for testing."""
        working_csv = tmp_path / "working.csv"
        return DatasetManager(source_csv, working_csv)

    def test_init(self, tmp_path: Path) -> None:
        """Test DatasetManager initialization."""
        source = tmp_path / "source.csv"
        working = tmp_path / "working.csv"
        manager = DatasetManager(source, working)
        assert manager.source_csv == source
        assert manager.working_csv == working

    def test_get_molecule_count(self, manager: DatasetManager) -> None:
        """Test get_molecule_count returns correct count."""
        count = manager.get_molecule_count()
        assert count == 4

    def test_get_molecule_count_missing_file(self, tmp_path: Path) -> None:
        """Test get_molecule_count raises FileNotFoundError for missing file."""
        manager = DatasetManager(
            tmp_path / "nonexistent.csv",
            tmp_path / "working.csv",
        )
        with pytest.raises(FileNotFoundError):
            manager.get_molecule_count()

    def test_sync_from_source(self, manager: DatasetManager) -> None:
        """Test sync_from_source returns correct DataFrame."""
        df = manager.sync_from_source(2)

        assert len(df) == 2
        assert "mol_id" in df.columns
        assert "status" in df.columns
        assert df["mol_id"].tolist() == ["mol_00001", "mol_00002"]
        assert df["status"].tolist() == ["PENDING", "PENDING"]
        assert df["smiles"].tolist() == ["CCO", "CC(C)O"]

    def test_sync_from_source_all_molecules(self, manager: DatasetManager) -> None:
        """Test sync_from_source with all molecules."""
        df = manager.sync_from_source(4)

        assert len(df) == 4
        assert df["mol_id"].tolist() == [
            "mol_00001",
            "mol_00002",
            "mol_00003",
            "mol_00004",
        ]

    def test_sync_from_source_exceeds_available(self, manager: DatasetManager) -> None:
        """Test sync_from_source raises ValueError when exceeding available."""
        with pytest.raises(ValueError, match="Requested 10 molecules"):
            manager.sync_from_source(10)

    def test_create_working_csv(self, manager: DatasetManager) -> None:
        """Test create_working_csv creates file with correct content."""
        manager.create_working_csv(3)

        assert manager.working_csv.exists()
        df = pd.read_csv(manager.working_csv)
        assert len(df) == 3
        assert df["mol_id"].tolist() == ["mol_00001", "mol_00002", "mol_00003"]
        assert all(df["status"] == "PENDING")

    def test_validate_sync_valid(self, manager: DatasetManager) -> None:
        """Test validate_sync returns True for valid sync."""
        manager.create_working_csv(3)
        assert manager.validate_sync() is True

    def test_validate_sync_missing_working(self, manager: DatasetManager) -> None:
        """Test validate_sync returns False for missing working CSV."""
        assert manager.validate_sync() is False

    def test_validate_sync_smiles_mismatch(
        self, manager: DatasetManager, tmp_path: Path
    ) -> None:
        """Test validate_sync returns False for SMILES mismatch."""
        manager.create_working_csv(2)

        # Corrupt the working CSV
        df = pd.read_csv(manager.working_csv)
        df.loc[0, "smiles"] = "CORRUPTED"
        df.to_csv(manager.working_csv, index=False)

        assert manager.validate_sync() is False

    def test_refresh_database_existing(self, manager: DatasetManager) -> None:
        """Test refresh_database resyncs existing working CSV."""
        manager.create_working_csv(2)

        # Corrupt the working CSV
        df = pd.read_csv(manager.working_csv)
        df.loc[0, "smiles"] = "CORRUPTED"
        df.to_csv(manager.working_csv, index=False)

        # Refresh should restore correct data
        manager.refresh_database()

        df = pd.read_csv(manager.working_csv)
        assert df["smiles"].tolist() == ["CCO", "CC(C)O"]

    def test_refresh_database_preserves_count(self, manager: DatasetManager) -> None:
        """Test refresh_database preserves molecule count."""
        manager.create_working_csv(3)
        manager.refresh_database()

        df = pd.read_csv(manager.working_csv)
        assert len(df) == 3

    def test_refresh_database_missing_source(self, tmp_path: Path) -> None:
        """Test refresh_database raises FileNotFoundError for missing source."""
        manager = DatasetManager(
            tmp_path / "nonexistent.csv",
            tmp_path / "working.csv",
        )
        with pytest.raises(FileNotFoundError):
            manager.refresh_database()

    def test_mol_id_format(self, manager: DatasetManager) -> None:
        """Test mol_id uses 5-digit zero-padded format."""
        df = manager.sync_from_source(4)

        for i, mol_id in enumerate(df["mol_id"], start=1):
            assert mol_id == f"mol_{i:05d}"
            assert len(mol_id) == 9  # "mol_" + 5 digits

    def test_column_order(self, manager: DatasetManager) -> None:
        """Test columns are ordered correctly (mol_id, status first)."""
        df = manager.sync_from_source(2)
        columns = df.columns.tolist()

        assert columns[0] == "mol_id"
        assert columns[1] == "status"
