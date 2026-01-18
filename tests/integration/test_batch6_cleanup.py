"""Integration tests for Batch 6: File Cleanup & Conformer Details Renaming."""

from pathlib import Path


class TestBatch6Cleanup:
    """Test suite verifying Batch 6 cleanup operations."""

    def test_csv_files_deleted(self):
        """Verify old CSV files are deleted."""
        assert not Path("data/test_batch_small.csv").exists()
        assert not Path("data/thermo_cbs_clean.csv").exists()
        assert not Path("data/thermo_cbs_opt.csv").exists()

    def test_required_csvs_exist(self):
        """Verify required CSV files still exist."""
        assert Path("data/thermo_cbs_chon.csv").exists(), "Primary dataset missing"
        assert Path("data/thermo_pm7.csv").exists(), "Secondary dataset (PM7) missing"

    def test_csv_count_exactly_two(self):
        """Verify exactly 2 CSV files in data/ (no orphaned CSVs)."""
        csv_files = list(Path("data").glob("*.csv"))
        expected = {"thermo_cbs_chon.csv", "thermo_pm7.csv"}
        actual = {f.name for f in csv_files}
        assert actual == expected, (
            f"Expected {expected}, got {actual}. "
            f"Found {len(csv_files)} CSV files: {csv_files}"
        )

    def test_conformer_details_renamed(self):
        """Verify all JSON files renamed to mol_XXXXX format."""
        conformer_dir = Path("data/molecules_pm7/conformer_details")

        # No numeric files should remain
        numeric_files = list(conformer_dir.glob("[0-9].json")) + list(
            conformer_dir.glob("[0-9][0-9].json")
        )
        assert len(numeric_files) == 0, f"Found numeric JSON files: {numeric_files}"

        # Should have exactly 3 mol_XXXXX files
        mol_files = list(conformer_dir.glob("mol_*.json"))
        assert (
            len(mol_files) == 3
        ), f"Expected 3 mol_XXXXX files, found {len(mol_files)}"

        # All should be in 5-digit format
        for f in mol_files:
            num_part = f.name.replace("mol_", "").replace(".json", "")
            assert len(num_part) == 5, f"Not 5-digit format: {f.name}"
            assert num_part.isdigit(), f"Not all digits: {f.name}"
