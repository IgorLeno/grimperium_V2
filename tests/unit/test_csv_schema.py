"""Unit tests for CSV schema in BatchCSVManager.

Tests the get_schema() method and column definitions.
"""

from pathlib import Path

import pytest

pytest.importorskip("rdkit", reason="rdkit not available")

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager


class TestCSVSchema:
    """Tests for CSV schema definition."""

    @pytest.fixture
    def manager(self, tmp_path: Path) -> BatchCSVManager:
        """Create BatchCSVManager with dummy path."""
        return BatchCSVManager(tmp_path / "test.csv")

    def test_schema_contains_crest_settings(self, manager: BatchCSVManager) -> None:
        """Verify CREST columns exist."""
        schema = manager.get_schema()

        crest_columns = [
            "v3",  # Renamed from crest_v3
            "qm",  # Renamed from crest_quick
            "nci",  # Renamed from crest_nci
            "c_method",  # Renamed from crest_gfnff
            "energy_window",  # Renamed from crest_ewin
            "rmsd_threshold",  # Renamed from crest_rthr
            "crest_optlev",  # Keep old name for now
            "threads",  # Renamed from crest_threads
            "xtb",  # Renamed from crest_xtb_preopt
        ]

        for col in crest_columns:
            assert col in schema, f"Missing CREST column: {col}"

    def test_schema_contains_mopac_settings(self, manager: BatchCSVManager) -> None:
        """Verify MOPAC columns exist (Phase A names)."""
        schema = manager.get_schema()

        mopac_columns = [
            "precise_scf",  # Renamed from mopac_precise
            "scf_threshold",  # Renamed from mopac_scfcrt
            # Note: mopac_itry, mopac_pulay, mopac_prtall, mopac_archive not in Phase A
        ]

        for col in mopac_columns:
            assert col in schema, f"Missing MOPAC column: {col}"

    def test_crest_xtb_preopt_in_schema(self, manager: BatchCSVManager) -> None:
        """Verify xtb column exists (renamed from crest_xtb_preopt)."""
        schema = manager.get_schema()
        assert "xtb" in schema

    def test_schema_length(self, manager: BatchCSVManager) -> None:
        """Verify schema has expected column count.

        Current: 55 columns (Phase A complete schema)
        - 44 original columns
        - 3 new Phase A columns (mopac_status, mopac_time, conformer_selected)
        - 8 reserved Phase B placeholder columns
        """
        schema = manager.get_schema()
        assert len(schema) == 55, f"Expected 55 columns, got {len(schema)}"

    def test_schema_order(self, manager: BatchCSVManager) -> None:
        """Verify schema column order."""
        schema = manager.get_schema()

        assert schema[0] == "mol_id"
        assert schema[1] == "smiles"
        assert "status" in schema
        assert "xtb" in schema  # Renamed from crest_xtb_preopt
        assert "retry_count" in schema
        assert "last_error_message" in schema

        retry_idx = schema.index("retry_count")
        last_error_idx = schema.index("last_error_message")
        assert retry_idx < last_error_idx

    def test_schema_no_duplicates(self, manager: BatchCSVManager) -> None:
        """Verify no duplicate columns."""
        schema = manager.get_schema()
        unique_columns = set(schema)
        assert len(schema) == len(unique_columns), "Duplicate columns found"

    def test_result_columns_in_schema(self, manager: BatchCSVManager) -> None:
        """Verify RESULT_COLUMNS are in schema."""
        schema = manager.get_schema()

        for col in manager.RESULT_COLUMNS:
            assert col in schema, f"RESULT_COLUMN {col} not in schema"

    def test_identity_columns(self, manager: BatchCSVManager) -> None:
        """Verify identity columns class attribute."""
        assert manager.IDENTITY_COLUMNS == ["mol_id"]

    def test_molecular_properties_columns(self, manager: BatchCSVManager) -> None:
        """Verify molecular properties columns class attribute."""
        expected = [
            "smiles",
            "nheavy",
            "nrotbonds",
            "tpsa",
            "aromatic_rings",
            "has_heteroatoms",
            "reference_hof",
        ]
        assert expected == manager.MOLECULAR_PROPERTIES_COLUMNS

    def test_batch_info_columns(self, manager: BatchCSVManager) -> None:
        """Verify batch info columns class attribute."""
        expected = ["status", "batch_id", "batch_order", "batch_failure_policy"]
        assert expected == manager.BATCH_INFO_COLUMNS

    def test_crest_config_columns(self, manager: BatchCSVManager) -> None:
        """Verify CREST config columns class attribute (Phase A names)."""
        expected = [
            "v3",  # Renamed from crest_v3
            "qm",  # Renamed from crest_quick
            "nci",  # Renamed from crest_nci
            "c_method",  # Renamed from crest_gfnff
            "energy_window",  # Renamed from crest_ewin
            "rmsd_threshold",  # Renamed from crest_rthr
            "crest_optlev",  # Keep old name for now
            "threads",  # Renamed from crest_threads
            "xtb",  # Renamed from crest_xtb_preopt
        ]
        assert expected == manager.CREST_CONFIG_COLUMNS

    def test_mopac_config_columns(self, manager: BatchCSVManager) -> None:
        """Verify MOPAC config columns class attribute (Phase A names)."""
        expected = [
            "precise_scf",  # Renamed from mopac_precise
            "scf_threshold",  # Renamed from mopac_scfcrt
        ]
        assert expected == manager.MOPAC_CONFIG_COLUMNS

    def test_retry_tracking_columns(self, manager: BatchCSVManager) -> None:
        """Verify retry tracking columns class attribute."""
        expected = ["retry_count", "last_error_message"]
        assert expected == manager.RETRY_TRACKING_COLUMNS
