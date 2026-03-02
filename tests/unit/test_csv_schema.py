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

    def test_schema_length(self, manager: BatchCSVManager) -> None:
        """Verify schema has expected 57-column count.

        - 3 identity columns
        - 7 molecular property columns
        - 4 batch info columns
        - 15 RDKit descriptor columns
        - 12 CREST columns
        - 13 MOPAC columns
        - 3 batch config columns (batch_order moved here)
        """
        schema = manager.get_schema()
        assert len(schema) == 58, f"Expected 58 columns, got {len(schema)}"

    def test_schema_exact_order(self, manager: BatchCSVManager) -> None:
        """Verify schema matches target 57-column order exactly."""
        schema = manager.get_schema()
        expected = [
            "mol_id", "status", "smiles",
            "multiplicity", "charge", "nheavy", "H298_cbs",
            "H298_pm7", "target_delta_kcalmol", "abs_diff_%",
            "batch_id", "timestamp", "total_time", "reruns",
            "rdkit_nrotbonds", "rdkit_tpsa", "rdkit_num_rings", "rdkit_fsp3",
            "rdkit_mol_weight", "rdkit_hbond_donors", "rdkit_hbond_acceptors",
            "rdkit_nC", "rdkit_nH", "rdkit_nO", "rdkit_nN",
            "rdkit_bonds_single", "rdkit_bonds_double", "rdkit_bonds_triple",
            "rdkit_bonds_aromatic",
            "crest_status", "xtb", "v3", "qm", "nci", "c_method",
            "energy_window", "rmsd_threshold", "opt_lvl",
            "crest_conformers_generated", "num_conformers_selected", "crest_time_s",
            "mopac_status", "k_selected_pm7",
            "mopac_dipole_debye", "mopac_ionization_potential_ev",
            "mopac_homo_ev", "mopac_lumo_ev", "mopac_gap_ev",
            "mopac_cosmo_area_a2", "mopac_cosmo_volume_a3",
            "mopac_gradient_norm", "mopac_num_scf_cycles", "mopac_point_group",
            "mopac_time_s",
            "batch_order", "batch_failure_policy",
            "assigned_crest_timeout", "assigned_mopac_timeout",
        ]
        assert schema == expected

    def test_schema_no_duplicates(self, manager: BatchCSVManager) -> None:
        """Verify no duplicate columns."""
        schema = manager.get_schema()
        unique_columns = set(schema)
        assert len(schema) == len(unique_columns), "Duplicate columns found"

    def test_schema_contains_crest_settings(self, manager: BatchCSVManager) -> None:
        """Verify CREST columns exist."""
        schema = manager.get_schema()

        crest_columns = [
            "crest_status", "xtb", "v3", "qm", "nci", "c_method",
            "energy_window", "rmsd_threshold", "opt_lvl",
            "crest_conformers_generated", "num_conformers_selected", "crest_time_s",
        ]

        for col in crest_columns:
            assert col in schema, f"Missing CREST column: {col}"

    def test_schema_contains_mopac_columns(self, manager: BatchCSVManager) -> None:
        """Verify MOPAC columns exist."""
        schema = manager.get_schema()

        mopac_columns = [
            "mopac_status", "k_selected_pm7",
            "mopac_dipole_debye", "mopac_ionization_potential_ev",
            "mopac_homo_ev", "mopac_lumo_ev", "mopac_gap_ev",
            "mopac_cosmo_area_a2", "mopac_cosmo_volume_a3",
            "mopac_gradient_norm", "mopac_num_scf_cycles", "mopac_point_group",
            "mopac_time_s",
        ]

        for col in mopac_columns:
            assert col in schema, f"Missing MOPAC column: {col}"

    def test_schema_contains_rdkit_columns(self, manager: BatchCSVManager) -> None:
        """Verify RDKit descriptor columns exist with rdkit_ prefix."""
        schema = manager.get_schema()

        rdkit_columns = [
            "rdkit_nrotbonds", "rdkit_tpsa", "rdkit_num_rings", "rdkit_fsp3",
            "rdkit_mol_weight", "rdkit_hbond_donors", "rdkit_hbond_acceptors",
            "rdkit_nC", "rdkit_nH", "rdkit_nO", "rdkit_nN",
            "rdkit_bonds_single", "rdkit_bonds_double", "rdkit_bonds_triple",
            "rdkit_bonds_aromatic",
        ]

        for col in rdkit_columns:
            assert col in schema, f"Missing RDKit column: {col}"

    def test_schema_removed_old_columns(self, manager: BatchCSVManager) -> None:
        """Verify old columns are no longer in schema."""
        schema = manager.get_schema()

        removed = [
            "reference_hof", "has_heteroatoms", "quality_grade",
            "success", "error_message", "total_execution_time",
            "crest_error", "mopac_time", "precise_scf", "scf_threshold",
            "crest_optlev", "threads", "retry_count", "last_error_message",
            "max_retries", "reserved_42", "reserved_43", "reserved_44",
        ]

        for col in removed:
            assert col not in schema, f"Old column still present: {col}"

    def test_result_columns_in_schema(self, manager: BatchCSVManager) -> None:
        """Verify RESULT_COLUMNS are in schema."""
        schema = manager.get_schema()

        for col in manager.RESULT_COLUMNS:
            assert col in schema, f"RESULT_COLUMN {col} not in schema"

    def test_identity_columns(self, manager: BatchCSVManager) -> None:
        """Verify identity columns class attribute."""
        assert ["mol_id", "status", "smiles"] == manager.IDENTITY_COLUMNS

    def test_molecular_properties_columns(self, manager: BatchCSVManager) -> None:
        """Verify molecular properties columns class attribute."""
        expected = [
            "multiplicity", "charge", "nheavy", "H298_cbs",
            "H298_pm7", "target_delta_kcalmol", "abs_diff_%",
        ]
        assert expected == manager.MOLECULAR_PROPERTIES_COLUMNS

    def test_batch_info_columns(self, manager: BatchCSVManager) -> None:
        """Verify batch info columns class attribute."""
        expected = ["batch_id", "timestamp", "total_time", "reruns"]
        assert expected == manager.BATCH_INFO_COLUMNS

    def test_batch_config_columns(self, manager: BatchCSVManager) -> None:
        """Verify batch config columns class attribute."""
        expected = [
            "batch_order", "batch_failure_policy",
            "assigned_crest_timeout", "assigned_mopac_timeout",
        ]
        assert expected == manager.BATCH_CONFIG_COLUMNS
