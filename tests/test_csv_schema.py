"""Tests for Phase A CSV schema (55 columns, derived from len(EXPECTED_COLUMNS))."""

from pathlib import Path

import pandas as pd
import pytest

EXPECTED_COLUMNS = [
    # Core molecule data (1-6)
    "mol_id",
    "status",
    "smiles",
    "multiplicity",
    "charge",
    "nheavy",
    # Thermodynamic data (7-11)
    "H298_cbs",
    "H298_pm7",
    "target_delta_kcalmol",
    "abs_diff",
    "abs_diff_%",
    # Batch metadata (12-15)
    "batch_id",
    "timestamp",
    "total_time",
    "reruns",
    # CREST settings (16-24)
    "crest_status",
    "xtb",
    "v3",
    "qm",
    "nci",
    "c_method",
    "energy_window",
    "rmsd_threshold",
    "opt_lvl",
    # RDKit descriptors (25-39)
    "rdkit_nrotbonds",
    "rdkit_tpsa",
    "rdkit_num_rings",
    "rdkit_fsp3",
    "rdkit_mol_weight",
    "rdkit_hbond_donors",
    "rdkit_hbond_acceptors",
    "rdkit_nC",
    "rdkit_nH",
    "rdkit_nO",
    "rdkit_nN",
    "rdkit_bonds_single",
    "rdkit_bonds_double",
    "rdkit_bonds_triple",
    "rdkit_bonds_aromatic",
    # CREST results (40-42)
    "crest_conformers_generated",
    "num_conformers_selected",
    "crest_time_s",
    # MOPAC status (43)
    "mopac_status",
    # PM7 selection & electronic descriptors (44-55)
    "k_selected_pm7",
    "mopac_dipole_debye",
    "mopac_ionization_potential_ev",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
    "mopac_cosmo_area_a2",
    "mopac_cosmo_volume_a3",
    "mopac_gradient_norm",
    "mopac_num_scf_cycles",
    "mopac_point_group",
    "mopac_time_s",
]


_CSV_PATH = Path("data/thermo_pm7.csv")


@pytest.mark.skipif(
    not _CSV_PATH.exists(),
    reason="CSV data file not yet generated; run the pipeline first",
)
def test_csv_has_expected_columns() -> None:
    """Test that CSV output has the expected column count."""
    df = pd.read_csv(_CSV_PATH)
    assert len(df.columns) == len(
        EXPECTED_COLUMNS
    ), f"Expected {len(EXPECTED_COLUMNS)} columns, got {len(df.columns)}"


def test_csv_column_order() -> None:
    """Test that CSV columns are in correct order."""
    csv_path = Path("data/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        # Skip test if CSV hasn't been regenerated yet
        if len(df.columns) < 50:
            pytest.skip("CSV not yet regenerated with new schema")
        assert list(df.columns) == EXPECTED_COLUMNS
    else:
        pytest.skip("CSV file does not exist yet")


def test_csv_no_old_columns() -> None:
    """Test that old column names are removed."""
    csv_path = Path("data/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        assert "crest_timeout" not in df.columns
        assert "mopac_timeout" not in df.columns
        assert "most_stable_hof" not in df.columns
        assert "delta_1" not in df.columns
        assert "delta_2" not in df.columns
        assert "delta_3" not in df.columns
        assert "conformer_selected" not in df.columns
    else:
        pytest.skip("CSV file does not exist yet")


def test_csv_has_new_columns() -> None:
    """Test that new columns are present."""
    csv_path = Path("data/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        # Skip test if CSV hasn't been regenerated yet
        if len(df.columns) < 50:
            pytest.skip("CSV not yet regenerated with new schema")
        # Metrics
        assert "abs_diff" in df.columns
        assert "abs_diff_%" in df.columns
        # PM7 selection & descriptors
        assert "target_delta_kcalmol" in df.columns
        assert "k_selected_pm7" in df.columns
        assert "mopac_dipole_debye" in df.columns
        assert "mopac_homo_ev" in df.columns
        assert "mopac_gap_ev" in df.columns
        # Settings
        assert "c_method" in df.columns
        assert "qm" in df.columns
        assert "opt_lvl" in df.columns
    else:
        pytest.skip("CSV file does not exist yet")
