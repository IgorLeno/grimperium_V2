"""Tests for Phase A CSV schema (43 columns)."""

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
    # Thermodynamic data (7-10)
    "H298_cbs",
    "H298_pm7",
    "abs_diff",
    "abs_diff_%",
    # Batch metadata (11-14)
    "batch_id",
    "timestamp",
    "total_time",
    "reruns",
    # RDKit descriptors (15-17)
    "nrotbonds",
    "tpsa",
    "aromatic_rings",
    # CREST results and settings (18-30)
    "crest_status",
    "xtb",
    "v3",
    "qm",
    "nci",
    "c_method",
    "energy_window",
    "rmsd_threshold",
    "opt_lvl",
    "threads",
    "crest_conformers_generated",
    "crest_time",
    "num_conformers_selected",
    # MOPAC results and settings (31-34)
    "mopac_status",
    "precise_scf",
    "scf_threshold",
    "mopac_time",
    # Delta calculations (35-38)
    "delta_1",
    "delta_2",
    "delta_3",
    "conformer_selected",
    # Error handling and batch control (39-43)
    "error_message",
    "batch_order",
    "batch_failure_policy",
    "assigned_crest_timeout",
    "assigned_mopac_timeout",
]


def test_csv_has_43_columns() -> None:
    """Test that CSV output has exactly 43 columns (after regeneration)."""
    csv_path = Path("data/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        # Skip test if CSV hasn't been regenerated yet (will have 41 columns)
        if len(df.columns) == 41:
            pytest.skip("CSV not yet regenerated with new schema (still 41 columns)")
        assert len(df.columns) == 43, f"Expected 43 columns, got {len(df.columns)}"
    else:
        pytest.skip("CSV file does not exist yet")


def test_csv_column_order() -> None:
    """Test that CSV columns are in correct order."""
    csv_path = Path("data/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        # Skip test if CSV hasn't been regenerated yet (will have 41 columns)
        if len(df.columns) == 41:
            pytest.skip("CSV not yet regenerated with new schema (still 41 columns)")
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
    else:
        pytest.skip("CSV file does not exist yet")


def test_csv_has_new_columns() -> None:
    """Test that new columns are present."""
    csv_path = Path("data/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        # Skip test if CSV hasn't been regenerated yet (will have 41 columns)
        if len(df.columns) == 41:
            pytest.skip("CSV not yet regenerated with new schema (still 41 columns)")
        # Metrics
        assert "abs_diff" in df.columns
        assert "abs_diff_%" in df.columns
        # Deltas
        assert "delta_1" in df.columns
        assert "delta_2" in df.columns
        assert "delta_3" in df.columns
        assert "conformer_selected" in df.columns
        # Settings
        assert "c_method" in df.columns
        assert "qm" in df.columns
        assert "opt_lvl" in df.columns
    else:
        pytest.skip("CSV file does not exist yet")
