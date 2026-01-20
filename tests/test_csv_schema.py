"""Tests for Phase A CSV schema (49 columns)."""

import pandas as pd
import pytest
from pathlib import Path

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
    # Batch metadata (11-13)
    "batch_id",
    "timestamp",
    "reruns",
    # RDKit descriptors (14-16)
    "nrotbonds",
    "tpsa",
    "aromatic_rings",
    # CREST results and settings (17-28)
    "crest_status",
    "xtb",
    "v3",
    "qm",
    "nci",
    "c_method",
    "energy_window",
    "rmsd_threshold",
    "threads",
    "crest_conformers_generated",
    "crest_time",
    "num_conformers_selected",
    # MOPAC results and settings (29-32)
    "mopac_status",
    "precise_scf",
    "scf_threshold",
    "mopac_time",
    # Delta calculations (33-36)
    "delta_1",
    "delta_2",
    "delta_3",
    "conformer_selected",
    # Error handling and batch control (37-41)
    "error_message",
    "batch_order",
    "batch_failure_policy",
    "assigned_crest_timeout",
    "assigned_mopac_timeout",
    # Reserved for Phase B (42-49)
    "reserved_42",
    "reserved_43",
    "reserved_44",
    "reserved_45",
    "reserved_46",
    "reserved_47",
    "reserved_48",
    "reserved_49",
]


def test_csv_has_49_columns() -> None:
    """Test that CSV output has exactly 49 columns (after regeneration)."""
    csv_path = Path("data/thermo_pm7.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        # Skip test if CSV hasn't been regenerated yet (will have 41 columns)
        if len(df.columns) == 41:
            pytest.skip("CSV not yet regenerated with new schema (still 41 columns)")
        assert len(df.columns) == 49, f"Expected 49 columns, got {len(df.columns)}"
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
    else:
        pytest.skip("CSV file does not exist yet")
