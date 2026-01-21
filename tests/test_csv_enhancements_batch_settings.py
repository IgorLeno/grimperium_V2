"""Tests for batch settings persistence in CSV.

This module tests that batch configuration settings (CREST/MOPAC parameters)
are correctly saved to CSV during batch execution.

Bug fixed: BATCH 13 - SPEC #1
Related files:
- src/grimperium/crest_pm7/csv_enhancements.py
- src/grimperium/crest_pm7/batch/csv_manager.py
"""

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.csv_enhancements import (
    CSVManagerExtensions,
)


def test_batch_settings_persist_to_csv(tmp_path):
    """Test that batch settings are saved to CSV via _update_extra_fields."""
    # Setup CSV with minimal schema
    csv_path = tmp_path / "test.csv"
    df = pd.DataFrame(
        {
            "mol_id": ["mol_001"],
            "smiles": ["CCO"],
            "nheavy": [2],
            "status": ["RUNNING"],
            "v3": [None],
            "qm": [None],
            "nci": [None],
            "c_method": [None],
            "energy_window": [None],
            "rmsd_threshold": [None],
            "threads": [None],
            "xtb": [None],
            "precise_scf": [None],
            "scf_threshold": [None],
            "delta_1": [None],
            "delta_2": [None],
            "delta_3": [None],
            "conformer_selected": [None],
            "abs_diff": [None],
            "abs_diff_%": [None],
        }
    )
    df.to_csv(csv_path, index=False)

    # Create manager
    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()

    # Mock batch settings (typical CREST/MOPAC config)
    batch_settings = {
        "v3": False,
        "qm": False,
        "nci": False,
        "c_method": "gfn2-xtb",
        "energy_window": 10.0,
        "rmsd_threshold": 0.125,
        "threads": 4,
        "xtb": True,
        "precise_scf": True,
        "scf_threshold": 1.0,
    }

    # Update molecule with settings
    success = CSVManagerExtensions.update_molecule_with_mopac_results(
        csv_manager=csv_manager,
        mol_id="mol_001",
        h298_cbs=-17.5,
        h298_pm7=-15.3,
        mopac_hof_values=[0.42, 0.87, 1.23],
        batch_settings=batch_settings,
    )

    assert success, "update_molecule_with_mopac_results should return True"

    # Verify CSV has settings
    df_updated = pd.read_csv(csv_path)

    # Check CREST settings (use == not is, pandas returns numpy bool types)
    assert df_updated.loc[0, "v3"] == False  # noqa: E712
    assert df_updated.loc[0, "qm"] == False  # noqa: E712
    assert df_updated.loc[0, "nci"] == False  # noqa: E712
    assert df_updated.loc[0, "c_method"] == "gfn2-xtb"
    assert df_updated.loc[0, "energy_window"] == 10.0
    assert df_updated.loc[0, "rmsd_threshold"] == 0.125
    assert df_updated.loc[0, "threads"] == 4
    assert df_updated.loc[0, "xtb"] == True  # noqa: E712

    # Check MOPAC settings
    assert df_updated.loc[0, "precise_scf"] == True  # noqa: E712
    assert df_updated.loc[0, "scf_threshold"] == 1.0

    # Check calculated fields also saved
    assert df_updated.loc[0, "delta_1"] == 0.0  # Best conformer delta
    assert df_updated.loc[0, "delta_2"] == pytest.approx(0.45, abs=0.01)
    assert df_updated.loc[0, "delta_3"] == pytest.approx(0.81, abs=0.01)
    assert df_updated.loc[0, "conformer_selected"] == 0


def test_batch_settings_with_nan_deltas(tmp_path):
    """Test that settings persist even when deltas are NaN (edge case)."""
    csv_path = tmp_path / "test.csv"
    df = pd.DataFrame(
        {
            "mol_id": ["mol_001"],
            "smiles": ["CCO"],
            "nheavy": [2],
            "status": ["RUNNING"],
            "v3": [None],
            "c_method": [None],
            "energy_window": [None],
            "delta_1": [None],
            "delta_2": [None],
        }
    )
    df.to_csv(csv_path, index=False)

    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()

    batch_settings = {
        "v3": True,
        "c_method": "gfn2-xtb",
        "energy_window": 5.0,
    }

    # Empty HOF values (should produce NaN deltas)
    success = CSVManagerExtensions.update_molecule_with_mopac_results(
        csv_manager=csv_manager,
        mol_id="mol_001",
        h298_cbs=None,
        h298_pm7=None,
        mopac_hof_values=[],  # Empty! Should give NaN deltas
        batch_settings=batch_settings,
    )

    assert success

    # Settings should still persist even with NaN deltas
    df_updated = pd.read_csv(csv_path)
    assert df_updated.loc[0, "v3"] == True  # noqa: E712
    assert df_updated.loc[0, "c_method"] == "gfn2-xtb"
    assert df_updated.loc[0, "energy_window"] == 5.0

    # Deltas should be NaN
    assert pd.isna(df_updated.loc[0, "delta_1"])
    assert pd.isna(df_updated.loc[0, "delta_2"])


def test_update_extra_fields_method(tmp_path):
    """Test the new _update_extra_fields() method directly."""
    csv_path = tmp_path / "test.csv"
    df = pd.DataFrame(
        {
            "mol_id": ["mol_001", "mol_002"],
            "smiles": ["CCO", "CC"],
            "nheavy": [2, 1],
            "status": ["OK", "RUNNING"],
            "v3": [None, None],
            "delta_1": [None, None],
        }
    )
    df.to_csv(csv_path, index=False)

    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()

    # Update extra fields for mol_001
    csv_manager._update_extra_fields(
        mol_id="mol_001",
        field_updates={
            "v3": False,
            "delta_1": 0.0,
        },
    )

    # Verify only mol_001 updated
    df_updated = pd.read_csv(csv_path)
    assert df_updated.loc[0, "v3"] == False  # noqa: E712
    assert df_updated.loc[0, "delta_1"] == 0.0
    assert pd.isna(df_updated.loc[1, "v3"])  # mol_002 unchanged
    assert pd.isna(df_updated.loc[1, "delta_1"])


def test_update_extra_fields_unknown_column_warning(tmp_path, caplog):
    """Test that unknown columns generate warnings but don't crash."""
    csv_path = tmp_path / "test.csv"
    df = pd.DataFrame(
        {
            "mol_id": ["mol_001"],
            "smiles": ["CCO"],
            "nheavy": [2],
            "status": ["OK"],
            "v3": [None],
        }
    )
    df.to_csv(csv_path, index=False)

    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()

    # Try to update a column that doesn't exist
    csv_manager._update_extra_fields(
        mol_id="mol_001",
        field_updates={
            "v3": True,
            "fake_column": 123,  # Doesn't exist!
        },
    )

    # Should log warning
    assert "fake_column" in caplog.text
    assert "not in CSV schema" in caplog.text

    # v3 should still be updated
    df_updated = pd.read_csv(csv_path)
    assert df_updated.loc[0, "v3"] == True  # noqa: E712


def test_multiple_molecules_batch_settings(tmp_path):
    """Test that each molecule gets its own settings (integration test)."""
    csv_path = tmp_path / "test.csv"
    df = pd.DataFrame(
        {
            "mol_id": ["mol_001", "mol_002", "mol_003"],
            "smiles": ["CCO", "CC", "CCC"],
            "nheavy": [2, 1, 3],
            "status": ["RUNNING", "RUNNING", "RUNNING"],
            "v3": [None, None, None],
            "c_method": [None, None, None],
            "threads": [None, None, None],
        }
    )
    df.to_csv(csv_path, index=False)

    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()

    # Update each molecule with same settings
    batch_settings = {
        "v3": False,
        "c_method": "gfn2-xtb",
        "threads": 4,
    }

    for mol_id in ["mol_001", "mol_002", "mol_003"]:
        CSVManagerExtensions.update_molecule_with_mopac_results(
            csv_manager=csv_manager,
            mol_id=mol_id,
            h298_cbs=-17.5,
            h298_pm7=-15.3,
            mopac_hof_values=[0.42],
            batch_settings=batch_settings,
        )

    # Verify all molecules have settings
    df_updated = pd.read_csv(csv_path)
    for i in range(3):
        assert df_updated.loc[i, "v3"] == False  # noqa: E712
        assert df_updated.loc[i, "c_method"] == "gfn2-xtb"
        assert df_updated.loc[i, "threads"] == 4
