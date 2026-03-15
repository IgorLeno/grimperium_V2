"""Fixtures for ML pipeline tests.

Provides a 10-row synthetic CSV that passes CSVDataLoader validation
(requires mol_id, smiles, nheavy, status with valid MoleculeStatus values).
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import pandas as pd
import pytest

# Features used by the ML pipeline
FEATURE_COLS = (
    "nheavy",
    "rdkit_nrotbonds",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
)

# 10-row synthetic CSV with valid MoleculeStatus values and realistic data.
# 2 rows have NaN in mopac columns to test imputation (~20%).
_SYNTHETIC_CSV = """\
mol_id,smiles,nheavy,status,charge,multiplicity,H298_cbs,H298_pm7,rdkit_nrotbonds,mopac_homo_ev,mopac_lumo_ev,mopac_gap_ev
mol_00001,C,1,OK,0,1,-17.9,-15.2,0,-10.5,1.2,11.7
mol_00002,CC,2,OK,0,1,-20.0,-17.8,0,-10.8,1.5,12.3
mol_00003,CCC,3,OK,0,1,-25.0,-22.5,1,-10.3,1.1,11.4
mol_00004,CCO,3,OK,0,1,-56.2,-52.0,0,-10.9,1.8,12.7
mol_00005,CC(=O)O,4,OK,0,1,-103.4,-98.0,1,-11.2,0.5,11.7
mol_00006,c1ccccc1,6,OK,0,1,19.8,22.5,0,-9.1,0.3,9.4
mol_00007,CC(C)C,4,OK,0,1,-32.1,-29.0,0,-10.6,1.4,12.0
mol_00008,CCCO,4,OK,0,1,-61.0,-57.0,1,-10.7,1.6,12.3
mol_00009,CC=O,3,OK,0,1,-39.7,-36.5,0,,,
mol_00010,COC,3,OK,0,1,-44.0,-41.0,0,,,
"""


@pytest.fixture
def synthetic_csv_path(tmp_path: Path) -> Path:
    """Write synthetic CSV to a temp file and return its path.

    Args:
        tmp_path: Temporary directory provided by pytest.

    Returns:
        Path: Path to the created synthetic CSV file.
    """
    csv_file = tmp_path / "thermo_pm7.csv"
    csv_file.write_text(_SYNTHETIC_CSV)
    return csv_file


@pytest.fixture
def synthetic_df() -> pd.DataFrame:
    """Load the synthetic CSV as a DataFrame (bypassing CSVDataLoader).

    Returns:
        pd.DataFrame: DataFrame loaded from _SYNTHETIC_CSV (bypassing CSVDataLoader).
    """
    return pd.read_csv(StringIO(_SYNTHETIC_CSV))


@pytest.fixture
def feature_cols() -> list[str]:
    """Feature column names used by the pipeline."""
    return list(FEATURE_COLS)
