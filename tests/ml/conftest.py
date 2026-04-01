"""Fixtures for ML pipeline tests.

Provides a 10-row synthetic CSV that passes CSVDataLoader validation
(requires mol_id, smiles, nheavy, status with valid MoleculeStatus values).
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

# Features used by the ML pipeline (31 validated columns)
FEATURE_COLS = (
    "nheavy",
    "multiplicity",
    "charge",
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
    "crest_conformers_generated",
    "num_conformers_selected",
    "crest_time_s",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
    "mopac_dipole_debye",
    "mopac_ionization_potential_ev",
    "mopac_cosmo_area_a2",
    "mopac_cosmo_volume_a3",
    "mopac_gradient_norm",
    "mopac_num_scf_cycles",
    "H298_pm7",
)

# 10-row synthetic CSV with valid MoleculeStatus values and realistic data.
# 2 rows have NaN in mopac columns to test imputation (~20%).
_SYNTHETIC_CSV = """\
mol_id,smiles,nheavy,status,cbs_quality_flag,charge,multiplicity,H298_cbs,H298_pm7,rdkit_nrotbonds,rdkit_tpsa,rdkit_num_rings,rdkit_fsp3,rdkit_mol_weight,rdkit_hbond_donors,rdkit_hbond_acceptors,rdkit_nC,rdkit_nH,rdkit_nO,rdkit_nN,rdkit_bonds_single,rdkit_bonds_double,rdkit_bonds_triple,rdkit_bonds_aromatic,crest_conformers_generated,num_conformers_selected,crest_time_s,mopac_homo_ev,mopac_lumo_ev,mopac_gap_ev,mopac_dipole_debye,mopac_ionization_potential_ev,mopac_cosmo_area_a2,mopac_cosmo_volume_a3,mopac_gradient_norm,mopac_num_scf_cycles
mol_00001,C,1,OK,OK,0,1,-17.9,-15.2,0,0.0,0,1.0,16.04,0,0,1,4,0,0,4,0,0,0,5,1,1.2,-10.5,1.2,11.7,0.0,12.5,42.1,28.3,0.01,8
mol_00002,CC,2,OK,OK,0,1,-20.0,-17.8,0,0.0,0,1.0,30.07,0,0,2,6,0,0,7,0,0,0,8,2,2.1,-10.8,1.5,12.3,0.0,12.8,58.2,45.1,0.02,10
mol_00003,CCC,3,OK,OK,0,1,-25.0,-22.5,1,0.0,0,1.0,44.10,0,0,3,8,0,0,10,0,0,0,12,3,3.5,-10.3,1.1,11.4,0.1,12.3,74.5,62.0,0.01,9
mol_00004,CCO,3,OK,OK,0,1,-56.2,-52.0,0,20.2,0,0.67,46.07,1,1,2,6,1,0,8,0,0,0,10,2,2.8,-10.9,1.8,12.7,1.7,13.0,68.9,55.2,0.03,11
mol_00005,CC(=O)O,4,OK,OK,0,1,-103.4,-98.0,1,37.3,0,0.25,60.05,1,2,2,4,2,0,6,1,0,0,7,1,4.2,-11.2,0.5,11.7,1.9,13.5,78.1,60.3,0.02,12
mol_00006,c1ccccc1,6,OK,OK,0,1,19.8,22.5,0,0.0,1,0.0,78.11,0,0,6,6,0,0,0,0,0,6,3,1,1.5,-9.1,0.3,9.4,0.0,11.0,95.3,80.1,0.01,7
mol_00007,CC(C)C,4,OK,OK,0,1,-32.1,-29.0,0,0.0,0,1.0,58.12,0,0,4,10,0,0,13,0,0,0,15,4,5.0,-10.6,1.4,12.0,0.0,12.6,88.7,75.4,0.02,9
mol_00008,CCCO,4,OK,OK,0,1,-61.0,-57.0,1,20.2,0,0.75,60.10,1,1,3,8,1,0,11,0,0,0,14,3,4.1,-10.7,1.6,12.3,1.5,12.9,85.2,70.8,0.01,10
mol_00009,CC=O,3,OK,OK,0,1,-39.7,-36.5,0,17.1,0,0.5,44.05,0,1,2,4,1,0,5,1,0,0,6,1,1.8,,,,,,,,,
mol_00010,COC,3,OK,OK,0,1,-44.0,-41.0,0,18.5,0,0.67,46.07,0,2,2,6,1,0,8,0,0,0,9,2,2.0,,,,,,,,,
mol_00011,CN,2,OK,SUSPECT,0,1,-99999.0,-20.0,0,0.0,0,1.0,29.04,1,0,1,5,0,1,4,0,0,0,6,1,1.0,-9.0,1.0,10.0,1.0,11.0,50.0,35.0,0.01,8
mol_00012,CO,2,PENDING,OK,0,1,-45.0,-42.0,0,18.5,0,0.5,32.04,1,1,1,4,1,0,5,0,0,0,6,1,1.5,-10.0,1.2,11.2,1.8,12.0,55.0,40.0,0.02,9
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
def synthetic_csv_path_mixed(tmp_path: Path) -> Path:
    """12-row CSV: 10 OK/OK, 1 OK/SUSPECT, 1 PENDING/OK."""
    csv_file = tmp_path / "thermo_pm7_mixed.csv"
    csv_file.write_text(_SYNTHETIC_CSV)
    return csv_file


@pytest.fixture(scope="module")
def _module_csv(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Write the synthetic CSV once per module for read-only model training tests."""
    base = tmp_path_factory.mktemp("ml_module")
    csv_file = base / "thermo_pm7.csv"
    csv_file.write_text(_SYNTHETIC_CSV)
    return csv_file


@pytest.fixture
def synthetic_df() -> pd.DataFrame:
    """Load the synthetic CSV as a DataFrame (bypassing CSVDataLoader).

    Returns:
        pd.DataFrame: DataFrame loaded from _SYNTHETIC_CSV
        (bypassing CSVDataLoader).
    """
    return pd.read_csv(StringIO(_SYNTHETIC_CSV))


@pytest.fixture
def feature_cols() -> list[str]:
    """Feature column names used by the pipeline.

    Returns:
        list[str]: Feature column names sourced from FEATURE_COLS.
    """
    return list(FEATURE_COLS)


@pytest.fixture
def csv_with_missing_target(tmp_path: Path) -> Path:
    """Create a CSV fixture with one row missing target value.

    Args:
        tmp_path: Temporary directory provided by pytest.

    Returns:
        Path: CSV path where one row has missing `H298_pm7`.
    """
    csv_content = (
        "mol_id,smiles,nheavy,status,cbs_quality_flag,charge,multiplicity,"
        "H298_cbs,H298_pm7,rdkit_nrotbonds,mopac_homo_ev,"
        "mopac_lumo_ev,mopac_gap_ev\n"
        "mol_01,C,1,OK,OK,0,1,-17.9,-15.2,0,-10.5,1.2,11.7\n"
        "mol_02,CC,2,OK,OK,0,1,-20.0,,0,-10.8,1.5,12.3\n"
    )
    csv_file = tmp_path / "missing_target.csv"
    csv_file.write_text(csv_content)
    return csv_file


@pytest.fixture
def csv_no_ok_status(tmp_path: Path) -> Path:
    """Create a CSV fixture with no rows in status OK.

    Args:
        tmp_path: Temporary directory provided by pytest.

    Returns:
        Path: CSV path where all rows are filtered out by status.
    """
    csv_content = (
        "mol_id,smiles,nheavy,status,cbs_quality_flag,charge,multiplicity,"
        "H298_cbs,H298_pm7,rdkit_nrotbonds,mopac_homo_ev,"
        "mopac_lumo_ev,mopac_gap_ev\n"
        "mol_01,C,1,Pending,OK,0,1,-17.9,-15.2,0,-10.5,1.2,11.7\n"
    )
    csv_file = tmp_path / "no_ok.csv"
    csv_file.write_text(csv_content)
    return csv_file


@pytest.fixture(scope="module")
def trained_model_fixture(_module_csv: Path) -> dict[str, Any]:
    """Train a model on synthetic data and return a save-ready bundle.

    Args:
        _module_csv: Path to the module-scoped synthetic CSV fixture.

    Returns:
        dict with keys ``learner``, ``pipeline``, and ``metrics``
        (``{"train": {...}, "test": {...}}``).
    """
    from grimperium.ml.trainer import train

    learner, train_metrics, test_metrics, pipeline = train(
        _module_csv,
        test_size=0.2,
        random_state=42,
        return_pipeline=True,
    )
    return {
        "learner": learner,
        "pipeline": pipeline,
        "metrics": {
            "train": train_metrics,
            "test": test_metrics,
        },
    }


@pytest.fixture
def csv_with_predictions(tmp_path: Path) -> Path:
    """Copy synthetic CSV (12 rows: 10 OK/OK + 2 non-eligible) for write testing.

    Args:
        tmp_path: Temporary directory provided by pytest.

    Returns:
        Path: Writable CSV path for batch prediction tests.
    """
    csv_file = tmp_path / "predict_target.csv"
    csv_file.write_text(_SYNTHETIC_CSV)
    return csv_file


@pytest.fixture
def csv_with_predictions_for_charts(
    tmp_path: Path, synthetic_csv_path_mixed: Path
) -> Path:
    """Train a model, predict, and return CSV path with H298_predicted column.

    This fixture trains on the synthetic mixed CSV, saves the model, runs
    ``predict_batch``, and returns the updated CSV path ready for chart tests.

    Args:
        tmp_path: Temporary directory provided by pytest.
        synthetic_csv_path_mixed: Path to the 12-row mixed-status CSV.

    Returns:
        Path: CSV path containing ``H298_predicted`` column with predictions.
    """
    from grimperium.ml.persistence import save_model
    from grimperium.ml.predictor import predict_batch
    from grimperium.ml.trainer import train

    # Train model on synthetic data
    learner, train_metrics, test_metrics, pipeline = train(
        synthetic_csv_path_mixed,
        test_size=0.2,
        random_state=42,
        return_pipeline=True,
    )

    # Save model to tmp_path
    model_path = tmp_path / "charts_model.joblib"
    bundle: dict[str, Any] = {
        "learner": learner,
        "pipeline": pipeline,
        "metrics": {"train": train_metrics, "test": test_metrics},
    }
    save_model(bundle, model_path)

    # Run batch prediction (overwrites CSV with H298_predicted column)
    predict_batch(synthetic_csv_path_mixed, model_path)

    return synthetic_csv_path_mixed
