from __future__ import annotations

from pathlib import Path

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus, WorkerStatus
from grimperium.crest_pm7.batch.output_contracts import BATCH_STATE_COLUMNS
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import PM7Config
from tests.test_csv_schema import EXPECTED_COLUMNS

LEGACY_STATE_ONLY_COLUMNS = ["assigned_worker", "worker_status", "assigned_at"]


def test_batch_csv_manager_creates_legacy_thermo_pm7_with_61_columns(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    manager = BatchCSVManager(csv_path)

    manager.create_empty_csv()

    df = pd.read_csv(csv_path)
    assert len(df.columns) == 61
    assert list(df.columns) == EXPECTED_COLUMNS


def test_legacy_thermo_pm7_excludes_state_only_columns(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    BatchCSVManager(csv_path).create_empty_csv()

    df = pd.read_csv(csv_path)

    for column in LEGACY_STATE_ONLY_COLUMNS:
        assert column not in df.columns


def test_state_only_columns_are_present_in_batch_state_csv(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"

    BatchStateManager(state_path, PM7Config()).update_molecule_status(
        "mol-001",
        MoleculeStatus.ASSIGNED.value,
        worker_id="worker-a",
        extra_fields={"worker_status": WorkerStatus.ONLINE.value},
    )

    state_df = pd.read_csv(state_path, keep_default_na=False)
    for column in LEGACY_STATE_ONLY_COLUMNS:
        assert column in state_df.columns
        assert column in BATCH_STATE_COLUMNS


def test_legacy_schema_matches_existing_61_column_contract() -> None:
    schema = BatchCSVManager(csv_path=None).get_schema()

    assert schema == EXPECTED_COLUMNS
    assert len(schema) == 61


def test_mark_success_persists_scientific_columns_without_state_columns(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    manager = BatchCSVManager(csv_path)
    schema = manager.get_schema() + LEGACY_STATE_ONLY_COLUMNS
    row = dict.fromkeys(schema)
    row.update(
        {
            "mol_id": "mol-001",
            "smiles": "CCO",
            "nheavy": 3,
            "status": MoleculeStatus.RUNNING.value,
            "assigned_worker": "stale-worker",
            "worker_status": WorkerStatus.OFFLINE.value,
            "assigned_at": "2026-06-24T00:00:00+00:00",
        }
    )
    pd.DataFrame([row], columns=schema).to_csv(csv_path, index=False)

    manager.load_csv()
    manager.mark_success(
        "mol-001",
        {
            "H298_pm7": -55.0,
            "H298_predicted": -54.5,
            "delta_correction": 0.5,
            "assigned_worker": "worker-a",
            "worker_status": WorkerStatus.ONLINE.value,
            "assigned_at": "2026-06-25T00:00:00+00:00",
        },
    )

    df = pd.read_csv(csv_path)
    assert df.loc[0, "status"] == MoleculeStatus.OK.value
    assert df.loc[0, "H298_pm7"] == -55.0
    assert df.loc[0, "H298_predicted"] == -54.5
    assert df.loc[0, "delta_correction"] == 0.5
    for column in LEGACY_STATE_ONLY_COLUMNS:
        assert column not in df.columns
