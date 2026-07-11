from pathlib import Path

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import BatchSortingStrategy, MoleculeStatus
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import PM7Config


def test_state_manager_marks_selected_without_worker_assignment(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    state_manager = BatchStateManager(state_path, PM7Config())
    state_manager.reconcile_molecules(
        [{"mol_id": "mol_001", "smiles": "CCO", "charge": -1, "multiplicity": 2}]
    )

    state_manager.mark_selected(["mol_001"], batch_id="batch_0001")

    row = pd.read_csv(state_path, keep_default_na=False).iloc[0]
    assert row["status"] == MoleculeStatus.SELECTED.value
    assert row["assigned_worker"] == ""
    assert row["assigned_at"] == ""
    assert int(row["charge"]) == -1
    assert int(row["multiplicity"]) == 2


def test_csv_selection_can_be_mirrored_into_batch_state(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        [
            {
                "mol_id": "mol_001",
                "smiles": "CCO",
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "charge": 1,
                "multiplicity": 2,
                "reruns": 0,
            }
        ]
    ).to_csv(csv_path, index=False)
    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()
    state_manager = BatchStateManager(tmp_path / "batch_state.csv", PM7Config())
    state_manager.reconcile_molecules(csv_manager.state_seed_rows())

    batch = csv_manager.select_batch(
        batch_id="batch_0001",
        batch_size=1,
        crest_timeout_minutes=30,
        mopac_timeout_minutes=60,
        strategy=BatchSortingStrategy.RERUN_FIRST_THEN_EASY,
    )
    state_manager.mark_selected_from_batch(batch)

    row = pd.read_csv(tmp_path / "batch_state.csv", keep_default_na=False).iloc[0]
    assert row["status"] == MoleculeStatus.SELECTED.value
    assert row["batch_id"] == "batch_0001"
    assert row["assigned_worker"] == ""
    assert batch.molecules[0].charge == 1
    assert batch.molecules[0].multiplicity == 2
