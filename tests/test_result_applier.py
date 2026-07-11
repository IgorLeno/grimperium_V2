from pathlib import Path

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus, WorkerStatus
from grimperium.crest_pm7.batch.result_applier import BatchResultApplier
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import PM7Config


def _write_scientific_csv(
    path: Path, status: str = MoleculeStatus.RUNNING.value
) -> None:
    pd.DataFrame(
        [
            {
                "mol_id": "mol_001",
                "smiles": "CCO",
                "nheavy": 3,
                "status": status,
                "reruns": 0,
                "crest_status": "NOT_ATTEMPTED",
                "mopac_status": "NOT_ATTEMPTED",
            }
        ]
    ).to_csv(path, index=False)


def _write_state_csv(path: Path, *, reruns: int = 0) -> None:
    pd.DataFrame(
        [
            {
                "mol_id": "mol_001",
                "smiles": "CCO",
                "status": MoleculeStatus.RUNNING.value,
                "reruns": reruns,
                "assigned_worker": "worker-1",
                "worker_status": WorkerStatus.ONLINE.value,
                "assigned_at": "2026-01-01T00:00:00+00:00",
            }
        ]
    ).to_csv(path, index=False)


def _applier(tmp_path: Path, *, reruns: int = 0) -> BatchResultApplier:
    csv_path = tmp_path / "thermo_pm7.csv"
    state_path = tmp_path / "batch_state.csv"
    _write_scientific_csv(csv_path)
    _write_state_csv(state_path, reruns=reruns)
    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()
    state_manager = BatchStateManager(state_path, PM7Config(), max_reruns=3)
    return BatchResultApplier(state_manager=state_manager, csv_manager=csv_manager)


def test_apply_success_updates_state_first_and_mirrors_scientific_csv(
    tmp_path: Path,
) -> None:
    applier = _applier(tmp_path)

    decision = applier.apply_success("mol_001", {"H298_pm7": -55.25})

    assert decision.final_status == MoleculeStatus.OK.value
    state = pd.read_csv(tmp_path / "batch_state.csv", keep_default_na=False)
    scientific = pd.read_csv(tmp_path / "thermo_pm7.csv")
    assert state.loc[0, "status"] == MoleculeStatus.OK.value
    assert state.loc[0, "assigned_worker"] == ""
    assert state.loc[0, "worker_status"] == WorkerStatus.UNASSIGNED.value
    assert state.loc[0, "assigned_at"] == ""
    assert scientific.loc[0, "status"] == MoleculeStatus.OK.value
    assert float(scientific.loc[0, "H298_pm7"]) == -55.25


def test_apply_failure_lets_state_decide_rerun_and_csv_only_mirrors(
    tmp_path: Path,
) -> None:
    applier = _applier(tmp_path)

    decision = applier.apply_failure("mol_001", "MOPAC failed")

    assert decision.final_status == MoleculeStatus.RERUN.value
    assert decision.reruns == 1
    state = pd.read_csv(tmp_path / "batch_state.csv", keep_default_na=False)
    scientific = pd.read_csv(tmp_path / "thermo_pm7.csv")
    assert state.loc[0, "status"] == MoleculeStatus.RERUN.value
    assert state.loc[0, "assigned_worker"] == ""
    assert scientific.loc[0, "status"] == MoleculeStatus.RERUN.value
    assert int(scientific.loc[0, "reruns"]) == 1


def test_apply_failure_mirrors_skip_after_retry_limit(tmp_path: Path) -> None:
    applier = _applier(tmp_path, reruns=2)

    decision = applier.apply_failure("mol_001", "CREST timeout")

    assert decision.final_status == MoleculeStatus.SKIP.value
    assert decision.reruns == 3
    scientific = pd.read_csv(tmp_path / "thermo_pm7.csv")
    assert scientific.loc[0, "status"] == MoleculeStatus.SKIP.value
    assert int(scientific.loc[0, "reruns"]) == 3
