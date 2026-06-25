from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus, WorkerStatus
from grimperium.crest_pm7.batch.output_contracts import BATCH_STATE_COLUMNS
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import PM7Config


def _state_row(mol_id: str, worker_id: str) -> dict[str, object]:
    row: dict[str, object] = dict.fromkeys(BATCH_STATE_COLUMNS, "")
    row.update(
        {
            "mol_id": mol_id,
            "status": MoleculeStatus.ASSIGNED.value,
            "smiles": "C",
            "nheavy": 1,
            "assigned_worker": worker_id,
            "worker_status": WorkerStatus.ONLINE.value,
            "assigned_at": (
                datetime.now(timezone.utc) - timedelta(minutes=90)
            ).isoformat(),
        }
    )
    return row


def test_batch_csv_manager_no_longer_exposes_distributed_state_api() -> None:
    manager = BatchCSVManager(csv_path=None)

    assert not hasattr(manager, "distribute_molecules")
    assert not hasattr(manager, "mark_worker_offline")
    assert not hasattr(manager, "reassign_offline_molecules")
    assert not hasattr(manager, "reset_stuck_assigned")


def test_batch_csv_manager_schema_excludes_distributed_state_columns() -> None:
    schema = BatchCSVManager(csv_path=None).get_schema()

    assert "assigned_worker" not in schema
    assert "worker_status" not in schema
    assert "assigned_at" not in schema


def test_batch_state_manager_owns_worker_reassignment(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    pd.DataFrame(
        [
            _state_row("mol-001", "worker-offline"),
            _state_row("mol-002", "worker-active"),
        ],
        columns=BATCH_STATE_COLUMNS,
    ).to_csv(state_path, index=False)
    manager = BatchStateManager(state_path, PM7Config())

    reassigned = manager.reassign_offline_molecules(
        active_worker_ids=["worker-active"],
        timeout_minutes=60,
    )

    assert reassigned == ["mol-001"]
    rows = pd.read_csv(state_path, keep_default_na=False).set_index("mol_id")
    assert rows.loc["mol-001", "status"] == MoleculeStatus.PENDING.value
    assert rows.loc["mol-001", "assigned_worker"] == ""
    assert rows.loc["mol-001", "worker_status"] == WorkerStatus.UNASSIGNED.value
    assert rows.loc["mol-002", "status"] == MoleculeStatus.ASSIGNED.value
    assert rows.loc["mol-002", "assigned_worker"] == "worker-active"
