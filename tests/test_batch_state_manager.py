from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
from threading import Thread

import pandas as pd

from grimperium.crest_pm7.batch.enums import MoleculeStatus, WorkerStatus
from grimperium.crest_pm7.batch.output_contracts import BATCH_STATE_COLUMNS
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import PM7Config


def _row(mol_id: str, status: str = MoleculeStatus.PENDING.value) -> dict[str, object]:
    row: dict[str, object] = dict.fromkeys(BATCH_STATE_COLUMNS, "")
    row.update(
        {
            "mol_id": mol_id,
            "status": status,
            "smiles": "C",
            "multiplicity": 1,
            "charge": 0,
            "nheavy": 1,
            "reruns": 0,
            "worker_status": WorkerStatus.UNASSIGNED.value,
        }
    )
    return row


def _write_state_csv(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows, columns=BATCH_STATE_COLUMNS).to_csv(path, index=False)


def _read_state_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, keep_default_na=False)


def _manager(path: Path) -> BatchStateManager:
    return BatchStateManager(path, PM7Config())


def _utc_minutes_ago(minutes: int) -> str:
    return (datetime.now(timezone.utc) - timedelta(minutes=minutes)).isoformat()


def test_claim_single_molecule_assigns_first_pending_molecule(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(
        state_path,
        [
            _row("mol-001"),
            _row("mol-002"),
        ],
    )

    claimed = _manager(state_path).claim_single_molecule("worker-a")

    assert claimed == "mol-001"
    rows = _read_state_csv(state_path)
    first = rows.loc[rows["mol_id"] == "mol-001"].iloc[0]
    second = rows.loc[rows["mol_id"] == "mol-002"].iloc[0]
    assert first["status"] == MoleculeStatus.ASSIGNED.value
    assert first["assigned_worker"] == "worker-a"
    assert first["worker_status"] == WorkerStatus.ONLINE.value
    assert first["assigned_at"]
    assert second["status"] == MoleculeStatus.PENDING.value


def test_claim_single_molecule_returns_none_without_pending_molecules(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(
        state_path,
        [
            _row("mol-001", MoleculeStatus.ASSIGNED.value),
            _row("mol-002", MoleculeStatus.OK.value),
        ],
    )

    assert _manager(state_path).claim_single_molecule("worker-a") is None


def test_distribute_molecules_assigns_round_robin_workers(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(
        state_path,
        [_row("mol-001"), _row("mol-002"), _row("mol-003"), _row("mol-004")],
    )

    assignments = _manager(state_path).distribute_molecules(
        ["mol-001", "mol-002", "mol-003", "mol-004"],
        ["worker-a", "worker-b"],
    )

    assert assignments == {
        "mol-001": "worker-a",
        "mol-002": "worker-b",
        "mol-003": "worker-a",
        "mol-004": "worker-b",
    }
    rows = _read_state_csv(state_path).set_index("mol_id")
    assert rows.loc["mol-001", "assigned_worker"] == "worker-a"
    assert rows.loc["mol-002", "assigned_worker"] == "worker-b"
    assert rows.loc["mol-003", "assigned_worker"] == "worker-a"
    assert rows.loc["mol-004", "assigned_worker"] == "worker-b"
    assert set(rows["status"]) == {MoleculeStatus.ASSIGNED.value}


def test_reassign_offline_molecules_returns_stale_inactive_workers_to_pending(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "batch_state.csv"
    stale_offline = _row("mol-001", MoleculeStatus.ASSIGNED.value)
    stale_offline.update(
        {
            "assigned_worker": "worker-offline",
            "worker_status": WorkerStatus.ONLINE.value,
            "assigned_at": _utc_minutes_ago(90),
        }
    )
    active_worker = _row("mol-002", MoleculeStatus.ASSIGNED.value)
    active_worker.update(
        {
            "assigned_worker": "worker-active",
            "worker_status": WorkerStatus.ONLINE.value,
            "assigned_at": _utc_minutes_ago(90),
        }
    )
    fresh_offline = _row("mol-003", MoleculeStatus.ASSIGNED.value)
    fresh_offline.update(
        {
            "assigned_worker": "worker-offline",
            "worker_status": WorkerStatus.ONLINE.value,
            "assigned_at": _utc_minutes_ago(5),
        }
    )
    _write_state_csv(state_path, [stale_offline, active_worker, fresh_offline])

    reassigned = _manager(state_path).reassign_offline_molecules(
        active_worker_ids=["worker-active"],
        timeout_minutes=30,
    )

    assert reassigned == ["mol-001"]
    rows = _read_state_csv(state_path).set_index("mol_id")
    assert rows.loc["mol-001", "status"] == MoleculeStatus.PENDING.value
    assert rows.loc["mol-001", "assigned_worker"] == ""
    assert rows.loc["mol-001", "worker_status"] == WorkerStatus.UNASSIGNED.value
    assert rows.loc["mol-001", "assigned_at"] == ""
    assert rows.loc["mol-002", "status"] == MoleculeStatus.ASSIGNED.value
    assert rows.loc["mol-003", "status"] == MoleculeStatus.ASSIGNED.value


def test_reset_stuck_assigned_clears_old_assignments(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    stuck = _row("mol-001", MoleculeStatus.ASSIGNED.value)
    stuck.update(
        {
            "assigned_worker": "worker-a",
            "worker_status": WorkerStatus.ONLINE.value,
            "assigned_at": _utc_minutes_ago(120),
        }
    )
    fresh = _row("mol-002", MoleculeStatus.ASSIGNED.value)
    fresh.update(
        {
            "assigned_worker": "worker-b",
            "worker_status": WorkerStatus.ONLINE.value,
            "assigned_at": _utc_minutes_ago(10),
        }
    )
    _write_state_csv(state_path, [stuck, fresh])

    reset = _manager(state_path).reset_stuck_assigned(timeout_minutes=60)

    assert reset == ["mol-001"]
    rows = _read_state_csv(state_path).set_index("mol_id")
    assert rows.loc["mol-001", "status"] == MoleculeStatus.PENDING.value
    assert rows.loc["mol-001", "assigned_worker"] == ""
    assert rows.loc["mol-001", "worker_status"] == WorkerStatus.UNASSIGNED.value
    assert rows.loc["mol-001", "assigned_at"] == ""
    assert rows.loc["mol-002", "status"] == MoleculeStatus.ASSIGNED.value


def test_state_manager_uses_batch_state_csv_not_thermo_pm7(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    legacy_path = tmp_path / "thermo_pm7.csv"
    _write_state_csv(state_path, [_row("state-mol")])
    _write_state_csv(legacy_path, [_row("legacy-mol")])
    before_legacy = legacy_path.read_text(encoding="utf-8")

    claimed = _manager(state_path).claim_single_molecule("worker-a")

    assert claimed == "state-mol"
    assert legacy_path.read_text(encoding="utf-8") == before_legacy
    assert "state-mol" in state_path.read_text(encoding="utf-8")


def test_get_pending_molecules_filters_by_status(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(
        state_path,
        [
            _row("mol-001", MoleculeStatus.PENDING.value),
            _row("mol-002", MoleculeStatus.ASSIGNED.value),
            _row("mol-003", MoleculeStatus.PENDING.value),
        ],
    )

    manager = _manager(state_path)

    assert manager.get_pending_molecules() == ["mol-001", "mol-003"]
    assert manager.get_molecules_by_status(MoleculeStatus.ASSIGNED.value) == ["mol-002"]


def test_update_molecule_status_writes_operational_fields(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(state_path, [_row("mol-001")])

    _manager(state_path).update_molecule_status(
        "mol-001",
        MoleculeStatus.ASSIGNED.value,
        worker_id="worker-a",
        extra_fields={"worker_status": WorkerStatus.ONLINE.value},
    )

    row = _read_state_csv(state_path).iloc[0]
    assert row["status"] == MoleculeStatus.ASSIGNED.value
    assert row["assigned_worker"] == "worker-a"
    assert row["worker_status"] == WorkerStatus.ONLINE.value


def test_concurrent_claim_allows_only_one_worker_per_pending_molecule(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(state_path, [_row("mol-001")])
    manager = _manager(state_path)
    results: list[str | None] = []

    def claim(worker_id: str) -> None:
        results.append(manager.claim_single_molecule(worker_id))

    threads = [Thread(target=claim, args=(worker_id,)) for worker_id in ["a", "b"]]
    _run_threads(threads)

    assert sorted(results, key=lambda value: value or "") == [None, "mol-001"]
    row = _read_state_csv(state_path).iloc[0]
    assert row["status"] == MoleculeStatus.ASSIGNED.value
    assert row["assigned_worker"] in {"a", "b"}


def _run_threads(threads: list[Thread]) -> None:
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
