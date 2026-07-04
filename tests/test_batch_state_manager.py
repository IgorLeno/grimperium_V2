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
    _row_with_smiles = {**_row("mol-001"), "smiles": "CCO"}
    _write_state_csv(
        state_path,
        [
            _row_with_smiles,
            _row("mol-002"),
        ],
    )

    result = _manager(state_path).claim_single_molecule("worker-a")

    assert result is not None
    mol_id, smiles = result
    assert mol_id == "mol-001"
    assert smiles == "CCO"
    rows = _read_state_csv(state_path)
    first = rows.loc[rows["mol_id"] == "mol-001"].iloc[0]
    second = rows.loc[rows["mol_id"] == "mol-002"].iloc[0]
    assert first["status"] == MoleculeStatus.ASSIGNED.value
    assert first["assigned_worker"] == "worker-a"
    assert first["worker_status"] == WorkerStatus.ONLINE.value
    assert first["assigned_at"]
    assert second["status"] == MoleculeStatus.PENDING.value


def test_claim_single_molecule_also_claims_rerun_molecules(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(
        state_path,
        [_row("mol-001", MoleculeStatus.RERUN.value)],
    )

    result = _manager(state_path).claim_single_molecule("worker-a")

    assert result is not None
    mol_id, _ = result
    assert mol_id == "mol-001"
    rows = _read_state_csv(state_path)
    assert rows.iloc[0]["status"] == MoleculeStatus.ASSIGNED.value


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

    result = _manager(state_path).claim_single_molecule("worker-a")

    assert result is not None
    mol_id, _ = result
    assert mol_id == "state-mol"
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
    results: list[tuple[str, str] | None] = []

    def claim(worker_id: str) -> None:
        results.append(manager.claim_single_molecule(worker_id))

    threads = [Thread(target=claim, args=(worker_id,)) for worker_id in ["a", "b"]]
    _run_threads(threads)

    # Exactly one worker gets the molecule; the other gets None.
    assert results.count(None) == 1
    claimed = [r for r in results if r is not None]
    assert len(claimed) == 1
    mol_id, _ = claimed[0]
    assert mol_id == "mol-001"
    row = _read_state_csv(state_path).iloc[0]
    assert row["status"] == MoleculeStatus.ASSIGNED.value
    assert row["assigned_worker"] in {"a", "b"}


def _run_threads(threads: list[Thread]) -> None:
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()


# ── reset_stuck_running ───────────────────────────────────────────────────────


def test_reset_stuck_running_resets_running_molecules_to_pending(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "batch_state.csv"
    running = _row("mol-001", MoleculeStatus.RUNNING.value)
    running.update(
        {"assigned_worker": "worker-a", "worker_status": WorkerStatus.ONLINE.value}
    )
    pending = _row("mol-002", MoleculeStatus.PENDING.value)
    ok_mol = _row("mol-003", MoleculeStatus.OK.value)
    _write_state_csv(state_path, [running, pending, ok_mol])

    count = _manager(state_path).reset_stuck_running()

    assert count == 1
    rows = _read_state_csv(state_path).set_index("mol_id")
    assert rows.loc["mol-001", "status"] == MoleculeStatus.PENDING.value
    assert rows.loc["mol-001", "assigned_worker"] == ""
    assert rows.loc["mol-001", "worker_status"] == WorkerStatus.UNASSIGNED.value
    assert rows.loc["mol-002", "status"] == MoleculeStatus.PENDING.value
    assert rows.loc["mol-003", "status"] == MoleculeStatus.OK.value


def test_reset_stuck_running_returns_zero_when_no_running_molecules(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(
        state_path, [_row("mol-001"), _row("mol-002", MoleculeStatus.OK.value)]
    )

    assert _manager(state_path).reset_stuck_running() == 0


# ── mark_worker_offline ───────────────────────────────────────────────────────


def test_mark_worker_offline_reclaims_assigned_and_running_molecules(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "batch_state.csv"
    assigned = _row("mol-001", MoleculeStatus.ASSIGNED.value)
    assigned.update(
        {"assigned_worker": "worker-x", "worker_status": WorkerStatus.ONLINE.value}
    )
    running = _row("mol-002", MoleculeStatus.RUNNING.value)
    running.update(
        {"assigned_worker": "worker-x", "worker_status": WorkerStatus.ONLINE.value}
    )
    other = _row("mol-003", MoleculeStatus.ASSIGNED.value)
    other.update(
        {"assigned_worker": "worker-y", "worker_status": WorkerStatus.ONLINE.value}
    )
    _write_state_csv(state_path, [assigned, running, other])

    count = _manager(state_path).mark_worker_offline("worker-x")

    assert count == 2
    rows = _read_state_csv(state_path).set_index("mol_id")
    assert rows.loc["mol-001", "status"] == MoleculeStatus.PENDING.value
    assert rows.loc["mol-001", "assigned_worker"] == ""
    assert rows.loc["mol-001", "worker_status"] == WorkerStatus.UNASSIGNED.value
    assert rows.loc["mol-002", "status"] == MoleculeStatus.PENDING.value
    assert rows.loc["mol-003", "status"] == MoleculeStatus.ASSIGNED.value  # untouched


def test_mark_worker_offline_returns_zero_for_unknown_worker(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(state_path, [_row("mol-001")])

    assert _manager(state_path).mark_worker_offline("nonexistent-worker") == 0


# ── seed_from_mol_list ────────────────────────────────────────────────────────


def test_seed_from_mol_list_creates_batch_state_csv(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"

    seeded = BatchStateManager(state_path, PM7Config()).seed_from_mol_list(
        [
            {"mol_id": "m1", "smiles": "C"},
            {"mol_id": "m2", "smiles": "CC"},
        ]
    )

    assert seeded == 2
    assert state_path.exists()
    rows = _read_state_csv(state_path)
    assert list(rows["mol_id"]) == ["m1", "m2"]
    assert list(rows["smiles"]) == ["C", "CC"]
    assert all(rows["status"] == MoleculeStatus.PENDING.value)


def test_seed_from_mol_list_skips_nonempty_file(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(state_path, [_row("existing-mol")])

    seeded = BatchStateManager(state_path, PM7Config()).seed_from_mol_list(
        [{"mol_id": "new-mol", "smiles": "CCC"}]
    )

    assert seeded == 0
    rows = _read_state_csv(state_path)
    assert list(rows["mol_id"]) == ["existing-mol"]


def test_seed_from_mol_list_preserves_rerun_status(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"

    BatchStateManager(state_path, PM7Config()).seed_from_mol_list(
        [
            {
                "mol_id": "m1",
                "smiles": "C",
                "status": MoleculeStatus.RERUN.value,
                "reruns": "2",
            }
        ]
    )

    rows = _read_state_csv(state_path)
    assert rows.iloc[0]["status"] == MoleculeStatus.RERUN.value
    assert int(rows.iloc[0]["reruns"]) == 2


def test_seed_from_mol_list_preserves_all_provided_fields(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"

    seeded = BatchStateManager(state_path, PM7Config()).seed_from_mol_list(
        [
            {
                "mol_id": "mol-001",
                "smiles": "CCO",
                "status": MoleculeStatus.PENDING.value,
                "charge": -1,
                "multiplicity": 2,
                "method_id": "semiempirical_am1_pm3_pm7",
                "method_version": "0.1.0",
                "batch_id": "batch-001",
                "not_a_column": "ignored",
            }
        ]
    )

    assert seeded == 1
    rows = _read_state_csv(state_path)
    row = rows.iloc[0]
    # Every provided field that is a real column is persisted.
    assert int(row["charge"]) == -1
    assert int(row["multiplicity"]) == 2
    assert str(row["method_id"]) == "semiempirical_am1_pm3_pm7"
    assert str(row["method_version"]) == "0.1.0"
    assert str(row["batch_id"]) == "batch-001"
    assert str(row["smiles"]) == "CCO"
    assert row["status"] == MoleculeStatus.PENDING.value
    # Defaults applied only where fields are absent.
    assert row["worker_status"] == WorkerStatus.UNASSIGNED.value
    assert int(row["reruns"]) == 0
    # Unknown keys are never written as columns.
    assert "not_a_column" not in rows.columns
    assert set(rows.columns) == set(BATCH_STATE_COLUMNS)


# ── reconcile_molecules ───────────────────────────────────────────────────────


def test_reconcile_molecules_creates_file_when_missing(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"

    added = _manager(state_path).reconcile_molecules(
        [{"mol_id": "m1", "smiles": "C"}, {"mol_id": "m2", "smiles": "CC"}]
    )

    assert added == 2
    assert state_path.exists()
    rows = _read_state_csv(state_path)
    assert list(rows["mol_id"]) == ["m1", "m2"]
    assert all(rows["status"] == MoleculeStatus.PENDING.value)


def test_reconcile_molecules_seeds_empty_file(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    _write_state_csv(state_path, [])

    added = _manager(state_path).reconcile_molecules([{"mol_id": "m1", "smiles": "C"}])

    assert added == 1
    assert list(_read_state_csv(state_path)["mol_id"]) == ["m1"]


def test_reconcile_molecules_adds_missing_and_preserves_existing(
    tmp_path: Path,
) -> None:
    state_path = tmp_path / "batch_state.csv"
    running = _row("mol-001", MoleculeStatus.RUNNING.value)
    running.update(
        {
            "assigned_worker": "worker-a",
            "worker_status": WorkerStatus.ONLINE.value,
            "reruns": 1,
            "method_id": "legacy",
            "method_version": "9.9.9",
        }
    )
    _write_state_csv(state_path, [running])

    added = _manager(state_path).reconcile_molecules(
        [
            {"mol_id": "mol-001", "smiles": "CCO", "status": "Pending", "reruns": 0},
            {"mol_id": "mol-002", "smiles": "CC"},
        ]
    )

    assert added == 1
    rows = _read_state_csv(state_path).set_index("mol_id")
    # Existing operational state untouched — never demoted from Running.
    assert rows.loc["mol-001", "status"] == MoleculeStatus.RUNNING.value
    assert rows.loc["mol-001", "assigned_worker"] == "worker-a"
    assert int(rows.loc["mol-001", "reruns"]) == 1
    assert rows.loc["mol-001", "method_id"] == "legacy"
    assert rows.loc["mol-001", "method_version"] == "9.9.9"
    # New molecule appended as Pending.
    assert rows.loc["mol-002", "status"] == MoleculeStatus.PENDING.value


def test_reconcile_molecules_is_idempotent(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    mols = [{"mol_id": "m1", "smiles": "C"}, {"mol_id": "m2", "smiles": "CC"}]

    assert _manager(state_path).reconcile_molecules(mols) == 2
    second = _manager(state_path).reconcile_molecules(mols)

    assert second == 0
    rows = _read_state_csv(state_path)
    assert list(rows["mol_id"]) == ["m1", "m2"]  # no duplicates


def test_reconcile_molecules_fills_only_missing_smiles(tmp_path: Path) -> None:
    state_path = tmp_path / "batch_state.csv"
    no_smiles = _row("mol-001", MoleculeStatus.OK.value)
    no_smiles["smiles"] = ""
    has_smiles = _row("mol-002", MoleculeStatus.OK.value)
    has_smiles["smiles"] = "CC"
    _write_state_csv(state_path, [no_smiles, has_smiles])

    _manager(state_path).reconcile_molecules(
        [
            {"mol_id": "mol-001", "smiles": "CCO"},
            {"mol_id": "mol-002", "smiles": "OVERWRITE"},
        ]
    )

    rows = _read_state_csv(state_path).set_index("mol_id")
    assert rows.loc["mol-001", "smiles"] == "CCO"  # empty identity filled
    assert rows.loc["mol-002", "smiles"] == "CC"  # existing identity preserved
    assert rows.loc["mol-001", "status"] == MoleculeStatus.OK.value  # state kept
