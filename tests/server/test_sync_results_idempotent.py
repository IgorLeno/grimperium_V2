from pathlib import Path

import pandas as pd
import pytest
from httpx import ASGITransport, AsyncClient

from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.batch.result_ledger import build_result_fingerprint
from grimperium.server.app import create_app
from grimperium.server.config import ServerConfig
from grimperium.server.models import SyncResult
from grimperium.server.sync_application import (
    SyncResultApplicationService,
    apply_worker_result,
    is_active_assignment_status,
    sync_result_payload,
)


def _write_csv(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "mol_id": "mol_001",
                "smiles": "CCO",
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            }
        ]
    ).to_csv(path, index=False)


async def _client(csv_path: Path) -> AsyncClient:
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    return AsyncClient(transport=ASGITransport(app=app), base_url="http://test")


async def _claim(client: AsyncClient) -> tuple[str, str]:
    await client.post("/register", json={"worker_id": "w1", "hostname": "lab"})
    await client.post("/dispatch/start")
    response = await client.post("/claim", json={"worker_id": "w1"})
    body = response.json()
    mol_id = body["mol_id"]
    attempt_id = body["attempt_id"]
    assert mol_id is not None
    assert attempt_id is not None
    return str(mol_id), str(attempt_id)


@pytest.mark.anyio
async def test_sync_results_duplicate_result_id_is_not_applied_twice(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "result_id": "result-1",
                    "mol_id": mol_id,
                    "attempt_id": attempt_id,
                    "success": False,
                    "result_update": None,
                    "error": "MOPAC timeout",
                    "completed_at": "2026-04-21T10:00:00Z",
                }
            ],
        }

        first = await client.post("/sync_results", json=payload)
        second = await client.post("/sync_results", json=payload)

    assert first.status_code == 200
    assert first.json()["accepted"] == 1
    assert second.status_code == 200
    assert second.json()["duplicate"] is True
    state = pd.read_csv(tmp_path / "batch_state.csv")
    assert int(state.loc[0, "reruns"]) == 1


@pytest.mark.anyio
async def test_sync_results_legacy_fallback_without_result_id_is_stable(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "mol_id": mol_id,
                    "attempt_id": attempt_id,
                    "success": False,
                    "result_update": None,
                    "error": "MOPAC timeout",
                    "completed_at": "2026-04-21T10:00:00Z",
                }
            ],
        }
        first = await client.post("/sync_results", json=payload)
        second = await client.post("/sync_results", json=payload)

    assert first.status_code == 200
    assert first.json()["accepted"] == 1
    assert second.status_code == 200
    assert second.json()["duplicate"] is True
    state = pd.read_csv(tmp_path / "batch_state.csv")
    assert int(state.loc[0, "reruns"]) == 1


@pytest.mark.anyio
async def test_sync_results_conflicting_result_id_returns_conflict_item(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        base_result = {
            "result_id": "result-1",
            "mol_id": mol_id,
            "attempt_id": attempt_id,
            "success": False,
            "result_update": None,
            "error": "MOPAC timeout",
            "completed_at": "2026-04-21T10:00:00Z",
        }
        first = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [base_result]}
        )
        conflicting = dict(base_result)
        conflicting["success"] = True
        conflicting["result_update"] = {"H298_pm7": -55.0}
        conflicting["error"] = None
        second = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [conflicting]}
        )

    assert first.status_code == 200
    assert second.status_code == 200
    assert second.json()["items"][0]["status"] == "conflict"
    assert second.json()["rejected"] == 1


@pytest.mark.anyio
async def test_sync_results_updates_worker_registry_once(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "result_id": "result-success",
                    "mol_id": mol_id,
                    "attempt_id": attempt_id,
                    "success": True,
                    "result_update": {"H298_pm7": -55.0},
                    "error": None,
                    "completed_at": "2026-04-21T10:00:00Z",
                }
            ],
        }
        first = await client.post("/sync_results", json=payload)
        second = await client.post("/sync_results", json=payload)
        status = await client.get("/workers/status")

    assert first.status_code == 200
    assert second.json()["duplicate"] is True
    workers = status.json()
    worker = next(w for w in workers if w["worker_id"] == "w1")
    assert worker["successful"] == 1
    assert worker["processed"] == 1
    assert worker["current_mol_id"] is None


@pytest.mark.anyio
async def test_sync_results_recovers_prepared_without_duplicating(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        mol_id, attempt_id = await _claim(client)
        result_id = "result-crash"
        fingerprint_payload = {
            "mol_id": mol_id,
            "attempt_id": attempt_id,
            "success": False,
            "result_update": None,
            "error": "MOPAC timeout",
            "completed_at": "2026-04-21T10:00:00Z",
        }
        fingerprint = build_result_fingerprint(fingerprint_payload)
        ledger = app.state.result_ledger
        ledger.prepare(
            result_id=result_id,
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=False,
            previous_status="Running",
            previous_reruns=0,
        )
        # Dual-write completed; commit never happened (crash window).
        app.state.state_manager.mark_rerun(mol_id, "MOPAC timeout")

        payload = {
            "worker_id": "w1",
            "results": [{"result_id": result_id, **fingerprint_payload}],
        }
        response = await client.post("/sync_results", json=payload)

    assert response.status_code == 200
    state = pd.read_csv(tmp_path / "batch_state.csv")
    assert int(state.loc[0, "reruns"]) == 1


@pytest.mark.anyio
async def test_sync_results_failure_and_skip_update_registry(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        [
            {
                "mol_id": "mol_001",
                "smiles": "CCO",
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            },
            {
                "mol_id": "mol_002",
                "smiles": "CCC",
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            },
        ]
    ).to_csv(csv_path, index=False)

    async with await _client(csv_path) as client:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab"})
        await client.post("/dispatch/start")
        first = await client.post("/claim", json={"worker_id": "w1"})
        mol_fail = first.json()["mol_id"]
        attempt_fail = first.json()["attempt_id"]
        await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "fail-1",
                        "mol_id": mol_fail,
                        "attempt_id": attempt_fail,
                        "success": False,
                        "result_update": None,
                        "error": "timeout",
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )
        second = await client.post("/claim", json={"worker_id": "w1"})
        mol_skip = second.json()["mol_id"]
        attempt_skip = second.json()["attempt_id"]
        # force skip via apply_failure semantics: success=False with special update
        await client.post(
            "/report/failure",
            json={
                "worker_id": "w1",
                "mol_id": mol_skip,
                "error": "skip",
                "force_skip": True,
                "attempt_id": attempt_skip,
            },
        )
        status = await client.get("/workers/status")

    worker = next(w for w in status.json() if w["worker_id"] == "w1")
    assert worker["failed"] == 1
    assert worker["skipped"] == 1
    assert worker["current_mol_id"] is None


@pytest.mark.anyio
async def test_sync_results_marks_failed_when_apply_raises(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as client:
        mol_id, attempt_id = await _claim(client)
        original = app.state.csv_manager

        class Boom:
            def __getattr__(self, name: str):  # noqa: ANN001
                raise RuntimeError("dual-write boom")

        app.state.csv_manager = Boom()  # type: ignore[assignment]
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "result_id": "prepared-fail",
                    "mol_id": mol_id,
                    "attempt_id": attempt_id,
                    "success": True,
                    "result_update": {"H298_pm7": -1.0},
                    "error": None,
                    "completed_at": "2026-04-21T10:00:00Z",
                }
            ],
        }
        response = await client.post("/sync_results", json=payload)
        app.state.csv_manager = original
        failed_entries = [
            entry
            for entry in app.state.result_ledger._journal.values()
            if entry.result_id == "prepared-fail"
        ]
        assert response.status_code == 200
        assert response.json()["rejected"] == 1
        assert failed_entries
        assert failed_entries[-1].txn_status.value == "failed"

        # Retry after restoring applier must recover safely.
        retry = await client.post("/sync_results", json=payload)
        assert retry.status_code == 200
        assert retry.json()["accepted"] == 1
        state = pd.read_csv(tmp_path / "batch_state.csv")
        assert int(state.loc[0, "reruns"]) == 0


@pytest.mark.anyio
async def test_second_result_same_attempt_is_stale_after_ok(tmp_path: Path) -> None:
    """R1 applies; R1 resend is duplicate; R2 same attempt is stale_attempt."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        science = {"H298_pm7": -55.0, "G298_pm7": -60.0, "gap": 5.0}
        r1 = {
            "result_id": "result-r1",
            "mol_id": mol_id,
            "attempt_id": attempt_id,
            "success": True,
            "result_update": science,
            "error": None,
            "completed_at": "2026-04-21T10:00:00Z",
        }
        first = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [r1]}
        )
        dup = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [r1]}
        )
        r2 = dict(r1)
        r2["result_id"] = "result-r2"
        r2["result_update"] = {"H298_pm7": -99.0}
        second = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [r2]}
        )

    assert first.json()["items"][0]["status"] == "applied"
    assert dup.json()["items"][0]["status"] == "duplicate"
    assert second.json()["items"][0]["status"] == "stale_attempt"
    state = pd.read_csv(tmp_path / "batch_state.csv", keep_default_na=False)
    scientific = pd.read_csv(csv_path)
    assert str(state.loc[0, "status"]) == MoleculeStatus.OK.value
    assert int(state.loc[0, "reruns"]) == 0
    assert float(scientific.loc[0, "H298_pm7"]) == -55.0


@pytest.mark.anyio
async def test_second_result_same_attempt_is_stale_after_skip(tmp_path: Path) -> None:
    """Terminal Skip: second success and second failure are stale_attempt."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        mol_id, attempt_id = await _claim(client)
        sync_service = SyncResultApplicationService(
            csv_manager=app.state.csv_manager,
            state_manager=app.state.state_manager,
            ledger=app.state.result_ledger,
            worker_registry=app.state.worker_registry,
            running_molecules=app.state.running_molecules,
        )
        skip_result = SyncResult(
            result_id="skip-r1",
            mol_id=mol_id,
            success=False,
            result_update=None,
            error="manual skip",
            completed_at="2026-04-21T10:00:00Z",
            attempt_id=attempt_id,
        )
        applied = sync_service.apply_force_skip("w1", skip_result, "manual skip")
        assert applied.item.status.value == "applied"
        assert applied.final_status == MoleculeStatus.SKIP.value
        reruns_after = app.state.state_manager.get_reruns(mol_id)

        success_retry = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "skip-r2-success",
                        "mol_id": mol_id,
                        "attempt_id": attempt_id,
                        "success": True,
                        "result_update": {"H298_pm7": -1.0},
                        "error": None,
                        "completed_at": "2026-04-21T10:01:00Z",
                    }
                ],
            },
        )
        failure_retry = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "skip-r2-failure",
                        "mol_id": mol_id,
                        "attempt_id": attempt_id,
                        "success": False,
                        "result_update": None,
                        "error": "late failure",
                        "completed_at": "2026-04-21T10:02:00Z",
                    }
                ],
            },
        )

    assert success_retry.json()["items"][0]["status"] == "stale_attempt"
    assert failure_retry.json()["items"][0]["status"] == "stale_attempt"
    assert app.state.state_manager.get_status(mol_id) == MoleculeStatus.SKIP.value
    assert app.state.state_manager.get_reruns(mol_id) == reruns_after


@pytest.mark.anyio
async def test_prepared_stale_after_reclaim_and_new_claim(tmp_path: Path) -> None:
    """PREPARED for A, reclaim, claim B, retry A → stale; B stays Assigned."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        mol_id, attempt_a = await _claim(client)
        result_id = "prepared-a"
        fingerprint_payload = {
            "mol_id": mol_id,
            "attempt_id": attempt_a,
            "success": True,
            "result_update": {"H298_pm7": -55.0},
            "error": None,
            "completed_at": "2026-04-21T10:00:00Z",
        }
        fingerprint = build_result_fingerprint(fingerprint_payload)
        app.state.result_ledger.prepare(
            result_id=result_id,
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=True,
            previous_status=MoleculeStatus.ASSIGNED.value,
            previous_reruns=0,
            expected_final_status=MoleculeStatus.OK.value,
            expected_reruns=0,
            attempt_id=attempt_a,
            worker_id="w1",
        )
        reclaimed = app.state.state_manager.mark_worker_offline("w1")
        assert reclaimed == 1
        app.state.running_molecules.pop(mol_id, None)

        await client.post("/register", json={"worker_id": "w2", "hostname": "lab2"})
        claimed = await client.post("/claim", json={"worker_id": "w2"})
        attempt_b = claimed.json()["attempt_id"]
        assert attempt_b != attempt_a
        assert (
            app.state.state_manager.get_status(mol_id) == MoleculeStatus.ASSIGNED.value
        )

        stale = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [{"result_id": result_id, **fingerprint_payload}],
            },
        )

    assert stale.json()["items"][0]["status"] == "stale_attempt"
    assert app.state.state_manager.get_attempt_id(mol_id) == attempt_b
    assert app.state.state_manager.get_status(mol_id) == MoleculeStatus.ASSIGNED.value


@pytest.mark.anyio
async def test_prepared_resume_same_lease_assigned_or_running(
    tmp_path: Path,
) -> None:
    """PREPARED with matching lease on Assigned/Running resumes dual-write."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        mol_id, attempt_id = await _claim(client)
        assert (
            app.state.state_manager.get_status(mol_id) == MoleculeStatus.ASSIGNED.value
        )
        result_id = "prepared-resume"
        fingerprint_payload = {
            "mol_id": mol_id,
            "attempt_id": attempt_id,
            "success": True,
            "result_update": {"H298_pm7": -42.0},
            "error": None,
            "completed_at": "2026-04-21T10:00:00Z",
        }
        fingerprint = build_result_fingerprint(fingerprint_payload)
        app.state.result_ledger.prepare(
            result_id=result_id,
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=True,
            previous_status=MoleculeStatus.ASSIGNED.value,
            previous_reruns=0,
            expected_final_status=MoleculeStatus.OK.value,
            expected_reruns=0,
            attempt_id=attempt_id,
            worker_id="w1",
        )

        response = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [{"result_id": result_id, **fingerprint_payload}],
            },
        )

    assert response.status_code == 200
    assert response.json()["items"][0]["status"] == "applied"
    assert app.state.state_manager.get_status(mol_id) == MoleculeStatus.OK.value


def test_active_assignment_status_uses_assigned_not_claimed() -> None:
    assert is_active_assignment_status(MoleculeStatus.ASSIGNED.value) is True
    assert is_active_assignment_status(MoleculeStatus.RUNNING.value) is True
    assert is_active_assignment_status(MoleculeStatus.SELECTED.value) is True
    assert is_active_assignment_status("Assigned") is True
    assert is_active_assignment_status("assigned") is True
    assert is_active_assignment_status("claimed") is False
    assert is_active_assignment_status("Claimed") is False
    assert is_active_assignment_status(MoleculeStatus.PENDING.value) is False
    assert is_active_assignment_status(MoleculeStatus.OK.value) is False


@pytest.mark.anyio
async def test_force_skip_prepared_recovery_without_metric_duplication(
    tmp_path: Path,
) -> None:
    """Crash after dual-write before commit: recover commits, no re-apply/metrics."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        mol_id, attempt_id = await _claim(client)
        sync_service = SyncResultApplicationService(
            csv_manager=app.state.csv_manager,
            state_manager=app.state.state_manager,
            ledger=app.state.result_ledger,
            worker_registry=app.state.worker_registry,
            running_molecules=app.state.running_molecules,
        )
        result = SyncResult(
            result_id="force-skip-1",
            mol_id=mol_id,
            success=False,
            result_update=None,
            error="manual skip",
            completed_at="2026-04-21T10:00:00Z",
            attempt_id=attempt_id,
        )
        fingerprint = build_result_fingerprint(sync_result_payload(result))
        previous_reruns = app.state.state_manager.get_reruns(mol_id)
        app.state.result_ledger.prepare(
            result_id="force-skip-1",
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=False,
            previous_status=MoleculeStatus.ASSIGNED.value,
            previous_reruns=previous_reruns,
            expected_final_status=MoleculeStatus.SKIP.value,
            expected_reruns=previous_reruns,
            attempt_id=attempt_id,
            worker_id="w1",
        )
        apply_worker_result(
            app.state.csv_manager,
            app.state.state_manager,
            mol_id=mol_id,
            success=False,
            error="manual skip",
            force_skip=True,
        )
        assert app.state.state_manager.get_status(mol_id) == MoleculeStatus.SKIP.value
        assert app.state.state_manager.get_reruns(mol_id) == previous_reruns

        recovered = sync_service.apply_force_skip("w1", result, "manual skip")
        assert recovered.item.status.value == "duplicate"
        worker = app.state.worker_registry.get_worker("w1")
        assert worker is not None
        assert worker.skipped == 0
        assert app.state.state_manager.get_reruns(mol_id) == previous_reruns


@pytest.mark.anyio
async def test_force_skip_keeps_reruns_and_metrics_once(tmp_path: Path) -> None:
    """force_skip → Skip with unchanged reruns; duplicate does not bump metrics."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    app = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab"})
        await client.post("/dispatch/start")
        claim = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = claim.json()["mol_id"]
        attempt_id = claim.json()["attempt_id"]
        before = app.state.state_manager.get_reruns(mol_id)
        payload = {
            "worker_id": "w1",
            "mol_id": mol_id,
            "error": "skip",
            "force_skip": True,
            "attempt_id": attempt_id,
            "result_id": "force-skip-c",
            "completed_at": "2026-04-21T10:00:00Z",
        }
        resp = await client.post("/report/failure", json=payload)
        assert resp.status_code == 200
        assert app.state.state_manager.get_status(mol_id) == MoleculeStatus.SKIP.value
        assert app.state.state_manager.get_reruns(mol_id) == before
        status = await client.get("/workers/status")
        worker = next(w for w in status.json() if w["worker_id"] == "w1")
        assert worker["skipped"] == 1

        resp2 = await client.post("/report/failure", json=payload)
        assert resp2.status_code == 200
        status2 = await client.get("/workers/status")
        worker2 = next(w for w in status2.json() if w["worker_id"] == "w1")
        assert worker2["skipped"] == 1
        assert app.state.state_manager.get_reruns(mol_id) == before
