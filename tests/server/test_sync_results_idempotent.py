from pathlib import Path

import pandas as pd
import pytest
from httpx import ASGITransport, AsyncClient

from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.batch.result_ledger import build_result_fingerprint
from grimperium.server.app import create_app
from grimperium.server.config import ServerConfig


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


async def _claim(client: AsyncClient) -> str:
    await client.post("/register", json={"worker_id": "w1", "hostname": "lab"})
    await client.post("/dispatch/start")
    response = await client.post("/claim", json={"worker_id": "w1"})
    mol_id = response.json()["mol_id"]
    assert mol_id is not None
    return str(mol_id)


@pytest.mark.anyio
async def test_sync_results_duplicate_result_id_is_not_applied_twice(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id = await _claim(client)
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "result_id": "result-1",
                    "mol_id": mol_id,
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
        mol_id = await _claim(client)
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "mol_id": mol_id,
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
async def test_sync_results_conflicting_result_id_returns_409(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id = await _claim(client)
        base_result = {
            "result_id": "result-1",
            "mol_id": mol_id,
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
    assert second.status_code == 409


@pytest.mark.anyio
async def test_sync_results_updates_worker_registry_once(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id = await _claim(client)
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "result_id": "result-success",
                    "mol_id": mol_id,
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
        mol_id = await _claim(client)
        result_id = "result-crash"
        fingerprint_payload = {
            "mol_id": mol_id,
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
        await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "fail-1",
                        "mol_id": mol_fail,
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
        # force skip via apply_failure semantics: success=False with special update
        await client.post(
            "/report/failure",
            json={
                "worker_id": "w1",
                "mol_id": mol_skip,
                "error": "skip",
                "force_skip": True,
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
        mol_id = await _claim(client)
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
