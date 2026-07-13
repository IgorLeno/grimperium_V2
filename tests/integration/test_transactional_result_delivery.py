"""Testes de entrega transacional única (online == offline via /sync_results)."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import pandas as pd
import pytest
from httpx import ASGITransport, AsyncClient

from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.batch.result_ledger import (
    JournalTxnStatus,
    ResultLedger,
    build_result_fingerprint,
)
from grimperium.server.app import create_app
from grimperium.server.config import ServerConfig
from grimperium.worker.client import ServerError, WorkerClient
from grimperium.worker.offline_queue import OfflineResultQueue
from grimperium.worker.runner import WorkerConfig, WorkerRunner


def _write_csv(path: Path, mols: list[tuple[str, str]] | None = None) -> None:
    rows = mols or [("mol_001", "CCO")]
    pd.DataFrame(
        [
            {
                "mol_id": mol_id,
                "smiles": smiles,
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            }
            for mol_id, smiles in rows
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


async def _claim(client: AsyncClient, worker_id: str = "w1") -> tuple[str, str]:
    await client.post("/register", json={"worker_id": worker_id, "hostname": "lab"})
    await client.post("/dispatch/start")
    response = await client.post("/claim", json={"worker_id": worker_id})
    body = response.json()
    mol_id = body["mol_id"]
    attempt_id = body["attempt_id"]
    assert mol_id is not None
    assert attempt_id is not None
    return str(mol_id), str(attempt_id)


@pytest.mark.anyio
async def test_lost_http_response_does_not_double_apply(tmp_path: Path) -> None:
    """Servidor aplica; cliente perde ACK; reenvio via sync não duplica efeito."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    queue = OfflineResultQueue(tmp_path / "offline.jsonl")

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        entry = queue.enqueue(
            mol_id=mol_id,
            success=False,
            result_update=None,
            error="MOPAC timeout",
            result_id="stable-result-1",
            attempt_id=attempt_id,
            completed_at="2026-04-21T10:00:00Z",
        )
        payload = {
            "worker_id": "w1",
            "results": [entry.to_sync_dict()],
        }
        first = await client.post("/sync_results", json=payload)
        assert first.status_code == 200
        assert first.json()["accepted"] == 1
        # Simula perda da resposta: fila não é confirmada.
        assert len(queue.pending()) == 1

        second = await client.post("/sync_results", json=payload)
        assert second.status_code == 200
        body = second.json()
        assert body["duplicate"] is True
        assert body["items"][0]["status"] == "duplicate"

        # Confirmação idempotente esvazia a fila.
        for item in body["items"]:
            if item["status"] in {"applied", "duplicate"}:
                queue.confirm(item["result_id"])

    state = pd.read_csv(tmp_path / "batch_state.csv")
    assert int(state.loc[0, "reruns"]) == 1
    assert len(queue.pending()) == 0
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    decision = ledger.check(
        "stable-result-1",
        build_result_fingerprint(
            {
                "mol_id": mol_id,
                "success": False,
                "result_update": None,
                "error": "MOPAC timeout",
                "completed_at": "2026-04-21T10:00:00Z",
                "attempt_id": attempt_id,
            }
        ),
    )
    assert decision.duplicate is True


@pytest.mark.anyio
async def test_sync_items_distinguish_partial_batch(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path, [("mol_001", "CCO"), ("mol_002", "CCC")])

    async with await _client(csv_path) as client:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab"})
        await client.post("/dispatch/start")
        first = await client.post("/claim", json={"worker_id": "w1"})
        mol_a = first.json()["mol_id"]
        attempt_a = first.json()["attempt_id"]
        second = await client.post("/claim", json={"worker_id": "w1"})
        mol_b = second.json()["mol_id"]
        attempt_b = second.json()["attempt_id"]

        # Pré-aplica um result_id para forçar duplicate no lote.
        await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "already-done",
                        "mol_id": mol_a,
                        "attempt_id": attempt_a,
                        "success": False,
                        "result_update": None,
                        "error": "timeout",
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )
        response = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "already-done",
                        "mol_id": mol_a,
                        "attempt_id": attempt_a,
                        "success": False,
                        "result_update": None,
                        "error": "timeout",
                        "completed_at": "2026-04-21T10:00:00Z",
                    },
                    {
                        "result_id": "new-one",
                        "mol_id": mol_b,
                        "attempt_id": attempt_b,
                        "success": True,
                        "result_update": {"H298_pm7": -55.0},
                        "error": None,
                        "completed_at": "2026-04-21T10:01:00Z",
                    },
                ],
            },
        )

    assert response.status_code == 200
    body = response.json()
    statuses = {item["result_id"]: item["status"] for item in body["items"]}
    assert statuses["already-done"] == "duplicate"
    assert statuses["new-one"] == "applied"
    assert body["accepted"] == 1


@pytest.mark.anyio
async def test_conflict_keeps_operational_state_unchanged(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        base = {
            "result_id": "same-id",
            "mol_id": mol_id,
            "attempt_id": attempt_id,
            "success": False,
            "result_update": None,
            "error": "timeout",
            "completed_at": "2026-04-21T10:00:00Z",
        }
        first = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [base]}
        )
        assert first.status_code == 200
        state_before = pd.read_csv(tmp_path / "batch_state.csv")
        conflict = dict(base)
        conflict["success"] = True
        conflict["result_update"] = {"H298_pm7": -1.0}
        conflict["error"] = None
        second = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [conflict]}
        )
        assert second.status_code == 200
        assert second.json()["items"][0]["status"] == "conflict"
        state_after = pd.read_csv(tmp_path / "batch_state.csv")

    assert int(state_before.loc[0, "reruns"]) == int(state_after.loc[0, "reruns"])
    assert str(state_before.loc[0, "status"]) == str(state_after.loc[0, "status"])


@pytest.mark.anyio
async def test_prepared_pending_not_committed_on_startup(tmp_path: Path) -> None:
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
        fingerprint = build_result_fingerprint(
            {
                "mol_id": mol_id,
                "success": True,
                "result_update": {"H298_pm7": -10.0},
                "error": None,
                "completed_at": "2026-04-21T10:00:00Z",
            }
        )
        app.state.result_ledger.prepare(
            result_id="prep-pending",
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=True,
            previous_status="Running",
            previous_reruns=0,
            expected_final_status=MoleculeStatus.OK.value,
            expected_reruns=0,
            expected_science_hash=build_result_fingerprint(
                {"H298_pm7": -10.0, "G298_pm7": None, "gap": None}
            ),
            worker_id="w1",
        )
        # Molécula volta a Pending sem dual-write.
        app.state.state_manager.update_molecule_status(
            mol_id, MoleculeStatus.PENDING.value
        )

    # Recria app (startup recovery).
    app2 = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    incomplete = app2.state.result_ledger.get_incomplete()
    assert any(entry.result_id == "prep-pending" for entry in incomplete)
    entry = next(e for e in incomplete if e.result_id == "prep-pending")
    assert entry.txn_status is JournalTxnStatus.PREPARED


@pytest.mark.anyio
async def test_prepared_success_verified_commits_on_startup(tmp_path: Path) -> None:
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
        science_hash = build_result_fingerprint(
            {"H298_pm7": -55.0, "G298_pm7": None, "gap": None}
        )
        fingerprint = build_result_fingerprint(
            {
                "mol_id": mol_id,
                "success": True,
                "result_update": {"H298_pm7": -55.0},
                "error": None,
                "completed_at": "2026-04-21T10:00:00Z",
            }
        )
        app.state.result_ledger.prepare(
            result_id="prep-ok",
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=True,
            previous_status="Running",
            previous_reruns=0,
            expected_final_status=MoleculeStatus.OK.value,
            expected_reruns=0,
            expected_science_hash=science_hash,
            worker_id="w1",
        )
        app.state.state_manager.mark_success(mol_id)
        app.state.csv_manager.mark_success(mol_id, {"H298_pm7": -55.0})

    app2 = create_app(
        ServerConfig(
            csv_path=str(csv_path),
            api_token="",
            startup_grace_s=0,
            watchdog_interval_s=999,
        )
    )
    assert not any(
        e.result_id == "prep-ok" for e in app2.state.result_ledger.get_incomplete()
    )
    committed = app2.state.result_ledger._journal["prep-ok"]
    assert committed.txn_status is JournalTxnStatus.COMMITTED


def test_worker_runner_persists_before_sync_and_confirms_on_ack(
    tmp_path: Path,
) -> None:
    queue_path = tmp_path / "offline.jsonl"
    client = MagicMock(spec=WorkerClient)
    client.claim.return_value = ("m1", "CCO", "att-1")
    client.sync_results.return_value = {
        "accepted": 1,
        "rejected": 0,
        "items": [
            {"result_id": "will-be-replaced", "mol_id": "m1", "status": "applied"}
        ],
    }

    def _sync(results: list[dict[str, Any]]) -> dict[str, Any]:
        result_id = results[0]["result_id"]
        return {
            "accepted": 1,
            "rejected": 0,
            "items": [{"result_id": result_id, "mol_id": "m1", "status": "applied"}],
        }

    client.sync_results.side_effect = _sync

    pipeline = MagicMock()
    result = MagicMock()
    result.mol_id = "m1"
    result.success = True
    result.error_message = None
    result.most_stable_hof = -42.0
    pipeline.process_molecule.return_value = result

    runner = WorkerRunner(
        WorkerConfig(
            server_url="http://test",
            worker_id="w1",
            offline_queue_path=str(queue_path),
            heartbeat_interval_s=999.0,
        ),
        pipeline=pipeline,
        client=client,
    )
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(
            "grimperium.worker.runner._pm7result_to_update",
            lambda *a, **k: {"H298_pm7": -42.0},
        )
        runner.run_one()

    assert client.sync_results.called
    assert client.report_success.call_count == 0
    assert len(OfflineResultQueue(queue_path).pending()) == 0


def test_worker_keeps_queue_on_temporary_failure(tmp_path: Path) -> None:
    queue_path = tmp_path / "offline.jsonl"
    client = MagicMock(spec=WorkerClient)
    client.claim.return_value = ("m1", "CCO", "att-1")
    client.sync_results.side_effect = ServerError("POST /sync_results → 503")

    pipeline = MagicMock()
    result = MagicMock()
    result.mol_id = "m1"
    result.success = False
    result.error_message = "timeout"
    result.most_stable_hof = None
    pipeline.process_molecule.return_value = result

    runner = WorkerRunner(
        WorkerConfig(
            server_url="http://test",
            worker_id="w1",
            offline_queue_path=str(queue_path),
            heartbeat_interval_s=999.0,
        ),
        pipeline=pipeline,
        client=client,
    )
    runner.run_one()
    pending = OfflineResultQueue(queue_path).pending()
    assert len(pending) == 1
    stable_id = pending[0].result_id

    # Retry posterior com o mesmo result_id.
    client.sync_results.side_effect = None
    client.sync_results.return_value = {
        "accepted": 1,
        "rejected": 0,
        "items": [{"result_id": stable_id, "mol_id": "m1", "status": "applied"}],
    }
    accepted, rejected = runner.flush_offline_queue()
    assert accepted == 1
    assert rejected == 0
    assert OfflineResultQueue(queue_path).pending() == []
    assert client.sync_results.call_args[0][0][0]["result_id"] == stable_id


@pytest.mark.anyio
async def test_concurrent_same_result_id_applies_once(tmp_path: Path) -> None:
    """Duas requisições concorrentes com o mesmo result_id → efeito único."""
    import asyncio

    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        payload = {
            "worker_id": "w1",
            "results": [
                {
                    "result_id": "concurrent-1",
                    "mol_id": mol_id,
                    "attempt_id": attempt_id,
                    "success": False,
                    "result_update": None,
                    "error": "timeout",
                    "completed_at": "2026-04-21T10:00:00Z",
                }
            ],
        }
        first, second = await asyncio.gather(
            client.post("/sync_results", json=payload),
            client.post("/sync_results", json=payload),
        )

    assert {first.status_code, second.status_code} <= {200}
    bodies = [first.json(), second.json()]
    applied = sum(1 for body in bodies if body.get("accepted", 0) == 1)
    duplicates = sum(1 for body in bodies if body.get("duplicate") is True)
    assert applied == 1
    assert duplicates >= 1
    state = pd.read_csv(tmp_path / "batch_state.csv")
    assert int(state.loc[0, "reruns"]) == 1


@pytest.mark.anyio
async def test_concurrent_same_id_different_fingerprint_conflicts(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        base = {
            "result_id": "conflict-race",
            "mol_id": mol_id,
            "attempt_id": attempt_id,
            "success": False,
            "result_update": None,
            "error": "timeout",
            "completed_at": "2026-04-21T10:00:00Z",
        }
        first = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [base]}
        )
        assert first.status_code == 200
        other = dict(base)
        other["success"] = True
        other["result_update"] = {"H298_pm7": -1.0}
        other["error"] = None
        second = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": [other]}
        )
        assert second.status_code == 200
        assert second.json()["items"][0]["status"] == "conflict"

    state = pd.read_csv(tmp_path / "batch_state.csv")
    assert int(state.loc[0, "reruns"]) == 1
    assert str(state.loc[0, "status"]) == MoleculeStatus.RERUN.value


@pytest.mark.anyio
async def test_stale_attempt_after_reclaim_does_not_affect_new_claim(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)

    async with await _client(csv_path) as client:
        mol_id, attempt_a = await _claim(client, worker_id="w1")
        # Watchdog/reclaim: molecule returns to pending without completing.
        state_path = tmp_path / "batch_state.csv"
        state = pd.read_csv(state_path)
        state.loc[0, "status"] = MoleculeStatus.PENDING.value
        state.loc[0, "assigned_worker"] = ""
        state.loc[0, "assigned_at"] = ""
        state.loc[0, "attempt_id"] = ""
        state.loc[0, "worker_status"] = "unassigned"
        state.to_csv(state_path, index=False)

        # Force server to reload state from disk on next claim.
        # Claim as worker-2 with a fresh attempt lease.
        await client.post("/register", json={"worker_id": "w2", "hostname": "lab2"})
        # Re-read via claim requires state_manager reload — recreate client/app.
    async with await _client(csv_path) as client:
        await client.post("/register", json={"worker_id": "w2", "hostname": "lab2"})
        await client.post("/dispatch/start")
        claimed = await client.post("/claim", json={"worker_id": "w2"})
        body = claimed.json()
        assert body["mol_id"] == mol_id
        attempt_b = body["attempt_id"]
        assert attempt_b != attempt_a

        stale = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "late-from-a",
                        "mol_id": mol_id,
                        "attempt_id": attempt_a,
                        "success": True,
                        "result_update": {"H298_pm7": -99.0},
                        "error": None,
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )
        assert stale.status_code == 200
        assert stale.json()["items"][0]["status"] == "stale_attempt"

        state_mid = pd.read_csv(tmp_path / "batch_state.csv")
        assert str(state_mid.loc[0, "attempt_id"]) == attempt_b
        assert str(state_mid.loc[0, "status"]) != MoleculeStatus.OK.value

        ok = await client.post(
            "/sync_results",
            json={
                "worker_id": "w2",
                "results": [
                    {
                        "result_id": "from-b",
                        "mol_id": mol_id,
                        "attempt_id": attempt_b,
                        "success": True,
                        "result_update": {"H298_pm7": -55.0},
                        "error": None,
                        "completed_at": "2026-04-21T10:01:00Z",
                    }
                ],
            },
        )
        assert ok.status_code == 200
        assert ok.json()["items"][0]["status"] == "applied"

    state_final = pd.read_csv(tmp_path / "batch_state.csv")
    assert str(state_final.loc[0, "status"]) == MoleculeStatus.OK.value


@pytest.mark.anyio
async def test_mixed_batch_with_conflict_returns_200_per_item(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path, [("mol_001", "CCO"), ("mol_002", "CCC")])
    queue = OfflineResultQueue(tmp_path / "offline.jsonl")

    async with await _client(csv_path) as client:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab"})
        await client.post("/dispatch/start")
        c1 = await client.post("/claim", json={"worker_id": "w1"})
        mol_a = c1.json()["mol_id"]
        att_a = c1.json()["attempt_id"]
        c2 = await client.post("/claim", json={"worker_id": "w1"})
        mol_b = c2.json()["mol_id"]
        att_b = c2.json()["attempt_id"]

        await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "dup-seed",
                        "mol_id": mol_a,
                        "attempt_id": att_a,
                        "success": False,
                        "result_update": None,
                        "error": "seed",
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )

        results = [
            {
                "result_id": "dup-seed",
                "mol_id": mol_a,
                "attempt_id": att_a,
                "success": False,
                "result_update": None,
                "error": "seed",
                "completed_at": "2026-04-21T10:00:00Z",
            },
            {
                "result_id": "apply-b",
                "mol_id": mol_b,
                "attempt_id": att_b,
                "success": True,
                "result_update": {"H298_pm7": -1.0},
                "error": None,
                "completed_at": "2026-04-21T10:02:00Z",
            },
            {
                "result_id": "conflict-payload",
                "mol_id": mol_a,
                "attempt_id": att_a,
                "success": True,
                "result_update": {"H298_pm7": -2.0},
                "error": None,
                "completed_at": "2026-04-21T10:00:00Z",
            },
        ]
        for item in results:
            queue.enqueue(
                mol_id=str(item["mol_id"]),
                success=bool(item["success"]),
                result_update=(
                    dict(item["result_update"])
                    if isinstance(item["result_update"], dict)
                    else None
                ),
                error=(str(item["error"]) if item["error"] is not None else None),
                result_id=str(item["result_id"]),
                attempt_id=str(item["attempt_id"]),
                completed_at=str(item["completed_at"]),
            )
        pending_payloads = [e.to_sync_dict() for e in queue.pending()]
        for payload in pending_payloads:
            if payload["result_id"] == "conflict-payload":
                payload["result_id"] = "dup-seed"

        response = await client.post(
            "/sync_results",
            json={"worker_id": "w1", "results": pending_payloads},
        )
        assert response.status_code == 200
        body = response.json()
        flat = [item["status"] for item in body["items"]]
        assert "duplicate" in flat
        assert "applied" in flat
        assert "conflict" in flat

        for item in body["items"]:
            if item["status"] in {
                "applied",
                "duplicate",
                "conflict",
                "stale_attempt",
            }:
                queue.confirm(item["result_id"])
                queue.confirm("conflict-payload")
        assert len(queue.pending()) == 0
