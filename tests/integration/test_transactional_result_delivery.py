"""Testes de entrega transacional única (online == offline via /sync_results)."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import MagicMock

import pandas as pd
import pytest
from httpx import ASGITransport, AsyncClient
from rich.console import Console

from grimperium.cli.session import SessionContext
from grimperium.cli.views.batch_view import BatchView
from grimperium.crest_pm7.batch.enums import BatchFailurePolicy, MoleculeStatus
from grimperium.crest_pm7.batch.result_ledger import (
    JournalTxnStatus,
    OperationKind,
    ResultLedger,
    build_result_fingerprint,
    with_operation_kind,
)
from grimperium.runs.models import RunStatus
from grimperium.runs.persistence import MANIFEST_FILENAME
from grimperium.runs.service import RunService
from grimperium.server.app import create_app
from grimperium.server.config import ServerConfig
from grimperium.server.models import SyncResult
from grimperium.server.sync_application import (
    SyncConflictError,
    SyncResultApplicationService,
    apply_worker_result,
    sync_result_payload,
)
from grimperium.worker.client import LeaseLostError, ServerError, WorkerClient
from grimperium.worker.offline_queue import (
    DeadLetterQueue,
    DeadLetterRecord,
    OfflineResultQueue,
    dead_letter_path_for,
)
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


# ---------------------------------------------------------------------------
# Cenários finais A–F (lease / recovery / dead-letter / Runs)
# ---------------------------------------------------------------------------


@pytest.mark.anyio
async def test_scenario_a_all_or_nothing_clears_attempt_rejects_late_result(
    tmp_path: Path,
) -> None:
    """AoN reset limpa attempt_id; resultado atrasado → stale_attempt; Pending intacto."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    science_before = csv_path.read_text(encoding="utf-8")
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
        sm = app.state.state_manager
        df = sm._ensure_loaded()
        row_idx = int(df.index[df["mol_id"] == mol_id][0])
        df.at[row_idx, "batch_id"] = "batch-aon"
        df.at[row_idx, "batch_failure_policy"] = BatchFailurePolicy.ALL_OR_NOTHING.value
        sm._save_csv()

        reset_count = sm.reset_all_or_nothing("batch-aon")
        assert reset_count == 1
        assert sm.get_attempt_id(mol_id) is None
        assert sm.get_status(mol_id) == MoleculeStatus.PENDING.value
        app.state.running_molecules.pop(mol_id, None)

        late = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "late-aon-a",
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
        assert late.status_code == 200
        assert late.json()["items"][0]["status"] == "stale_attempt"

    state = pd.read_csv(tmp_path / "batch_state.csv")
    row = state.loc[state["mol_id"] == mol_id].iloc[0]
    assert str(row["status"]) == MoleculeStatus.PENDING.value
    assert str(row["attempt_id"]) in {"", "nan"}
    assert csv_path.read_text(encoding="utf-8") == science_before


@pytest.mark.anyio
async def test_scenario_b_second_result_id_same_attempt_rejected(
    tmp_path: Path,
) -> None:
    """Uma attempt_id produz no máximo um resultado terminal."""
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
        r1 = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "r1-terminal",
                        "mol_id": mol_id,
                        "attempt_id": attempt_id,
                        "success": True,
                        "result_update": {"H298_pm7": -55.0},
                        "error": None,
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )
        assert r1.json()["items"][0]["status"] == "applied"
        science_after_r1 = pd.read_csv(csv_path)
        reruns_after_r1 = app.state.state_manager.get_reruns(mol_id)

        dup = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "r1-terminal",
                        "mol_id": mol_id,
                        "attempt_id": attempt_id,
                        "success": True,
                        "result_update": {"H298_pm7": -55.0},
                        "error": None,
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )
        assert dup.json()["items"][0]["status"] == "duplicate"

        r2 = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "r2-same-attempt",
                        "mol_id": mol_id,
                        "attempt_id": attempt_id,
                        "success": False,
                        "result_update": None,
                        "error": "late failure",
                        "completed_at": "2026-04-21T10:01:00Z",
                    }
                ],
            },
        )
        assert r2.json()["items"][0]["status"] == "stale_attempt"

    assert app.state.state_manager.get_status(mol_id) == MoleculeStatus.OK.value
    assert app.state.state_manager.get_reruns(mol_id) == reruns_after_r1
    science_final = pd.read_csv(csv_path)
    assert float(science_final.loc[0, "H298_pm7"]) == float(
        science_after_r1.loc[0, "H298_pm7"]
    )
    worker = app.state.worker_registry.get_worker("w1")
    assert worker is not None
    assert worker.successful == 1


@pytest.mark.anyio
async def test_scenario_c_prepared_old_lease_vs_new_claim(tmp_path: Path) -> None:
    """PREPARED da tentativa A não retoma após reclaim + claim B."""
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
        mol_id, attempt_a = await _claim(client, worker_id="w1")
        fingerprint = build_result_fingerprint(
            {
                "mol_id": mol_id,
                "success": True,
                "result_update": {"H298_pm7": -11.0},
                "error": None,
                "completed_at": "2026-04-21T10:00:00Z",
            }
        )
        app.state.result_ledger.prepare(
            result_id="prep-a",
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=True,
            previous_status=MoleculeStatus.ASSIGNED.value,
            previous_reruns=0,
            expected_final_status=MoleculeStatus.OK.value,
            expected_reruns=0,
            expected_science_hash=build_result_fingerprint(
                {"H298_pm7": -11.0, "G298_pm7": None, "gap": None}
            ),
            worker_id="w1",
            attempt_id=attempt_a,
        )
        app.state.state_manager.mark_worker_offline("w1")
        app.state.running_molecules.pop(mol_id, None)

        await client.post("/register", json={"worker_id": "w2", "hostname": "lab2"})
        claimed = await client.post("/claim", json={"worker_id": "w2"})
        body = claimed.json()
        assert body["mol_id"] == mol_id
        attempt_b = body["attempt_id"]
        assert attempt_b != attempt_a
        assert app.state.state_manager.get_attempt_id(mol_id) == attempt_b

        retry_a = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "prep-a",
                        "mol_id": mol_id,
                        "attempt_id": attempt_a,
                        "success": True,
                        "result_update": {"H298_pm7": -11.0},
                        "error": None,
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )
        assert retry_a.status_code == 200
        assert retry_a.json()["items"][0]["status"] == "stale_attempt"
        assert app.state.state_manager.get_attempt_id(mol_id) == attempt_b
        assert (
            app.state.state_manager.get_status(mol_id) == MoleculeStatus.ASSIGNED.value
        )


@pytest.mark.anyio
async def test_scenario_d_force_skip_crash_recovery(tmp_path: Path) -> None:
    """force_skip: dual-write + crash antes do commit → retry sem reaplicar."""
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
            result_id="force-skip-crash",
            mol_id=mol_id,
            success=False,
            result_update=None,
            error="manual skip",
            completed_at="2026-04-21T10:00:00Z",
            attempt_id=attempt_id,
        )
        fingerprint = build_result_fingerprint(
            with_operation_kind(sync_result_payload(result), OperationKind.FORCE_SKIP)
        )
        previous_reruns = app.state.state_manager.get_reruns(mol_id)
        app.state.result_ledger.prepare(
            result_id="force-skip-crash",
            mol_id=mol_id,
            fingerprint=fingerprint,
            desired_success=False,
            previous_status=MoleculeStatus.ASSIGNED.value,
            previous_reruns=previous_reruns,
            expected_final_status=MoleculeStatus.SKIP.value,
            expected_reruns=previous_reruns,
            attempt_id=attempt_id,
            worker_id="w1",
            operation_kind=OperationKind.FORCE_SKIP,
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
        assert app.state.state_manager.get_reruns(mol_id) == previous_reruns
        entry = app.state.result_ledger._journal["force-skip-crash"]
        assert entry.txn_status is JournalTxnStatus.COMMITTED
        worker = app.state.worker_registry.get_worker("w1")
        assert worker is not None
        assert worker.skipped == 0


def test_scenario_e_dead_letter_flush_mixed_batch(tmp_path: Path) -> None:
    """applied/duplicate confirmados; conflict/stale arquivados; fila não bloqueada."""
    queue_path = tmp_path / "offline.jsonl"
    queue = OfflineResultQueue(queue_path)
    for result_id, status_hint in (
        ("ok-1", "applied"),
        ("dup-1", "duplicate"),
        ("conflict-1", "conflict"),
        ("stale-1", "stale_attempt"),
    ):
        queue.enqueue(
            mol_id="m1",
            success=True,
            result_update={"H298_pm7": -1.0},
            error=None,
            result_id=result_id,
            attempt_id=f"att-{result_id}",
            completed_at="2026-04-21T10:00:00Z",
        )
        del status_hint

    client = MagicMock(spec=WorkerClient)

    def _sync(results: list[dict[str, Any]]) -> dict[str, Any]:
        status_by_id = {
            "ok-1": "applied",
            "dup-1": "duplicate",
            "conflict-1": "conflict",
            "stale-1": "stale_attempt",
        }
        items = [
            {
                "result_id": item["result_id"],
                "mol_id": item["mol_id"],
                "status": status_by_id[str(item["result_id"])],
                "detail": "test",
            }
            for item in results
        ]
        return {"accepted": 2, "rejected": 0, "duplicate": True, "items": items}

    client.sync_results.side_effect = _sync
    runner = WorkerRunner(
        WorkerConfig(
            server_url="http://test",
            worker_id="w1",
            offline_queue_path=str(queue_path),
            heartbeat_interval_s=999.0,
        ),
        pipeline=MagicMock(),
        client=client,
    )
    runner.flush_offline_queue()

    assert len(OfflineResultQueue(queue_path).pending()) == 0
    dl = DeadLetterQueue(dead_letter_path_for(queue_path)).entries()
    statuses = {entry.returned_status for entry in dl}
    assert statuses == {"conflict", "stale_attempt"}
    assert {entry.result_id for entry in dl} == {"conflict-1", "stale-1"}
    for entry in dl:
        assert entry.original_payload["result_id"] == entry.result_id
        assert entry.worker_id == "w1"
        assert entry.rejection_origin == "sync_results"

    # Restart: dead-letter persiste; fila principal permanece vazia.
    assert len(OfflineResultQueue(queue_path).pending()) == 0
    assert len(DeadLetterQueue(dead_letter_path_for(queue_path)).entries()) == 2


def test_scenario_f_batch_view_portable_runs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """BatchView cria Run, grava sob runs/<run_id>/ e manifest recarrega outputs."""
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    service = RunService(runs_root)
    controller = SimpleNamespace(
        console=Console(record=True),
        settings_manager=None,
        session=SessionContext(),
        run_service=service,
    )
    view = BatchView(controller)  # type: ignore[arg-type]
    view.csv_path = csv_path
    view.detail_dir = tmp_path / "details"

    exec_manager, batch, manifest = view._prepare_batch()
    assert batch.batch_id
    assert exec_manager.canonical_run_id == manifest.run_id
    layout = exec_manager._output_layout
    assert layout is not None
    assert layout.output_dir == service.run_dir(manifest.run_id)
    assert not (tmp_path / "batch_outputs").exists()

    layout.output_dir.mkdir(parents=True, exist_ok=True)
    layout.calculation_results_csv.write_text(
        f"mol_id,H298_pm7,run_id\nmol_001,-10.0,{manifest.run_id}\n",
        encoding="utf-8",
    )
    service.start_run(manifest.run_id)
    view._attach_existing_outputs(manifest, layout)
    service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    raw = json.loads(
        (runs_root / manifest.run_id / MANIFEST_FILENAME).read_text(encoding="utf-8")
    )
    assert raw["output_paths"]["calculation_results_csv"] == (
        f"{manifest.run_id}/calculation_results.csv"
    )
    reloaded = RunService(runs_root).get_run(manifest.run_id)
    assert reloaded.status is RunStatus.COMPLETED
    assert reloaded.output_paths["calculation_results_csv"].exists()
    # Results consegue localizar o output via path portátil do manifest.
    assert "mol_001" in reloaded.output_paths["calculation_results_csv"].read_text(
        encoding="utf-8"
    )


# ── Cenários da estabilização terminal (prompt §14) ───────────────────────────


@pytest.mark.anyio
async def test_cenario_1_legacy_on_terminal_ok(tmp_path: Path) -> None:
    """R1 OK → R2 legado sem attempt_id → stale; CSV/reruns intactos."""
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    async with await _client(csv_path) as client:
        mol_id, attempt_id = await _claim(client)
        r1 = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "c1-r1",
                        "mol_id": mol_id,
                        "attempt_id": attempt_id,
                        "success": True,
                        "result_update": {"H298_pm7": -12.0},
                        "error": None,
                        "completed_at": "2026-07-13T00:00:00Z",
                    }
                ],
            },
        )
        assert r1.json()["items"][0]["status"] == "applied"
        r2 = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "result_id": "c1-r2",
                        "mol_id": mol_id,
                        "success": True,
                        "result_update": {"H298_pm7": -99.0},
                        "error": None,
                        "completed_at": "2026-07-13T00:01:00Z",
                    }
                ],
            },
        )
    assert r2.json()["items"][0]["status"] == "stale_attempt"
    scientific = pd.read_csv(csv_path)
    state = pd.read_csv(tmp_path / "batch_state.csv", keep_default_na=False)
    assert float(scientific.loc[0, "H298_pm7"]) == -12.0
    assert int(state.loc[0, "reruns"]) == 0


@pytest.mark.anyio
async def test_cenario_2_operation_kind_conflict(tmp_path: Path) -> None:
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
        service = SyncResultApplicationService(
            csv_manager=app.state.csv_manager,
            state_manager=app.state.state_manager,
            ledger=app.state.result_ledger,
            worker_registry=app.state.worker_registry,
            running_molecules=app.state.running_molecules,
        )
        result = SyncResult(
            result_id="c2-r1",
            mol_id=mol_id,
            success=False,
            result_update=None,
            error="timeout",
            completed_at="2026-07-13T00:00:00Z",
            attempt_id=attempt_id,
        )
        fp_normal = build_result_fingerprint(
            with_operation_kind(
                sync_result_payload(result), OperationKind.NORMAL_RESULT
            )
        )
        app.state.result_ledger.prepare(
            result_id="c2-r1",
            mol_id=mol_id,
            fingerprint=fp_normal,
            desired_success=False,
            attempt_id=attempt_id,
            worker_id="w1",
            operation_kind=OperationKind.NORMAL_RESULT,
        )
        with pytest.raises(SyncConflictError):
            service.apply_force_skip("w1", result, "admin")

        # Inverso: force_skip preparado → normal conflita
        result2 = SyncResult(
            result_id="c2-r2",
            mol_id=mol_id,
            success=False,
            result_update=None,
            error="skip",
            completed_at="2026-07-13T00:00:00Z",
            attempt_id=attempt_id,
        )
        fp_force = build_result_fingerprint(
            with_operation_kind(sync_result_payload(result2), OperationKind.FORCE_SKIP)
        )
        app.state.result_ledger.prepare(
            result_id="c2-r2",
            mol_id=mol_id,
            fingerprint=fp_force,
            desired_success=False,
            attempt_id=attempt_id,
            worker_id="w1",
            operation_kind=OperationKind.FORCE_SKIP,
        )
        with pytest.raises(SyncConflictError):
            service.apply_one("w1", result2)


def test_cenario_3_finalize_attach_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    service = RunService(runs_root)
    controller = SimpleNamespace(
        console=Console(record=True),
        settings_manager=None,
        session=SessionContext(),
        run_service=service,
    )
    view = BatchView(controller)  # type: ignore[arg-type]
    view.csv_path = csv_path
    view.detail_dir = tmp_path / "details"
    exec_manager, batch, manifest = view._prepare_batch()
    layout = exec_manager._output_layout
    assert layout is not None
    started = service.start_run(manifest.run_id)
    artifact = layout.output_dir / "calculation_results.csv"
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_text("mol_id\nmol_001\n", encoding="utf-8")
    result = SimpleNamespace(
        batch_id=batch.batch_id,
        total_count=1,
        success_count=1,
        failed_count=0,
        rerun_count=0,
        skip_count=0,
        total_time=1.0,
        success_rate=100.0,
        min_hof=None,
        max_hof=None,
        min_hof_mol_id=None,
        max_hof_mol_id=None,
        invalidated=False,
    )
    monkeypatch.setattr(view, "_display_batch_result", lambda *_a, **_k: None)
    monkeypatch.setattr(
        view,
        "_attach_existing_outputs",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("attach boom")),
    )
    with pytest.raises(RuntimeError, match="attach boom"):
        view._finalize_batch_run_safely(started, layout, result)
    assert service.get_run(manifest.run_id).status is RunStatus.FAILED
    assert artifact.exists()


def test_cenario_4_selection_compensated_on_create_run_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    csv_path = tmp_path / "thermo_pm7.csv"
    _write_csv(csv_path)
    service = RunService(runs_root)
    controller = SimpleNamespace(
        console=Console(record=True),
        settings_manager=None,
        session=SessionContext(),
        run_service=service,
    )
    view = BatchView(controller)  # type: ignore[arg-type]
    view.csv_path = csv_path
    view.detail_dir = tmp_path / "details"
    monkeypatch.setattr(
        service,
        "create_run",
        lambda **_k: (_ for _ in ()).throw(RuntimeError("create_run boom")),
    )
    with pytest.raises(RuntimeError, match="create_run boom"):
        view._prepare_batch()
    state = pd.read_csv(csv_path)
    assert str(state.loc[0, "status"]) == MoleculeStatus.PENDING.value


def test_cenario_5_dead_letter_crash_window(tmp_path: Path) -> None:
    queue_path = tmp_path / "offline.jsonl"
    client = MagicMock(spec=WorkerClient)
    client.sync_results.return_value = {
        "accepted": 0,
        "items": [
            {
                "result_id": "c5-r1",
                "mol_id": "mol_001",
                "status": "conflict",
                "detail": "fp",
            }
        ],
    }
    runner = WorkerRunner(
        WorkerConfig(
            server_url="http://test",
            worker_id="w1",
            offline_queue_path=str(queue_path),
            heartbeat_interval_s=999.0,
        ),
        pipeline=MagicMock(),
        client=client,
    )
    entry = runner._offline_queue.enqueue(
        mol_id="mol_001",
        success=False,
        error="x",
        result_id="c5-r1",
        attempt_id="att",
    )
    runner._dead_letter.append(
        DeadLetterRecord(
            result_id=entry.result_id,
            mol_id=entry.mol_id,
            attempt_id=entry.attempt_id,
            original_payload=entry.to_sync_dict(),
            returned_status="conflict",
            detail="fp",
            worker_id="w1",
            rejected_at="2026-07-13T00:00:00Z",
            rejection_origin="sync_results",
        )
    )
    restarted = WorkerRunner(
        WorkerConfig(
            server_url="http://test",
            worker_id="w1",
            offline_queue_path=str(queue_path),
            heartbeat_interval_s=999.0,
        ),
        pipeline=MagicMock(),
        client=client,
    )
    restarted.flush_offline_queue()
    assert len(restarted._dead_letter.entries()) == 1
    assert restarted._offline_queue.pending() == []


def test_cenario_6_lease_lost_does_not_publish_valid_result(tmp_path: Path) -> None:
    import threading
    import time
    from unittest.mock import patch

    client = MagicMock(spec=WorkerClient)
    client.claim.return_value = ("mol_001", "CCO", "att-1")
    client.sync_results.side_effect = lambda results: {
        "accepted": len(results),
        "items": [
            {
                "result_id": str(item.get("result_id", "")),
                "mol_id": str(item.get("mol_id", "")),
                "status": "applied",
            }
            for item in results
        ],
    }
    pipeline = MagicMock()
    release = threading.Event()

    def slow(_mol: str, _smiles: str) -> MagicMock:
        release.wait(timeout=2.0)
        out = MagicMock()
        out.mol_id = "mol_001"
        out.success = True
        out.error_message = None
        out.most_stable_hof = -1.0
        return out

    pipeline.process_molecule.side_effect = slow
    client.heartbeat.side_effect = LeaseLostError("mol_001", 409, "lost")
    runner = WorkerRunner(
        WorkerConfig(
            server_url="http://test",
            worker_id="w1",
            offline_queue_path=str(tmp_path / "offline.jsonl"),
            heartbeat_interval_s=0.05,
        ),
        pipeline=pipeline,
        client=client,
    )
    with patch(
        "grimperium.worker.runner._pm7result_to_update",
        return_value={"H298_pm7": -1.0},
    ):
        thread = threading.Thread(target=runner.run_one, daemon=True)
        thread.start()
        time.sleep(0.2)
        release.set()
        thread.join(timeout=5.0)
    client.sync_results.assert_not_called()
    assert any(
        e.rejection_origin == "lease_lost" for e in runner._dead_letter.entries()
    )
