"""Tests for server/app.py — FastAPI endpoints via httpx.AsyncClient."""

from pathlib import Path
from typing import Any, AsyncGenerator

import pytest
from httpx import ASGITransport, AsyncClient

from grimperium.server.app import create_app
from grimperium.server.config import ServerConfig

# ── Minimal CSV fixture ───────────────────────────────────────────────────────

_MINIMAL_CSV = """\
mol_id,smiles,nheavy,status,reruns
mol_001,CCO,3,Pending,0
mol_002,CCCO,4,Pending,0
mol_003,CCC,3,Rerun,1
mol_004,CCCC,4,OK,0
"""


@pytest.fixture
def csv_path(tmp_path: Path) -> Path:
    p = tmp_path / "test.csv"
    p.write_text(_MINIMAL_CSV)
    return p


@pytest.fixture
def server_config(csv_path: Path) -> ServerConfig:
    return ServerConfig(
        csv_path=str(csv_path),
        api_token="",  # auth disabled
        startup_grace_s=0,
        watchdog_interval_s=999,  # prevent watchdog from firing in tests
    )


@pytest.fixture
def auth_config(csv_path: Path) -> ServerConfig:
    return ServerConfig(
        csv_path=str(csv_path),
        api_token="test-token",
        startup_grace_s=0,
        watchdog_interval_s=999,
    )


@pytest.fixture
async def client(server_config: ServerConfig) -> AsyncGenerator[AsyncClient, None]:
    app = create_app(server_config)
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as c:
        yield c


@pytest.fixture
async def auth_client(auth_config: ServerConfig) -> AsyncGenerator[AsyncClient, None]:
    app = create_app(auth_config)
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as c:
        yield c


# ── /register ─────────────────────────────────────────────────────────────────


class TestRegister:
    async def test_register_returns_200(self, client: AsyncClient) -> None:
        r = await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        assert r.status_code == 200

    async def test_register_missing_field_returns_422(self, client: AsyncClient) -> None:
        r = await client.post("/register", json={"worker_id": "w1"})
        assert r.status_code == 422

    async def test_duplicate_register_ok(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        assert r.status_code == 200


# ── /assign ───────────────────────────────────────────────────────────────────


class TestAssign:
    async def test_assign_returns_200(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post(
            "/assign",
            json={
                "worker_id": "w1",
                "server_url": "http://localhost:8000",
                "crest_timeout_minutes": 30,
                "mopac_timeout_minutes": 60,
                "molecules": [{"mol_id": "mol_001", "smiles": "CCO"}],
            },
        )
        assert r.status_code == 200

    async def test_assign_empty_molecules_ok(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post(
            "/assign",
            json={
                "worker_id": "w1",
                "server_url": "http://localhost:8000",
                "crest_timeout_minutes": 30,
                "mopac_timeout_minutes": 60,
                "molecules": [],
            },
        )
        assert r.status_code == 200


# ── /claim ────────────────────────────────────────────────────────────────────


class TestClaim:
    async def test_claim_returns_molecule(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        assert r.status_code == 200
        data = r.json()
        assert data["mol_id"] is not None
        assert data["smiles"] is not None

    async def test_claim_marks_running(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        # Claim again — same mol_id should not be returned twice
        r2 = await client.post("/claim", json={"worker_id": "w1"})
        data2 = r2.json()
        assert data2.get("mol_id") != mol_id or data2["mol_id"] is None

    async def test_claim_empty_queue_returns_null(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        # Drain all pending molecules
        for _ in range(5):
            await client.post("/claim", json={"worker_id": "w1"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        assert r.status_code == 200
        assert r.json()["mol_id"] is None


# ── /heartbeat/{mol_id} ───────────────────────────────────────────────────────


class TestHeartbeat:
    async def test_heartbeat_updates_registry(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        assert mol_id is not None
        hb = await client.put(f"/heartbeat/{mol_id}", json={"worker_id": "w1"})
        assert hb.status_code == 200

    async def test_heartbeat_unknown_mol_returns_404(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.put("/heartbeat/ghost_mol", json={"worker_id": "w1"})
        assert r.status_code == 404


# ── /report/success ───────────────────────────────────────────────────────────


class TestReportSuccess:
    async def test_success_returns_200(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        assert mol_id is not None
        rep = await client.post(
            "/report/success",
            json={"worker_id": "w1", "mol_id": mol_id, "result_update": {}},
        )
        assert rep.status_code == 200

    async def test_success_unknown_mol_returns_404(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        rep = await client.post(
            "/report/success",
            json={"worker_id": "w1", "mol_id": "ghost", "result_update": {}},
        )
        assert rep.status_code == 404


# ── /report/failure ───────────────────────────────────────────────────────────


class TestReportFailure:
    async def test_failure_rerun_returns_200(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        assert mol_id is not None
        rep = await client.post(
            "/report/failure",
            json={"worker_id": "w1", "mol_id": mol_id, "error": "CREST failed"},
        )
        assert rep.status_code == 200

    async def test_failure_force_skip_returns_200(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        assert mol_id is not None
        rep = await client.post(
            "/report/failure",
            json={
                "worker_id": "w1",
                "mol_id": mol_id,
                "error": "permanent failure",
                "force_skip": True,
            },
        )
        assert rep.status_code == 200


# ── /sync_results ─────────────────────────────────────────────────────────────


class TestSyncResults:
    async def test_sync_success_results(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        assert mol_id is not None
        resp = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "mol_id": mol_id,
                        "success": True,
                        "result_update": {},
                        "error": None,
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["accepted"] == 1
        assert data["rejected"] == 0

    async def test_sync_unknown_mol_is_rejected(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        resp = await client.post(
            "/sync_results",
            json={
                "worker_id": "w1",
                "results": [
                    {
                        "mol_id": "ghost_mol",
                        "success": True,
                        "result_update": {},
                        "error": None,
                        "completed_at": "2026-04-21T10:00:00Z",
                    }
                ],
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["rejected"] == 1

    async def test_sync_empty_results(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        resp = await client.post(
            "/sync_results", json={"worker_id": "w1", "results": []}
        )
        assert resp.status_code == 200
        assert resp.json() == {"accepted": 0, "rejected": 0}


# ── /status ───────────────────────────────────────────────────────────────────


class TestStatus:
    async def test_status_returns_counts(self, client: AsyncClient) -> None:
        r = await client.get("/status")
        assert r.status_code == 200
        data = r.json()
        assert "counts" in data
        assert "workers" in data

    async def test_status_shows_registered_worker(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.get("/status")
        workers = r.json()["workers"]
        ids = [w["worker_id"] for w in workers]
        assert "w1" in ids


# ── Authentication ────────────────────────────────────────────────────────────


class TestAuthentication:
    async def test_no_token_header_rejected(self, auth_client: AsyncClient) -> None:
        r = await auth_client.get("/status")
        assert r.status_code == 401

    async def test_wrong_token_rejected(self, auth_client: AsyncClient) -> None:
        r = await auth_client.get(
            "/status", headers={"X-Token": "wrong-token"}
        )
        assert r.status_code == 401

    async def test_correct_token_accepted(self, auth_client: AsyncClient) -> None:
        r = await auth_client.get(
            "/status", headers={"X-Token": "test-token"}
        )
        assert r.status_code == 200

    async def test_auth_disabled_no_token_needed(self, client: AsyncClient) -> None:
        r = await client.get("/status")
        assert r.status_code == 200
