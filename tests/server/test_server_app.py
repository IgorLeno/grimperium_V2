"""Tests for server/app.py — FastAPI endpoints via httpx.AsyncClient."""

from collections.abc import AsyncGenerator
from pathlib import Path

import pandas as pd
import pytest
from httpx import ASGITransport, AsyncClient

from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.batch.output_contracts import BATCH_STATE_COLUMNS
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
        r = await client.post(
            "/register", json={"worker_id": "w1", "hostname": "lab-01"}
        )
        assert r.status_code == 200

    async def test_register_returns_config_payload(self, client: AsyncClient) -> None:
        r = await client.post(
            "/register", json={"worker_id": "w1", "hostname": "lab-01"}
        )
        data = r.json()
        assert data["worker_id"] == "w1"
        assert data["hostname"] == "lab-01"
        assert "crest_timeout_minutes" in data
        assert "mopac_timeout_minutes" in data
        assert "batch_size" in data
        assert "profile_name" in data

    async def test_register_missing_field_returns_422(
        self, client: AsyncClient
    ) -> None:
        r = await client.post("/register", json={"worker_id": "w1"})
        assert r.status_code == 422

    async def test_duplicate_register_ok(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post(
            "/register", json={"worker_id": "w1", "hostname": "lab-01"}
        )
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
    def test_create_app_reconciles_all_molecules_into_batch_state_csv(
        self, csv_path: Path, server_config: ServerConfig
    ) -> None:
        state_path = csv_path.parent / "batch_state.csv"
        state_path.write_text(",".join(BATCH_STATE_COLUMNS) + "\n", encoding="utf-8")

        create_app(server_config)

        rows = pd.read_csv(state_path, keep_default_na=False).set_index("mol_id")
        # Reconciliation copies ALL molecules (not just Pending/Rerun) so status
        # counts stay complete; legacy status is preserved on first creation.
        assert list(rows.index) == ["mol_001", "mol_002", "mol_003", "mol_004"]
        assert rows.loc["mol_001", "status"] == MoleculeStatus.PENDING.value
        assert rows.loc["mol_003", "status"] == MoleculeStatus.RERUN.value
        assert rows.loc["mol_004", "status"] == MoleculeStatus.OK.value

    async def test_claim_returns_molecule(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post("/dispatch/start")
        r = await client.post("/claim", json={"worker_id": "w1"})
        assert r.status_code == 200
        data = r.json()
        assert data["mol_id"] is not None
        assert data["smiles"] is not None

    async def test_claim_returns_null_when_dispatch_disabled(
        self, client: AsyncClient
    ) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        assert r.status_code == 200
        assert r.json()["mol_id"] is None

    async def test_claim_marks_running(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post("/dispatch/start")
        r = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        # Claim again — same mol_id should not be returned twice
        r2 = await client.post("/claim", json={"worker_id": "w1"})
        data2 = r2.json()
        assert data2.get("mol_id") != mol_id or data2["mol_id"] is None

    async def test_claim_empty_queue_returns_null(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post("/dispatch/start")
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
        await client.post("/dispatch/start")
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
        await client.post("/dispatch/start")
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
        await client.post("/dispatch/start")
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
        await client.post("/dispatch/start")
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
        await client.post("/dispatch/start")
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
        assert resp.json() == {
            "accepted": 0,
            "rejected": 0,
            "duplicate": False,
            "items": [],
        }


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
        r = await auth_client.get("/status", headers={"X-Token": "wrong-token"})
        assert r.status_code == 401

    async def test_correct_token_accepted(self, auth_client: AsyncClient) -> None:
        r = await auth_client.get("/status", headers={"X-Token": "test-token"})
        assert r.status_code == 200

    async def test_auth_disabled_no_token_needed(self, client: AsyncClient) -> None:
        r = await client.get("/status")
        assert r.status_code == 200


# ── /workers ──────────────────────────────────────────────────────────────────


class TestWorkers:
    async def test_workers_empty_initially(self, client: AsyncClient) -> None:
        r = await client.get("/workers")
        assert r.status_code == 200
        assert r.json() == []

    async def test_workers_shows_registered(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.get("/workers")
        assert r.status_code == 200
        data = r.json()
        assert len(data) == 1
        assert data[0]["worker_id"] == "w1"
        assert data[0]["hostname"] == "lab-01"

    async def test_workers_status_empty_initially(self, client: AsyncClient) -> None:
        r = await client.get("/workers/status")
        assert r.status_code == 200
        assert r.json() == []

    async def test_workers_status_includes_metrics(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.get("/workers/status")
        assert r.status_code == 200
        data = r.json()
        assert len(data) == 1
        w = data[0]
        assert w["worker_id"] == "w1"
        assert w["processed"] == 0
        assert w["successful"] == 0
        assert w["failed"] == 0
        assert w["skipped"] == 0
        assert w["shutdown_requested"] is False
        assert "registered_at" in w

    async def test_workers_status_updates_after_success(
        self, client: AsyncClient
    ) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post("/dispatch/start")
        r = await client.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        assert mol_id is not None
        await client.post(
            "/report/success",
            json={"worker_id": "w1", "mol_id": mol_id, "result_update": {}},
        )
        r2 = await client.get("/workers/status")
        w = r2.json()[0]
        assert w["processed"] == 1
        assert w["successful"] == 1


# ── /configure/{worker_id} ────────────────────────────────────────────────────


class TestConfigure:
    async def test_configure_known_worker(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post(
            "/configure/w1",
            json={
                "batch_size": 5,
                "crest_timeout_minutes": 45,
                "mopac_timeout_minutes": 20,
                "profile_name": "fast",
            },
        )
        assert r.status_code == 200
        assert r.json()["status"] == "configured"

    async def test_configure_unknown_worker_returns_404(
        self, client: AsyncClient
    ) -> None:
        r = await client.post(
            "/configure/ghost",
            json={
                "batch_size": 5,
                "crest_timeout_minutes": 45,
                "mopac_timeout_minutes": 20,
                "profile_name": "fast",
            },
        )
        assert r.status_code == 404

    async def test_get_config_returns_defaults_when_not_set(
        self, client: AsyncClient
    ) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.get("/configure/w1")
        assert r.status_code == 200
        data = r.json()
        assert "batch_size" in data
        assert "crest_timeout_minutes" in data

    async def test_get_config_returns_override_after_configure(
        self, client: AsyncClient
    ) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post(
            "/configure/w1",
            json={
                "batch_size": 3,
                "crest_timeout_minutes": 25,
                "mopac_timeout_minutes": 15,
                "profile_name": "custom",
            },
        )
        r = await client.get("/configure/w1")
        assert r.status_code == 200
        data = r.json()
        assert data["batch_size"] == 3
        assert data["profile_name"] == "custom"

    async def test_configure_reflected_in_register_response(
        self, client: AsyncClient
    ) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post(
            "/configure/w1",
            json={
                "batch_size": 7,
                "crest_timeout_minutes": 90,
                "mopac_timeout_minutes": 45,
                "profile_name": "heavy",
            },
        )
        r = await client.post(
            "/register", json={"worker_id": "w1", "hostname": "lab-01"}
        )
        data = r.json()
        assert data["batch_size"] == 7
        assert data["profile_name"] == "heavy"


# ── /shutdown ─────────────────────────────────────────────────────────────────


class TestShutdown:
    async def test_shutdown_known_worker(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/shutdown/w1")
        assert r.status_code == 200
        data = r.json()
        assert data["status"] == "shutdown_requested"
        assert data["worker_id"] == "w1"

    async def test_shutdown_unknown_worker_returns_404(
        self, client: AsyncClient
    ) -> None:
        r = await client.post("/shutdown/ghost")
        assert r.status_code == 404

    async def test_shutdown_all(self, client: AsyncClient) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post("/register", json={"worker_id": "w2", "hostname": "lab-02"})
        r = await client.post("/shutdown/all")
        assert r.status_code == 200
        data = r.json()
        assert data["status"] == "shutdown_requested"
        assert set(data["workers_signalled"]) == {"w1", "w2"}

    async def test_shutdown_reflected_in_workers_status(
        self, client: AsyncClient
    ) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post("/shutdown/w1")
        r = await client.get("/workers/status")
        w = r.json()[0]
        assert w["shutdown_requested"] is True


# ── /dispatch/start ───────────────────────────────────────────────────────────


class TestDispatch:
    async def test_dispatch_start_returns_200(self, client: AsyncClient) -> None:
        r = await client.post("/dispatch/start")
        assert r.status_code == 200
        assert r.json()["status"] == "dispatch_enabled"

    async def test_claim_returns_null_before_dispatch(
        self, client: AsyncClient
    ) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        r = await client.post("/claim", json={"worker_id": "w1"})
        assert r.json()["mol_id"] is None

    async def test_claim_returns_molecule_after_dispatch(
        self, client: AsyncClient
    ) -> None:
        await client.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await client.post("/dispatch/start")
        r = await client.post("/claim", json={"worker_id": "w1"})
        assert r.json()["mol_id"] is not None


# ── Operational-state authority + dual-write consistency ──────────────────────


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, keep_default_na=False).set_index("mol_id")


def _make_config(csv_path: Path, *, max_reruns: int = 3) -> ServerConfig:
    return ServerConfig(
        csv_path=str(csv_path),
        api_token="",
        startup_grace_s=0,
        watchdog_interval_s=999,
        max_reruns=max_reruns,
    )


async def _client_for(app: object) -> AsyncClient:
    return AsyncClient(transport=ASGITransport(app=app), base_url="http://test")


class TestOperationalStateAuthority:
    async def _claim_one(self, c: AsyncClient) -> str:
        await c.post("/register", json={"worker_id": "w1", "hostname": "lab-01"})
        await c.post("/dispatch/start")
        r = await c.post("/claim", json={"worker_id": "w1"})
        mol_id = r.json()["mol_id"]
        assert mol_id is not None
        return str(mol_id)

    async def test_status_counts_come_from_operational_state(
        self, csv_path: Path
    ) -> None:
        app = create_app(_make_config(csv_path))
        async with await _client_for(app) as c:
            mol_id = await self._claim_one(c)
            counts = (await c.get("/status")).json()["counts"]

        # Claim moved the molecule to Assigned in batch_state.csv; /status must
        # reflect that even though the scientific CSV still shows it Pending.
        assert counts.get("Assigned", 0) >= 1
        sci = _read_csv(csv_path)
        assert sci.loc[mol_id, "status"] == MoleculeStatus.PENDING.value
        assert "Assigned" not in set(sci["status"])

    async def test_failure_mirrors_single_rerun_in_both_files(
        self, csv_path: Path
    ) -> None:
        app = create_app(_make_config(csv_path))
        state_path = csv_path.parent / "batch_state.csv"
        async with await _client_for(app) as c:
            mol_id = await self._claim_one(c)
            rep = await c.post(
                "/report/failure",
                json={"worker_id": "w1", "mol_id": mol_id, "error": "boom"},
            )
            assert rep.status_code == 200

        state = _read_csv(state_path)
        sci = _read_csv(csv_path)
        assert state.loc[mol_id, "status"] == MoleculeStatus.RERUN.value
        assert sci.loc[mol_id, "status"] == MoleculeStatus.RERUN.value
        # Single decision: reruns incremented exactly once, identical both sides.
        assert int(state.loc[mol_id, "reruns"]) == 1
        assert int(sci.loc[mol_id, "reruns"]) == 1

    async def test_failure_at_limit_mirrors_skip_in_both_files(
        self, csv_path: Path
    ) -> None:
        app = create_app(_make_config(csv_path, max_reruns=1))
        state_path = csv_path.parent / "batch_state.csv"
        async with await _client_for(app) as c:
            mol_id = await self._claim_one(c)
            await c.post(
                "/report/failure",
                json={"worker_id": "w1", "mol_id": mol_id, "error": "boom"},
            )

        state = _read_csv(state_path)
        sci = _read_csv(csv_path)
        assert state.loc[mol_id, "status"] == MoleculeStatus.SKIP.value
        assert sci.loc[mol_id, "status"] == MoleculeStatus.SKIP.value

    async def test_success_updates_both_files_and_clears_running(
        self, csv_path: Path
    ) -> None:
        app = create_app(_make_config(csv_path))
        state_path = csv_path.parent / "batch_state.csv"
        async with await _client_for(app) as c:
            mol_id = await self._claim_one(c)
            await c.post(
                "/report/success",
                json={"worker_id": "w1", "mol_id": mol_id, "result_update": {}},
            )
            # Accepted molecule is no longer running.
            assert mol_id not in app.state.running_molecules  # type: ignore[attr-defined]

        assert _read_csv(state_path).loc[mol_id, "status"] == MoleculeStatus.OK.value
        assert _read_csv(csv_path).loc[mol_id, "status"] == MoleculeStatus.OK.value

    async def test_sync_second_persistence_failure_rolls_back_and_rejects(
        self, csv_path: Path
    ) -> None:
        app = create_app(_make_config(csv_path))
        state_path = csv_path.parent / "batch_state.csv"

        # Force the operational (second) persistence to fail, exercising the
        # compensating rollback of the scientific-CSV write.
        def _boom(mol_id: str) -> None:
            raise RuntimeError("state write failed")

        app.state.state_manager.mark_success = _boom  # type: ignore[attr-defined]

        async with await _client_for(app) as c:
            mol_id = await self._claim_one(c)
            resp = await c.post(
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

        data = resp.json()
        assert data["accepted"] == 0
        assert data["rejected"] == 1
        # Scientific CSV write was rolled back — molecule is not falsely OK.
        assert _read_csv(csv_path).loc[mol_id, "status"] != MoleculeStatus.OK.value
        # Operational state untouched (still Assigned from the claim).
        assert (
            _read_csv(state_path).loc[mol_id, "status"] == MoleculeStatus.ASSIGNED.value
        )
