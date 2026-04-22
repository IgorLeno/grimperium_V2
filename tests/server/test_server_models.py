"""Tests for server/models.py — Pydantic request/response schemas."""

import pytest
from pydantic import ValidationError

from grimperium.server.models import (
    AssignPayload,
    ClaimResponse,
    HeartbeatRequest,
    MoleculeRef,
    RegisterRequest,
    ReportFailureRequest,
    ReportSuccessRequest,
    StatusResponse,
    SyncResponse,
    SyncResult,
    WorkerInfo,
)


# ── RegisterRequest ───────────────────────────────────────────────────────────


class TestRegisterRequest:
    def test_valid(self) -> None:
        r = RegisterRequest(worker_id="w1", hostname="lab-01")
        assert r.worker_id == "w1"
        assert r.hostname == "lab-01"

    def test_missing_worker_id(self) -> None:
        with pytest.raises(ValidationError):
            RegisterRequest(hostname="lab-01")  # type: ignore[call-arg]

    def test_missing_hostname(self) -> None:
        with pytest.raises(ValidationError):
            RegisterRequest(worker_id="w1")  # type: ignore[call-arg]


# ── MoleculeRef ───────────────────────────────────────────────────────────────


class TestMoleculeRef:
    def test_valid(self) -> None:
        m = MoleculeRef(mol_id="mol_001", smiles="CCO")
        assert m.mol_id == "mol_001"
        assert m.smiles == "CCO"


# ── AssignPayload ─────────────────────────────────────────────────────────────


class TestAssignPayload:
    def test_valid(self) -> None:
        p = AssignPayload(
            worker_id="w1",
            server_url="http://localhost:8000",
            crest_timeout_minutes=30,
            mopac_timeout_minutes=60,
            molecules=[MoleculeRef(mol_id="mol_001", smiles="CCO")],
        )
        assert p.worker_id == "w1"
        assert len(p.molecules) == 1

    def test_empty_molecules_list_allowed(self) -> None:
        p = AssignPayload(
            worker_id="w1",
            server_url="http://localhost:8000",
            crest_timeout_minutes=30,
            mopac_timeout_minutes=60,
            molecules=[],
        )
        assert p.molecules == []


# ── HeartbeatRequest ──────────────────────────────────────────────────────────


class TestHeartbeatRequest:
    def test_valid(self) -> None:
        r = HeartbeatRequest(worker_id="w1")
        assert r.worker_id == "w1"


# ── ReportSuccessRequest ──────────────────────────────────────────────────────


class TestReportSuccessRequest:
    def test_valid(self) -> None:
        r = ReportSuccessRequest(
            worker_id="w1",
            mol_id="mol_001",
            result_update={"H298_pm7": -45.3, "mopac_status": "OK"},
        )
        assert r.mol_id == "mol_001"
        assert r.result_update["H298_pm7"] == -45.3

    def test_empty_result_update_allowed(self) -> None:
        r = ReportSuccessRequest(worker_id="w1", mol_id="mol_001", result_update={})
        assert r.result_update == {}


# ── ReportFailureRequest ──────────────────────────────────────────────────────


class TestReportFailureRequest:
    def test_valid_defaults(self) -> None:
        r = ReportFailureRequest(worker_id="w1", mol_id="mol_001", error="CREST failed")
        assert r.force_skip is False
        assert r.error == "CREST failed"

    def test_force_skip_true(self) -> None:
        r = ReportFailureRequest(
            worker_id="w1", mol_id="mol_001", error="crash", force_skip=True
        )
        assert r.force_skip is True


# ── SyncResult ────────────────────────────────────────────────────────────────


class TestSyncResult:
    def test_success_result(self) -> None:
        r = SyncResult(
            mol_id="mol_001",
            success=True,
            result_update={"H298_pm7": -45.3},
            error=None,
            completed_at="2026-04-21T11:45:00Z",
        )
        assert r.success is True
        assert r.error is None

    def test_failure_result(self) -> None:
        r = SyncResult(
            mol_id="mol_002",
            success=False,
            result_update=None,
            error="MOPAC timeout",
            completed_at="2026-04-21T11:50:00Z",
        )
        assert r.success is False
        assert r.result_update is None


# ── ClaimResponse ─────────────────────────────────────────────────────────────


class TestClaimResponse:
    def test_with_molecule(self) -> None:
        r = ClaimResponse(mol_id="mol_001", smiles="CCO")
        assert r.mol_id == "mol_001"

    def test_empty_queue(self) -> None:
        r = ClaimResponse(mol_id=None, smiles=None)
        assert r.mol_id is None
        assert r.smiles is None


# ── WorkerInfo / StatusResponse ───────────────────────────────────────────────


class TestStatusResponse:
    def test_valid(self) -> None:
        info = WorkerInfo(
            worker_id="w1",
            hostname="lab-01",
            last_seen="2026-04-21T10:00:00Z",
            mol_id_current="mol_005",
        )
        resp = StatusResponse(
            counts={"PENDING": 100, "RUNNING": 2, "OK": 50},
            workers=[info],
        )
        assert resp.counts["PENDING"] == 100
        assert resp.workers[0].worker_id == "w1"

    def test_no_workers(self) -> None:
        resp = StatusResponse(counts={}, workers=[])
        assert resp.workers == []


# ── SyncResponse ──────────────────────────────────────────────────────────────


class TestSyncResponse:
    def test_valid(self) -> None:
        r = SyncResponse(accepted=5, rejected=1)
        assert r.accepted == 5
        assert r.rejected == 1
