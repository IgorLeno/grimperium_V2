"""API the cells call: register, claim, heartbeat and result, by bearer token."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from grimperium.webapi.cells import create_cell_app
from grimperium.webapi.config import DashboardPaths
from semi_imperium.training.cells import (
    ADMISSION_BLOCKED,
    ADMISSION_PASSED,
    ADMISSION_REFERENCES,
    AdmissionTest,
    CellRegistry,
)
from semi_imperium.training.db import database_path
from semi_imperium.training.geometry_store import GeometryStore
from semi_imperium.training.job_state import JobState
from semi_imperium.training.molecule_queue import TrainingQueue
from semi_imperium.training.protocol import canonical_settings, canonical_signature
from tests.training.conftest import accept_references, reference_rows, row

PROGRAMS = {
    "crest": "3.0.2",
    "crest_sha256": "aaa",
    "crest_path": "/central/bin/crest",
    "xtb": "6.7.1",
    "xtb_sha256": "bbb",
    "xtb_path": "/central/bin/xtb",
    "mopac": "23.2.2",
    "mopac_sha256": "ccc",
    "mopac_path": "/central/bin/mopac",
}
SIGNATURE = canonical_signature(canonical_settings())


def _paths(tmp_path: Path, rows: list[dict[str, Any]] | None = None) -> DashboardPaths:
    if rows is None:
        rows = [row("a", nheavy=1), row("b", nheavy=2), row("c", nheavy=3)]
    pd.DataFrame(rows).to_csv(tmp_path / "cbs.csv", index=False)
    return DashboardPaths(
        semi_imperium_store=tmp_path / "store",
        exp_csv=tmp_path / "exp.csv",
        exp_vs_cbs_csv=tmp_path / "vs_cbs.csv",
        cbs_csv=tmp_path / "cbs.csv",
    )


class Env:
    def __init__(
        self, tmp_path: Path, rows: list[dict[str, Any]] | None = None
    ) -> None:
        self.paths = _paths(tmp_path, rows)
        self.db = database_path(self.paths.semi_imperium_store)
        self.queue = TrainingQueue(self.db, self.paths.cbs_csv)
        self.queue.ensure_seeded()
        self.cells = CellRegistry(self.db)
        self.probes = 0

        def central_programs() -> dict[str, str]:
            self.probes += 1
            return dict(PROGRAMS)

        self.client = TestClient(
            create_cell_app(self.paths, central_programs=central_programs)
        )

    def new_cell(self, *, admitted: bool = True) -> tuple[str, dict[str, str]]:
        cell, token = self.cells.create("Lab")
        if admitted:
            self.cells.set_admission(cell.cell_id, ADMISSION_PASSED)
        return cell.cell_id, {"Authorization": f"Bearer {token}"}

    def register(self, headers: dict[str, str], **changes: Any) -> Any:
        body: dict[str, Any] = {
            "host": "lab-pc",
            "cores": 16,
            "workers": 2,
            "crest_threads": 8,
            "program_versions": {**PROGRAMS, "crest_path": "/home/aluno/crest"},
            "protocol_signature": SIGNATURE,
        }
        body.update(changes)
        return self.client.post("/api/cell/register", json=body, headers=headers)

    def ready_cell(self) -> tuple[str, dict[str, str]]:
        cell_id, headers = self.new_cell()
        assert self.register(headers).status_code == 200
        return cell_id, headers

    def claim(self, headers: dict[str, str], count: int = 1) -> Any:
        return self.client.post(
            "/api/cell/claim", json={"count": count}, headers=headers
        )


@pytest.fixture
def env(tmp_path: Path) -> Env:
    return Env(tmp_path)


def _result(mol_id: str, **changes: Any) -> dict[str, Any]:
    body: dict[str, Any] = {
        "mol_id": mol_id,
        "smiles": f"smiles-{mol_id}",
        "signature": SIGNATURE,
        "calculation_id": "calc-1",
        "conformer_index": 0,
        "elements": ["C", "H"],
        "coordinates": [[0.0, 0.0, 0.0], [0.0, 0.0, 1.09]],
        "hof_kcal_mol": -12.5,
        "lowest_frequency_cm1": 1300.0,
        "conformers_evaluated": 3,
        "elapsed_s": 42.0,
    }
    body.update(changes)
    return body


# ---------------------------------------------------------------- auth


def test_every_endpoint_needs_a_valid_token(env: Env) -> None:
    for path in ("register", "claim", "heartbeat", "result"):
        for headers in ({}, {"Authorization": "Bearer nope"}, {"Authorization": "x"}):
            response = env.client.post(f"/api/cell/{path}", json={}, headers=headers)
            assert response.status_code == 401, (path, headers)
            assert response.headers["www-authenticate"] == "Bearer"


def test_a_revoked_token_is_refused(env: Env) -> None:
    cell_id, headers = env.ready_cell()
    env.cells.revoke(cell_id)
    assert env.claim(headers).status_code == 401


def test_the_cell_app_has_no_dashboard_routes(env: Env) -> None:
    assert env.client.post("/api/training/start", json={}).status_code == 404
    assert env.client.get("/api/cells").status_code == 404


# ------------------------------------------------------------ register


def test_register_with_the_central_programs_is_accepted(env: Env) -> None:
    cell_id, headers = env.new_cell()
    response = env.register(headers)
    assert response.status_code == 200
    body = response.json()
    assert body["cell_id"] == cell_id and body["admission"] == ADMISSION_PASSED
    assert body["released"] == []
    assert body["lease_s"] > 0 and 0 < body["heartbeat_s"] < 90
    stored = env.cells.get(cell_id)
    assert stored.registered_at is not None and stored.refusal is None
    assert stored.host == "lab-pc" and stored.workers == 2


def test_register_with_other_programs_or_protocol_is_refused(env: Env) -> None:
    cell_id, headers = env.new_cell()
    response = env.register(
        headers,
        program_versions={**PROGRAMS, "xtb": "6.6.1"},
        protocol_signature="other",
    )
    assert response.status_code == 409
    detail = response.json()["detail"]
    assert "xtb: 6.6.1 (central 6.7.1)" in detail
    assert "protocol_signature: other" in detail
    stored = env.cells.get(cell_id)
    assert stored.registered_at is None and stored.refusal is not None
    assert env.claim(headers).status_code == 409


def test_the_central_programs_are_probed_once(env: Env) -> None:
    for _ in range(2):
        env.ready_cell()
    assert env.probes == 1


def test_registering_again_abandons_what_the_cell_held(env: Env) -> None:
    # A cell registers after a restart; its interrupted attempts start over.
    _, headers = env.ready_cell()
    claimed = env.claim(headers, 2).json()["molecules"]
    response = env.register(headers)
    assert sorted(response.json()["released"]) == sorted(m["mol_id"] for m in claimed)
    assert env.queue.held_attempts(env.cells.list()[0].cell_id) == set()


# --------------------------------------------------------------- claim


def test_claim_hands_out_molecules_up_to_the_free_workers(env: Env) -> None:
    cell_id, headers = env.ready_cell()  # 2 workers
    first = env.claim(headers, 5).json()
    assert [m["mol_id"] for m in first["molecules"]] == ["a", "b"]
    assert set(first["molecules"][0]) == {
        "attempt_id",
        "mol_id",
        "smiles",
        "charge",
        "multiplicity",
        "nheavy",
    }
    assert first["lease_s"] > 0
    assert env.claim(headers, 1).json()["molecules"] == []
    assert env.queue.held_attempts(cell_id) == {
        m["attempt_id"] for m in first["molecules"]
    }


def test_claim_needs_registration_and_admission(env: Env) -> None:
    _, unregistered = env.new_cell()
    assert env.claim(unregistered).status_code == 409
    pending_id, pending = env.new_cell(admitted=False)
    env.register(pending)
    # Pending cells get the admission test, which needs the central's values.
    response = env.claim(pending)
    assert response.status_code == 503 and "cbs_09698" in response.json()["detail"]
    env.cells.set_admission(pending_id, ADMISSION_BLOCKED)
    response = env.claim(pending)
    assert response.status_code == 403 and "bloqueada" in response.json()["detail"]


def test_claim_count_is_validated(env: Env) -> None:
    _, headers = env.ready_cell()
    assert env.claim(headers, 0).status_code == 422


def test_a_cell_never_gets_two_spellings_of_one_molecule(tmp_path: Path) -> None:
    # Same MolecularIdentity: they would share the CREST directory and the
    # calculation_id with the central or between this cell's own workers.
    env = Env(
        tmp_path,
        [
            row("cbs_52946", smiles="C=C[C]=C", multiplicity=2, H298_cbs=78.63363),
            row("cbs_52947", smiles="C=[C]C=C", multiplicity=2, H298_cbs=86.79913),
        ],
    )
    _, headers = env.ready_cell()
    claimed = env.claim(headers, 2).json()["molecules"]
    assert [m["mol_id"] for m in claimed] == ["cbs_52946"]
    assert env.queue.claim("central", 2) == []


def test_a_cell_does_not_retry_its_own_failure(env: Env) -> None:
    _, headers = env.ready_cell()
    (first, _) = env.claim(headers, 2).json()["molecules"]
    env.client.post(
        "/api/cell/result",
        json={"attempt_id": first["attempt_id"], "failure": "crest timeout"},
        headers=headers,
    )
    assert [m["mol_id"] for m in env.claim(headers, 1).json()["molecules"]] == ["c"]


# ----------------------------------------------------------- heartbeat


def test_heartbeat_renews_and_reports_attempts_no_longer_held(env: Env) -> None:
    cell_id, headers = env.ready_cell()
    env.cells.touch(cell_id)
    before = env.cells.get(cell_id).last_seen_at
    held = [m["attempt_id"] for m in env.claim(headers, 1).json()["molecules"]]
    response = env.client.post(
        "/api/cell/heartbeat",
        json={"attempt_ids": [*held, "gone"]},
        headers=headers,
    )
    assert response.status_code == 200
    assert response.json() == {"renewed": 1, "not_held": ["gone"]}
    after = env.cells.get(cell_id).last_seen_at
    assert after is not None and before is not None and after > before


# -------------------------------------------------------------- result


def test_a_valid_result_is_accepted_with_its_geometry(env: Env) -> None:
    cell_id, headers = env.ready_cell()
    (molecule,) = env.claim(headers, 1).json()["molecules"]
    response = env.client.post(
        "/api/cell/result",
        json={"attempt_id": molecule["attempt_id"], "result": _result("a")},
        headers=headers,
    )
    assert response.status_code == 200 and response.json() == {"outcome": "accepted"}
    geometry = GeometryStore(env.db).get("a")
    assert geometry is not None and geometry.hof_kcal_mol == -12.5
    assert geometry.elements == ("C", "H")
    again = env.client.post(
        "/api/cell/result",
        json={"attempt_id": molecule["attempt_id"], "result": _result("a")},
        headers=headers,
    )
    assert again.json() == {"outcome": "duplicate"}
    assert env.queue.held_attempts(cell_id) == set()


def test_a_failure_spends_a_rerun_and_is_listed_with_the_cell(env: Env) -> None:
    cell_id, headers = env.ready_cell()
    (molecule,) = env.claim(headers, 1).json()["molecules"]
    response = env.client.post(
        "/api/cell/result",
        json={"attempt_id": molecule["attempt_id"], "failure": "saddle point"},
        headers=headers,
    )
    assert response.json() == {"outcome": "rerun"}
    ((mol_id, reason, _),) = JobState(env.db).recent_failures()
    assert mol_id == "a" and reason == f"[{cell_id}] saddle point"


@pytest.mark.parametrize(
    ("body", "status"),
    [
        ({"attempt_id": "nope", "failure": "x"}, 404),
        ({"attempt_id": "ATTEMPT"}, 422),
        ({"attempt_id": "ATTEMPT", "failure": "x", "result": _result("a")}, 422),
        ({"attempt_id": "ATTEMPT", "result": _result("a", signature="x")}, 422),
        ({"attempt_id": "ATTEMPT", "result": _result("b")}, 422),
        ({"attempt_id": "ATTEMPT", "result": _result("a", elements=["C"])}, 422),
        ({"attempt_id": "ATTEMPT", "result": _result("a", elements=[])}, 422),
    ],
)
def test_inconsistent_results_are_refused(
    env: Env, body: dict[str, Any], status: int
) -> None:
    _, headers = env.ready_cell()
    (molecule,) = env.claim(headers, 1).json()["molecules"]
    if body["attempt_id"] == "ATTEMPT":
        body = {**body, "attempt_id": molecule["attempt_id"]}
    response = env.client.post("/api/cell/result", json=body, headers=headers)
    assert response.status_code == status
    assert env.queue.held_attempts(env.cells.list()[0].cell_id) == {
        molecule["attempt_id"]
    }


def test_non_finite_numbers_are_refused(env: Env) -> None:
    _, headers = env.ready_cell()
    (molecule,) = env.claim(headers, 1).json()["molecules"]
    result = _result("a", hof_kcal_mol=float("nan"))
    # json.dumps writes NaN (not valid JSON, but what a careless cell sends).
    payload = json.dumps({"attempt_id": molecule["attempt_id"], "result": result})
    response = env.client.post(
        "/api/cell/result",
        content=payload,
        headers={**headers, "Content-Type": "application/json"},
    )
    assert response.status_code == 422


def test_a_cell_cannot_report_another_cells_attempt(env: Env) -> None:
    _, owner = env.ready_cell()
    _, intruder = env.ready_cell()
    (molecule,) = env.claim(owner, 1).json()["molecules"]
    response = env.client.post(
        "/api/cell/result",
        json={"attempt_id": molecule["attempt_id"], "failure": "x"},
        headers=intruder,
    )
    assert response.status_code == 403


def test_a_result_after_release_is_still_accepted_while_open(env: Env) -> None:
    cell_id, headers = env.ready_cell()
    (molecule,) = env.claim(headers, 1).json()["molecules"]
    env.queue.release(cell_id)
    response = env.client.post(
        "/api/cell/result",
        json={"attempt_id": molecule["attempt_id"], "result": _result("a")},
        headers=headers,
    )
    assert response.json() == {"outcome": "accepted"}


# ------------------------------------------------------- admission test


REFERENCE_HOF = -12.5


@pytest.fixture
def admission_env(tmp_path: Path) -> Env:
    env = Env(tmp_path, [*reference_rows(), row("a", nheavy=1), row("b", nheavy=2)])
    accept_references(env.db, hof=REFERENCE_HOF)
    return env


def _pending_cell(env: Env) -> tuple[str, dict[str, str]]:
    cell_id, headers = env.new_cell(admitted=False)
    assert env.register(headers, workers=5).status_code == 200
    return cell_id, headers


def _send(
    env: Env, headers: dict[str, str], molecule: dict[str, Any], **kw: Any
) -> Any:
    body: dict[str, Any] = {"attempt_id": molecule["attempt_id"]}
    if "failure" in kw:
        body["failure"] = kw["failure"]
    else:
        body["result"] = _result(molecule["mol_id"], **kw)
    return env.client.post("/api/cell/result", json=body, headers=headers)


def test_a_pending_cell_computes_the_references_then_gets_work(
    admission_env: Env,
) -> None:
    env = admission_env
    cell_id, headers = _pending_cell(env)
    response = env.claim(headers, 5)
    assert response.status_code == 200
    references = response.json()["molecules"]
    assert [m["mol_id"] for m in references] == list(ADMISSION_REFERENCES)
    assert all(m["attempt_id"].startswith("adm-") for m in references)
    assert references[1]["multiplicity"] == 2
    # Every worker is busy with a reference: nothing more, not even them again.
    assert env.claim(headers, 5).json()["molecules"] == []
    before = env.queue.held_attempts(cell_id), GeometryStore(env.db).count()

    for n, molecule in enumerate(references):
        reply = _send(env, headers, molecule, hof_kcal_mol=REFERENCE_HOF + 0.05)
        assert reply.status_code == 200 and reply.json() == {"outcome": "accepted"}
        admission = env.cells.get(cell_id).admission
        assert admission == (ADMISSION_PASSED if n == 4 else "pending")

    # The references never touched the queue nor the central's geometries.
    assert (env.queue.held_attempts(cell_id), GeometryStore(env.db).count()) == before
    assert GeometryStore(env.db).get("cbs_09698").hof_kcal_mol == REFERENCE_HOF  # type: ignore[union-attr]
    work = env.claim(headers, 5).json()["molecules"]
    assert [m["mol_id"] for m in work] == ["a", "b"]


def test_a_reference_outside_tolerance_blocks_the_cell(admission_env: Env) -> None:
    env = admission_env
    cell_id, headers = _pending_cell(env)
    first, second, *_ = env.claim(headers, 5).json()["molecules"]
    assert _send(env, headers, first, hof_kcal_mol=REFERENCE_HOF + 0.2).json() == {
        "outcome": "accepted"
    }
    assert env.cells.get(cell_id).admission == ADMISSION_BLOCKED
    assert env.claim(headers).status_code == 403
    # Late references are acknowledged (the outbox lets go) but not judged.
    late = _send(env, headers, second, failure="crest timeout")
    assert late.status_code == 200 and late.json() == {"outcome": "discarded"}
    note = AdmissionTest(env.db).note(cell_id, ADMISSION_BLOCKED)
    assert note is not None and note.startswith("cbs_09698 (CO): HOF -12.30000")


def test_a_failed_reference_blocks_the_cell(admission_env: Env) -> None:
    env = admission_env
    cell_id, headers = _pending_cell(env)
    (first, *_) = env.claim(headers, 5).json()["molecules"]
    assert _send(env, headers, first, failure="mopac error").status_code == 200
    assert env.cells.get(cell_id).admission == ADMISSION_BLOCKED
    # The failure of a reference is not a molecule failure at the central.
    assert JobState(env.db).recent_failures() == []


def test_reference_results_are_checked_like_any_result(admission_env: Env) -> None:
    env = admission_env
    _, headers = _pending_cell(env)
    first, second, *_ = env.claim(headers, 5).json()["molecules"]
    _other_id, other = _pending_cell(env)
    assert _send(env, other, first).status_code == 403
    assert _send(env, headers, first, signature="other").status_code == 422
    wrong = {**second, "mol_id": first["mol_id"]}
    assert _send(env, headers, wrong).status_code == 422
    unknown = {"attempt_id": "adm-nope", "mol_id": "cbs_09698"}
    assert _send(env, headers, unknown).status_code == 404


def test_the_heartbeat_keeps_reference_attempts(admission_env: Env) -> None:
    env = admission_env
    _, headers = _pending_cell(env)
    ids = [m["attempt_id"] for m in env.claim(headers, 2).json()["molecules"]]
    reply = env.client.post(
        "/api/cell/heartbeat", json={"attempt_ids": [*ids, "gone"]}, headers=headers
    ).json()
    assert reply == {"renewed": 2, "not_held": ["gone"]}


def test_a_restart_hands_unfinished_references_again(admission_env: Env) -> None:
    env = admission_env
    _, headers = _pending_cell(env)
    first, *_ = env.claim(headers, 5).json()["molecules"]
    _send(env, headers, first)
    assert env.register(headers, workers=5).status_code == 200
    again = env.claim(headers, 5).json()["molecules"]
    assert [m["mol_id"] for m in again] == list(ADMISSION_REFERENCES[1:])
