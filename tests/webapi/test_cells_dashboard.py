"""Dashboard endpoints for the cells: list, create, rename, release, revoke."""

from __future__ import annotations

from contextlib import closing
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from grimperium.webapi.app import create_app
from grimperium.webapi.config import DashboardPaths
from semi_imperium.training.cells import ADMISSION_PASSED, AdmissionTest, CellRegistry
from semi_imperium.training.db import connect, database_path
from semi_imperium.training.molecule_queue import TrainingQueue
from tests.training.conftest import accept_references, reference_rows, row

LOCAL = {"Origin": "http://127.0.0.1:5173"}


class Env:
    def __init__(self, tmp_path: Path) -> None:
        pd.DataFrame([row("a", nheavy=1), row("b", nheavy=2)]).to_csv(
            tmp_path / "cbs.csv", index=False
        )
        self.paths = DashboardPaths(
            semi_imperium_store=tmp_path / "store",
            exp_csv=tmp_path / "exp.csv",
            exp_vs_cbs_csv=tmp_path / "vs_cbs.csv",
            cbs_csv=tmp_path / "cbs.csv",
        )
        self.db = database_path(self.paths.semi_imperium_store)
        self.queue = TrainingQueue(self.db, self.paths.cbs_csv)
        self.cells = CellRegistry(self.db)
        self.client = TestClient(create_app(self.paths, launcher=lambda _: None))

    def create(self, name: str = "Lab 1") -> dict[str, object]:
        response = self.client.post("/api/cells", json={"name": name}, headers=LOCAL)
        assert response.status_code == 200
        body: dict[str, object] = response.json()
        return body


@pytest.fixture
def env(tmp_path: Path) -> Env:
    return Env(tmp_path)


def test_listing_without_a_database_is_empty_and_writes_nothing(env: Env) -> None:
    assert env.client.get("/api/cells").json() == []
    assert not env.db.exists()


def test_create_shows_the_token_once(env: Env) -> None:
    body = env.create("PC do Prof. Silva")
    token = body["token"]
    assert isinstance(token, str) and len(token) >= 40
    cell = body["cell"]
    assert isinstance(cell, dict)
    assert cell["name"] == "PC do Prof. Silva" and cell["state"] == "never_seen"
    assert env.cells.authenticate(token) is not None
    listed = env.client.get("/api/cells").json()
    assert [c["cell_id"] for c in listed] == [cell["cell_id"]]
    assert "token" not in listed[0] and "token_hash" not in listed[0]


def test_list_shows_state_and_molecules_held(env: Env) -> None:
    cell_id = env.create()["cell"]["cell_id"]  # type: ignore[index]
    env.cells.record_registration(
        cell_id,
        host="lab-pc",
        cores=8,
        workers=2,
        crest_threads=4,
        program_versions={"crest": "3.0.2"},
        protocol_signature="sig",
        refusal=None,
    )
    env.cells.set_admission(cell_id, ADMISSION_PASSED)
    env.queue.ensure_seeded()
    env.queue.claim(cell_id, 2)
    (cell,) = env.client.get("/api/cells").json()
    assert cell["state"] == "online" and cell["held"] == 2
    assert cell["host"] == "lab-pc" and cell["admission"] == ADMISSION_PASSED
    assert cell["program_versions"] == {"crest": "3.0.2"}
    # Staleness is judged by the central's clock at request time.
    old = (datetime.now(timezone.utc) - timedelta(minutes=10)).isoformat()
    with closing(connect(env.db)) as connection:
        connection.execute("UPDATE cells SET last_seen_at = ?", (old,))
    assert env.client.get("/api/cells").json()[0]["state"] == "silent"


def test_rename(env: Env) -> None:
    cell_id = env.create()["cell"]["cell_id"]  # type: ignore[index]
    response = env.client.post(
        f"/api/cells/{cell_id}/rename", json={"name": "  Lab 2 "}, headers=LOCAL
    )
    assert response.status_code == 200 and response.json()["name"] == "Lab 2"
    blank = env.client.post(
        f"/api/cells/{cell_id}/rename", json={"name": "  "}, headers=LOCAL
    )
    assert blank.status_code == 422


def test_release_returns_the_cells_molecules_to_the_queue(env: Env) -> None:
    cell_id = env.create()["cell"]["cell_id"]  # type: ignore[index]
    env.queue.ensure_seeded()
    env.queue.claim(cell_id, 2)
    response = env.client.post(f"/api/cells/{cell_id}/release", headers=LOCAL)
    assert response.status_code == 200 and response.json() == {"released": ["a", "b"]}
    assert env.queue.held_attempts(cell_id) == set()


def test_revoke_invalidates_the_token(env: Env) -> None:
    body = env.create()
    cell_id = body["cell"]["cell_id"]  # type: ignore[index]
    response = env.client.post(f"/api/cells/{cell_id}/revoke", headers=LOCAL)
    assert response.status_code == 200 and response.json()["state"] == "revoked"
    assert env.cells.authenticate(str(body["token"])) is None


@pytest.mark.parametrize("action", ["rename", "release", "revoke"])
def test_unknown_cells_are_404_and_central_is_not_a_cell(env: Env, action: str) -> None:
    env.create()
    for cell_id in ("nope", "central"):
        response = env.client.post(
            f"/api/cells/{cell_id}/{action}", json={"name": "x"}, headers=LOCAL
        )
        assert response.status_code == 404


@pytest.mark.parametrize(
    ("path", "body"),
    [
        ("/api/cells", {"name": "x"}),
        ("/api/cells/ID/rename", {"name": "x"}),
        ("/api/cells/ID/release", None),
        ("/api/cells/ID/revoke", None),
    ],
)
def test_actions_refuse_foreign_origins(
    env: Env, path: str, body: dict[str, str] | None
) -> None:
    cell_id = env.create()["cell"]["cell_id"]  # type: ignore[index]
    response = env.client.post(
        path.replace("ID", str(cell_id)),
        json=body,
        headers={"Origin": "http://evil.example"},
    )
    assert response.status_code == 403


def test_the_list_says_why_the_admission_test_blocked(env: Env, tmp_path: Path) -> None:
    pd.DataFrame(reference_rows()).to_csv(tmp_path / "cbs.csv", index=False)
    env.queue.ensure_seeded()
    accept_references(env.db, hof=-48.9388)
    cell_id = str(env.create()["cell"]["cell_id"])  # type: ignore[index]
    test = AdmissionTest(env.db)
    (first,) = test.claim(cell_id, 1)
    test.record(first.attempt_id, hof_kcal_mol=-48.5, failure=None)
    (cell,) = env.client.get("/api/cells").json()
    assert cell["admission"] == "blocked" and cell["state"] == "blocked"
    assert cell["admission_note"].startswith("cbs_09698 (CO): HOF -48.50000")

    # "Refazer admissão": back to pending, the earlier results set aside.
    response = env.client.post(
        f"/api/cells/{cell_id}/admission", json={"admission": "pending"}, headers=LOCAL
    )
    assert response.status_code == 200
    assert response.json()["admission"] == "pending"
    assert response.json()["admission_note"] is None
    assert [m.mol_id for m in test.claim(cell_id, 1)] == ["cbs_09698"]
