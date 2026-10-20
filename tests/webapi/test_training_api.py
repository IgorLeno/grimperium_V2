from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from grimperium.webapi.app import create_app
from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.training import TrainingLaunch
from semi_imperium.training.db import database_path
from semi_imperium.training.job_state import JobState
from semi_imperium.training.molecule_queue import TrainingQueue
from tests.training.conftest import row

CBS_ROWS = [row("a", nheavy=1), row("b", nheavy=2), row("c", nheavy=3)]


def _paths(
    tmp_path: Path, rows: list[dict[str, object]] | None = None
) -> DashboardPaths:
    pd.DataFrame(rows or CBS_ROWS).to_csv(tmp_path / "cbs.csv", index=False)
    for name in ("exp.csv", "vs_cbs.csv"):
        (tmp_path / name).write_text("x\n")
    return DashboardPaths(
        semi_imperium_store=tmp_path / "store",
        exp_csv=tmp_path / "exp.csv",
        exp_vs_cbs_csv=tmp_path / "vs_cbs.csv",
        cbs_csv=tmp_path / "cbs.csv",
    )


class Recorder:
    def __init__(self) -> None:
        self.launches: list[TrainingLaunch] = []

    def __call__(self, launch: TrainingLaunch) -> None:
        self.launches.append(launch)


def _client(paths: DashboardPaths, recorder: Recorder) -> TestClient:
    return TestClient(create_app(paths, launcher=recorder))


def test_status_reports_progress_and_an_idle_job(tmp_path: Path) -> None:
    client = _client(_paths(tmp_path), Recorder())
    body = client.get("/api/training/status").json()
    assert body["progress"] == {
        "total": 3,
        "ok": 0,
        "pending": 3,
        "rerun": 0,
        "running": 0,
        "skipped": 0,
    }
    assert body["job"]["state"] == "idle"
    assert body["geometries"] == 0 and body["failures"] == []
    assert not database_path(tmp_path / "store").exists()  # GET never creates it


def test_status_counts_results_from_the_database(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    queue = TrainingQueue(database_path(paths.semi_imperium_store), paths.cbs_csv)
    queue.ensure_seeded()
    queue.record_success(
        "a", hof_kcal_mol=-1.0, conformers_evaluated=1, elapsed_s=1.0, batch_id="b"
    )
    body = _client(paths, Recorder()).get("/api/training/status").json()
    assert (body["progress"]["total"], body["progress"]["ok"]) == (3, 1)
    assert "protocol1" not in body["progress"] and "migrated" not in body["progress"]


def test_status_with_a_missing_cbs_is_503(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    paths.cbs_csv.unlink()
    response = _client(paths, Recorder()).get("/api/training/status")
    assert response.status_code == 503


def test_start_launches_with_the_requested_limits(tmp_path: Path) -> None:
    recorder = Recorder()
    response = _client(_paths(tmp_path), recorder).post(
        "/api/training/start", json={"workers": 2, "crest_threads": 6, "limit": 5}
    )
    assert response.status_code == 200
    (launch,) = recorder.launches
    assert (launch.workers, launch.crest_threads, launch.limit) == (2, 6, 5)
    assert launch.cbs_csv == tmp_path / "cbs.csv"


def test_start_is_refused_while_a_job_is_alive(tmp_path: Path) -> None:
    recorder = Recorder()
    paths = _paths(tmp_path)
    JobState(database_path(paths.semi_imperium_store)).begin(
        pid=os.getpid(), workers=1, limit=None
    )
    response = _client(paths, recorder).post("/api/training/start", json={})
    assert response.status_code == 409
    assert recorder.launches == []


def test_a_stale_job_does_not_block_start(tmp_path: Path) -> None:
    recorder = Recorder()
    paths = _paths(tmp_path)
    JobState(database_path(paths.semi_imperium_store)).begin(
        pid=2_000_000_000, workers=1, limit=None
    )
    status = _client(paths, recorder).get("/api/training/status").json()
    assert status["job"]["state"] == "stale"
    started = _client(paths, recorder).post("/api/training/start", json={})
    assert started.status_code == 200
    assert len(recorder.launches) == 1


def test_start_defaults_to_four_crest_threads(tmp_path: Path) -> None:
    recorder = Recorder()
    _client(_paths(tmp_path), recorder).post("/api/training/start", json={})
    assert recorder.launches[0].crest_threads == 4


def test_status_reports_activity_and_cpu_count(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    queue = TrainingQueue(database_path(paths.semi_imperium_store), paths.cbs_csv)
    queue.ensure_seeded()
    JobState(database_path(paths.semi_imperium_store)).begin(
        pid=os.getpid(), workers=2, limit=None
    )
    queue.mark_running("b")
    queue.mark_running("a")
    queue.record_success(
        "a", hof_kcal_mol=-20.0, conformers_evaluated=3, elapsed_s=12.0, batch_id="x"
    )
    body = _client(paths, Recorder()).get("/api/training/status").json()
    activity = body["activity"]
    assert [m["mol_id"] for m in activity["in_flight"]] == ["b"]
    assert activity["in_flight"][0]["smiles"] == "smiles-b"
    (done,) = activity["recent"]
    assert (done["mol_id"], done["h298_pm7"], done["elapsed_s"]) == ("a", -20.0, 12.0)
    assert done["delta_kcal_mol"] == pytest.approx(-17.8 + 20.0)
    assert activity["mean_elapsed_s"] == 12.0
    assert body["cpu_count"] >= 1


def test_status_activity_is_empty_without_a_database(tmp_path: Path) -> None:
    body = _client(_paths(tmp_path), Recorder()).get("/api/training/status").json()
    assert body["activity"] == {"in_flight": [], "recent": [], "mean_elapsed_s": None}


@pytest.mark.parametrize(
    "body",
    [
        {"workers": 0},
        {"workers": 9},
        {"limit": 0},
        {"crest_threads": 0},
        {"crest_threads": 65},
    ],
)
def test_start_validates_bounds(tmp_path: Path, body: dict[str, int]) -> None:
    recorder = Recorder()
    response = _client(_paths(tmp_path), recorder).post(
        "/api/training/start", json=body
    )
    assert response.status_code == 422
    assert recorder.launches == []


def test_stop_flags_the_running_job(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    state = JobState(database_path(paths.semi_imperium_store))
    state.begin(pid=os.getpid(), workers=1, limit=None)
    response = _client(paths, Recorder()).post("/api/training/stop")
    assert response.status_code == 200
    assert state.stop_requested() is True


def test_stop_without_a_job_is_409(tmp_path: Path) -> None:
    response = _client(_paths(tmp_path), Recorder()).post("/api/training/stop")
    assert response.status_code == 409


def test_posts_from_foreign_origins_are_rejected(tmp_path: Path) -> None:
    recorder = Recorder()
    client = _client(_paths(tmp_path), recorder)
    foreign = client.post(
        "/api/training/start", json={}, headers={"Origin": "https://evil.example"}
    )
    assert foreign.status_code == 403
    local = client.post(
        "/api/training/start", json={}, headers={"Origin": "http://localhost:5173"}
    )
    assert local.status_code == 200
    assert len(recorder.launches) == 1


def test_health_says_the_api_can_dispatch(tmp_path: Path) -> None:
    body = _client(_paths(tmp_path), Recorder()).get("/api/health").json()
    assert body == {"status": "ok", "can_dispatch": True}
