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
from tests.training.conftest import row

SCHEMA_ROWS = [
    row("a", "OK", protocol="1"),
    row("b", "Pending", protocol=""),
    row("c", "Pending", protocol=""),
]


def _paths(
    tmp_path: Path, rows: list[dict[str, object]] | None = None
) -> DashboardPaths:
    from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager

    schema = BatchCSVManager(csv_path=None).get_schema()
    pd.DataFrame(rows or SCHEMA_ROWS).reindex(columns=schema).to_csv(
        tmp_path / "pm7.csv", index=False
    )
    for name in ("exp.csv", "vs_cbs.csv"):
        (tmp_path / name).write_text("x\n")
    return DashboardPaths(
        pm7_csv=tmp_path / "pm7.csv",
        semi_imperium_store=tmp_path / "store",
        exp_csv=tmp_path / "exp.csv",
        exp_vs_cbs_csv=tmp_path / "vs_cbs.csv",
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
    assert body["progress"]["total"] == 3
    assert (body["progress"]["protocol1"], body["progress"]["pending"]) == (1, 2)
    assert body["progress"]["migrated"] is True
    assert body["job"]["state"] == "idle"
    assert body["geometries"] == 0 and body["failures"] == []
    assert not database_path(tmp_path / "store").exists()  # GET never creates it


def test_status_of_an_unmigrated_csv_is_not_an_error(tmp_path: Path) -> None:
    paths = _paths(tmp_path, [row("a", "OK"), row("b", "Pending")])
    body = _client(paths, Recorder()).get("/api/training/status").json()
    assert body["progress"]["migrated"] is False
    assert body["progress"]["ok_unmarked"] == 1


def test_status_with_a_missing_csv_is_503(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    paths.pm7_csv.unlink()
    response = _client(paths, Recorder()).get("/api/training/status")
    assert response.status_code == 503


def test_start_launches_with_the_requested_limits(tmp_path: Path) -> None:
    recorder = Recorder()
    response = _client(_paths(tmp_path), recorder).post(
        "/api/training/start", json={"workers": 2, "limit": 5}
    )
    assert response.status_code == 200
    (launch,) = recorder.launches
    assert (launch.workers, launch.limit) == (2, 5)
    assert launch.csv_path == tmp_path / "pm7.csv"


def test_start_is_refused_when_the_csv_is_not_migrated(tmp_path: Path) -> None:
    recorder = Recorder()
    paths = _paths(tmp_path, [row("a", "OK"), row("b", "Pending")])
    response = _client(paths, recorder).post("/api/training/start", json={})
    assert response.status_code == 409
    assert "migrate --apply" in response.json()["detail"]
    assert recorder.launches == []


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


@pytest.mark.parametrize("body", [{"workers": 0}, {"workers": 9}, {"limit": 0}])
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
