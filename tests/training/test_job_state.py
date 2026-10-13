from __future__ import annotations

import os
from pathlib import Path

import pytest

from semi_imperium.training.job_state import (
    JobAlreadyRunning,
    JobState,
)


def _state(tmp_path: Path) -> JobState:
    return JobState(tmp_path / "training.sqlite3")


def test_missing_database_is_idle_and_is_not_created(tmp_path: Path) -> None:
    state = _state(tmp_path)
    snapshot = state.snapshot()
    assert snapshot.state == "idle" and not snapshot.is_active
    assert state.recent_failures() == []
    assert state.stop_requested() is False
    assert not (tmp_path / "training.sqlite3").exists()


def test_begin_heartbeat_finish_cycle(tmp_path: Path) -> None:
    state = _state(tmp_path)
    state.begin(pid=os.getpid(), workers=3, limit=10)
    snapshot = state.snapshot()
    assert snapshot.state == "running" and snapshot.is_active
    assert (snapshot.workers, snapshot.limit, snapshot.pid) == (3, 10, os.getpid())
    state.heartbeat(done=4, failed=1)
    assert (state.snapshot().done, state.snapshot().failed) == (4, 1)
    state.finish("finished", "ok")
    assert state.snapshot().state == "finished" and not state.snapshot().is_active


def test_second_begin_is_refused_while_alive(tmp_path: Path) -> None:
    state = _state(tmp_path)
    state.begin(pid=os.getpid(), workers=1, limit=None)
    with pytest.raises(JobAlreadyRunning):
        state.begin(pid=os.getpid(), workers=1, limit=None)


def test_dead_pid_is_reported_stale_and_can_be_replaced(tmp_path: Path) -> None:
    state = _state(tmp_path)
    state.begin(pid=2_000_000_000, workers=1, limit=None)  # no such process
    assert state.snapshot().state == "stale"
    assert state.snapshot().is_active is False
    state.begin(pid=os.getpid(), workers=2, limit=None)
    assert state.snapshot().state == "running"


def test_old_heartbeat_with_live_pid_is_stale(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import semi_imperium.training.job_state as module

    state = _state(tmp_path)
    state.begin(pid=os.getpid(), workers=1, limit=None)
    monkeypatch.setattr(module, "STALE_AFTER_SECONDS", -1.0)
    assert state.snapshot().state == "stale"


def test_request_stop_only_when_active(tmp_path: Path) -> None:
    state = _state(tmp_path)
    assert state.request_stop() is False
    state.begin(pid=os.getpid(), workers=1, limit=None)
    assert state.request_stop() is True
    assert state.stop_requested() is True
    state.mark_stopping()
    assert state.snapshot().state == "stopping"
    assert state.snapshot().is_active is True


def test_failures_are_kept_newest_first(tmp_path: Path) -> None:
    state = _state(tmp_path)
    state.record_failure("a", "saddle point")
    state.record_failure("b", "CREST timeout")
    failures = state.recent_failures(limit=1)
    assert [(mol, reason) for mol, reason, _ in failures] == [("b", "CREST timeout")]
    assert len(state.recent_failures()) == 2
