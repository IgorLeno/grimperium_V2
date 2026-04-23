"""Tests for session_store module."""

import json

import pytest

from grimperium.cli.session_store import (
    SessionState,
    WorkerSessionInfo,
    delete_session,
    load_session,
    save_session,
)


class TestWorkerSessionInfo:
    def test_defaults(self) -> None:
        w = WorkerSessionInfo(worker_id="IQ", hostname="iq-machine")
        assert w.batch_size == 10
        assert w.profile_name == "default"
        assert w.crest_timeout_minutes == 60
        assert w.mopac_timeout_minutes == 30

    def test_to_dict_round_trip(self) -> None:
        w = WorkerSessionInfo(worker_id="A", hostname="h1", batch_size=5)
        restored = WorkerSessionInfo.from_dict(w.to_dict())
        assert restored.worker_id == "A"
        assert restored.batch_size == 5

    def test_from_dict_ignores_unknown_keys(self) -> None:
        w = WorkerSessionInfo.from_dict(
            {"worker_id": "X", "hostname": "h", "extra": "ignored"}
        )
        assert w.worker_id == "X"


class TestSessionState:
    def test_empty_workers(self) -> None:
        s = SessionState(started_at="2026-04-23T14:00:00")
        assert s.workers == []
        assert s.server_url == "http://localhost:8000"

    def test_to_dict_round_trip(self) -> None:
        workers = [WorkerSessionInfo(worker_id="IQ", hostname="iq")]
        s = SessionState(
            started_at="2026-04-23T14:00:00",
            server_url="http://10.0.0.1:8000",
            workers=workers,
        )
        d = s.to_dict()
        restored = SessionState.from_dict(d)
        assert restored.started_at == "2026-04-23T14:00:00"
        assert restored.server_url == "http://10.0.0.1:8000"
        assert len(restored.workers) == 1
        assert restored.workers[0].worker_id == "IQ"

    def test_from_dict_missing_workers_defaults_empty(self) -> None:
        s = SessionState.from_dict({"started_at": "2026-04-23T00:00:00"})
        assert s.workers == []


class TestSessionPersistence:
    def _patch_session_path(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: pytest.TempPathFactory
    ) -> None:
        import grimperium.cli.session_store as ss

        monkeypatch.setattr(ss, "_session_path", lambda: tmp_path / "session.json")  # type: ignore[arg-type]

    def test_load_session_returns_none_when_missing(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: pytest.TempPathFactory
    ) -> None:
        self._patch_session_path(monkeypatch, tmp_path)
        assert load_session() is None

    def test_save_and_load_round_trip(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: pytest.TempPathFactory
    ) -> None:
        self._patch_session_path(monkeypatch, tmp_path)
        state = SessionState(
            started_at="2026-04-23T14:00:00",
            workers=[WorkerSessionInfo(worker_id="IQ", hostname="iq")],
        )
        assert save_session(state) is True
        loaded = load_session()
        assert loaded is not None
        assert loaded.started_at == "2026-04-23T14:00:00"
        assert loaded.workers[0].worker_id == "IQ"

    def test_delete_session_removes_file(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: pytest.TempPathFactory
    ) -> None:
        self._patch_session_path(monkeypatch, tmp_path)
        save_session(SessionState(started_at="2026-04-23T00:00:00"))
        assert delete_session() is True
        assert load_session() is None

    def test_delete_session_returns_false_when_missing(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: pytest.TempPathFactory
    ) -> None:
        self._patch_session_path(monkeypatch, tmp_path)
        assert delete_session() is False

    def test_load_session_returns_none_on_corrupt_json(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: pytest.TempPathFactory
    ) -> None:
        self._patch_session_path(monkeypatch, tmp_path)
        (tmp_path / "session.json").write_text("not json", encoding="utf-8")  # type: ignore[operator]
        assert load_session() is None

    def test_saved_json_is_readable(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: pytest.TempPathFactory
    ) -> None:
        self._patch_session_path(monkeypatch, tmp_path)
        state = SessionState(started_at="2026-04-23T14:00:00")
        save_session(state)
        raw = json.loads((tmp_path / "session.json").read_text(encoding="utf-8"))  # type: ignore[operator]
        assert raw["started_at"] == "2026-04-23T14:00:00"
        assert "workers" in raw
