"""Tests for server/config.py — ServerConfig via pydantic-settings."""

import pytest
from pydantic import ValidationError

from grimperium.server.config import ServerConfig


class TestServerConfigDefaults:
    def test_required_csv_path_raises_without_env(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("GRIMPERIUM_CSV_PATH", raising=False)
        with pytest.raises(ValidationError):
            ServerConfig(_env_file=None, **{"csv_path": None})  # type: ignore[arg-type]

    def test_defaults_all_optional_fields(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("GRIMPERIUM_CSV_PATH", "/tmp/test.csv")
        cfg = ServerConfig()
        assert cfg.api_token == ""
        assert cfg.heartbeat_timeout_s == 300
        assert cfg.watchdog_interval_s == 30
        assert cfg.startup_grace_s == 120
        assert cfg.stuck_assigned_threshold_h == 2.0
        assert cfg.host == "0.0.0.0"
        assert cfg.port == 8000
        assert cfg.max_reruns == 3

    def test_csv_path_from_env(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("GRIMPERIUM_CSV_PATH", "/data/molecules.csv")
        cfg = ServerConfig()
        assert cfg.csv_path == "/data/molecules.csv"

    def test_api_token_from_env(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("GRIMPERIUM_CSV_PATH", "/tmp/test.csv")
        monkeypatch.setenv("GRIMPERIUM_API_TOKEN", "secret-token-123")
        cfg = ServerConfig()
        assert cfg.api_token == "secret-token-123"

    def test_heartbeat_timeout_overridable(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("GRIMPERIUM_CSV_PATH", "/tmp/test.csv")
        monkeypatch.setenv("GRIMPERIUM_HEARTBEAT_TIMEOUT_S", "600")
        cfg = ServerConfig()
        assert cfg.heartbeat_timeout_s == 600

    def test_port_overridable(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("GRIMPERIUM_CSV_PATH", "/tmp/test.csv")
        monkeypatch.setenv("GRIMPERIUM_PORT", "9000")
        cfg = ServerConfig()
        assert cfg.port == 9000

    def test_stuck_threshold_overridable(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("GRIMPERIUM_CSV_PATH", "/tmp/test.csv")
        monkeypatch.setenv("GRIMPERIUM_STUCK_ASSIGNED_THRESHOLD_H", "4.5")
        cfg = ServerConfig()
        assert cfg.stuck_assigned_threshold_h == 4.5

    def test_constructed_directly(self) -> None:
        cfg = ServerConfig(csv_path="/tmp/direct.csv", port=7777)
        assert cfg.csv_path == "/tmp/direct.csv"
        assert cfg.port == 7777
        assert cfg.host == "0.0.0.0"
