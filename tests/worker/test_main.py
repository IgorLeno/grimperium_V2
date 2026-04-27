"""Tests for worker/__main__.py — grimperium-worker CLI entry point."""

import socket
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from grimperium.worker.__main__ import (
    _connect_with_retry,
    _try_register,
    build_parser,
    main,
)
from grimperium.worker.client import ServerError
from grimperium.worker.runner import WorkerConfig

_SERVER_CFG: dict[str, Any] = {
    "worker_id": "w1",
    "hostname": "lab-01",
    "crest_timeout_minutes": 60,
    "mopac_timeout_minutes": 30,
    "batch_size": 10,
    "profile_name": "default",
}


# ── Parser ────────────────────────────────────────────────────────────────────


class TestParser:
    def test_defaults(self) -> None:
        args = build_parser().parse_args(["--server-url", "http://x:8000"])
        assert args.server_url == "http://x:8000"
        assert args.worker_id == socket.gethostname()
        assert args.api_token == ""
        assert args.max_molecules == 0
        assert args.idle_stop == 300

    def test_no_crest_timeout_flag(self) -> None:
        # Verify the removed flags are gone
        with pytest.raises(SystemExit):
            build_parser().parse_args(
                ["--server-url", "http://x", "--crest-timeout", "3600"]
            )

    def test_no_mopac_timeout_flag(self) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(
                ["--server-url", "http://x", "--mopac-timeout", "1800"]
            )

    def test_custom_args(self) -> None:
        args = build_parser().parse_args(
            [
                "--server-url",
                "http://lab:9000",
                "--worker-id",
                "rig-07",
                "--api-token",
                "secret",  # noqa: S106
                "--max-molecules",
                "50",
                "--idle-stop",
                "120",
            ]
        )
        assert args.server_url == "http://lab:9000"
        assert args.worker_id == "rig-07"
        assert args.api_token == "secret"  # noqa: S105
        assert args.max_molecules == 50
        assert args.idle_stop == 120

    def test_missing_server_url_exits_2(self) -> None:
        with pytest.raises(SystemExit) as exc:
            build_parser().parse_args([])
        assert exc.value.code == 2


# ── _try_register ──────────────────────────────────────────────────────────────


class TestTryRegister:
    def test_returns_dict_on_success(self) -> None:
        client = MagicMock()
        client.register.return_value = _SERVER_CFG.copy()
        result = _try_register(client)
        assert result == _SERVER_CFG

    def test_returns_none_on_server_error(self) -> None:
        client = MagicMock()
        client.register.side_effect = ServerError("connection refused")
        assert _try_register(client) is None

    def test_returns_none_on_any_exception(self) -> None:
        client = MagicMock()
        client.register.side_effect = RuntimeError("unexpected")
        assert _try_register(client) is None


# ── _connect_with_retry ────────────────────────────────────────────────────────


class TestConnectWithRetry:
    def test_returns_on_first_success(self) -> None:
        client = MagicMock()
        client.register.return_value = _SERVER_CFG.copy()
        result = _connect_with_retry(client, "http://x:8000")
        assert result["batch_size"] == 10
        assert client.register.call_count == 1

    def test_exits_1_on_failure_in_non_tty(self) -> None:
        client = MagicMock()
        client.register.side_effect = ServerError("refused")
        with (
            patch("grimperium.worker.__main__.sys.stdin.isatty", return_value=False),
            pytest.raises(SystemExit) as exc,
        ):
            _connect_with_retry(client, "http://x:8000")
        assert exc.value.code == 1

    def test_retries_after_user_chooses_retry(self) -> None:
        client = MagicMock()
        # First call fails, second succeeds
        client.register.side_effect = [ServerError("refused"), _SERVER_CFG.copy()]
        with (
            patch("grimperium.worker.__main__.sys.stdin.isatty", return_value=True),
            patch("grimperium.worker.__main__._show_retry_menu", return_value=True),
        ):
            result = _connect_with_retry(client, "http://x:8000")
        assert result["batch_size"] == 10
        assert client.register.call_count == 2

    def test_exits_0_when_user_chooses_encerrar(self) -> None:
        client = MagicMock()
        client.register.side_effect = ServerError("refused")
        with (
            patch("grimperium.worker.__main__.sys.stdin.isatty", return_value=True),
            patch("grimperium.worker.__main__._show_retry_menu", return_value=False),
            pytest.raises(SystemExit) as exc,
        ):
            _connect_with_retry(client, "http://x:8000")
        assert exc.value.code == 0


# ── main() ────────────────────────────────────────────────────────────────────


def _patched_main(argv: list[str]) -> tuple[int, MagicMock, MagicMock, MagicMock]:
    """Run main() with WorkerRunner/WorkerClient/connect patched.

    The mock client's register() returns _SERVER_CFG so config propagation is testable.
    """
    with (
        patch("grimperium.worker.__main__.WorkerRunner") as runner_cls,
        patch("grimperium.worker.__main__.WorkerClient") as client_cls,
        patch("grimperium.worker.__main__.setup_logging"),
        patch(
            "grimperium.worker.__main__._connect_with_retry",
            return_value=_SERVER_CFG.copy(),
        ),
    ):
        runner_instance = MagicMock()
        runner_instance.run.return_value = 0
        runner_cls.return_value = runner_instance
        code = main(argv)
        return code, runner_cls, client_cls, runner_instance


class TestMainWiring:
    def test_config_uses_server_timeouts(self) -> None:
        code, runner_cls, _, _ = _patched_main(["--server-url", "http://x"])
        assert code == 0
        config = runner_cls.call_args.kwargs["config"]
        assert isinstance(config, WorkerConfig)
        # Timeouts come from _SERVER_CFG, not from CLI flags
        assert config.crest_timeout_minutes == 60
        assert config.mopac_timeout_minutes == 30
        assert config.batch_size == 10

    def test_config_poll_interval_and_idle_stop(self) -> None:
        _, runner_cls, _, _ = _patched_main(
            ["--server-url", "http://x", "--idle-stop", "150"]
        )
        config = runner_cls.call_args.kwargs["config"]
        assert config.poll_interval_s == 5.0
        assert config.max_idle_polls == 30  # ceil(150 / 5)

    def test_config_server_url_and_worker_id(self) -> None:
        _, runner_cls, client_cls, _ = _patched_main(
            ["--server-url", "http://lab:9000", "--worker-id", "rig-07"]
        )
        config = runner_cls.call_args.kwargs["config"]
        assert config.server_url == "http://lab:9000"
        assert config.worker_id == "rig-07"

    def test_max_molecules_zero_becomes_none(self) -> None:
        _, _, _, runner_instance = _patched_main(
            ["--server-url", "http://x", "--max-molecules", "0"]
        )
        runner_instance.run.assert_called_once_with(max_molecules=None)

    def test_max_molecules_positive_passed_through(self) -> None:
        _, _, _, runner_instance = _patched_main(
            ["--server-url", "http://x", "--max-molecules", "42"]
        )
        runner_instance.run.assert_called_once_with(max_molecules=42)

    def test_returns_zero_on_success(self) -> None:
        code, _, _, _ = _patched_main(["--server-url", "http://x"])
        assert code == 0


class TestMainErrorPaths:
    def test_keyboard_interrupt_calls_stop_and_returns_zero(self) -> None:
        with (
            patch("grimperium.worker.__main__.WorkerRunner") as runner_cls,
            patch("grimperium.worker.__main__.WorkerClient"),
            patch("grimperium.worker.__main__.setup_logging"),
            patch(
                "grimperium.worker.__main__._connect_with_retry",
                return_value=_SERVER_CFG.copy(),
            ),
        ):
            runner_instance = MagicMock()
            runner_instance.run.side_effect = KeyboardInterrupt
            runner_cls.return_value = runner_instance

            code = main(["--server-url", "http://x"])

        assert code == 0
        runner_instance.stop.assert_called_once_with()

    def test_unexpected_exception_returns_one(self) -> None:
        with (
            patch("grimperium.worker.__main__.WorkerRunner") as runner_cls,
            patch("grimperium.worker.__main__.WorkerClient"),
            patch("grimperium.worker.__main__.setup_logging"),
            patch(
                "grimperium.worker.__main__._connect_with_retry",
                return_value=_SERVER_CFG.copy(),
            ),
        ):
            runner_instance = MagicMock()
            runner_instance.run.side_effect = RuntimeError("boom")
            runner_cls.return_value = runner_instance

            code = main(["--server-url", "http://x"])

        assert code == 1
