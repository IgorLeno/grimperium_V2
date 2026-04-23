"""Tests for worker/__main__.py — grimperium-worker CLI entry point."""

import socket
from unittest.mock import MagicMock, patch

import pytest

from grimperium.worker.__main__ import build_parser, main
from grimperium.worker.client import WorkerClientConfig
from grimperium.worker.runner import WorkerConfig

# ── Parser ────────────────────────────────────────────────────────────────────


class TestParser:
    def test_defaults(self) -> None:
        args = build_parser().parse_args(["--server-url", "http://x:8000"])
        assert args.server_url == "http://x:8000"
        assert args.worker_id == socket.gethostname()
        assert args.api_token == ""
        assert args.max_molecules == 0
        assert args.idle_stop == 300
        assert args.crest_timeout == 3600
        assert args.mopac_timeout == 1800

    def test_custom_args(self) -> None:
        args = build_parser().parse_args(
            [
                "--server-url",
                "http://lab:9000",
                "--worker-id",
                "rig-07",
                "--api-token",
                "secret",  # noqa: S106 — synthetic test token, not a real credential
                "--max-molecules",
                "50",
                "--idle-stop",
                "120",
                "--crest-timeout",
                "7200",
                "--mopac-timeout",
                "900",
            ]
        )
        assert args.server_url == "http://lab:9000"
        assert args.worker_id == "rig-07"
        assert args.api_token == "secret"  # noqa: S105
        assert args.max_molecules == 50
        assert args.idle_stop == 120
        assert args.crest_timeout == 7200
        assert args.mopac_timeout == 900

    def test_missing_server_url_exits_2(self) -> None:
        with pytest.raises(SystemExit) as exc:
            build_parser().parse_args([])
        assert exc.value.code == 2


# ── main() ────────────────────────────────────────────────────────────────────


def _patched_main(argv: list[str]) -> tuple[int, MagicMock, MagicMock, MagicMock]:
    """Run main() with WorkerRunner/WorkerClient/setup_logging patched.

    Returns (exit_code, runner_cls_mock, client_cls_mock, runner_instance_mock).
    """
    with (
        patch("grimperium.worker.__main__.WorkerRunner") as runner_cls,
        patch("grimperium.worker.__main__.WorkerClient") as client_cls,
        patch("grimperium.worker.__main__.setup_logging"),
    ):
        runner_instance = MagicMock()
        runner_instance.run.return_value = 0
        runner_cls.return_value = runner_instance
        code = main(argv)
        return code, runner_cls, client_cls, runner_instance


class TestMainWiring:
    def test_builds_config_with_unit_conversion(self) -> None:
        code, runner_cls, client_cls, _ = _patched_main(
            [
                "--server-url",
                "http://x",
                "--worker-id",
                "w1",
                "--api-token",
                "t",
                "--crest-timeout",
                "7200",
                "--mopac-timeout",
                "900",
                "--idle-stop",
                "150",
            ]
        )
        assert code == 0

        # WorkerConfig passed to WorkerRunner
        config = runner_cls.call_args.kwargs["config"]
        assert isinstance(config, WorkerConfig)
        assert config.server_url == "http://x"
        assert config.worker_id == "w1"
        assert config.api_token == "t"  # noqa: S105
        assert config.crest_timeout_minutes == 120  # 7200 / 60
        assert config.mopac_timeout_minutes == 15  # 900 / 60
        assert config.poll_interval_s == 5.0
        assert config.max_idle_polls == 30  # ceil(150 / 5)

        # WorkerClientConfig passed to WorkerClient
        client_config = client_cls.call_args.args[0]
        assert isinstance(client_config, WorkerClientConfig)
        assert client_config.server_url == "http://x"
        assert client_config.worker_id == "w1"
        assert client_config.api_token == "t"  # noqa: S105

    def test_sub_minute_timeouts_clamp_to_one(self) -> None:
        # 30 seconds must not collapse to 0 minutes (would disable the timeout).
        _, runner_cls, _, _ = _patched_main(
            [
                "--server-url",
                "http://x",
                "--crest-timeout",
                "30",
                "--mopac-timeout",
                "30",
                "--idle-stop",
                "1",
            ]
        )
        config = runner_cls.call_args.kwargs["config"]
        assert config.crest_timeout_minutes == 1
        assert config.mopac_timeout_minutes == 1
        assert config.max_idle_polls == 1

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


class TestMainErrorPaths:
    def test_keyboard_interrupt_calls_stop_and_returns_zero(self) -> None:
        with (
            patch("grimperium.worker.__main__.WorkerRunner") as runner_cls,
            patch("grimperium.worker.__main__.WorkerClient"),
            patch("grimperium.worker.__main__.setup_logging"),
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
        ):
            runner_instance = MagicMock()
            runner_instance.run.side_effect = RuntimeError("boom")
            runner_cls.return_value = runner_instance

            code = main(["--server-url", "http://x"])

        assert code == 1
