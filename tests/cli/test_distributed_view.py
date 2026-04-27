"""Tests for DatabasesView distributed mode — state machine refactor."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
from rich.console import Console

from grimperium.cli.views.databases_view import DatabasesView

# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture
def view(tmp_path: Path) -> DatabasesView:
    controller = MagicMock()
    controller.console = Console(force_terminal=False, no_color=True)
    v = DatabasesView(controller)
    db = MagicMock()
    db.csv_path = "thermo_pm7.csv"
    v.selected_db = db
    return v


def _worker(
    worker_id: str, hostname: str = "host", mol_id_current: str | None = None
) -> dict[str, Any]:
    return {
        "worker_id": worker_id,
        "hostname": hostname,
        "last_seen": "2026-04-22T10:00:00Z",
        "mol_id_current": mol_id_current,
    }


def _worker_extended(
    worker_id: str,
    hostname: str = "host",
    processed: int = 0,
    successful: int = 0,
    failed: int = 0,
    skipped: int = 0,
    current_mol_id: str | None = None,
) -> dict[str, Any]:
    return {
        "worker_id": worker_id,
        "hostname": hostname,
        "processed": processed,
        "successful": successful,
        "failed": failed,
        "skipped": skipped,
        "current_mol_id": current_mol_id,
    }


# ── _handle_run_calculation_menu ──────────────────────────────────────────────


class TestHandleRunCalculationMenu:
    @patch("grimperium.cli.views.databases_view.show_back_menu", return_value="local")
    def test_local_mode_calls_handle_calculate_pm7(
        self, _mock_menu: MagicMock, view: DatabasesView
    ) -> None:
        view.handle_calculate_pm7 = MagicMock()  # type: ignore[method-assign]
        view._handle_run_calculation_menu()
        view.handle_calculate_pm7.assert_called_once()

    @patch(
        "grimperium.cli.views.databases_view.show_back_menu", return_value="distributed"
    )
    def test_distributed_mode_calls_distributed_handler(
        self, _mock_menu: MagicMock, view: DatabasesView
    ) -> None:
        view._handle_distributed_mode = MagicMock()  # type: ignore[method-assign]
        view._handle_run_calculation_menu()
        view._handle_distributed_mode.assert_called_once()

    @patch("grimperium.cli.views.databases_view.show_back_menu", return_value=None)
    def test_cancel_does_nothing(
        self, _mock_menu: MagicMock, view: DatabasesView
    ) -> None:
        view.handle_calculate_pm7 = MagicMock()  # type: ignore[method-assign]
        view._handle_distributed_mode = MagicMock()  # type: ignore[method-assign]
        view._handle_run_calculation_menu()
        view.handle_calculate_pm7.assert_not_called()
        view._handle_distributed_mode.assert_not_called()

    @patch("grimperium.cli.views.databases_view.show_back_menu", return_value="local")
    def test_local_mode_does_not_call_distributed(
        self, _mock_menu: MagicMock, view: DatabasesView
    ) -> None:
        view.handle_calculate_pm7 = MagicMock()  # type: ignore[method-assign]
        view._handle_distributed_mode = MagicMock()  # type: ignore[method-assign]
        view._handle_run_calculation_menu()
        view._handle_distributed_mode.assert_not_called()


# ── _render_workers_table ─────────────────────────────────────────────────────


class TestRenderWorkersTable:
    def test_renders_empty_table(self, view: DatabasesView) -> None:
        table = view._render_workers_table([])
        assert len(table.columns) == 4

    def test_renders_workers(self, view: DatabasesView) -> None:
        workers = [_worker("w1", "lab-01"), _worker("w2", "lab-02")]
        table = view._render_workers_table(workers)
        assert len(table.columns) == 4
        assert table.row_count == 2

    def test_renders_mol_id_in_processing_column(self, view: DatabasesView) -> None:
        workers = [_worker("w1", "lab-01", mol_id_current="mol_042")]
        table = view._render_workers_table(workers)
        assert table.row_count == 1


# ── _state_check_port ─────────────────────────────────────────────────────────


class TestStateCheckPort:
    def test_free_port_returns_check_session(self, view: DatabasesView) -> None:
        with patch.object(DatabasesView, "_is_port_free", return_value=True):
            assert view._state_check_port() == "check_session"

    def test_busy_port_with_server_responding_returns_check_session(
        self, view: DatabasesView
    ) -> None:
        with (
            patch.object(DatabasesView, "_is_port_free", return_value=False),
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
        ):
            assert view._state_check_port() == "check_session"

    def test_busy_port_user_retries_then_free(self, view: DatabasesView) -> None:
        with (
            patch.object(DatabasesView, "_is_port_free", side_effect=[False, True]),
            patch.object(DatabasesView, "_server_is_responding", return_value=False),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu",
                return_value="retry",
            ),
        ):
            assert view._state_check_port() == "check_session"

    def test_busy_port_user_exits_returns_none(self, view: DatabasesView) -> None:
        with (
            patch.object(DatabasesView, "_is_port_free", return_value=False),
            patch.object(DatabasesView, "_server_is_responding", return_value=False),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu",
                return_value=None,
            ),
        ):
            assert view._state_check_port() is None


# ── _state_check_session ──────────────────────────────────────────────────────


class TestStateCheckSession:
    def test_no_session_returns_config_menu(self, view: DatabasesView) -> None:
        with patch(
            "grimperium.cli.views.databases_view.load_session", return_value=None
        ):
            assert view._state_check_session() == "config_menu"

    def test_stale_session_deletes_and_returns_config_menu(
        self, view: DatabasesView
    ) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(DatabasesView, "_server_is_responding", return_value=False),
            patch("grimperium.cli.views.databases_view.delete_session") as mock_delete,
        ):
            result = view._state_check_session()
        assert result == "config_menu"
        mock_delete.assert_called_once()

    def test_join_starts_local_worker_and_returns_monitoring(
        self, view: DatabasesView
    ) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu",
                return_value="join",
            ),
            patch.object(DatabasesView, "_start_local_worker") as mock_start,
        ):
            result = view._state_check_session()
        assert result == "monitoring"
        mock_start.assert_called_once_with("http://localhost:8000")

    def test_monitor_returns_monitoring(self, view: DatabasesView) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu",
                return_value="monitor",
            ),
        ):
            assert view._state_check_session() == "monitoring"

    def test_new_shuts_down_deletes_and_returns_config_menu(
        self, view: DatabasesView
    ) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu", return_value="new"
            ),
            patch.object(DatabasesView, "_shutdown_all_workers") as mock_shutdown,
            patch("grimperium.cli.views.databases_view.delete_session") as mock_delete,
        ):
            result = view._state_check_session()
        assert result == "config_menu"
        mock_shutdown.assert_called_once()
        mock_delete.assert_called_once()

    def test_end_shuts_down_deletes_and_returns_none(self, view: DatabasesView) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu", return_value="end"
            ),
            patch.object(DatabasesView, "_shutdown_all_workers") as mock_shutdown,
            patch("grimperium.cli.views.databases_view.delete_session") as mock_delete,
        ):
            result = view._state_check_session()
        assert result is None
        mock_shutdown.assert_called_once()
        mock_delete.assert_called_once()

    def test_cancel_returns_none(self, view: DatabasesView) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu", return_value=None
            ),
        ):
            assert view._state_check_session() is None


# ── _state_config_menu ────────────────────────────────────────────────────────


class TestStateConfigMenu:
    def test_back_returns_none(self, view: DatabasesView) -> None:
        with (
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch.object(DatabasesView, "_fetch_worker_status", return_value=[]),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu", return_value=None
            ),
        ):
            assert view._state_config_menu() is None

    def test_run_with_no_workers_warns_and_continues(self, view: DatabasesView) -> None:
        with (
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch.object(DatabasesView, "_fetch_worker_status", return_value=[]),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu",
                side_effect=["run", None],
            ),
        ):
            assert view._state_config_menu() is None

    def test_run_with_workers_configures_dispatch_and_returns_monitoring(
        self, view: DatabasesView
    ) -> None:
        workers = [_worker("w1", "lab-01")]
        from grimperium.cli.settings_manager import DistributedDefaults

        with (
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch.object(DatabasesView, "_fetch_worker_status", return_value=workers),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu", return_value="run"
            ),
            patch.object(DatabasesView, "_configure_worker") as mock_cfg,
            patch("grimperium.cli.views.databases_view.save_session"),
            patch(
                "grimperium.cli.views.databases_view.SettingsManager.load_distributed_defaults",
                return_value=DistributedDefaults(),
            ),
            patch("httpx.post"),
        ):
            result = view._state_config_menu()

        assert result == "monitoring"
        mock_cfg.assert_called_once_with(
            "http://localhost:8000", "w1", DistributedDefaults()
        )

    def test_refresh_then_back_fetches_twice(self, view: DatabasesView) -> None:
        with (
            patch.object(DatabasesView, "_server_is_responding", return_value=True),
            patch.object(
                DatabasesView, "_fetch_worker_status", return_value=[]
            ) as mock_fetch,
            patch(
                "grimperium.cli.views.databases_view.show_back_menu",
                side_effect=["refresh", None],
            ),
        ):
            view._state_config_menu()
        assert mock_fetch.call_count == 2

    def test_starts_server_when_not_responding(self, view: DatabasesView) -> None:
        with (
            patch.object(
                DatabasesView,
                "_server_is_responding",
                side_effect=[False, True, True, True],
            ),
            patch.object(DatabasesView, "_start_server_in_background") as mock_start,
            patch.object(DatabasesView, "_fetch_worker_status", return_value=[]),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu", return_value=None
            ),
        ):
            view._state_config_menu()
        mock_start.assert_called_once()


# ── _state_monitoring ─────────────────────────────────────────────────────────


class TestStateMonitoring:
    def test_quit_returns_none(self, view: DatabasesView) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(DatabasesView, "_fetch_workers_extended", return_value=[]),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu",
                return_value="quit",
            ),
        ):
            assert view._state_monitoring() is None

    def test_cancel_returns_none(self, view: DatabasesView) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(DatabasesView, "_fetch_workers_extended", return_value=[]),
            patch(
                "grimperium.cli.views.databases_view.show_back_menu", return_value=None
            ),
        ):
            assert view._state_monitoring() is None

    def test_refresh_fetches_again(self, view: DatabasesView) -> None:
        session = MagicMock()
        session.server_url = "http://localhost:8000"
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=session
            ),
            patch.object(
                DatabasesView, "_fetch_workers_extended", return_value=[]
            ) as mock_fetch,
            patch(
                "grimperium.cli.views.databases_view.show_back_menu",
                side_effect=["refresh", "quit"],
            ),
        ):
            view._state_monitoring()
        assert mock_fetch.call_count == 2

    def test_no_session_falls_back_to_localhost(self, view: DatabasesView) -> None:
        with (
            patch(
                "grimperium.cli.views.databases_view.load_session", return_value=None
            ),
            patch.object(
                DatabasesView, "_fetch_workers_extended", return_value=[]
            ) as mock_fetch,
            patch(
                "grimperium.cli.views.databases_view.show_back_menu", return_value=None
            ),
        ):
            view._state_monitoring()
        mock_fetch.assert_called_once_with("http://localhost:8000")


# ── _render_monitoring_table ──────────────────────────────────────────────────


class TestRenderMonitoringTable:
    def test_empty_table_has_seven_columns(self, view: DatabasesView) -> None:
        table = view._render_monitoring_table([])
        assert len(table.columns) == 7

    def test_workers_shown(self, view: DatabasesView) -> None:
        workers = [_worker_extended("w1", "lab-01", processed=5, successful=4)]
        table = view._render_monitoring_table(workers)
        assert table.row_count == 1

    def test_multiple_workers(self, view: DatabasesView) -> None:
        workers = [_worker_extended("w1"), _worker_extended("w2")]
        table = view._render_monitoring_table(workers)
        assert table.row_count == 2


# ── _start_local_worker ───────────────────────────────────────────────────────


class TestStartLocalWorker:
    def test_returns_daemon_thread(self, view: DatabasesView) -> None:
        with (
            patch(
                "grimperium.cli.views.databases_view.WorkerClient"
            ) as mock_client_cls,
            patch(
                "grimperium.cli.views.databases_view.WorkerRunner"
            ) as mock_runner_cls,
        ):
            mock_client = MagicMock()
            mock_client.register.return_value = {
                "crest_timeout_minutes": 60,
                "mopac_timeout_minutes": 30,
                "batch_size": 10,
            }
            mock_client_cls.return_value = mock_client
            mock_runner = MagicMock()
            mock_runner_cls.return_value = mock_runner

            t = view._start_local_worker("http://localhost:8000")

        assert t.daemon is True
        mock_runner.run.assert_called_once()

    def test_register_failure_uses_defaults(self, view: DatabasesView) -> None:
        with (
            patch(
                "grimperium.cli.views.databases_view.WorkerClient"
            ) as mock_client_cls,
            patch("grimperium.cli.views.databases_view.WorkerRunner"),
            patch("grimperium.cli.views.databases_view.WorkerConfig") as mock_cfg_cls,
        ):
            mock_client = MagicMock()
            mock_client.register.side_effect = Exception("connection refused")
            mock_client_cls.return_value = mock_client
            mock_cfg_cls.return_value = MagicMock()

            view._start_local_worker("http://localhost:8000")

        call_kwargs = mock_cfg_cls.call_args.kwargs
        assert call_kwargs["crest_timeout_minutes"] == 60
        assert call_kwargs["mopac_timeout_minutes"] == 30
        assert call_kwargs["batch_size"] == 10


# ── _screen_check_offline_workers ─────────────────────────────────────────────


class TestScreenCheckOfflineWorkers:
    @patch.object(DatabasesView, "wait_for_enter")
    def test_no_offline_molecules_no_prompt(
        self, _mock_wait: MagicMock, view: DatabasesView, tmp_path: Path
    ) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.touch()

        mock_mgr = MagicMock()
        mock_mgr.reassign_offline_molecules.return_value = 0

        import pandas as pd

        mock_mgr.df = pd.DataFrame(
            {
                "status": ["OK", "Pending"],
                "worker_status": ["unassigned", "unassigned"],
            }
        )

        with (
            patch(
                "grimperium.cli.views.databases_view.BatchCSVManager",
                return_value=mock_mgr,
            ),
            patch("grimperium.cli.views.databases_view.menu_confirm") as mock_confirm,
        ):
            view._screen_check_offline_workers(csv_path)
            mock_confirm.assert_not_called()
            mock_mgr.reassign_offline_molecules.assert_not_called()

    @patch.object(DatabasesView, "wait_for_enter")
    def test_offline_molecules_yes_reassigns(
        self, _mock_wait: MagicMock, view: DatabasesView, tmp_path: Path
    ) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.touch()

        import pandas as pd

        mock_mgr = MagicMock()
        mock_mgr.df = pd.DataFrame(
            {
                "status": ["Assigned", "Assigned"],
                "worker_status": ["offline", "offline"],
            }
        )
        mock_mgr.reassign_offline_molecules.return_value = 2

        with (
            patch(
                "grimperium.cli.views.databases_view.BatchCSVManager",
                return_value=mock_mgr,
            ),
            patch(
                "grimperium.cli.views.databases_view.menu_confirm", return_value=True
            ),
        ):
            view._screen_check_offline_workers(csv_path)
            mock_mgr.reassign_offline_molecules.assert_called_once()

    @patch.object(DatabasesView, "wait_for_enter")
    def test_offline_molecules_no_skips_reassign(
        self, _mock_wait: MagicMock, view: DatabasesView, tmp_path: Path
    ) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.touch()

        import pandas as pd

        mock_mgr = MagicMock()
        mock_mgr.df = pd.DataFrame(
            {
                "status": ["Assigned"],
                "worker_status": ["offline"],
            }
        )
        mock_mgr.reassign_offline_molecules.return_value = 1

        with (
            patch(
                "grimperium.cli.views.databases_view.BatchCSVManager",
                return_value=mock_mgr,
            ),
            patch(
                "grimperium.cli.views.databases_view.menu_confirm", return_value=False
            ),
        ):
            view._screen_check_offline_workers(csv_path)
            mock_mgr.reassign_offline_molecules.assert_not_called()

    def test_no_worker_status_column_returns_silently(
        self, view: DatabasesView, tmp_path: Path
    ) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.touch()

        import pandas as pd

        mock_mgr = MagicMock()
        mock_mgr.df = pd.DataFrame({"status": ["OK"]})

        with patch(
            "grimperium.cli.views.databases_view.BatchCSVManager",
            return_value=mock_mgr,
        ):
            view._screen_check_offline_workers(csv_path)
            mock_mgr.reassign_offline_molecules.assert_not_called()


# ── handle_action routing ─────────────────────────────────────────────────────


class TestHandleActionRouting:
    def test_calculate_run_calls_new_menu(self, view: DatabasesView) -> None:
        view._handle_run_calculation_menu = MagicMock()  # type: ignore[method-assign]
        result = view.handle_action("calculate_run")
        view._handle_run_calculation_menu.assert_called_once()
        assert result is None
