"""Tests for DatabasesView distributed mode — Phase 4 CLI integration."""

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
        # Rich Table is returned with 4 columns
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


# ── _screen_connect_workers ───────────────────────────────────────────────────


class TestScreenConnectWorkers:
    @patch("grimperium.cli.views.databases_view.show_back_menu", return_value="proceed")
    @patch.object(DatabasesView, "_fetch_worker_status", return_value=[])
    def test_proceed_returns_true(
        self,
        _mock_fetch: MagicMock,
        _mock_menu: MagicMock,
        view: DatabasesView,
    ) -> None:
        assert view._screen_connect_workers("http://localhost:8000") is True

    @patch("grimperium.cli.views.databases_view.show_back_menu", return_value=None)
    @patch.object(DatabasesView, "_fetch_worker_status", return_value=[])
    def test_cancel_returns_false(
        self,
        _mock_fetch: MagicMock,
        _mock_menu: MagicMock,
        view: DatabasesView,
    ) -> None:
        assert view._screen_connect_workers("http://localhost:8000") is False

    @patch(
        "grimperium.cli.views.databases_view.show_back_menu",
        side_effect=["refresh", "proceed"],
    )
    @patch.object(DatabasesView, "_fetch_worker_status", return_value=[])
    def test_refresh_loops_and_proceeds(
        self,
        _mock_fetch: MagicMock,
        mock_menu: MagicMock,
        view: DatabasesView,
    ) -> None:
        assert view._screen_connect_workers("http://localhost:8000") is True
        assert mock_menu.call_count == 2

    @patch("grimperium.cli.views.databases_view.show_back_menu", return_value="proceed")
    @patch.object(
        DatabasesView,
        "_fetch_worker_status",
        return_value=[_worker("w1", "lab-01")],
    )
    def test_shows_connected_workers(
        self,
        _mock_fetch: MagicMock,
        _mock_menu: MagicMock,
        view: DatabasesView,
    ) -> None:
        result = view._screen_connect_workers("http://localhost:8000")
        assert result is True
        _mock_fetch.assert_called_once()


# ── _screen_assign_workload ───────────────────────────────────────────────────


class TestScreenAssignWorkload:
    def _setup_csv_manager(
        self,
        pending: int = 5,
        rerun: int = 2,
    ) -> MagicMock:
        mgr = MagicMock()
        mgr.get_status_counts.return_value = {"Pending": pending, "Rerun": rerun}
        mgr.distribute_molecules.return_value = {"central": ["m1"]}
        return mgr

    @patch.object(DatabasesView, "_fetch_worker_status", return_value=[])
    @patch.object(DatabasesView, "_prompt_positive_int", return_value=3)
    @patch.object(DatabasesView, "wait_for_enter")
    def test_calls_distribute_molecules_with_correct_args(
        self,
        _mock_wait: MagicMock,
        _mock_prompt: MagicMock,
        _mock_fetch: MagicMock,
        view: DatabasesView,
        tmp_path: Path,
    ) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.touch()

        mock_mgr = self._setup_csv_manager()
        with patch(
            "grimperium.cli.views.databases_view.BatchCSVManager",
            return_value=mock_mgr,
        ):
            result = view._screen_assign_workload(
                "http://localhost:8000",
                crest_timeout=30,
                mopac_timeout=60,
                csv_path=csv_path,
            )

        assert result is True
        mock_mgr.distribute_molecules.assert_called_once()
        call_kwargs = mock_mgr.distribute_molecules.call_args
        assert call_kwargs.kwargs["crest_timeout_minutes"] == 30
        assert call_kwargs.kwargs["mopac_timeout_minutes"] == 60

    @patch.object(DatabasesView, "_fetch_worker_status", return_value=[_worker("w1")])
    @patch.object(DatabasesView, "_prompt_positive_int", return_value=2)
    @patch.object(DatabasesView, "wait_for_enter")
    def test_includes_central_and_remote_workers(
        self,
        _mock_wait: MagicMock,
        _mock_prompt: MagicMock,
        _mock_fetch: MagicMock,
        view: DatabasesView,
        tmp_path: Path,
    ) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.touch()

        mock_mgr = self._setup_csv_manager()
        with patch(
            "grimperium.cli.views.databases_view.BatchCSVManager",
            return_value=mock_mgr,
        ):
            view._screen_assign_workload(
                "http://localhost:8000",
                crest_timeout=30,
                mopac_timeout=60,
                csv_path=csv_path,
            )

        allocations = mock_mgr.distribute_molecules.call_args.kwargs["allocations"]
        assert "central" in allocations
        assert "w1" in allocations

    @patch.object(DatabasesView, "_fetch_worker_status", return_value=[])
    @patch.object(DatabasesView, "_prompt_positive_int", return_value=None)
    def test_cancel_on_prompt_returns_false(
        self,
        _mock_prompt: MagicMock,
        _mock_fetch: MagicMock,
        view: DatabasesView,
        tmp_path: Path,
    ) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.touch()

        mock_mgr = self._setup_csv_manager()
        with patch(
            "grimperium.cli.views.databases_view.BatchCSVManager",
            return_value=mock_mgr,
        ):
            result = view._screen_assign_workload(
                "http://localhost:8000",
                crest_timeout=30,
                mopac_timeout=60,
                csv_path=csv_path,
            )

        assert result is False
        mock_mgr.distribute_molecules.assert_not_called()

    @patch.object(DatabasesView, "_fetch_worker_status", return_value=[])
    @patch.object(DatabasesView, "_prompt_positive_int", return_value=0)
    @patch.object(DatabasesView, "wait_for_enter")
    def test_all_zero_allocations_returns_false(
        self,
        _mock_wait: MagicMock,
        _mock_prompt: MagicMock,
        _mock_fetch: MagicMock,
        view: DatabasesView,
        tmp_path: Path,
    ) -> None:
        csv_path = tmp_path / "test.csv"
        csv_path.touch()

        mock_mgr = self._setup_csv_manager()
        with patch(
            "grimperium.cli.views.databases_view.BatchCSVManager",
            return_value=mock_mgr,
        ):
            result = view._screen_assign_workload(
                "http://localhost:8000",
                crest_timeout=30,
                mopac_timeout=60,
                csv_path=csv_path,
            )

        assert result is False
        mock_mgr.distribute_molecules.assert_not_called()


# ── _screen_check_offline_workers ─────────────────────────────────────────────


class TestScreenCheckOfflineWorkers:
    def _make_mgr(self, offline_count: int) -> MagicMock:
        mgr = MagicMock()
        mgr.df = MagicMock()

        import pandas as pd

        if offline_count > 0:
            mgr.df.__getitem__ = MagicMock(
                side_effect=lambda key: (
                    pd.Series(["Assigned"] * offline_count)
                    if key == "status"
                    else pd.Series(["offline"] * offline_count)
                )
            )
            # Make boolean mask work
            mock_df = MagicMock()
            mock_df.columns = ["status", "worker_status"]
            mock_df.__getitem__ = MagicMock(
                side_effect=lambda k: pd.Series(["Assigned"] * offline_count)
            )
            mgr.df = mock_df
        else:
            mock_df = MagicMock()
            mock_df.columns = ["status", "worker_status"]
            mgr.df = mock_df

        mgr.reassign_offline_molecules.return_value = offline_count
        return mgr

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
