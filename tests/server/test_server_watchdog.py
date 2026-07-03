"""Tests for server/watchdog.py — asyncio watchdog task logic."""

import asyncio
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus, WorkerStatus
from grimperium.crest_pm7.batch.output_contracts import BATCH_STATE_COLUMNS
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import PM7Config
from grimperium.server.config import ServerConfig
from grimperium.server.watchdog import (
    HeartbeatRegistry,
    check_offline_workers,
    make_heartbeat_registry,
    run_startup_recovery,
)

# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture
def config() -> ServerConfig:
    return ServerConfig(
        csv_path="/tmp/test.csv",
        heartbeat_timeout_s=300,
        watchdog_interval_s=30,
        startup_grace_s=0,  # no grace in tests
        stuck_assigned_threshold_h=2.0,
    )


@pytest.fixture
def minimal_csv(tmp_path: Path) -> Path:
    p = tmp_path / "test.csv"
    p.write_text("mol_id,smiles,nheavy,status\nmol_001,CCO,3,Pending\n")
    return p


@pytest.fixture
def manager(minimal_csv: Path) -> BatchCSVManager:
    m = BatchCSVManager(minimal_csv)
    m.load_csv()
    return m


# ── HeartbeatRegistry ─────────────────────────────────────────────────────────


class TestMakeHeartbeatRegistry:
    def test_returns_empty_dict(self) -> None:
        reg = make_heartbeat_registry()
        assert reg == {}

    def test_type_is_dict(self) -> None:
        reg = make_heartbeat_registry()
        assert isinstance(reg, dict)


# ── check_offline_workers ─────────────────────────────────────────────────────


class TestCheckOfflineWorkers:
    @pytest.mark.asyncio
    async def test_no_workers_no_action(
        self, manager: BatchCSVManager, config: ServerConfig
    ) -> None:
        registry: HeartbeatRegistry = {}
        lock = asyncio.Lock()
        removed = await check_offline_workers(registry, manager, config, lock)
        assert removed == []

    @pytest.mark.asyncio
    async def test_fresh_heartbeat_not_removed(
        self, manager: BatchCSVManager, config: ServerConfig
    ) -> None:
        now = datetime.now(timezone.utc)
        registry: HeartbeatRegistry = {"w1": ("lab-01", now)}
        lock = asyncio.Lock()
        removed = await check_offline_workers(registry, manager, config, lock)
        assert removed == []
        assert "w1" in registry

    @pytest.mark.asyncio
    async def test_stale_heartbeat_removed_from_registry(
        self, manager: BatchCSVManager, config: ServerConfig
    ) -> None:
        old_time = datetime.now(timezone.utc) - timedelta(seconds=400)
        registry: HeartbeatRegistry = {"w1": ("lab-01", old_time)}
        lock = asyncio.Lock()

        with patch(
            "grimperium.server.watchdog.asyncio.to_thread",
            new_callable=AsyncMock,
            return_value=0,
        ):
            removed = await check_offline_workers(registry, manager, config, lock)

        assert "w1" in removed
        assert "w1" not in registry

    @pytest.mark.asyncio
    async def test_stale_worker_calls_mark_offline(
        self, manager: BatchCSVManager, config: ServerConfig
    ) -> None:
        old_time = datetime.now(timezone.utc) - timedelta(seconds=400)
        registry: HeartbeatRegistry = {"w1": ("lab-01", old_time)}
        lock = asyncio.Lock()

        called_with: list[str] = []

        async def fake_to_thread(fn: object, *args: object, **kwargs: object) -> int:
            if callable(fn):
                called_with.append(str(args[0]) if args else "no-arg")
            return 0

        with patch(
            "grimperium.server.watchdog.asyncio.to_thread", side_effect=fake_to_thread
        ):
            await check_offline_workers(registry, manager, config, lock)

        assert "w1" in called_with

    @pytest.mark.asyncio
    async def test_only_stale_workers_removed(
        self, manager: BatchCSVManager, config: ServerConfig
    ) -> None:
        now = datetime.now(timezone.utc)
        old_time = now - timedelta(seconds=400)
        registry: HeartbeatRegistry = {
            "fresh": ("lab-01", now),
            "stale": ("lab-02", old_time),
        }
        lock = asyncio.Lock()

        with patch(
            "grimperium.server.watchdog.asyncio.to_thread",
            new_callable=AsyncMock,
            return_value=0,
        ):
            removed = await check_offline_workers(registry, manager, config, lock)

        assert "stale" in removed
        assert "fresh" not in removed
        assert "fresh" in registry

    @pytest.mark.asyncio
    async def test_multiple_stale_workers_all_removed(
        self, manager: BatchCSVManager, config: ServerConfig
    ) -> None:
        old_time = datetime.now(timezone.utc) - timedelta(seconds=400)
        registry: HeartbeatRegistry = {
            "w1": ("lab-01", old_time),
            "w2": ("lab-02", old_time),
        }
        lock = asyncio.Lock()

        with patch(
            "grimperium.server.watchdog.asyncio.to_thread",
            new_callable=AsyncMock,
            return_value=0,
        ):
            removed = await check_offline_workers(registry, manager, config, lock)

        assert set(removed) == {"w1", "w2"}
        assert registry == {}


# ── run_watchdog startup recovery regression tests ────────────────────────────


def _state_csv(tmp_path: Path, rows: list[dict[str, object]]) -> Path:
    p = tmp_path / "batch_state.csv"
    pd.DataFrame(rows, columns=BATCH_STATE_COLUMNS).to_csv(p, index=False)
    return p


def _base_state_row(
    mol_id: str,
    status: str = MoleculeStatus.PENDING.value,
    worker: str = "",
) -> dict[str, object]:
    row: dict[str, object] = dict.fromkeys(BATCH_STATE_COLUMNS, "")
    row.update(
        {
            "mol_id": mol_id,
            "status": status,
            "smiles": "C",
            "reruns": 0,
            "worker_status": WorkerStatus.UNASSIGNED.value,
            "assigned_worker": worker,
        }
    )
    return row


async def _run_without_thread(fn: object, *args: object, **kwargs: object) -> object:
    if not callable(fn):
        raise TypeError("expected callable")
    return fn(*args, **kwargs)


class TestRunWatchdogStartupRecovery:
    """Regression tests for the startup recovery bug (broken csv_manager calls).

    Before the fix, run_watchdog called csv_manager.reset_stuck_assigned(...),
    which had been removed from BatchCSVManager; this raised AttributeError on
    every server startup.  Now startup recovery goes through BatchStateManager.
    """

    @pytest.mark.asyncio
    async def test_startup_recovery_calls_state_manager_reset_stuck_running(
        self, tmp_path: Path, config: ServerConfig, minimal_csv: Path
    ) -> None:
        state_path = _state_csv(
            tmp_path,
            [_base_state_row("mol-001", MoleculeStatus.RUNNING.value, "w1")],
        )
        sm = BatchStateManager(state_path, PM7Config())
        csv_mgr = BatchCSVManager(minimal_csv)
        csv_mgr.load_csv()
        lock = asyncio.Lock()

        with patch(
            "grimperium.server.watchdog.asyncio.to_thread",
            side_effect=_run_without_thread,
        ):
            await run_startup_recovery(csv_mgr, config, lock, sm)

        rows = pd.read_csv(state_path, keep_default_na=False)
        # RUNNING molecule must have been reset to PENDING by recovery.
        assert rows.loc[rows["mol_id"] == "mol-001", "status"].values[0] == (
            MoleculeStatus.PENDING.value
        )

    @pytest.mark.asyncio
    async def test_startup_recovery_calls_state_manager_reset_stuck_assigned(
        self, tmp_path: Path, config: ServerConfig, minimal_csv: Path
    ) -> None:
        """reset_stuck_assigned on csv_manager was removed; must use state_manager."""
        state_path = _state_csv(
            tmp_path,
            [
                {
                    **_base_state_row("mol-001", MoleculeStatus.ASSIGNED.value, "w1"),
                    # assigned_at far in the past so it exceeds the threshold
                    "assigned_at": "2000-01-01T00:00:00+00:00",
                    "worker_status": WorkerStatus.ONLINE.value,
                }
            ],
        )
        sm = BatchStateManager(state_path, PM7Config())
        csv_mgr = BatchCSVManager(minimal_csv)
        csv_mgr.load_csv()
        lock = asyncio.Lock()

        with patch(
            "grimperium.server.watchdog.asyncio.to_thread",
            side_effect=_run_without_thread,
        ):
            await run_startup_recovery(csv_mgr, config, lock, sm)

        rows = pd.read_csv(state_path, keep_default_na=False)
        assert rows.loc[rows["mol_id"] == "mol-001", "status"].values[0] == (
            MoleculeStatus.PENDING.value
        )

    @pytest.mark.asyncio
    async def test_check_offline_workers_marks_offline_via_state_manager(
        self, tmp_path: Path, config: ServerConfig, minimal_csv: Path
    ) -> None:
        state_path = _state_csv(
            tmp_path,
            [
                {
                    **_base_state_row("mol-001", MoleculeStatus.RUNNING.value, "w1"),
                    "worker_status": WorkerStatus.ONLINE.value,
                }
            ],
        )
        sm = BatchStateManager(state_path, PM7Config())
        csv_mgr = BatchCSVManager(minimal_csv)
        csv_mgr.load_csv()

        old_time = datetime.now(timezone.utc) - timedelta(seconds=400)
        registry: HeartbeatRegistry = {"w1": ("lab-01", old_time)}
        lock = asyncio.Lock()

        with patch(
            "grimperium.server.watchdog.asyncio.to_thread",
            side_effect=_run_without_thread,
        ):
            removed = await check_offline_workers(registry, csv_mgr, config, lock, sm)

        assert "w1" in removed
        rows = pd.read_csv(state_path, keep_default_na=False)
        # Molecule must have been reclaimed by mark_worker_offline.
        assert rows.loc[rows["mol_id"] == "mol-001", "status"].values[0] == (
            MoleculeStatus.PENDING.value
        )
