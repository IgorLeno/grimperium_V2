"""Tests for server/watchdog.py — asyncio watchdog task logic."""

import asyncio
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.server.config import ServerConfig
from grimperium.server.watchdog import (
    HeartbeatRegistry,
    check_offline_workers,
    make_heartbeat_registry,
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
