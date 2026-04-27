"""Asyncio watchdog for monitoring worker heartbeats and recovering stuck molecules."""

import asyncio
import logging
from datetime import datetime, timezone

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.server.config import ServerConfig

LOG = logging.getLogger(__name__)

# {worker_id: (hostname, last_heartbeat_utc)}
HeartbeatRegistry = dict[str, tuple[str, datetime]]


def make_heartbeat_registry() -> HeartbeatRegistry:
    return {}


async def check_offline_workers(
    registry: HeartbeatRegistry,
    csv_manager: BatchCSVManager,
    config: ServerConfig,
    lock: asyncio.Lock,
) -> list[str]:
    """Identify stale workers, mark them offline in CSV, and remove from registry.

    Returns list of worker_ids that were evicted.
    """
    now = datetime.now(timezone.utc)
    stale = [
        worker_id
        for worker_id, (_, last_seen) in registry.items()
        if (now - last_seen).total_seconds() > config.heartbeat_timeout_s
    ]

    for worker_id in stale:
        LOG.warning(
            "Worker %r heartbeat expired — marking offline and removing from registry",
            worker_id,
        )
        async with lock:
            await asyncio.to_thread(csv_manager.mark_worker_offline, worker_id)
        del registry[worker_id]

    return stale


async def run_watchdog(
    csv_manager: BatchCSVManager,
    heartbeat_registry: HeartbeatRegistry,
    config: ServerConfig,
    lock: asyncio.Lock,
) -> None:
    """Asyncio task: periodic heartbeat check with startup recovery."""
    LOG.info(
        "Watchdog starting — grace period %ds before first cycle",
        config.startup_grace_s,
    )
    await asyncio.sleep(config.startup_grace_s)

    # First cycle: run startup recovery
    async with lock:
        reset_running = await asyncio.to_thread(csv_manager.reset_stuck_running)
        auto_reset, warned = await asyncio.to_thread(
            csv_manager.reset_stuck_assigned, config.stuck_assigned_threshold_h
        )
    LOG.info(
        "Startup recovery: reset_stuck_running=%d, reset_stuck_assigned=%d, warned=%d",
        reset_running,
        auto_reset,
        warned,
    )

    while True:
        await asyncio.sleep(config.watchdog_interval_s)
        try:
            await check_offline_workers(heartbeat_registry, csv_manager, config, lock)
        except Exception:
            LOG.exception("Watchdog cycle failed — continuing")
