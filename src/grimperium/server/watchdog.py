"""Asyncio watchdog for monitoring worker heartbeats and recovering stuck molecules."""

from __future__ import annotations

import asyncio
import logging
from datetime import datetime, timezone

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.server.config import ServerConfig

LOG = logging.getLogger(__name__)

# {worker_id: (hostname, last_heartbeat_utc)}
HeartbeatRegistry = dict[str, tuple[str, datetime]]


def make_heartbeat_registry() -> HeartbeatRegistry:
    return {}


def _mark_worker_offline(
    worker_id: str,
    state_manager: BatchStateManager | None,
) -> int:
    """Return a worker's molecules to pending in the operational state CSV.

    Args:
        worker_id: Worker whose heartbeat expired.
        state_manager: Operational state manager; when ``None``, the call is a
            no-op for backward compatibility with tests that do not supply one.

    Returns:
        Number of molecules reclaimed (0 when ``state_manager`` is ``None``).
    """
    if state_manager is not None:
        return state_manager.mark_worker_offline(worker_id)
    return 0


async def check_offline_workers(
    registry: HeartbeatRegistry,
    csv_manager: BatchCSVManager,
    config: ServerConfig,
    lock: asyncio.Lock,
    state_manager: BatchStateManager | None = None,
) -> list[str]:
    """Identify stale workers, mark them offline in the state CSV, and evict.

    Args:
        registry: Heartbeat registry keyed by worker_id.
        csv_manager: Legacy scientific CSV manager (retained for signature
            compatibility; operational marking now goes through state_manager).
        config: Server configuration with heartbeat threshold.
        lock: Asyncio lock protecting CSV access.
        state_manager: Operational state manager for worker-offline recovery.
            When ``None``, the offline-marking step is a no-op.

    Returns:
        List of worker_ids that were evicted from the registry.
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
            await asyncio.to_thread(_mark_worker_offline, worker_id, state_manager)
        del registry[worker_id]

    return stale


async def run_watchdog(
    csv_manager: BatchCSVManager,
    heartbeat_registry: HeartbeatRegistry,
    config: ServerConfig,
    lock: asyncio.Lock,
    state_manager: BatchStateManager | None = None,
) -> None:
    """Asyncio task: periodic heartbeat check with startup recovery.

    Args:
        csv_manager: Legacy scientific CSV manager (retained for backward
            compatibility; startup recovery now delegates to state_manager).
        heartbeat_registry: Worker heartbeat timestamps.
        config: Server configuration.
        lock: Asyncio lock protecting CSV access.
        state_manager: Operational state manager for ``batch_state.csv``.
            Startup recovery and offline-marking use this when provided.
    """
    LOG.info(
        "Watchdog starting — grace period %ds before first cycle",
        config.startup_grace_s,
    )
    await asyncio.sleep(config.startup_grace_s)

    await run_startup_recovery(csv_manager, config, lock, state_manager)

    while True:
        await asyncio.sleep(config.watchdog_interval_s)
        try:
            await check_offline_workers(
                heartbeat_registry, csv_manager, config, lock, state_manager
            )
        except Exception:
            LOG.exception("Watchdog cycle failed — continuing")


async def run_startup_recovery(
    csv_manager: BatchCSVManager,
    config: ServerConfig,
    lock: asyncio.Lock,
    state_manager: BatchStateManager | None = None,
) -> None:
    """Recover stuck molecules once before the watchdog loop starts."""
    if state_manager is not None:
        timeout_minutes = int(config.stuck_assigned_threshold_h * 60)
        async with lock:
            reset_running = await asyncio.to_thread(state_manager.reset_stuck_running)
            reset_stuck = await asyncio.to_thread(
                state_manager.reset_stuck_assigned, timeout_minutes
            )
        LOG.info(
            "Startup recovery: reset_stuck_running=%d, reset_stuck_assigned=%d",
            reset_running,
            len(reset_stuck),
        )
        return

    # Legacy fallback: only reset RUNNING in the scientific CSV.
    async with lock:
        reset_running = await asyncio.to_thread(csv_manager.reset_stuck_running)
    LOG.info("Startup recovery (legacy): reset_stuck_running=%d", reset_running)
