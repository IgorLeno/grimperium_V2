"""Training endpoints: read status, launch and stop the routine subprocess.

The API never runs the science. ``start`` spawns the same command a person
would type (``python -m semi_imperium.training run``), detached, and returns;
the dashboard then polls ``status``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from fastapi import HTTPException

from grimperium.webapi import readers
from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.models import (
    TrainingActivity,
    TrainingCompleted,
    TrainingFailure,
    TrainingInFlight,
    TrainingJob,
    TrainingStartRequest,
    TrainingStatus,
)


@dataclass(frozen=True)
class TrainingLaunch:
    """What the launcher needs to start one run."""

    cbs_csv: Path
    store_root: Path
    workers: int
    crest_threads: int
    limit: int | None


Launcher = Callable[[TrainingLaunch], None]


def default_launcher(launch: TrainingLaunch) -> None:
    """Start the routine detached; its output goes to ``training.log``."""
    command = [
        sys.executable,
        "-m",
        "semi_imperium.training",
        "--csv",
        str(launch.cbs_csv),
        "--store",
        str(launch.store_root),
        "run",
        "--workers",
        str(launch.workers),
        "--crest-threads",
        str(launch.crest_threads),
    ]
    if launch.limit is not None:
        command += ["--limit", str(launch.limit)]
    launch.store_root.mkdir(parents=True, exist_ok=True)
    log = (launch.store_root / "training.log").open("ab")
    subprocess.Popen(  # noqa: S603 - fixed argv, no shell, local paths only
        command,
        stdin=subprocess.DEVNULL,
        stdout=log,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )


def build_status(paths: DashboardPaths) -> TrainingStatus:
    """Queue progress, job and failures from SQLite (the CBS only seeds the queue)."""
    from semi_imperium.training.db import database_path
    from semi_imperium.training.geometry_store import GeometryStore
    from semi_imperium.training.job_state import JobState
    from semi_imperium.training.molecule_queue import read_activity

    progress = readers.queue_progress(paths.semi_imperium_store, paths.cbs_csv)
    db = database_path(paths.semi_imperium_store)
    state = JobState(db)
    snapshot = state.snapshot()
    activity = read_activity(db, since=snapshot.started_at)
    return TrainingStatus(
        progress=progress,
        job=TrainingJob(
            state=snapshot.state,
            pid=snapshot.pid,
            started_at=snapshot.started_at,
            heartbeat_at=snapshot.heartbeat_at,
            finished_at=snapshot.finished_at,
            workers=snapshot.workers,
            limit=snapshot.limit,
            done=snapshot.done,
            failed=snapshot.failed,
            message=snapshot.message,
            stop_requested=snapshot.stop_requested,
        ),
        geometries=GeometryStore(db).count(),
        failures=[
            TrainingFailure(mol_id=mol, reason=reason, at=at)
            for mol, reason, at in state.recent_failures()
        ],
        activity=TrainingActivity(
            in_flight=[TrainingInFlight(**vars(m)) for m in activity.in_flight],
            recent=[TrainingCompleted(**vars(m)) for m in activity.recent],
            mean_elapsed_s=activity.mean_elapsed_s,
        ),
        cpu_count=os.cpu_count() or 1,
    )


def start(
    paths: DashboardPaths, body: TrainingStartRequest, launcher: Launcher
) -> TrainingStatus:
    """Launch the routine unless one is alive."""
    status = build_status(paths)
    if status.job.state in ("running", "stopping"):
        raise HTTPException(409, "Já existe um treinamento em execução.")
    launcher(
        TrainingLaunch(
            cbs_csv=paths.cbs_csv,
            store_root=paths.semi_imperium_store,
            workers=body.workers,
            crest_threads=body.crest_threads,
            limit=body.limit,
        )
    )
    return build_status(paths)


def stop(paths: DashboardPaths) -> TrainingStatus:
    """Ask the running routine to stop dispatching new molecules."""
    from semi_imperium.training.db import database_path
    from semi_imperium.training.job_state import JobState

    if not JobState(database_path(paths.semi_imperium_store)).request_stop():
        raise HTTPException(409, "Nenhum treinamento em execução.")
    return build_status(paths)
