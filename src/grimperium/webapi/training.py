"""Training endpoints: read status, launch and stop the routine subprocess.

The API never runs the science. ``start`` spawns the same command a person
would type (``python -m semi_imperium.training run``), detached, and returns;
the dashboard then polls ``status``.
"""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from fastapi import HTTPException

from grimperium.webapi import readers
from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.models import (
    TrainingFailure,
    TrainingJob,
    TrainingProgress,
    TrainingStartRequest,
    TrainingStatus,
)

_MIGRATE_HINT = (
    "O CSV ainda não foi migrado para o protocolo 2. Rode no terminal: "
    "poetry run python -m semi_imperium.training migrate --apply"
)


@dataclass(frozen=True)
class TrainingLaunch:
    """What the launcher needs to start one run."""

    csv_path: Path
    store_root: Path
    workers: int
    limit: int | None


Launcher = Callable[[TrainingLaunch], None]


def default_launcher(launch: TrainingLaunch) -> None:
    """Start the routine detached; its output goes to ``training.log``."""
    command = [
        sys.executable,
        "-m",
        "semi_imperium.training",
        "--csv",
        str(launch.csv_path),
        "--store",
        str(launch.store_root),
        "run",
        "--workers",
        str(launch.workers),
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
    """Progress from the CSV, job and failures from SQLite."""
    from semi_imperium.training.db import database_path
    from semi_imperium.training.geometry_store import GeometryStore
    from semi_imperium.training.job_state import JobState
    from semi_imperium.training.ledger import read_progress

    try:
        progress = read_progress(paths.pm7_csv)
    except FileNotFoundError as exc:
        raise readers.DataUnavailable(paths.pm7_csv, "file not found") from exc
    db = database_path(paths.semi_imperium_store)
    state = JobState(db)
    snapshot = state.snapshot()
    return TrainingStatus(
        progress=TrainingProgress(**progress.__dict__),
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
    )


def start(
    paths: DashboardPaths, body: TrainingStartRequest, launcher: Launcher
) -> TrainingStatus:
    """Launch the routine unless one is alive or the CSV is not migrated."""
    status = build_status(paths)
    if status.job.state in ("running", "stopping"):
        raise HTTPException(409, "Já existe um treinamento em execução.")
    if not status.progress.migrated:
        raise HTTPException(409, _MIGRATE_HINT)
    launcher(
        TrainingLaunch(
            csv_path=paths.pm7_csv,
            store_root=paths.semi_imperium_store,
            workers=body.workers,
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
