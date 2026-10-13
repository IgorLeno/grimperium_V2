"""Coordinator of the training routine: the single writer.

It reads the queue, dispatches one molecule per free worker, and records each
finished molecule: first the geometry (SQLite), then the CSV row, so a row
marked ``OK`` always has its geometry. Stopping is cooperative: no new
molecule is dispatched, in-flight ones finish, untouched rows stay ``Pending``.
"""

from __future__ import annotations

import multiprocessing
import os
import signal
import threading
import uuid
from collections import deque
from collections.abc import Callable
from concurrent.futures import (
    FIRST_COMPLETED,
    Executor,
    Future,
    ProcessPoolExecutor,
    wait,
)
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from semi_imperium.domain import RunRecord, RunState
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.training.db import database_path
from semi_imperium.training.geometry_store import GeometryRecord, GeometryStore
from semi_imperium.training.job_state import JobState
from semi_imperium.training.ledger import QueuedMolecule, TrainingLedger
from semi_imperium.training.protocol import (
    CANONICAL_PROTOCOL,
    HAMILTONIAN,
    canonical_configuration,
    canonical_settings,
)
from semi_imperium.training.worker import WorkerResult, WorkerTask, run_task
from semi_imperium.workflows.calculation import build_provenance

PoolFactory = Callable[[int], Executor]
TaskFn = Callable[[WorkerTask], WorkerResult]


@dataclass(frozen=True)
class RunSettings:
    """Inputs of one training run (execution knobs, not scientific choices)."""

    csv_path: Path
    store_root: Path
    work_dir: Path
    workers: int = 3
    limit: int | None = None
    crest_executable: str = "crest"
    xtb_executable: str = "xtb"
    mopac_executable: str = "mopac"
    crest_threads: int = 4
    crest_timeout_s: float = 3600.0
    mopac_timeout_s: float = 900.0
    max_reruns: int = 3
    keep_artifacts: bool = False

    def __post_init__(self) -> None:
        if self.workers < 1:
            raise ValueError(f"RunSettings.workers must be >= 1, got {self.workers}")
        if self.limit is not None and self.limit < 1:
            raise ValueError(f"RunSettings.limit must be >= 1, got {self.limit}")


@dataclass(frozen=True)
class RunSummary:
    """Outcome of one run."""

    attempted: int
    ok: int
    failed: int
    stopped: bool


def _spawn_pool(workers: int) -> Executor:
    return ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context("spawn")
    )


def _run_id() -> str:
    return f"run-{datetime.now().strftime('%Y%m%dT%H%M%S')}-{uuid.uuid4().hex[:6]}"


def run_training(
    settings: RunSettings,
    *,
    pool_factory: PoolFactory = _spawn_pool,
    task_fn: TaskFn = run_task,
    install_signals: bool = True,
) -> RunSummary:
    """Run the queue to exhaustion (or to ``limit``, or until stopped)."""
    ledger = TrainingLedger(settings.csv_path, max_reruns=settings.max_reruns)
    ledger.load()  # raises NotMigratedError before any state is touched
    queue = ledger.queue(settings.limit)

    db_path = database_path(settings.store_root)
    job = JobState(db_path)
    geometries = GeometryStore(db_path)
    job.begin(pid=os.getpid(), workers=settings.workers, limit=settings.limit)

    stop = threading.Event()
    if install_signals:
        signal.signal(signal.SIGTERM, lambda *_: stop.set())
        signal.signal(signal.SIGINT, lambda *_: stop.set())

    run_id = _run_id()
    batch_id = f"protocol{CANONICAL_PROTOCOL}-{run_id}"
    sem_settings = canonical_settings()
    store = SemiImperiumStore(settings.store_root)
    run = RunRecord(
        run_id=run_id,
        configuration=canonical_configuration(sem_settings),
        provenance=build_provenance(sem_settings),
        label="delta-entalpia-treinamento",
    )
    store.save_run(run)
    running = run.transition_to(RunState.RUNNING)
    store.save_run(running)

    ok = failed = 0
    attempted = 0
    stopped = False
    try:
        pending = deque(queue)
        in_flight: dict[Future[WorkerResult], QueuedMolecule] = {}
        with pool_factory(settings.workers) as pool:
            while pending or in_flight:
                halt = stop.is_set() or job.stop_requested()
                if halt and not stopped:
                    stopped = True
                    job.mark_stopping()
                while pending and len(in_flight) < settings.workers and not halt:
                    molecule = pending.popleft()
                    ledger.mark_running(molecule.mol_id)
                    in_flight[
                        pool.submit(task_fn, _task(settings, run_id, molecule))
                    ] = molecule
                    attempted += 1
                if not in_flight:
                    break
                done, _ = wait(in_flight, timeout=30.0, return_when=FIRST_COMPLETED)
                for future in done:
                    molecule = in_flight.pop(future)
                    if _record(ledger, geometries, job, molecule, future, batch_id):
                        ok += 1
                    else:
                        failed += 1
                job.heartbeat(done=ok, failed=failed)
        run_state = (
            RunState.COMPLETED
            if failed == 0
            else RunState.FAILED if ok == 0 and attempted > 0 else RunState.PARTIAL
        )
        store.save_run(
            running.transition_to(
                run_state,
                error_message=(
                    "every molecule failed" if run_state is RunState.FAILED else None
                ),
            )
        )
        job.finish("finished", "parado a pedido" if stopped else "concluído")
    except BaseException as exc:
        job.finish("failed", f"{type(exc).__name__}: {exc}")
        raise
    return RunSummary(attempted=attempted, ok=ok, failed=failed, stopped=stopped)


def _task(settings: RunSettings, run_id: str, molecule: QueuedMolecule) -> WorkerTask:
    return WorkerTask(
        mol_id=molecule.mol_id,
        smiles=molecule.smiles,
        charge=molecule.charge,
        multiplicity=molecule.multiplicity,
        run_id=run_id,
        store_root=str(settings.store_root),
        work_dir=str(settings.work_dir),
        crest_executable=settings.crest_executable,
        xtb_executable=settings.xtb_executable,
        mopac_executable=settings.mopac_executable,
        crest_threads=settings.crest_threads,
        crest_timeout_s=settings.crest_timeout_s,
        mopac_timeout_s=settings.mopac_timeout_s,
        keep_artifacts=settings.keep_artifacts,
    )


def _record(
    ledger: TrainingLedger,
    geometries: GeometryStore,
    job: JobState,
    molecule: QueuedMolecule,
    future: Future[WorkerResult],
    batch_id: str,
) -> bool:
    """Persist one finished molecule; return whether it was accepted."""
    try:
        result = future.result()
    except Exception as exc:  # worker process died or was killed
        reason = f"{type(exc).__name__}: {exc}"
        ledger.record_failure(molecule.mol_id, reason)
        job.record_failure(molecule.mol_id, reason)
        return False
    if not result.ok or result.hof_kcal_mol is None:
        ledger.record_failure(molecule.mol_id, result.error or result.state)
        job.record_failure(molecule.mol_id, result.error or result.state)
        return False
    geometries.save(
        GeometryRecord(
            mol_id=result.mol_id,
            smiles=result.smiles,
            hamiltonian=HAMILTONIAN,
            protocol=CANONICAL_PROTOCOL,
            signature=result.signature,
            calculation_id=result.calculation_id,
            conformer_index=result.conformer_index,
            elements=result.elements,
            coordinates=result.coordinates,
            hof_kcal_mol=result.hof_kcal_mol,
            lowest_frequency_cm1=result.lowest_frequency_cm1,
            minimum_confirmed=True,
        )
    )
    ledger.record_success(
        result.mol_id,
        hof_kcal_mol=result.hof_kcal_mol,
        conformers_evaluated=result.conformers_evaluated,
        elapsed_s=result.elapsed_s,
        batch_id=batch_id,
    )
    return True
