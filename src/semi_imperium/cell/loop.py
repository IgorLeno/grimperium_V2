"""The cell's loop: register, claim, compute, send, with one heartbeat for all.

A cell is a remote ``central``: it claims one molecule per free worker,
computes it with the same ``run_task`` the central uses, and reports each
result by ``attempt_id``. What it adds is distance:

* every result goes to the on-disk outbox first, and is resent with a
  growing wait while the central cannot be reached (``outbox``);
* a start is a restart: the work folder is wiped (an interrupted CREST
  starts over), finished results still on disk are sent *before*
  registering, so the central does not hand those molecules to someone
  else, and the registration itself returns whatever was still held;
* one heartbeat every ``heartbeat_s`` names the attempt of every busy
  worker, which keeps their leases and the cell ``online``.

Stopping is cooperative: the first SIGTERM/SIGINT stops claiming and lets
the running molecules finish and be sent; the second aborts at once.
"""

from __future__ import annotations

import json
import logging
import math
import multiprocessing
import os
import shutil
import signal
import socket
import threading
import time
import uuid
from collections.abc import Callable, Mapping
from concurrent.futures import (
    FIRST_COMPLETED,
    Executor,
    Future,
    ProcessPoolExecutor,
    wait,
)
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from semi_imperium.cell.client import CentralClient
from semi_imperium.cell.config import CellConfig, cpu_cores, worker_count
from semi_imperium.cell.outbox import (
    Backoff,
    Outbox,
    Reply,
    Verdict,
    atomic_write_json,
    classify,
    flush,
)
from semi_imperium.domain import RunRecord, RunState
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.training.program_info import ProgramProbe, probe_programs
from semi_imperium.training.protocol import (
    canonical_configuration,
    canonical_settings,
    canonical_signature,
)
from semi_imperium.training.worker import WorkerResult, WorkerTask, run_task
from semi_imperium.workflows.calculation import build_provenance

log = logging.getLogger(__name__)

EXIT_OK = 0
EXIT_TOKEN = 4  # token missing, wrong or revoked
EXIT_REFUSED = 5  # registration refused: programs or protocol differ

STATE_FILENAME = "state.json"
_FAILURE_MAX = 2000  # the central's limit for a failure reason

PoolFactory = Callable[[int], Executor]
TaskFn = Callable[[WorkerTask], WorkerResult]


@dataclass(frozen=True)
class Timings:
    #: Wait before asking again when the queue is empty.
    idle_s: float = 60.0
    #: Wait before asking again while the cell's admission is pending.
    admission_s: float = 60.0
    #: Override of the central's ``heartbeat_s`` (tests only).
    heartbeat_s: float | None = None
    backoff: Backoff = field(default_factory=Backoff)


def _spawn_pool(workers: int) -> Executor:
    return ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context("spawn")
    )


def _abort_process_group() -> None:
    """Second signal: stop now, CREST and MOPAC included.

    Started under ``setsid`` (the scheduled task) the cell leads its process
    group and takes every child with it; from a terminal, Ctrl-C already
    reached the whole foreground group.
    """
    if os.getpgrp() == os.getpid():
        os.killpg(os.getpgrp(), signal.SIGKILL)
    os._exit(130)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_state(root: Path) -> dict[str, Any]:
    """The last state the cell wrote (``{}`` before its first start)."""
    try:
        data = json.loads((root / STATE_FILENAME).read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    return data if isinstance(data, dict) else {}


def result_request(attempt_id: str, result: WorkerResult) -> dict[str, Any]:
    """The ``/api/cell/result`` body for one finished attempt."""
    hof = result.hof_kcal_mol
    if result.ok and hof is not None and math.isfinite(hof):
        return {
            "attempt_id": attempt_id,
            "result": {
                "mol_id": result.mol_id,
                "smiles": result.smiles,
                "signature": result.signature,
                "calculation_id": result.calculation_id,
                "conformer_index": result.conformer_index,
                "elements": list(result.elements),
                "coordinates": [list(xyz) for xyz in result.coordinates],
                "hof_kcal_mol": hof,
                "lowest_frequency_cm1": result.lowest_frequency_cm1,
                "conformers_evaluated": result.conformers_evaluated,
                "elapsed_s": result.elapsed_s,
            },
        }
    reason = result.error or result.state
    if result.ok:  # verified, yet no usable energy: never sent as a success
        reason = f"HOF não finito: {hof}"
    return failure_request(attempt_id, reason)


def failure_request(attempt_id: str, reason: str) -> dict[str, Any]:
    return {"attempt_id": attempt_id, "failure": (reason or "falha")[:_FAILURE_MAX]}


@dataclass
class _Busy:
    attempt_id: str
    mol_id: str
    smiles: str
    started_at: str


class CellLoop:
    def __init__(
        self,
        config: CellConfig,
        client: CentralClient,
        *,
        cores: int | None = None,
        executables: Mapping[str, str] | None = None,
        probe: ProgramProbe = probe_programs,
        pool_factory: PoolFactory = _spawn_pool,
        task_fn: TaskFn = run_task,
        timings: Timings | None = None,
        stop: threading.Event | None = None,
        install_signals: bool = True,
        hostname: str | None = None,
        abort: Callable[[], None] = _abort_process_group,
    ) -> None:
        self.config = config
        self.client = client
        self.cores = cores or cpu_cores()
        self.workers = worker_count(config, self.cores)
        self._executables = executables
        self._probe = probe
        self._pool_factory = pool_factory
        self._task_fn = task_fn
        self.timings = timings or Timings()
        self.stop = stop or threading.Event()
        self._install_signals = install_signals
        self.hostname = hostname or socket.gethostname()
        self._abort = abort
        self.outbox = Outbox(config.root)
        self._versions: dict[str, str] = {}
        self._heartbeat_s = 30.0
        self._state: dict[str, Any] = {}
        self._busy: dict[Future[WorkerResult], _Busy] = {}

    def __repr__(self) -> str:
        return f"CellLoop({self.config!r})"

    # ------------------------------------------------------------- state
    def _set(self, **changes: Any) -> None:
        self._state.update(changes)
        self._state["busy"] = [vars(b) for b in self._busy.values()]
        self._state["updated_at"] = _now()
        atomic_write_json(self.config.root / STATE_FILENAME, self._state)

    def _contact(self) -> None:
        self._state["last_contact_at"] = _now()

    def _end(self, code: int, phase: str, error: str | None = None) -> int:
        if error:
            log.error("%s: %s", phase, error)
            self._set(phase=phase, last_error=error)
        else:
            self._set(phase=phase)
        return code

    # ----------------------------------------------------------- results
    def queue_result(self, attempt_id: str, mol_id: str, result: WorkerResult) -> None:
        """Put one finished attempt on disk; it is sent from there."""
        self.outbox.put(attempt_id, mol_id, result_request(attempt_id, result))

    def _settle(self, future: Future[WorkerResult]) -> None:
        busy = self._busy.pop(future)
        try:
            result = future.result()
        except Exception as exc:  # worker process died or was killed
            request = failure_request(busy.attempt_id, f"{type(exc).__name__}: {exc}")
            self.outbox.put(busy.attempt_id, busy.mol_id, request)
            log.warning("%s: worker falhou: %s", busy.mol_id, exc)
            return
        self.queue_result(busy.attempt_id, busy.mol_id, result)
        log.info("%s: %s", busy.mol_id, "ok" if result.ok else result.error)

    # -------------------------------------------------------------- run
    def run(self) -> int:
        """Run until stopped or refused; return the process exit code."""
        if self._install_signals:
            self._signals()
        root = self.config.root
        self.config.store_root.mkdir(parents=True, exist_ok=True)
        # Whatever is in the work folder was interrupted: it starts over.
        shutil.rmtree(self.config.work_dir, ignore_errors=True)
        self.config.work_dir.mkdir(parents=True, exist_ok=True)
        self._state = {
            "pid": os.getpid(),
            "started_at": _now(),
            "central_url": self.config.central_url,
            "host": self.hostname,
            "cores": self.cores,
            "workers": self.workers,
            "crest_threads": self.config.crest_threads,
            "root": str(root),
            "last_error": None,
        }
        self._set(phase="iniciando")
        executables = self._executables or self._default_executables()
        self._versions = self._probe(executables)
        self._executables = executables

        self._set(phase="enviando pendentes")
        report = flush(self.outbox, self.client.result)
        if report.stop:
            return self._end(EXIT_TOKEN, "token recusado", report.last_error)
        if report.sent or report.dead:
            self._contact()
            log.info("pendentes: %d enviados, %d recusados", report.sent, report.dead)

        code = self._register()
        if code is not None:
            return code
        with self._pool_factory(self.workers) as pool:
            return self._compute(pool)

    def _default_executables(self) -> dict[str, str]:
        from semi_imperium.training.__main__ import _executables

        return _executables(self.config.bin_dir)

    def _signals(self) -> None:
        def handler(_signum: int, _frame: object) -> None:
            if self.stop.is_set():
                self._abort()
            self.stop.set()

        signal.signal(signal.SIGTERM, handler)
        signal.signal(signal.SIGINT, handler)

    def _register(self) -> int | None:
        """Register (again); ``None`` when registered, else the exit code."""
        body = {
            "host": self.hostname,
            "cores": self.cores,
            "workers": self.workers,
            "crest_threads": self.config.crest_threads,
            "program_versions": self._versions,
            "protocol_signature": canonical_signature(canonical_settings()),
        }
        failures = 0
        while not self.stop.is_set():
            self._set(phase="registrando")
            reply = self.client.register(body)
            verdict = classify(reply.status)
            if verdict is Verdict.SENT:
                self._contact()
                self._heartbeat_s = float(reply.body.get("heartbeat_s", 30.0))
                released = reply.body.get("released") or []
                if released:
                    log.info("a central devolveu à fila: %s", ", ".join(released))
                self._set(
                    phase="registrada",
                    cell_id=reply.body.get("cell_id"),
                    name=reply.body.get("name"),
                    admission=reply.body.get("admission"),
                    last_error=None,
                )
                return None
            if verdict is Verdict.STOP:
                return self._end(EXIT_TOKEN, "token recusado", reply.detail)
            if verdict is Verdict.DEAD:
                return self._end(EXIT_REFUSED, "recusada", reply.detail)
            failures += 1
            self._set(phase="sem conexão", last_error=self._describe(reply))
            self.stop.wait(self.timings.backoff.delay(failures))
        return self._end(EXIT_OK, "parado")

    @staticmethod
    def _describe(reply: Reply) -> str:
        return f"{reply.status or 'sem resposta'}: {reply.detail}"

    def _task(self, molecule: Mapping[str, Any], run_id: str) -> WorkerTask:
        executables = self._executables or {}
        return WorkerTask(
            mol_id=str(molecule["mol_id"]),
            smiles=str(molecule["smiles"]),
            charge=int(molecule["charge"]),
            multiplicity=int(molecule["multiplicity"]),
            run_id=run_id,
            store_root=str(self.config.store_root),
            work_dir=str(self.config.work_dir),
            crest_executable=executables["crest"],
            xtb_executable=executables["xtb"],
            mopac_executable=executables["mopac"],
            crest_threads=self.config.crest_threads,
            crest_timeout_s=self.config.crest_timeout_s,
            mopac_timeout_s=self.config.mopac_timeout_s,
            keep_artifacts=self.config.keep_artifacts,
            program_versions=tuple(sorted(self._versions.items())),
            hostname=self.hostname,
        )

    def _start_run(self) -> tuple[SemiImperiumStore, RunRecord]:
        stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
        run_id = f"cell-{self._state.get('cell_id')}-{stamp}-{uuid.uuid4().hex[:6]}"
        settings = canonical_settings()
        store = SemiImperiumStore(self.config.store_root)
        run = RunRecord(
            run_id=run_id,
            configuration=canonical_configuration(settings),
            provenance=build_provenance(
                settings, hostname=self.hostname, program_versions=self._versions
            ),
            label="delta-entalpia-celula",
        ).transition_to(RunState.RUNNING)
        store.save_run(run)
        return store, run

    def _compute(self, pool: Executor) -> int:
        timings = self.timings
        beat_every = (
            self._heartbeat_s if timings.heartbeat_s is None else timings.heartbeat_s
        )
        store, run = self._start_run()
        code, final = EXIT_OK, "parado"
        draining = token_refused = False
        next_claim = next_flush = next_beat = 0.0
        claim_failures = flush_failures = 0
        while True:
            now = time.monotonic()
            if self.stop.is_set() and not draining:
                draining = True
                log.info("parando: termina o que está em andamento")
                self._set(phase="parando")

            pending, _dead = self.outbox.counts()
            if pending and now >= next_flush and not token_refused:
                report = flush(self.outbox, self.client.result)
                if report.sent or report.dead:
                    self._contact()
                if report.stop:
                    token_refused, draining = True, True
                    code, final = EXIT_TOKEN, "token recusado"
                    self._set(last_error=report.last_error)
                elif report.retry:
                    flush_failures += 1
                    next_flush = now + timings.backoff.delay(flush_failures)
                    self._set(last_error=report.last_error)
                else:
                    flush_failures = 0

            free = self.workers - len(self._busy)
            if not draining and free > 0 and now >= next_claim:
                reply = self.client.claim(free)
                verdict = classify(reply.status)
                if verdict is Verdict.SENT:
                    self._contact()
                    claim_failures = 0
                    molecules = reply.body.get("molecules") or []
                    for molecule in molecules:
                        future = pool.submit(
                            self._task_fn, self._task(molecule, run.run_id)
                        )
                        self._busy[future] = _Busy(
                            attempt_id=str(molecule["attempt_id"]),
                            mol_id=str(molecule["mol_id"]),
                            smiles=str(molecule["smiles"]),
                            started_at=_now(),
                        )
                    if not molecules:
                        next_claim = now + timings.idle_s
                    self._set(
                        phase="calculando" if self._busy else "ociosa", last_error=None
                    )
                elif verdict is Verdict.STOP:
                    token_refused, draining = True, True
                    code, final = EXIT_TOKEN, "token recusado"
                    self._set(last_error=reply.detail)
                elif reply.status == 403:  # admission pending or blocked
                    self._contact()
                    next_claim = now + timings.admission_s
                    self._set(phase="em admissão", last_error=reply.detail)
                elif reply.status == 409:  # the central forgot or refused us
                    log.warning(
                        "claim recusado (%s); registrando de novo", reply.detail
                    )
                    again = self._register()
                    if again is not None:
                        draining, code = True, again
                        final = str(self._state.get("phase", final))
                else:
                    claim_failures += 1
                    next_claim = now + timings.backoff.delay(claim_failures)
                    self._set(last_error=self._describe(reply))

            if now >= next_beat and not token_refused:
                ids = [b.attempt_id for b in self._busy.values()]
                reply = self.client.heartbeat(ids)
                next_beat = now + beat_every
                if reply.status == 200:
                    self._contact()
                    lost = reply.body.get("not_held") or []
                    if lost:
                        # Released at the central; the result is still sent
                        # and accepted if nobody finished the molecule first.
                        log.warning("a central não reserva mais: %s", ", ".join(lost))
                elif classify(reply.status) is Verdict.STOP:
                    token_refused, draining = True, True
                    code, final = EXIT_TOKEN, "token recusado"
                    self._set(last_error=reply.detail)
                else:
                    self._set(last_error=self._describe(reply))

            if draining and not self._busy:
                break

            deadlines = [next_beat]
            if not draining and self.workers > len(self._busy):
                deadlines.append(next_claim)
            if self.outbox.counts()[0]:
                deadlines.append(next_flush)
            timeout = min(max(min(deadlines) - time.monotonic(), 0.0), 1.0)
            if self._busy:
                done, _ = wait(self._busy, timeout=timeout, return_when=FIRST_COMPLETED)
                for future in done:
                    self._settle(future)
                    next_flush = 0.0  # send it now
                if done:
                    self._set(phase="parando" if draining else "calculando")
            else:
                self.stop.wait(timeout)

        if not token_refused and self.outbox.counts()[0]:
            report = flush(self.outbox, self.client.result)
            if report.retry or report.stop:
                log.warning("resultados ficam no disco para o próximo início")
        store.save_run(run.transition_to(RunState.COMPLETED))
        return self._end(code, final)
