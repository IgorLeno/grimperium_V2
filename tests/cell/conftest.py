"""A central (cell API over a temporary training database) for cell tests."""

from __future__ import annotations

import threading
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from grimperium.webapi.cells import create_cell_app
from grimperium.webapi.config import DashboardPaths
from semi_imperium.cell.client import CentralClient
from semi_imperium.cell.config import CellConfig
from semi_imperium.cell.outbox import Reply
from semi_imperium.training.cells import ADMISSION_PASSED, CellRegistry
from semi_imperium.training.db import database_path
from semi_imperium.training.molecule_queue import TrainingQueue, read_progress
from semi_imperium.training.protocol import canonical_settings, canonical_signature
from semi_imperium.training.worker import WorkerResult, WorkerTask
from tests.training.conftest import row

PROGRAMS = {
    "crest": "3.0.2",
    "crest_sha256": "aaa",
    "crest_path": "/central/bin/crest",
    "xtb": "6.7.1",
    "xtb_sha256": "bbb",
    "xtb_path": "/central/bin/xtb",
    "mopac": "23.2.2",
    "mopac_sha256": "ccc",
    "mopac_path": "/central/bin/mopac",
    "environment_sha256": "ddd",
}
#: The cell's own install: other folders, same binaries.
CELL_PROGRAMS = {
    **PROGRAMS,
    "crest_path": "/home/aluno/grimperium-celula/env/bin/crest",
}
SIGNATURE = canonical_signature(canonical_settings())
EXECUTABLES = {"crest": "crest", "xtb": "xtb", "mopac": "mopac"}


class SpyClient(CentralClient):
    """The real client over the in-process central, recording each call."""

    def __init__(self, central: Central, token: str) -> None:
        super().__init__("http://testserver", token, http=central.http)
        self.calls: list[tuple[str, Any]] = []
        self.fail_results = 0  # the next N results get no answer
        self.on_call: Callable[[str, Reply], None] | None = None

    def _record(self, name: str, body: Any, reply: Reply) -> Reply:
        self.calls.append((name, body))
        if self.on_call is not None:
            self.on_call(name, reply)
        return reply

    def names(self) -> list[str]:
        return [name for name, _ in self.calls]

    def register(self, body: dict[str, Any]) -> Reply:
        return self._record("register", body, super().register(body))

    def claim(self, count: int) -> Reply:
        return self._record("claim", count, super().claim(count))

    def heartbeat(self, attempt_ids: list[str]) -> Reply:
        return self._record("heartbeat", attempt_ids, super().heartbeat(attempt_ids))

    def result(self, request: dict[str, Any]) -> Reply:
        if self.fail_results > 0:
            self.fail_results -= 1
            return self._record("result", request, Reply(None, {}, "ConnectError"))
        return self._record("result", request, super().result(request))


class Central:
    def __init__(self, tmp_path: Path, molecules: int = 3) -> None:
        rows = [row(f"m{n}", nheavy=n + 1) for n in range(molecules)]
        pd.DataFrame(rows).to_csv(tmp_path / "cbs.csv", index=False)
        self.paths = DashboardPaths(
            semi_imperium_store=tmp_path / "central",
            exp_csv=tmp_path / "exp.csv",
            exp_vs_cbs_csv=tmp_path / "vs_cbs.csv",
            cbs_csv=tmp_path / "cbs.csv",
        )
        self.db = database_path(self.paths.semi_imperium_store)
        self.queue = TrainingQueue(self.db, self.paths.cbs_csv)
        self.queue.ensure_seeded()
        self.cells = CellRegistry(self.db)
        self.http = TestClient(
            create_cell_app(self.paths, central_programs=lambda: dict(PROGRAMS))
        )

    def new_cell(self, *, admitted: bool = True) -> tuple[str, str]:
        cell, token = self.cells.create("Lab")
        if admitted:
            self.cells.set_admission(cell.cell_id, ADMISSION_PASSED)
        return cell.cell_id, token

    def progress(self) -> Any:
        return read_progress(self.db, self.paths.cbs_csv)

    def status(self, mol_id: str) -> str:
        import sqlite3

        with sqlite3.connect(self.db) as conn:
            (status,) = conn.execute(
                "SELECT status FROM training_queue WHERE mol_id = ?", (mol_id,)
            ).fetchone()
        return str(status)


@pytest.fixture
def central(tmp_path: Path) -> Central:
    return Central(tmp_path)


def cell_config(tmp_path: Path, token: str, **changes: Any) -> CellConfig:
    values: dict[str, Any] = {
        "central_url": "http://testserver",
        "token": token,
        "root": tmp_path / "celula",
        "crest_threads": 2,
    }
    values.update(changes)
    return CellConfig(**values)


def computed(task: WorkerTask, *, ok: bool = True) -> WorkerResult:
    """What ``run_task`` returns for ``task``, without CREST or MOPAC."""
    return WorkerResult(
        mol_id=task.mol_id,
        smiles=task.smiles,
        ok=ok,
        state="verified" if ok else "failed",
        hof_kcal_mol=-12.5 if ok else None,
        lowest_frequency_cm1=101.0 if ok else None,
        conformers_evaluated=3 if ok else 0,
        conformer_index=1 if ok else None,
        elements=("C", "O") if ok else (),
        coordinates=((0.0, 0.0, 0.0), (0.0, 0.0, 1.2)) if ok else (),
        calculation_id=f"calc-{task.mol_id}",
        signature=SIGNATURE,
        elapsed_s=1.5,
        error=None if ok else "failed: crest timeout",
    )


class Tasks:
    """A fake ``run_task`` that stops the cell after ``stop_after`` results."""

    def __init__(self, stop: threading.Event, stop_after: int) -> None:
        self.stop = stop
        self.stop_after = stop_after
        self.seen: list[WorkerTask] = []
        self.fail: set[str] = set()
        self.gate: threading.Event | None = None
        self._lock = threading.Lock()

    def __call__(self, task: WorkerTask) -> WorkerResult:
        if self.gate is not None and not self.gate.wait(5):
            self.stop.set()  # a failed test must end, not idle forever
            raise AssertionError("gate never opened")
        with self._lock:
            self.seen.append(task)
            if len(self.seen) >= self.stop_after:
                self.stop.set()
        return computed(task, ok=task.mol_id not in self.fail)
