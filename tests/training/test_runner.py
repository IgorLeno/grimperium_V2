from __future__ import annotations

import os
from collections.abc import Callable
from concurrent.futures import Executor, ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import pytest

from semi_imperium.training.geometry_store import GeometryStore
from semi_imperium.training.job_state import JobAlreadyRunning, JobState
from semi_imperium.training.ledger import NotMigratedError
from semi_imperium.training.runner import RunSettings, run_training
from semi_imperium.training.worker import WorkerResult, WorkerTask
from tests.training.conftest import CsvFactory, row


def _threads(workers: int) -> Executor:
    return ThreadPoolExecutor(max_workers=workers)


def _fake(fail: set[str] | None = None) -> Callable[[WorkerTask], WorkerResult]:
    fail = fail or set()

    def task_fn(task: WorkerTask) -> WorkerResult:
        bad = task.mol_id in fail
        return WorkerResult(
            mol_id=task.mol_id,
            smiles=task.smiles,
            ok=not bad,
            state="saddle" if bad else "verified",
            hof_kcal_mol=None if bad else -1.5,
            lowest_frequency_cm1=None if bad else 100.0,
            conformers_evaluated=0 if bad else 2,
            conformer_index=None if bad else 1,
            elements=() if bad else ("C",),
            coordinates=() if bad else ((0.0, 0.0, 0.0),),
            calculation_id=f"calc-{task.mol_id}",
            signature="sig",
            elapsed_s=6.0,
            error="saddle: imaginary mode" if bad else None,
        )

    return task_fn


def _settings(tmp_path: Path, csv: Path, **overrides: object) -> RunSettings:
    values: dict[str, object] = {
        "csv_path": csv,
        "store_root": tmp_path / "store",
        "work_dir": tmp_path / "work",
        "workers": 2,
    }
    values.update(overrides)
    return RunSettings(**values)  # type: ignore[arg-type]


def _csv(make_csv: CsvFactory) -> Path:
    return make_csv(
        [
            row("a", "Pending", nheavy=1, protocol=""),
            row("b", "Pending", nheavy=2, protocol=""),
            row("c", "Pending", nheavy=3, protocol=""),
            row("old", "OK", nheavy=1, protocol="1", H298_pm7=-3.83),
        ]
    )


def test_run_records_results_in_csv_geometry_store_and_job(
    tmp_path: Path, make_csv: CsvFactory
) -> None:
    csv = _csv(make_csv)
    summary = run_training(
        _settings(tmp_path, csv),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert (summary.attempted, summary.ok, summary.failed, summary.stopped) == (
        3,
        3,
        0,
        False,
    )
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert [saved.loc[m, "protocol"] for m in ("a", "b", "c")] == ["2", "2", "2"]
    assert saved.loc["old", "protocol"] == "1"
    store = GeometryStore(tmp_path / "store" / "training.sqlite3")
    assert store.mol_ids() == {"a", "b", "c"}
    job = JobState(tmp_path / "store" / "training.sqlite3").snapshot()
    assert job.state == "finished" and job.done == 3 and job.failed == 0


def test_limit_caps_the_number_of_molecules(
    tmp_path: Path, make_csv: CsvFactory
) -> None:
    csv = _csv(make_csv)
    summary = run_training(
        _settings(tmp_path, csv, limit=1),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert summary.attempted == 1
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["a", "status"] == "OK"
    assert saved.loc["b", "status"] == "Pending"


def test_failure_goes_to_rerun_and_its_reason_is_kept(
    tmp_path: Path, make_csv: CsvFactory
) -> None:
    csv = _csv(make_csv)
    summary = run_training(
        _settings(tmp_path, csv),
        pool_factory=_threads,
        task_fn=_fake(fail={"b"}),
        install_signals=False,
    )
    assert (summary.ok, summary.failed) == (2, 1)
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert (saved.loc["b", "status"], saved.loc["b", "protocol"]) == ("Rerun", "")
    store = GeometryStore(tmp_path / "store" / "training.sqlite3")
    assert "b" not in store.mol_ids()
    failures = JobState(tmp_path / "store" / "training.sqlite3").recent_failures()
    assert failures[0][0] == "b" and "saddle" in failures[0][1]


def test_orphan_running_rows_are_requeued(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = make_csv([row("orphan", "Running", nheavy=1, protocol="")])
    summary = run_training(
        _settings(tmp_path, csv),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert summary.ok == 1


def test_stop_request_prevents_new_dispatch(
    tmp_path: Path, make_csv: CsvFactory
) -> None:
    csv = _csv(make_csv)
    state = JobState(tmp_path / "store" / "training.sqlite3")

    def stopping_task(task: WorkerTask) -> WorkerResult:
        state.request_stop()  # the first molecule asks the job to stop
        return _fake()(task)

    summary = run_training(
        _settings(tmp_path, csv, workers=1),
        pool_factory=_threads,
        task_fn=stopping_task,
        install_signals=False,
    )
    assert summary.stopped is True and summary.attempted == 1
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["b", "status"] == "Pending"


def test_a_live_job_blocks_a_second_run(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = _csv(make_csv)
    JobState(tmp_path / "store" / "training.sqlite3").begin(
        pid=os.getpid(), workers=1, limit=None
    )
    with pytest.raises(JobAlreadyRunning):
        run_training(
            _settings(tmp_path, csv),
            pool_factory=_threads,
            task_fn=_fake(),
            install_signals=False,
        )


def test_unmigrated_csv_is_refused_before_any_work(
    tmp_path: Path, make_csv: CsvFactory
) -> None:
    csv = make_csv([row("a", "OK"), row("b", "Pending")])
    with pytest.raises(NotMigratedError):
        run_training(
            _settings(tmp_path, csv),
            pool_factory=_threads,
            task_fn=_fake(),
            install_signals=False,
        )
    assert JobState(tmp_path / "store" / "training.sqlite3").snapshot().state == "idle"


def test_empty_queue_finishes_cleanly(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = make_csv([row("old", "OK", protocol="1")])
    summary = run_training(
        _settings(tmp_path, csv),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert summary.attempted == 0
