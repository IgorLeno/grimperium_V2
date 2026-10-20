from __future__ import annotations

import os
from collections.abc import Callable, Mapping
from concurrent.futures import Executor, ThreadPoolExecutor
from pathlib import Path

import pytest

from semi_imperium.training.geometry_store import GeometryStore
from semi_imperium.training.job_state import JobAlreadyRunning, JobState
from semi_imperium.training.molecule_queue import CbsReadError, read_progress
from semi_imperium.training.runner import RunSettings, run_training
from semi_imperium.training.worker import WorkerResult, WorkerTask
from tests.training.conftest import CbsFactory, row


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


def _settings(tmp_path: Path, cbs: Path, **overrides: object) -> RunSettings:
    values: dict[str, object] = {
        "cbs_csv": cbs,
        "store_root": tmp_path / "store",
        "work_dir": tmp_path / "work",
        "workers": 2,
    }
    values.update(overrides)
    return RunSettings(**values)  # type: ignore[arg-type]


def _cbs(make_cbs: CbsFactory) -> Path:
    return make_cbs([row("a", nheavy=1), row("b", nheavy=2), row("c", nheavy=3)])


def _db(tmp_path: Path) -> Path:
    return tmp_path / "store" / "training.sqlite3"


def test_run_records_results_in_queue_geometry_store_and_job(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = _cbs(make_cbs)
    summary = run_training(
        _settings(tmp_path, cbs),
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
    progress = read_progress(_db(tmp_path), cbs)
    assert (progress.total, progress.ok, progress.pending) == (3, 3, 0)
    store = GeometryStore(_db(tmp_path))
    assert store.mol_ids() == {"a", "b", "c"}
    job = JobState(_db(tmp_path)).snapshot()
    assert job.state == "finished" and job.done == 3 and job.failed == 0


def test_the_cbs_csv_is_never_written(tmp_path: Path, make_cbs: CbsFactory) -> None:
    cbs = _cbs(make_cbs)
    before = (cbs.read_bytes(), cbs.stat().st_mtime_ns)
    run_training(
        _settings(tmp_path, cbs),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert (cbs.read_bytes(), cbs.stat().st_mtime_ns) == before


def test_limit_caps_the_number_of_molecules(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = _cbs(make_cbs)
    summary = run_training(
        _settings(tmp_path, cbs, limit=1),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert summary.attempted == 1
    progress = read_progress(_db(tmp_path), cbs)
    assert (progress.ok, progress.pending) == (1, 2)
    assert GeometryStore(_db(tmp_path)).mol_ids() == {"a"}


def test_a_second_run_continues_where_the_first_stopped(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = _cbs(make_cbs)
    for expected_ok in (1, 1, 1, 0):
        summary = run_training(
            _settings(tmp_path, cbs, limit=1),
            pool_factory=_threads,
            task_fn=_fake(),
            install_signals=False,
        )
        assert summary.ok == expected_ok
    assert GeometryStore(_db(tmp_path)).mol_ids() == {"a", "b", "c"}


def test_failure_goes_to_rerun_and_its_reason_is_kept(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = _cbs(make_cbs)
    summary = run_training(
        _settings(tmp_path, cbs),
        pool_factory=_threads,
        task_fn=_fake(fail={"b"}),
        install_signals=False,
    )
    assert (summary.ok, summary.failed) == (2, 1)
    progress = read_progress(_db(tmp_path), cbs)
    assert (progress.ok, progress.rerun) == (2, 1)
    assert "b" not in GeometryStore(_db(tmp_path)).mol_ids()
    failures = JobState(_db(tmp_path)).recent_failures()
    assert failures[0][0] == "b" and "saddle" in failures[0][1]


def test_repeated_failures_end_in_skip(tmp_path: Path, make_cbs: CbsFactory) -> None:
    cbs = make_cbs([row("a")])
    for _ in range(3):
        run_training(
            _settings(tmp_path, cbs, max_reruns=3),
            pool_factory=_threads,
            task_fn=_fake(fail={"a"}),
            install_signals=False,
        )
    progress = read_progress(_db(tmp_path), cbs)
    assert (progress.skipped, progress.rerun) == (1, 0)
    again = run_training(
        _settings(tmp_path, cbs),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert again.attempted == 0


def test_orphan_running_molecules_are_requeued(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    from semi_imperium.training.molecule_queue import TrainingQueue

    cbs = make_cbs([row("orphan")])
    seeded = TrainingQueue(_db(tmp_path), cbs)
    seeded.ensure_seeded()
    seeded.mark_running("orphan")  # a previous run died while it was in flight
    summary = run_training(
        _settings(tmp_path, cbs),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert summary.ok == 1


def test_stop_request_prevents_new_dispatch(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = _cbs(make_cbs)
    state = JobState(_db(tmp_path))

    def stopping_task(task: WorkerTask) -> WorkerResult:
        state.request_stop()  # the first molecule asks the job to stop
        return _fake()(task)

    summary = run_training(
        _settings(tmp_path, cbs, workers=1),
        pool_factory=_threads,
        task_fn=stopping_task,
        install_signals=False,
    )
    assert summary.stopped is True and summary.attempted == 1
    assert read_progress(_db(tmp_path), cbs).pending == 2


def test_a_live_job_blocks_a_second_run(tmp_path: Path, make_cbs: CbsFactory) -> None:
    cbs = _cbs(make_cbs)
    JobState(_db(tmp_path)).begin(pid=os.getpid(), workers=1, limit=None)
    with pytest.raises(JobAlreadyRunning):
        run_training(
            _settings(tmp_path, cbs),
            pool_factory=_threads,
            task_fn=_fake(),
            install_signals=False,
        )


def test_an_unreadable_cbs_is_refused_before_any_job_state(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = make_cbs([{"mol_id": "a", "smiles": "C"}])  # no CBS columns
    with pytest.raises(CbsReadError):
        run_training(
            _settings(tmp_path, cbs),
            pool_factory=_threads,
            task_fn=_fake(),
            install_signals=False,
        )
    assert JobState(_db(tmp_path)).snapshot().state == "idle"


def test_the_run_and_every_task_carry_program_versions_and_hostname(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    import socket

    from semi_imperium.persistence import SemiImperiumStore

    probed: list[dict[str, str]] = []

    def probe(executables: Mapping[str, str]) -> dict[str, str]:
        probed.append(dict(executables))
        return {"crest": "3.0.2", "crest_sha256": "abc"}

    tasks: list[WorkerTask] = []
    fake = _fake()

    def task_fn(task: WorkerTask) -> WorkerResult:
        tasks.append(task)
        return fake(task)

    run_training(
        _settings(tmp_path, _cbs(make_cbs), crest_executable="/opt/env/bin/crest"),
        pool_factory=_threads,
        task_fn=task_fn,
        install_signals=False,
        probe=probe,
    )

    assert probed == [{"crest": "/opt/env/bin/crest", "xtb": "xtb", "mopac": "mopac"}]
    (run,) = SemiImperiumStore(tmp_path / "store").list_runs()
    assert run.provenance.program_versions == {"crest": "3.0.2", "crest_sha256": "abc"}
    assert run.provenance.hostname == socket.gethostname()
    assert {t.program_versions for t in tasks} == {
        (("crest", "3.0.2"), ("crest_sha256", "abc"))
    }
    assert {t.hostname for t in tasks} == {socket.gethostname()}
