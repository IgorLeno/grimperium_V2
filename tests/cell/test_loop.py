"""The cell's loop against the real cell API: register, claim, compute, send."""

from __future__ import annotations

import json
import logging
import threading
from concurrent.futures import Executor, ThreadPoolExecutor
from pathlib import Path
from typing import Any

import pytest

from semi_imperium.cell.loop import (
    EXIT_OK,
    EXIT_REFUSED,
    EXIT_TOKEN,
    CellLoop,
    Timings,
    read_state,
)
from semi_imperium.cell.outbox import Backoff, Outbox
from semi_imperium.training.worker import WorkerTask
from tests.cell.conftest import (
    CELL_PROGRAMS,
    EXECUTABLES,
    SIGNATURE,
    Central,
    SpyClient,
    Tasks,
    cell_config,
    computed,
)

FAST = Timings(
    idle_s=0.01, admission_s=0.01, heartbeat_s=0.0, backoff=Backoff(0.01, 0.05)
)


def _threads(workers: int) -> Executor:
    return ThreadPoolExecutor(max_workers=workers)


def _loop(
    tmp_path: Path,
    central: Central,
    token: str,
    tasks: Tasks,
    *,
    programs: dict[str, str] | None = None,
    **config: Any,
) -> tuple[CellLoop, SpyClient]:
    client = SpyClient(central, token)
    loop = CellLoop(
        cell_config(tmp_path, token, **config),
        client,
        cores=4,
        executables=EXECUTABLES,
        probe=lambda _exe: dict(programs or CELL_PROGRAMS),
        pool_factory=_threads,
        task_fn=tasks,
        timings=FAST,
        stop=tasks.stop,
        install_signals=False,
    )
    return loop, client


def test_a_cell_computes_the_queue_and_reports_each_result(
    tmp_path: Path, central: Central
) -> None:
    cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=3)
    loop, client = _loop(tmp_path, central, token, tasks)

    assert loop.run() == EXIT_OK

    assert central.progress().ok == 3
    register = client.calls[client.names().index("register")][1]
    assert register["workers"] == 2 and register["crest_threads"] == 2
    assert register["program_versions"] == CELL_PROGRAMS
    # No claim ever asks for more than the free workers.
    assert all(count <= 2 for name, count in client.calls if name == "claim")
    assert "heartbeat" in client.names()
    assert Outbox(tmp_path / "celula").counts() == (0, 0)
    state = read_state(tmp_path / "celula")
    assert state["cell_id"] == cell_id and state["phase"] == "parado"
    assert state["busy"] == []
    # Each task runs under the cell's own folders and threads.
    task = tasks.seen[0]
    assert task.store_root == str(tmp_path / "celula" / "store")
    assert task.work_dir == str(tmp_path / "celula" / "work")
    assert task.crest_threads == 2
    assert dict(task.program_versions) == CELL_PROGRAMS


def test_a_failure_is_reported_and_not_claimed_again_by_the_cell(
    tmp_path: Path, central: Central
) -> None:
    _cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=3)
    tasks.fail = {"m0"}
    loop, client = _loop(tmp_path, central, token, tasks)

    assert loop.run() == EXIT_OK

    assert central.status("m0") == "rerun"
    assert [t.mol_id for t in tasks.seen].count("m0") == 1
    failures = [b for n, b in client.calls if n == "result" and "failure" in b]
    assert [f["failure"] for f in failures] == ["failed: crest timeout"]


def test_results_wait_on_disk_while_the_central_is_unreachable(
    tmp_path: Path, central: Central
) -> None:
    _cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=3)
    loop, client = _loop(tmp_path, central, token, tasks)
    client.fail_results = 4

    assert loop.run() == EXIT_OK

    assert central.progress().ok == 3
    assert Outbox(tmp_path / "celula").counts() == (0, 0)


def test_a_restart_sends_finished_results_before_registering(
    tmp_path: Path, central: Central
) -> None:
    cell_id, token = central.new_cell()
    before = SpyClient(central, token)
    assert before.register(REGISTRATION).status == 200
    (held,) = central.queue.claim(cell_id, 1)
    tasks = Tasks(threading.Event(), stop_after=2)
    loop, client = _loop(tmp_path, central, token, tasks)
    # Power cut after the result reached the disk but before it was sent.
    loop.queue_result(held.attempt_id, held.mol_id, computed(_task_for(held)))

    assert loop.run() == EXIT_OK

    assert client.names()[:2] == ["result", "register"]
    assert central.status(held.mol_id) == "ok"
    assert held.mol_id not in [t.mol_id for t in tasks.seen]


REGISTRATION = {
    "host": "lab",
    "cores": 4,
    "workers": 2,
    "crest_threads": 2,
    "program_versions": CELL_PROGRAMS,
    "protocol_signature": SIGNATURE,
}


def _task_for(held: Any) -> WorkerTask:
    return WorkerTask(
        mol_id=held.mol_id,
        smiles=held.smiles,
        charge=held.charge,
        multiplicity=held.multiplicity,
        run_id="r",
        store_root="s",
        work_dir="w",
        crest_executable="crest",
        xtb_executable="xtb",
        mopac_executable="mopac",
        crest_threads=2,
        crest_timeout_s=1.0,
        mopac_timeout_s=1.0,
        keep_artifacts=False,
    )


def test_programs_other_than_the_central_s_stop_the_cell(
    tmp_path: Path, central: Central
) -> None:
    _cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=99)
    other = {**CELL_PROGRAMS, "xtb": "6.6.1"}
    loop, client = _loop(tmp_path, central, token, tasks, programs=other)

    assert loop.run() == EXIT_REFUSED

    assert "claim" not in client.names() and tasks.seen == []
    state = read_state(tmp_path / "celula")
    assert "xtb: 6.6.1" in state["last_error"] and state["phase"] == "recusada"


def test_a_revoked_token_stops_the_cell(tmp_path: Path, central: Central) -> None:
    cell_id, token = central.new_cell()
    central.cells.revoke(cell_id)
    tasks = Tasks(threading.Event(), stop_after=99)
    loop, _client = _loop(tmp_path, central, token, tasks)

    assert loop.run() == EXIT_TOKEN
    assert read_state(tmp_path / "celula")["phase"] == "token recusado"


def test_a_new_cell_passes_the_admission_test_then_gets_work(
    tmp_path: Path,
) -> None:
    from semi_imperium.training.cells import ADMISSION_PASSED, ADMISSION_REFERENCES

    central = Central(tmp_path, references=True)
    cell_id, token = central.new_cell(admitted=False)
    tasks = Tasks(threading.Event(), stop_after=8)
    loop, client = _loop(tmp_path, central, token, tasks)

    assert loop.run() == EXIT_OK
    computed = [t.mol_id for t in tasks.seen]
    assert sorted(computed[:5]) == sorted(ADMISSION_REFERENCES)
    assert sorted(computed[5:]) == ["m0", "m1", "m2"]
    assert central.cells.get(cell_id).admission == ADMISSION_PASSED
    assert central.progress().ok == 3 + len(ADMISSION_REFERENCES)


def test_a_blocked_cell_waits_and_starts_once_admitted(
    tmp_path: Path, central: Central
) -> None:
    from semi_imperium.training.cells import ADMISSION_BLOCKED, ADMISSION_PASSED

    cell_id, token = central.new_cell(admitted=False)
    central.cells.set_admission(cell_id, ADMISSION_BLOCKED)
    tasks = Tasks(threading.Event(), stop_after=3)
    loop, client = _loop(tmp_path, central, token, tasks)
    phases: list[str] = []

    def admit_after_refusal(name: str, reply: Any) -> None:
        if name == "claim" and reply.status == 403:
            phases.append(read_state(tmp_path / "celula")["phase"])
            central.cells.set_admission(cell_id, ADMISSION_PASSED)

    client.on_call = admit_after_refusal

    assert loop.run() == EXIT_OK
    assert central.progress().ok == 3
    assert phases == ["registrada"]  # state written before the refused claim


def test_the_heartbeat_names_every_molecule_in_flight(
    tmp_path: Path, central: Central
) -> None:
    _cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=2)
    tasks.gate = threading.Event()
    loop, client = _loop(tmp_path, central, token, tasks)
    beats: list[list[str]] = []
    states: list[dict[str, Any]] = []

    def open_gate_after_a_full_beat(name: str, reply: Any) -> None:
        if name == "heartbeat" and len(client.calls[-1][1]) == 2:
            beats.append(client.calls[-1][1])
            states.append(read_state(tmp_path / "celula"))
            tasks.gate.set()  # type: ignore[union-attr]

    client.on_call = open_gate_after_a_full_beat

    assert loop.run() == EXIT_OK
    held = {
        b["attempt_id"]
        for n, b in client.calls
        if n == "result" and b.get("result") is not None
    }
    assert len(beats[0]) == 2 and set(beats[0]) <= held
    assert {w["mol_id"] for w in states[0]["busy"]} == {"m0", "m1"}
    assert states[0]["phase"] == "calculando"


def test_stopping_finishes_what_is_running_and_claims_nothing_more(
    tmp_path: Path, central: Central
) -> None:
    _cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=1)  # stop after the first result
    loop, client = _loop(tmp_path, central, token, tasks)

    assert loop.run() == EXIT_OK

    claims = client.names().count("claim")
    assert claims == 1  # the first claim took both workers; nothing after stop
    assert central.progress().ok == 2 and central.progress().pending == 1
    assert central.queue.held_attempts(_cell_id) == set()


def test_leftovers_of_an_interrupted_run_are_wiped(
    tmp_path: Path, central: Central
) -> None:
    _cell_id, token = central.new_cell()
    leftover = tmp_path / "celula" / "work" / "crest" / "m9" / "crest_rotamers.xyz"
    leftover.parent.mkdir(parents=True)
    leftover.write_text("partial")
    tasks = Tasks(threading.Event(), stop_after=3)
    loop, _client = _loop(tmp_path, central, token, tasks)

    assert loop.run() == EXIT_OK
    assert not leftover.exists()


def test_the_token_is_never_written_anywhere(
    tmp_path: Path, central: Central, caplog: pytest.LogCaptureFixture
) -> None:
    _cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=3)
    tasks.fail = {"m1"}
    loop, client = _loop(tmp_path, central, token, tasks)
    client.fail_results = 1
    caplog.set_level(logging.DEBUG)

    assert loop.run() == EXIT_OK

    assert token not in caplog.text
    assert token not in repr(loop) and token not in repr(client)
    for path in (tmp_path / "celula").rglob("*"):
        if path.is_file():
            assert token not in path.read_text(errors="replace"), path
    state = json.loads((tmp_path / "celula" / "state.json").read_text())
    assert "token" not in state


def test_each_heartbeat_refreshes_the_last_contact_on_disk(
    tmp_path: Path, central: Central, monkeypatch: pytest.MonkeyPatch
) -> None:
    import itertools

    from semi_imperium.cell import loop as loop_module

    ticks = itertools.count()
    monkeypatch.setattr(loop_module, "_now", lambda: f"t{next(ticks):06d}")
    _cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=2)
    tasks.gate = threading.Event()  # both workers stay busy: no claims
    loop, client = _loop(tmp_path, central, token, tasks)
    contacts: list[str] = []

    def watch(name: str, _reply: Any) -> None:
        if name == "heartbeat" and len(client.calls[-1][1]) == 2:
            contacts.append(read_state(tmp_path / "celula")["last_contact_at"])
            if len(contacts) == 3:
                tasks.gate.set()  # type: ignore[union-attr]

    client.on_call = watch

    assert loop.run() == EXIT_OK
    assert len(set(contacts)) == 3  # written after every heartbeat, not only claims
