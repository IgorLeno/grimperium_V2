"""``celula painel``: what the panel shows for each state of the cell."""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import pytest
from rich.console import Console

from semi_imperium.cell.config import CellConfig
from semi_imperium.cell.outbox import Outbox
from semi_imperium.cell.panel import (
    STALE_S,
    Snapshot,
    duration,
    phase_of,
    read_snapshot,
    render,
    tail,
)

NOW = datetime(2026, 10, 23, 17, 0, 0, tzinfo=timezone.utc)


def _ago(seconds: float) -> str:
    return (NOW - timedelta(seconds=seconds)).isoformat(timespec="seconds")


def _state(**changes: Any) -> dict[str, Any]:
    state: dict[str, Any] = {
        "phase": "calculando",
        "name": "pc_elias",
        "cell_id": "cell-2f0068a3",
        "host": "IQ",
        "admission": "passed",
        "workers": 2,
        "crest_threads": 4,
        "started_at": _ago(3 * 3600 + 41 * 60),
        "updated_at": _ago(5),
        "last_contact_at": _ago(12),
        "last_error": None,
        "busy": [
            {
                "attempt_id": "a1",
                "mol_id": "cbs_04263",
                "smiles": "CC(C)OC(=O)c1ccccc1N",
                "started_at": _ago(23 * 60 + 41),
            },
            {
                "attempt_id": "adm-x",
                "mol_id": "cbs_08982",
                "smiles": "OCCO",
                "started_at": _ago(420),
            },
        ],
        "totals": {
            "date": "2026-10-23",
            "session_ok": 7,
            "session_failed": 0,
            "day_ok": 12,
            "day_failed": 1,
        },
        "recent": [
            {
                "at": _ago(60),
                "attempt_id": "a9",
                "mol_id": "cbs_31877",
                "smiles": "CCO",
                "ok": True,
                "hof_kcal_mol": -87.42,
                "elapsed_s": 1102.0,
                "error": None,
            },
            {
                "at": _ago(600),
                "attempt_id": "a8",
                "mol_id": "cbs_77801",
                "smiles": "CCN",
                "ok": False,
                "hof_kcal_mol": None,
                "elapsed_s": None,
                "error": "CREST timeout (3600 s)",
            },
        ],
    }
    state.update(changes)
    return state


def _snap(state: dict[str, Any] | None = None, **changes: Any) -> Snapshot:
    values: dict[str, Any] = {
        "now": NOW,
        "running": True,
        "state": _state() if state is None else state,
        "central_url": "http://100.108.94.19:8100",
        "limit_s": 4500.0,
        "unsent": frozenset({"a9"}),
        "dead": 0,
        "log": (
            "2026-10-23 13:59:58,120 INFO semi_imperium.cell.loop: cbs_31877: ok",
            "2026-10-23 14:00:02,004 WARNING semi_imperium.cell.loop: sem resposta",
        ),
    }
    values.update(changes)
    return Snapshot(**values)


def _text(snap: Snapshot, width: int = 120, height: int = 40) -> str:
    console = Console(record=True, width=width, height=height, color_system=None)
    console.print(render(snap, width, height))
    return console.export_text()


def test_a_computing_cell_shows_workers_results_and_log() -> None:
    text = _text(_snap())

    assert "CALCULANDO" in text and "pc_elias" in text and "cell-2f0068a3" in text
    assert "✓ admitida" in text and "há 12 s" in text and "3 h 41 min" in text
    assert "cbs_04263" in text and "CC(C)OC(=O)c1ccccc1N" in text
    assert "00:23:41" in text
    assert "REF" in text  # the reference molecule of the admission test
    assert "Workers (2 × 4 threads)" in text
    assert "7 ✓" in text and "12 ✓" in text and "1 ✗" in text
    assert "ΔHf   -87.42 kcal/mol" in text and "18 min 22 s" in text
    assert "CREST timeout (3600 s)" in text
    assert "sem resposta" in text and "13:59:58" in text
    assert "http://100.108.94.19:8100" in text


def test_unsent_results_are_marked() -> None:
    line = next(line for line in _text(_snap()).splitlines() if "cbs_31877" in line)
    assert "⇡" in line
    sent = _text(_snap(unsent=frozenset()))
    line = next(line for line in sent.splitlines() if "cbs_31877" in line)
    assert "⇡" not in line


@pytest.mark.parametrize(
    ("phase", "label"),
    [
        ("calculando", "CALCULANDO"),
        ("ociosa", "OCIOSA"),
        ("sem conexão", "SEM CONEXÃO"),
        ("bloqueada", "BLOQUEADA"),
        ("parando", "PARANDO"),
        ("recusada", "RECUSADA"),
        ("token recusado", "TOKEN RECUSADO"),
        ("registrando", "REGISTRANDO"),
    ],
)
def test_each_phase_has_its_label(phase: str, label: str) -> None:
    snap = _snap(_state(phase=phase))
    assert phase_of(snap).label == label
    assert label in _text(snap)


def test_a_stopped_cell_says_so_and_what_it_was_doing() -> None:
    phase = phase_of(_snap(_state(phase="token recusado"), running=False))
    assert phase.label == "PARADA" and "token recusado" in phase.detail
    assert phase_of(_snap({}, running=False)).label == "NUNCA INICIADA"


def test_a_running_cell_that_stopped_writing_is_flagged() -> None:
    stale = _snap(_state(updated_at=_ago(STALE_S + 60)))
    assert phase_of(stale).label == "SEM NOTÍCIAS"
    fresh = _snap(_state(updated_at=_ago(STALE_S - 60)))
    assert phase_of(fresh).label == "CALCULANDO"


def test_admission_in_progress_is_named() -> None:
    snap = _snap(_state(admission="pending"))
    assert phase_of(snap).detail == "teste de admissão"
    assert "em admissão" in _text(snap)


def test_an_old_state_file_without_totals_still_renders() -> None:
    state = _state()
    del state["totals"], state["recent"]
    text = _text(_snap(state))
    assert "nenhuma ainda" in text and "CALCULANDO" in text


def test_a_new_cell_before_its_first_start_renders() -> None:
    text = _text(_snap({}, running=False, log=()), width=70)
    assert "NUNCA INICIADA" in text and "celula start" in text


@pytest.mark.parametrize("width", [70, 99, 100, 160])
def test_no_line_is_wider_than_the_terminal(width: int) -> None:
    text = _text(_snap(), width=width)
    assert max(len(line) for line in text.splitlines()) <= width


def test_the_narrow_panel_keeps_the_essentials() -> None:
    text = _text(_snap(), width=60, height=60)
    for expected in ("CALCULANDO", "cbs_04263", "00:23:41", "cbs_31877", "-87.42"):
        assert expected in text


def test_a_short_terminal_drops_the_log() -> None:
    text = _text(_snap(), height=16)
    assert "sem resposta" not in text and "Ctrl+C" in text


def test_tail_reads_only_the_end_and_survives_bad_bytes(tmp_path: Path) -> None:
    log = tmp_path / "celula.log"
    assert tail(log, 5) == ()
    lines = [f"linha {n}" for n in range(5000)]
    log.write_bytes(("\n".join(lines) + "\n").encode() + b"\xff quebrada\n")
    assert tail(log, 3) == ("linha 4998", "linha 4999", "� quebrada")


def test_read_snapshot_gathers_the_local_files(tmp_path: Path) -> None:
    config = CellConfig(
        central_url="http://100.64.0.1:8100", token="t", root=tmp_path / "c"
    )
    config.root.mkdir()
    (config.root / "state.json").write_text(json.dumps(_state()))
    outbox = Outbox(config.root)
    outbox.put("a9", "cbs_31877", {"attempt_id": "a9"})
    (config.root / "celula.log").write_text("um\ndois\n")

    snap = read_snapshot(config, running=True, now=NOW)

    assert snap.state["name"] == "pc_elias"
    assert snap.unsent == frozenset({"a9"}) and snap.dead == 0
    assert snap.log == ("um", "dois")
    assert snap.limit_s == config.crest_timeout_s + config.mopac_timeout_s


def test_a_broken_state_file_is_reported_not_raised(tmp_path: Path) -> None:
    config = CellConfig(central_url="http://h:1", token="t", root=tmp_path)
    (tmp_path / "state.json").write_text("{meio")

    snap = read_snapshot(config, running=True, now=NOW)

    assert snap.state == {} and "state.json ilegível" in snap.state_error
    assert "state.json ilegível" in _text(snap)


def test_durations_read_naturally() -> None:
    assert duration(None) == "-"
    assert duration(45) == "45 s"
    assert duration(1102) == "18 min 22 s"
    assert duration(3 * 3600 + 41 * 60) == "3 h 41 min"
    assert duration(2 * 86400 + 4 * 3600) == "2 d 4 h"


def test_the_panel_reads_what_a_real_loop_wrote(tmp_path: Path) -> None:
    import threading
    from concurrent.futures import ThreadPoolExecutor

    from semi_imperium.cell.loop import CellLoop, Timings
    from semi_imperium.cell.outbox import Backoff
    from tests.cell.conftest import (
        CELL_PROGRAMS,
        EXECUTABLES,
        Central,
        SpyClient,
        Tasks,
        cell_config,
    )

    central = Central(tmp_path)
    _cell_id, token = central.new_cell()
    tasks = Tasks(threading.Event(), stop_after=3)
    tasks.fail = {"m2"}
    config = cell_config(tmp_path, token)
    loop = CellLoop(
        config,
        SpyClient(central, token),
        cores=4,
        executables=EXECUTABLES,
        probe=lambda _exe: dict(CELL_PROGRAMS),
        pool_factory=lambda workers: ThreadPoolExecutor(max_workers=workers),
        task_fn=tasks,
        timings=Timings(
            idle_s=0.01,
            admission_s=0.01,
            heartbeat_s=0.0,
            backoff=Backoff(0.01, 0.05),
        ),
        stop=tasks.stop,
        install_signals=False,
    )
    assert loop.run() == 0

    text = _text(read_snapshot(config, running=False))

    assert "PARADA" in text and "última fase: parado" in text
    assert "2 ✓" in text and "1 ✗" in text
    assert "m0" in text and "-12.50" in text and "failed: crest timeout" in text
    assert "a enviar" in text and "✓ admitida" in text
