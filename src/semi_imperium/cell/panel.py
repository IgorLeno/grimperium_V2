"""``celula painel``: the cell at a glance, live, from local files only.

Like ``status`` it never talks to the central: it reads ``state.json``, the
outbox, the lock and the end of ``celula.log``, so it works with the network
down. It is a process of its own: closing its window (or Ctrl+C) never
touches the cell, which the scheduled task runs in another process group.

``read_snapshot`` gathers the files, ``render`` turns a snapshot into a Rich
renderable (pure, so tests render it to text) and ``run_panel`` redraws it
every second.
"""

from __future__ import annotations

import json
import os
import re
import time
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from rich import box
from rich.console import Console, Group, RenderableType
from rich.live import Live
from rich.panel import Panel
from rich.progress_bar import ProgressBar
from rich.spinner import Spinner
from rich.table import Table
from rich.text import Text

from semi_imperium.cell.config import CellConfig
from semi_imperium.cell.loop import STATE_FILENAME
from semi_imperium.cell.outbox import Outbox
from semi_imperium.training.cells import (
    ADMISSION_BLOCKED,
    ADMISSION_PASSED,
    ADMISSION_PENDING,
    ADMISSION_PREFIX,
)

LOG_FILENAME = "celula.log"

#: Below this width the panel stacks into one column.
WIDE_MIN = 100
#: Longest silence of a healthy cell: a registration retry waits up to 300 s
#: (``Backoff.cap_s``); otherwise the state is rewritten every heartbeat.
STALE_S = 400.0
_LOG_TAIL_BYTES = 16_384
_LOG_LINE = re.compile(
    r"^\d{4}-\d\d-\d\d (?P<time>\d\d:\d\d:\d\d),\d+ (?P<level>\w+) [\w.]+: (?P<msg>.*)$"
)


@dataclass(frozen=True)
class Snapshot:
    """Everything the panel shows, read at one moment."""

    now: datetime
    running: bool
    state: Mapping[str, Any]
    central_url: str
    #: CREST + MOPAC timeouts: past this a molecule is about to fail.
    limit_s: float
    #: attempt ids computed but not yet accepted by the central.
    unsent: frozenset[str] = frozenset()
    dead: int = 0
    log: tuple[str, ...] = ()
    state_error: str = ""


def tail(path: Path, lines: int, chunk: int = _LOG_TAIL_BYTES) -> tuple[str, ...]:
    """The last ``lines`` lines of ``path``; empty when it cannot be read."""
    try:
        with path.open("rb") as handle:
            size = handle.seek(0, os.SEEK_END)
            handle.seek(max(0, size - chunk))
            data = handle.read()
    except OSError:
        return ()
    text = data.decode("utf-8", errors="replace").splitlines()
    if size > chunk:
        text = text[1:]  # the first line was cut in the middle
    return tuple(text[-lines:]) if lines > 0 else ()


def read_snapshot(
    config: CellConfig,
    *,
    running: bool,
    now: datetime | None = None,
    log_lines: int = 10,
) -> Snapshot:
    root = config.root
    state: dict[str, Any] = {}
    state_error = ""
    try:
        data = json.loads((root / STATE_FILENAME).read_text(encoding="utf-8"))
        if isinstance(data, dict):
            state = data
    except FileNotFoundError:
        pass
    except (OSError, ValueError) as exc:
        state_error = f"state.json ilegível: {exc}"
    outbox = Outbox(root)
    return Snapshot(
        now=now or datetime.now(timezone.utc),
        running=running,
        state=state,
        central_url=config.central_url,
        limit_s=config.crest_timeout_s + config.mopac_timeout_s,
        unsent=frozenset(p.stem for p in outbox.pending_dir.glob("*.json")),
        dead=outbox.counts()[1],
        log=tail(root / LOG_FILENAME, log_lines),
        state_error=state_error,
    )


# ------------------------------------------------------------------ format
def _parse(stamp: Any) -> datetime | None:
    try:
        moment = datetime.fromisoformat(str(stamp))
    except (TypeError, ValueError):
        return None
    return moment if moment.tzinfo else moment.replace(tzinfo=timezone.utc)


def _seconds_since(stamp: Any, now: datetime) -> float | None:
    moment = _parse(stamp)
    return None if moment is None else max(0.0, (now - moment).total_seconds())


def duration(seconds: float | None) -> str:
    """``45 s``, ``18 min 22 s``, ``3 h 41 min``, ``2 d 4 h``."""
    if seconds is None:
        return "-"
    total = int(seconds)
    days, rest = divmod(total, 86_400)
    hours, rest = divmod(rest, 3_600)
    minutes, secs = divmod(rest, 60)
    if days:
        return f"{days} d {hours} h"
    if hours:
        return f"{hours} h {minutes:02d} min"
    if minutes:
        return f"{minutes} min {secs:02d} s"
    return f"{secs} s"


def clock(seconds: float | None) -> str:
    if seconds is None:
        return "--:--:--"
    hours, rest = divmod(int(seconds), 3_600)
    return f"{hours:02d}:{rest // 60:02d}:{rest % 60:02d}"


def _local_time(stamp: Any, pattern: str = "%H:%M") -> str:
    moment = _parse(stamp)
    return "--:--" if moment is None else moment.astimezone().strftime(pattern)


def _is_reference(attempt_id: Any) -> bool:
    return str(attempt_id or "").startswith(ADMISSION_PREFIX)


# ------------------------------------------------------------------- phase
@dataclass(frozen=True)
class Phase:
    label: str
    style: str
    detail: str = ""
    #: Rich spinner name while the cell is doing something; ``None``: a dot.
    spinner: str | None = None


_PHASES: dict[str, Phase] = {
    "calculando": Phase("CALCULANDO", "bold green", "", "dots"),
    "ociosa": Phase("OCIOSA", "bold blue", "fila vazia; pergunta de novo", "point"),
    "sem conexão": Phase("SEM CONEXÃO", "bold yellow", "tentando de novo", "dots"),
    "bloqueada": Phase(
        "BLOQUEADA", "bold red", "aguardando decisão no dashboard", "point"
    ),
    "parando": Phase("PARANDO", "bold magenta", "termina o que está rodando", "dots"),
    "recusada": Phase("RECUSADA", "bold red", "programas diferentes da central"),
    "token recusado": Phase("TOKEN RECUSADO", "bold red", "instale um pacote novo"),
}


def phase_of(snap: Snapshot) -> Phase:
    """What the cell is doing, including what the state file cannot say."""
    state = snap.state
    if not snap.running:
        if not state:
            return Phase("NUNCA INICIADA", "bold bright_black", "rode: celula start")
        last = str(state.get("phase") or "?")
        return Phase("PARADA", "bold bright_black", f"última fase: {last}")
    silence = _seconds_since(state.get("updated_at"), snap.now)
    if silence is None or silence > STALE_S:
        return Phase(
            "SEM NOTÍCIAS",
            "bold yellow blink",
            f"estado sem mudança há {duration(silence)}",
        )
    phase = str(state.get("phase") or "")
    known = _PHASES.get(phase)
    if known is not None:
        if phase == "calculando" and state.get("admission") == ADMISSION_PENDING:
            return Phase(known.label, known.style, "teste de admissão", known.spinner)
        return known
    return Phase(phase.upper() or "?", "bold cyan", "", "dots")


def _admission(value: Any) -> Text:
    if value == ADMISSION_PASSED:
        return Text("✓ admitida", style="green")
    if value == ADMISSION_BLOCKED:
        return Text("✗ bloqueada", style="red")
    if value == ADMISSION_PENDING:
        return Text("em admissão", style="yellow")
    return Text("?", style="bright_black")


def _contact(snap: Snapshot) -> Text:
    since = _seconds_since(snap.state.get("last_contact_at"), snap.now)
    if since is None:
        return Text("○ nunca", style="bright_black")
    style = "green" if since < 120 else "yellow" if since < STALE_S else "red"
    return Text.assemble(("● ", style), f"há {duration(since)}")


# ------------------------------------------------------------------ render
@dataclass
class Spinners:
    """Spinners kept between frames, so redrawing does not restart them."""

    cache: dict[str, Spinner] = field(default_factory=dict)

    def get(self, name: str, text: Text, style: str) -> Spinner:
        spinner = self.cache.get(name)
        if spinner is None:
            spinner = self.cache[name] = Spinner(name, text=text, style=style)
        else:
            spinner.update(text=text, style=style)
        return spinner


def _phase_badge(phase: Phase, spinners: Spinners) -> RenderableType:
    text = Text(f" {phase.label}", style=phase.style)
    if phase.spinner is None:
        return Text.assemble(("●", phase.style), text)
    return spinners.get(phase.spinner, text, phase.style)


def _header(snap: Snapshot, wide: bool, spinners: Spinners) -> Panel:
    state = snap.state
    phase = phase_of(snap)
    name = state.get("name") or "?"
    contact = Text.assemble("central ", _contact(snap))
    admission = Text.assemble("admissão ", _admission(state.get("admission")))
    up = (
        _seconds_since(state.get("started_at"), snap.now)
        if snap.running and state
        else None
    )
    grid = Table.grid(padding=(0, 3), expand=True)
    detail = Text(phase.detail, style="italic bright_black")
    if wide:
        grid.add_column(no_wrap=True)
        grid.add_column(ratio=1)
        grid.add_column(no_wrap=True)
        grid.add_column(no_wrap=True)
        grid.add_column(no_wrap=True, justify="right")
        grid.add_row(
            _phase_badge(phase, spinners),
            detail,
            admission,
            contact,
            Text(f"no ar há {duration(up)}" if up is not None else ""),
        )
    else:
        grid.add_column()
        grid.add_row(_phase_badge(phase, spinners))
        if phase.detail:
            grid.add_row(detail)
        grid.add_row(Text.assemble(admission, " · ", contact))
    rows: list[RenderableType] = [grid]
    if state.get("last_error"):
        rows.append(Text(f"último erro: {state['last_error']}", style="yellow"))
    if snap.state_error:
        rows.append(Text(snap.state_error, style="red"))
    subtitle = " · ".join(
        str(part) for part in (state.get("cell_id"), state.get("host")) if part
    )
    return Panel(
        Group(*rows),
        title=f"[bold]Grimperium[/] · célula [bold cyan]{name}[/]",
        title_align="left",
        subtitle=subtitle or None,
        subtitle_align="right",
        border_style=phase.style.replace("bold ", "").replace(" blink", ""),
        box=box.ROUNDED,
    )


def _bar(snap: Snapshot, elapsed: float | None, width: int) -> Table:
    fraction = (elapsed or 0.0) / snap.limit_s if snap.limit_s > 0 else 0.0
    style = "green" if fraction < 0.75 else "yellow" if fraction < 1.0 else "red"
    bar = ProgressBar(
        total=1.0,
        completed=min(fraction, 1.0),
        width=width,
        complete_style=style,
        finished_style="red",
    )
    grid = Table.grid(padding=(0, 1))
    grid.add_row(bar, Text(f"{fraction * 100:3.0f} %", style=style))
    return grid


def _workers(snap: Snapshot, wide: bool) -> Panel:
    state = snap.state
    busy = [b for b in state.get("busy") or [] if isinstance(b, dict)]
    slots = max(int(state.get("workers") or 0), len(busy))
    table = Table.grid(padding=(0, 2), expand=True)
    table.add_column(justify="right", style="bright_black", no_wrap=True)
    table.add_column(no_wrap=True)
    if wide:
        table.add_column(ratio=1, no_wrap=True, overflow="ellipsis")
    table.add_column(no_wrap=True)
    table.add_column(no_wrap=True)
    for index in range(1, slots + 1):
        if index > len(busy):
            free = Text("livre", style="bright_black")
            cells: list[RenderableType] = [str(index), free, Text("")]
            table.add_row(*cells, *([Text("")] if wide else []))
            continue
        worker = busy[index - 1]
        elapsed = _seconds_since(worker.get("started_at"), snap.now)
        mol = Text(str(worker.get("mol_id") or "?"), style="bold")
        if _is_reference(worker.get("attempt_id")):
            mol.append(" REF", style="black on yellow")
        smiles = Text(str(worker.get("smiles") or ""), style="cyan")
        timer = Text(clock(elapsed), style="bold")
        if wide:
            table.add_row(str(index), mol, smiles, timer, _bar(snap, elapsed, width=24))
        else:
            table.add_row(str(index), mol, timer, _bar(snap, elapsed, width=10))
            smiles.overflow = "ellipsis"
            smiles.no_wrap = True
            table.add_row("", smiles, "", "")
    if not slots:
        table.add_row("", Text("nenhum worker ainda", style="bright_black"), "", "")
    threads = state.get("crest_threads")
    title = f"Workers ({slots} × {threads} threads)" if threads else "Workers"
    return Panel(table, title=title, title_align="left", box=box.ROUNDED)


def _totals(snap: Snapshot) -> Table:
    totals = snap.state.get("totals")
    grid = Table.grid(padding=(0, 2))
    grid.add_column(style="bright_black")
    grid.add_column(justify="right", style="bold green")
    grid.add_column(justify="right", style="bold red")
    if isinstance(totals, dict):
        grid.add_row(
            "sessão",
            f"{totals.get('session_ok', 0)} ✓",
            f"{totals.get('session_failed', 0)} ✗",
        )
        grid.add_row(
            "hoje", f"{totals.get('day_ok', 0)} ✓", f"{totals.get('day_failed', 0)} ✗"
        )
    else:
        grid.add_row("sessão", "-", "")
        grid.add_row("hoje", "-", "")
    grid.add_row("", "", "")
    unsent = len(snap.unsent)
    grid.add_row("a enviar", Text(str(unsent), style="yellow" if unsent else "dim"), "")
    grid.add_row(
        "dead-letter", Text(str(snap.dead), style="red" if snap.dead else "dim"), ""
    )
    return grid


def _recent(snap: Snapshot, wide: bool, limit: int) -> Table:
    table = Table.grid(padding=(0, 2), expand=True)
    table.add_column(no_wrap=True)
    table.add_column(no_wrap=True, style="bright_black")
    table.add_column(no_wrap=True, style="bold")
    table.add_column(ratio=1, no_wrap=True, overflow="ellipsis")
    if wide:
        table.add_column(no_wrap=True, justify="right")
    table.add_column(no_wrap=True)
    recent = [r for r in snap.state.get("recent") or [] if isinstance(r, dict)]
    for entry in recent[:limit]:
        ok = bool(entry.get("ok"))
        mark = Text("✓", style="bold green") if ok else Text("✗", style="bold red")
        hof = entry.get("hof_kcal_mol")
        if ok and isinstance(hof, int | float):
            what = Text(
                f"ΔHf {hof:8.2f} kcal/mol" if wide else f"{hof:8.2f}", style="cyan"
            )
        else:
            what = Text(str(entry.get("error") or "falha"), style="red")
        flags = Text()
        if str(entry.get("attempt_id")) in snap.unsent:
            flags.append("⇡ ", style="yellow")
        if _is_reference(entry.get("attempt_id")):
            flags.append("REF", style="black on yellow")
        took = duration(entry.get("elapsed_s")) if ok else ""
        row: list[RenderableType] = [
            mark,
            _local_time(entry.get("at")),
            str(entry.get("mol_id") or "?"),
            what,
        ]
        if wide:
            row.append(took)
        row.append(flags)
        table.add_row(*row)
    if not recent:
        table.add_row("", "", Text("nenhuma ainda", style="bright_black"), "", "")
    return table


def _log(snap: Snapshot, lines: int) -> Text:
    text = Text(no_wrap=True, overflow="ellipsis")
    for raw in snap.log[-lines:]:
        match = _LOG_LINE.match(raw)
        if match is None:
            text.append(raw + "\n", style="bright_black")
            continue
        level = match["level"]
        style = {"WARNING": "yellow", "ERROR": "bold red"}.get(level, "")
        text.append(match["time"] + " ", style="bright_black")
        if level != "INFO":
            text.append(level + " ", style=style)
        text.append(match["msg"] + "\n", style=style)
    text.rstrip()
    return text


def render(
    snap: Snapshot, width: int, height: int, spinners: Spinners | None = None
) -> RenderableType:
    """The whole panel for a terminal of ``width`` × ``height``."""
    spinners = spinners or Spinners()
    wide = width >= WIDE_MIN
    state = snap.state
    slots = max(int(state.get("workers") or 0), len(state.get("busy") or []), 1)
    updated = _local_time(snap.now.isoformat(), "%H:%M:%S")
    footer = (
        f"{snap.central_url} · atualizado {updated} · Ctrl+C fecha o painel"
        if wide
        else f"{snap.central_url} · {updated} · Ctrl+C fecha"
    )
    parts: list[RenderableType] = [
        _header(snap, wide, spinners),
        _workers(snap, wide),
    ]
    if wide:
        recent_rows = 5
        middle = Table.grid(padding=(0, 1), expand=True)
        middle.add_column(ratio=1)
        middle.add_column(ratio=2)
        # Same height side by side: the totals take five lines.
        middle.add_row(
            Panel(_totals(snap), title="Concluídas", title_align="left", height=7),
            Panel(
                _recent(snap, wide, recent_rows),
                title="Últimos resultados",
                title_align="left",
                height=7,
            ),
        )
        parts.append(middle)
        used = 3 + (2 + slots) + (2 + recent_rows)
    else:
        recent_rows = 5
        parts.append(
            Panel(
                Group(_totals(snap), Text(""), _recent(snap, wide, recent_rows)),
                title="Concluídas",
                title_align="left",
            )
        )
        used = 5 + (2 + 2 * slots) + (2 + 6 + recent_rows)
    log_lines = min(10, height - used - 2)
    if log_lines >= 2 and snap.log:
        parts.append(
            Panel(
                _log(snap, log_lines),
                title=LOG_FILENAME,
                title_align="left",
                subtitle=footer,
                subtitle_align="right",
                border_style="bright_black",
            )
        )
    else:
        parts.append(Text(footer, style="bright_black", justify="right"))
    return Group(*parts)


# --------------------------------------------------------------------- run
def run_panel(
    config: CellConfig,
    is_running: Callable[[], bool],
    *,
    once: bool = False,
    console: Console | None = None,
) -> int:
    """Draw the panel until Ctrl+C (or once); return the exit code."""
    console = console or Console()
    spinners = Spinners()

    def frame() -> RenderableType:
        width, height = console.size
        snap = read_snapshot(config, running=is_running())
        return render(snap, width, height, spinners)

    # Rich draws nothing live on a pipe or a TERM=dumb terminal: one frame.
    if once or not console.is_interactive:
        console.print(frame())
        return 0
    try:
        with Live(frame(), console=console, screen=True, refresh_per_second=8) as live:
            while True:
                time.sleep(1.0)
                live.update(frame())
    except KeyboardInterrupt:
        pass
    return 0
