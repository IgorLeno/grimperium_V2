"""Command line of a cell: ``celula start | status | painel | update``.

``start`` runs in the foreground (the scheduled task of the Windows host
calls it inside WSL); ``status`` and ``painel`` read only local files, so
they answer even with the central unreachable (``painel`` redraws itself
every second until Ctrl+C, which never reaches the cell); ``update``
installs a new wheel while the cell is stopped. One cell per folder:
``start`` holds ``celula.lock``.
"""

from __future__ import annotations

import argparse
import fcntl
import logging
import subprocess
import sys
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import IO

from semi_imperium.cell.client import CentralClient
from semi_imperium.cell.config import (
    CONFIG_FILENAME,
    CellConfig,
    ConfigError,
    default_home,
    load_config,
)
from semi_imperium.cell.loop import CellLoop, read_state
from semi_imperium.cell.outbox import Outbox
from semi_imperium.cell.panel import LOG_FILENAME, run_panel

LOCK_FILENAME = "celula.lock"

EXIT_CONFIG = 2
EXIT_RUNNING = 3

Pip = Callable[[list[str]], int]


class AlreadyRunning(RuntimeError):
    pass


@contextmanager
def _lock(root: Path) -> Iterator[IO[str]]:
    """Hold the folder's lock; ``AlreadyRunning`` if a cell holds it."""
    root.mkdir(parents=True, exist_ok=True)
    with (root / LOCK_FILENAME).open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise AlreadyRunning(f"uma célula já está rodando em {root}") from None
        try:
            yield handle
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def _is_running(root: Path) -> bool:
    try:
        with _lock(root):
            return False
    except AlreadyRunning:
        return True


def _logging(root: Path) -> None:
    handler = RotatingFileHandler(
        root / LOG_FILENAME, maxBytes=5_000_000, backupCount=3, encoding="utf-8"
    )
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    )
    console = logging.StreamHandler()
    console.setFormatter(logging.Formatter("%(asctime)s %(message)s", "%H:%M:%S"))
    logger = logging.getLogger("semi_imperium.cell")
    logger.setLevel(logging.INFO)
    for existing in list(logger.handlers):
        logger.removeHandler(existing)
        existing.close()
    logger.addHandler(handler)
    logger.addHandler(console)
    # One line per request every 30 s says nothing the state file does not.
    logging.getLogger("httpx").setLevel(logging.WARNING)


def _start(config: CellConfig) -> int:
    try:
        with _lock(config.root):
            _logging(config.root)
            client = CentralClient(config.central_url, config.token)
            try:
                return CellLoop(config, client).run()
            finally:
                client.close()
    except AlreadyRunning as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return EXIT_RUNNING


def _status(config: CellConfig) -> int:
    root = config.root
    state = read_state(root)
    pending, dead = Outbox(root).counts()
    print(f"Célula em {root}")
    print(f"  central: {config.central_url}")
    print("  token: presente")
    print(f"  rodando: {'sim' if _is_running(root) else 'não'}")
    if not state:
        print("  estado: nunca iniciada")
    else:
        name = state.get("name") or "?"
        print(f"  estado: {state.get('phase')} (atualizado {state.get('updated_at')})")
        print(
            f"  registro: {name} [{state.get('cell_id') or '?'}], "
            f"admissão {state.get('admission') or '?'}"
        )
        print(
            f"  workers: {state.get('workers')} × {state.get('crest_threads')} threads"
        )
        print(f"  último contato com a central: {state.get('last_contact_at') or '-'}")
        totals = state.get("totals")
        if isinstance(totals, dict):
            print(
                f"  concluídas: sessão {totals.get('session_ok', 0)} ok, "
                f"{totals.get('session_failed', 0)} falhas; "
                f"hoje {totals.get('day_ok', 0)} ok, "
                f"{totals.get('day_failed', 0)} falhas"
            )
        if state.get("last_error"):
            print(f"  último erro: {state['last_error']}")
        for busy in state.get("busy") or []:
            print(
                f"    {busy.get('mol_id')} {busy.get('smiles')} "
                f"desde {busy.get('started_at')}"
            )
    print(f"  resultados a enviar: {pending}")
    print(f"  dead-letter: {dead}")
    for letter in Outbox(root).dead_letters():
        print(f"    {letter.mol_id} ({letter.status}): {letter.detail}")
    return 0


def _pip(args: list[str]) -> int:
    return subprocess.run(args, check=False).returncode  # noqa: S603


def _update(config: CellConfig, wheel: Path, pip: Pip) -> int:
    if wheel.suffix != ".whl" or not wheel.is_file():
        print(f"erro: {wheel} não é um arquivo .whl", file=sys.stderr)
        return EXIT_CONFIG
    try:
        with _lock(config.root):
            # Only the package: its dependencies come with the environment
            # of the cell package (Lote 4), not from the internet here.
            code = pip(
                [
                    sys.executable,
                    "-m",
                    "pip",
                    "install",
                    "--no-deps",
                    "--force-reinstall",
                    str(wheel),
                ]
            )
    except AlreadyRunning as exc:
        print(f"erro: {exc}; pare a célula antes de atualizar", file=sys.stderr)
        return EXIT_RUNNING
    if code != 0:
        print(f"erro: pip terminou com código {code}", file=sys.stderr)
        return 1
    print(f"Atualizado com {wheel.name}. Inicie de novo: celula start")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="celula",
        description="Célula do Grimperium: calcula moléculas do protocolo 2 "
        "para a central.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_home() / CONFIG_FILENAME,
        help="arquivo de configuração (padrão: ~/grimperium-celula/config.json)",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("start", help="registra na central e calcula até parar")
    commands.add_parser("status", help="mostra o estado local (sem rede)")
    painel = commands.add_parser(
        "painel", help="painel ao vivo da célula (sem rede; Ctrl+C fecha só ele)"
    )
    painel.add_argument(
        "--uma-vez",
        action="store_true",
        dest="once",
        help="desenha uma vez e sai",
    )
    update = commands.add_parser("update", help="instala um wheel novo")
    update.add_argument("wheel", type=Path)
    return parser


def main(argv: Sequence[str] | None = None, *, pip: Pip = _pip) -> int:
    """Entry point; returns the process exit code."""
    args = _build_parser().parse_args(argv)
    try:
        config = load_config(args.config)
    except ConfigError as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return EXIT_CONFIG
    if args.command == "start":
        return _start(config)
    if args.command == "status":
        return _status(config)
    if args.command == "painel":
        return run_panel(config, lambda: _is_running(config.root), once=args.once)
    return _update(config, args.wheel, pip)


if __name__ == "__main__":
    raise SystemExit(main())
