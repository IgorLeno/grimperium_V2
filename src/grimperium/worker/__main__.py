"""CLI entry point for the Grimperium distributed worker.

Invoked as either ``python -m grimperium.worker`` or ``grimperium-worker``
(after ``poetry install``).

Connection behaviour (Bloco 1 of distributed-mode refactor):
  - Calls POST /register exactly once per startup attempt with a 10-second timeout.
  - On success: reads crest/mopac timeouts and batch_size from the server response
    (config comes from the server, not from CLI flags).
  - On failure: if stdin is a TTY, shows a Rich retry menu; otherwise exits with
    code 1 immediately (background/piped mode).
"""

from __future__ import annotations

import argparse
import logging
import math
import socket
import sys
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.panel import Panel

from grimperium.utils.logging import setup_logging
from grimperium.worker.client import WorkerClient, WorkerClientConfig
from grimperium.worker.runner import WorkerConfig, WorkerRunner

LOG = logging.getLogger(__name__)
CONSOLE = Console()

_REPO_ROOT = Path(__file__).resolve().parents[3]
_POLL_INTERVAL_S = 5.0
_REGISTER_TIMEOUT_S = 10.0


def build_parser() -> argparse.ArgumentParser:
    """Build the argparse parser. Isolated for unit testing."""
    parser = argparse.ArgumentParser(
        prog="grimperium-worker",
        description="Grimperium distributed worker — claims molecules from the server, "
        "runs CREST+MOPAC locally, reports results.",
    )
    parser.add_argument(
        "--server-url",
        required=True,
        help="Base URL of the Grimperium server, e.g. http://192.168.31.186:8000",
    )
    parser.add_argument(
        "--worker-id",
        default=socket.gethostname(),
        help="Identifier reported to the server (default: this host's hostname).",
    )
    parser.add_argument(
        "--api-token",
        default="",
        help="Optional X-Token auth header value.",
    )
    parser.add_argument(
        "--max-molecules",
        type=int,
        default=0,
        help="Maximum molecules to process before exiting (0 = unlimited).",
    )
    parser.add_argument(
        "--idle-stop",
        type=int,
        default=300,
        help="Seconds of empty-queue polling before the worker exits (default: 300).",
    )
    return parser


def _try_register(client: WorkerClient) -> dict[str, Any] | None:
    """Attempt a single /register call.

    Returns the server response dict on success, or None on any failure
    (ServerError, transport error, timeout, or any other exception).
    """
    try:
        return client.register()
    except Exception:
        return None


def _show_retry_menu(server_url: str) -> bool:
    """Show an interactive retry-or-exit menu. Returns True to retry.

    Only called when stdin is a TTY (interactive session).
    """
    try:
        import questionary
    except ImportError:
        CONSOLE.print("[red]questionary not installed — exiting[/red]")
        return False

    CONSOLE.print()
    CONSOLE.print(
        Panel(
            f"[red]✗ Não foi possível conectar ao servidor[/red]\n"
            f"  [dim]{server_url}[/dim]",
            border_style="red",
            padding=(1, 2),
        )
    )

    choice: str | None = questionary.select(
        "O que deseja fazer?",
        choices=["Tentar novamente", "Encerrar"],
        pointer="❯",
    ).ask()

    return choice == "Tentar novamente"


def _connect_with_retry(client: WorkerClient, server_url: str) -> dict[str, Any]:
    """Register once; on failure show interactive menu (TTY) or exit (non-TTY).

    Returns the server response dict on success.
    Calls sys.exit(1) on permanent failure.
    """
    while True:
        result = _try_register(client)
        if result is not None:
            return result

        interactive = sys.stdin.isatty()
        if not interactive:
            LOG.error(
                "Cannot connect to server %s — exiting (non-interactive mode)",
                server_url,
            )
            sys.exit(1)

        retry = _show_retry_menu(server_url)
        if not retry:
            sys.exit(0)


def main(argv: list[str] | None = None) -> int:
    """Entry point. Returns a shell exit code."""
    args = build_parser().parse_args(argv)

    setup_logging(
        level="INFO",
        log_file=_REPO_ROOT / "logs" / "worker.log",
        console_level="INFO",
    )

    client_cfg = WorkerClientConfig(
        server_url=args.server_url,
        worker_id=args.worker_id,
        api_token=args.api_token,
    )
    client = WorkerClient(client_cfg)

    # Single register attempt; interactive retry on failure
    server_cfg = _connect_with_retry(client, args.server_url)

    # Config comes from server — not from CLI flags
    crest_minutes = int(server_cfg.get("crest_timeout_minutes", 60))
    mopac_minutes = int(server_cfg.get("mopac_timeout_minutes", 30))
    batch_size = int(server_cfg.get("batch_size", 10))

    config = WorkerConfig(
        server_url=args.server_url,
        worker_id=args.worker_id,
        api_token=args.api_token,
        poll_interval_s=_POLL_INTERVAL_S,
        max_idle_polls=max(1, math.ceil(args.idle_stop / _POLL_INTERVAL_S)),
        crest_timeout_minutes=crest_minutes,
        mopac_timeout_minutes=mopac_minutes,
        batch_size=batch_size,
    )
    runner = WorkerRunner(config=config, client=client)

    max_mols = args.max_molecules or None
    LOG.info(
        "Worker %s connected to %s (crest=%dm mopac=%dm batch=%d idle_stop=%ss)",
        args.worker_id,
        args.server_url,
        crest_minutes,
        mopac_minutes,
        batch_size,
        args.idle_stop,
    )
    CONSOLE.print(
        f"[green]✓ Conectado a {args.server_url}[/green] — "
        f"CREST {crest_minutes}m / MOPAC {mopac_minutes}m / batch {batch_size}"
    )

    try:
        runner.run(max_molecules=max_mols)
        return 0
    except KeyboardInterrupt:
        LOG.info("KeyboardInterrupt received — stopping worker gracefully.")
        runner.stop()
        return 0
    except Exception:
        LOG.exception("Worker exited with an unexpected error.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
