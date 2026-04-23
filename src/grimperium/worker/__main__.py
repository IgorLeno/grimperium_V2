"""CLI entry point for the Grimperium distributed worker.

Invoked as either ``python -m grimperium.worker`` or ``grimperium-worker``
(after ``poetry install``). Parses arguments, wires them to ``WorkerConfig``
and ``WorkerClient``, and runs ``WorkerRunner.run()`` until the queue is
idle, ``--max-molecules`` is reached, or the user sends SIGINT.
"""

from __future__ import annotations

import argparse
import logging
import math
import socket
import sys
from pathlib import Path

from grimperium.utils.logging import setup_logging
from grimperium.worker.client import WorkerClient, WorkerClientConfig
from grimperium.worker.runner import WorkerConfig, WorkerRunner

LOG = logging.getLogger(__name__)

# Repo root is 3 levels up from this file:
# src/grimperium/worker/__main__.py -> parents[3] == repo root.
# Resolved directly from __file__ so the worker runs identically on the
# primary host and on remote worker machines (e.g. ~/grimperium_V2) without
# depending on grimperium.cli.constants.
_REPO_ROOT = Path(__file__).resolve().parents[3]
_POLL_INTERVAL_S = 5.0


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
    parser.add_argument(
        "--crest-timeout",
        type=int,
        default=3600,
        help="Per-molecule CREST timeout in seconds (default: 3600).",
    )
    parser.add_argument(
        "--mopac-timeout",
        type=int,
        default=1800,
        help="Per-molecule MOPAC timeout in seconds (default: 1800).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point. Returns a shell exit code."""
    args = build_parser().parse_args(argv)

    setup_logging(
        level="INFO",
        log_file=_REPO_ROOT / "logs" / "worker.log",
        console_level="INFO",
    )

    config = WorkerConfig(
        server_url=args.server_url,
        worker_id=args.worker_id,
        api_token=args.api_token,
        poll_interval_s=_POLL_INTERVAL_S,
        max_idle_polls=max(1, math.ceil(args.idle_stop / _POLL_INTERVAL_S)),
        crest_timeout_minutes=max(1, args.crest_timeout // 60),
        mopac_timeout_minutes=max(1, args.mopac_timeout // 60),
    )
    client = WorkerClient(
        WorkerClientConfig(
            server_url=args.server_url,
            worker_id=args.worker_id,
            api_token=args.api_token,
        )
    )
    runner = WorkerRunner(config=config, client=client)

    max_mols = args.max_molecules or None
    LOG.info(
        "Starting worker %s against %s (max_molecules=%s, idle_stop=%ss)",
        args.worker_id,
        args.server_url,
        max_mols,
        args.idle_stop,
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
