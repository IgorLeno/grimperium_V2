"""Command line of the training routine: ``python -m semi_imperium.training``."""

from __future__ import annotations

import argparse
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from semi_imperium.training.db import database_path
from semi_imperium.training.job_state import JobAlreadyRunning, JobState
from semi_imperium.training.ledger import NotMigratedError, read_progress
from semi_imperium.training.migration import migrate
from semi_imperium.training.runner import RunSettings, run_training

# src/semi_imperium/training/__main__.py -> repository root
_ROOT = Path(__file__).resolve().parents[3]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m semi_imperium.training",
        description="Delta Entalpia: rotina de treinamento (Semi-Imperium, protocolo 2).",
    )
    parser.add_argument("--csv", type=Path, default=_ROOT / "data" / "thermo_pm7.csv")
    parser.add_argument("--store", type=Path, default=_ROOT / "data" / "semi_imperium")
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "semi_imperium_training",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("status", help="mostra o progresso e o estado do job")
    mig = commands.add_parser(
        "migrate", help="marca o protocolo 1 e devolve Selected à fila"
    )
    mig.add_argument("--apply", action="store_true", help="grava (padrão: ensaio)")
    run = commands.add_parser("run", help="processa a fila com o protocolo 2")
    run.add_argument("--workers", type=int, default=3)
    run.add_argument("--limit", type=int, default=None)
    run.add_argument("--crest-threads", type=int, default=4)
    run.add_argument("--crest-timeout", type=float, default=3600.0)
    run.add_argument("--mopac-timeout", type=float, default=900.0)
    run.add_argument("--keep-artifacts", action="store_true")
    commands.add_parser("stop", help="pede ao job em execução que pare")
    return parser


def _status(args: argparse.Namespace) -> int:
    progress = read_progress(args.csv)
    total = progress.total

    def line(label: str, count: int) -> str:
        share = 0.0 if total == 0 else 100.0 * count / total
        return f"  {label}: {count} de {total} ({share:.1f}%)"

    print(f"Progresso de {args.csv}")
    print(line("protocolo 2 (canônico)", progress.protocol2))
    print(line("protocolo 1 (preliminar)", progress.protocol1))
    print(line("na fila (Pending)", progress.pending))
    print(line("selecionadas (Selected)", progress.selected))
    print(line("para repetir (Rerun)", progress.rerun))
    print(line("em execução (Running)", progress.running))
    print(line("descartadas (Skip)", progress.skipped))
    if not progress.migrated:
        print("  CSV ainda não migrado: rode 'migrate --apply'.")
    job = JobState(database_path(args.store)).snapshot()
    label = {"idle": "ocioso"}.get(job.state, job.state)
    print(
        f"Job: {label} (feitas {job.done}, falhas {job.failed}) {job.message}".rstrip()
    )
    return 0


def _migrate(args: argparse.Namespace) -> int:
    report = migrate(args.csv, apply=args.apply)
    mode = "aplicado" if report.applied else "ensaio (nada foi gravado)"
    print(f"Migração {mode}: {report.rows} linhas")
    print(f"  marcadas como protocolo 1: {report.legacy_marked}")
    print(f"  Selected devolvidas a Pending: {report.selected_reset}")
    print(f"  já marcadas: {report.already_marked}")
    if report.backup_path is not None:
        print(f"  cópia de segurança: {report.backup_path}")
    return 0


def _run(args: argparse.Namespace) -> int:
    settings = RunSettings(
        csv_path=args.csv,
        store_root=args.store,
        work_dir=args.work_dir,
        workers=args.workers,
        limit=args.limit,
        crest_threads=args.crest_threads,
        crest_timeout_s=args.crest_timeout,
        mopac_timeout_s=args.mopac_timeout,
        keep_artifacts=args.keep_artifacts,
    )
    try:
        summary = run_training(settings)
    except NotMigratedError as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return 2
    except JobAlreadyRunning as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return 3
    print(
        f"Concluído: {summary.attempted} tentadas, {summary.ok} ok, "
        f"{summary.failed} falhas{' (parado a pedido)' if summary.stopped else ''}"
    )
    return 0


def _stop(args: argparse.Namespace) -> int:
    if JobState(database_path(args.store)).request_stop():
        print("Parada solicitada: as moléculas em andamento terminam e o job encerra.")
        return 0
    print("Nenhum job em execução.", file=sys.stderr)
    return 1


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point; returns the process exit code."""
    args = _build_parser().parse_args(argv)
    handlers = {"status": _status, "migrate": _migrate, "run": _run, "stop": _stop}
    return handlers[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
