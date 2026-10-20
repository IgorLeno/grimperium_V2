"""Command line of the training routine: ``python -m semi_imperium.training``."""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from semi_imperium.training.db import database_path
from semi_imperium.training.job_state import JobAlreadyRunning, JobState
from semi_imperium.training.molecule_queue import CbsReadError, read_progress
from semi_imperium.training.runner import RunSettings, run_training

# src/semi_imperium/training/__main__.py -> repository root
_ROOT = Path(__file__).resolve().parents[3]

# Fixed environment from environment.yml / conda-lock.yml. Using it by default
# keeps a `dnf upgrade` from swapping CREST, xtb or MOPAC under a running queue.
FIXED_ENV = "grimperium-quimica"
_PROGRAMS = ("crest", "xtb", "mopac")


def _fixed_env_bin() -> Path:
    root = os.environ.get("MAMBA_ROOT_PREFIX") or str(Path.home() / "micromamba")
    return Path(root) / "envs" / FIXED_ENV / "bin"


def _executables(bin_dir: Path | None) -> dict[str, str]:
    """Programs of ``bin_dir``, else of the fixed environment, else of PATH.

    An explicit ``bin_dir`` missing a program is an error; only the implicit
    default may fall back to PATH, and it says so on stderr.
    """
    if bin_dir is not None:
        missing = [name for name in _PROGRAMS if not (bin_dir / name).is_file()]
        if missing:
            raise FileNotFoundError(f"{bin_dir} não tem: {', '.join(missing)}")
        return {name: str(bin_dir / name) for name in _PROGRAMS}
    fixed = _fixed_env_bin()
    if all((fixed / name).is_file() for name in _PROGRAMS):
        return {name: str(fixed / name) for name in _PROGRAMS}
    print(
        f"aviso: ambiente {FIXED_ENV} não encontrado em {fixed}; "
        "usando crest/xtb/mopac do PATH (versões não fixadas)",
        file=sys.stderr,
    )
    return {name: name for name in _PROGRAMS}


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m semi_imperium.training",
        description="Delta Entalpia: rotina de treinamento (Semi-Imperium, protocolo 2).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=_ROOT / "data" / "cbs_referencia.csv",
        help="referência CBS (somente leitura) de onde a fila é semeada",
    )
    parser.add_argument("--store", type=Path, default=_ROOT / "data" / "semi_imperium")
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "semi_imperium_training",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("status", help="mostra o progresso e o estado do job")
    run = commands.add_parser("run", help="processa a fila com o protocolo 2")
    run.add_argument("--workers", type=int, default=3)
    run.add_argument("--limit", type=int, default=None)
    run.add_argument("--crest-threads", type=int, default=4)
    run.add_argument("--crest-timeout", type=float, default=3600.0)
    run.add_argument("--mopac-timeout", type=float, default=900.0)
    run.add_argument("--keep-artifacts", action="store_true")
    run.add_argument(
        "--bin-dir",
        type=Path,
        default=None,
        help=f"pasta com crest, xtb e mopac (padrão: ambiente {FIXED_ENV})",
    )
    commands.add_parser("stop", help="pede ao job em execução que pare")
    return parser


def _status(args: argparse.Namespace) -> int:
    try:
        progress = read_progress(database_path(args.store), args.csv)
    except (FileNotFoundError, CbsReadError) as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return 2
    total = progress.total

    def line(label: str, count: int) -> str:
        share = 0.0 if total == 0 else 100.0 * count / total
        return f"  {label}: {count} de {total} ({share:.1f}%)"

    print(f"Fila de treinamento (protocolo 2) a partir de {args.csv}")
    print(line("concluídas (ok)", progress.ok))
    print(line("na fila (pending)", progress.pending))
    print(line("para repetir (rerun)", progress.rerun))
    print(line("em execução (running)", progress.running))
    print(line("descartadas (skip)", progress.skipped))
    job = JobState(database_path(args.store)).snapshot()
    label = {"idle": "ocioso"}.get(job.state, job.state)
    print(
        f"Job: {label} (feitas {job.done}, falhas {job.failed}) {job.message}".rstrip()
    )
    return 0


def _run(args: argparse.Namespace) -> int:
    try:
        executables = _executables(args.bin_dir)
    except FileNotFoundError as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return 2
    settings = RunSettings(
        cbs_csv=args.csv,
        store_root=args.store,
        work_dir=args.work_dir,
        workers=args.workers,
        limit=args.limit,
        crest_threads=args.crest_threads,
        crest_timeout_s=args.crest_timeout,
        mopac_timeout_s=args.mopac_timeout,
        keep_artifacts=args.keep_artifacts,
        crest_executable=executables["crest"],
        xtb_executable=executables["xtb"],
        mopac_executable=executables["mopac"],
    )
    try:
        summary = run_training(settings)
    except (FileNotFoundError, CbsReadError) as exc:
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
    handlers = {"status": _status, "run": _run, "stop": _stop}
    return handlers[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
