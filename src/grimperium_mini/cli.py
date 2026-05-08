"""Command line interface for grimperium_mini."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import MiniConfig
from .io import read_csv_rows, validate_workbook, write_simple_xlsx
from .logging_utils import configure_logging
from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="grimperium_mini")
    parser.add_argument("--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate-data")
    validate.add_argument("--xlsx", type=Path, required=True)

    run = sub.add_parser("run")
    run.add_argument("--xlsx", type=Path, required=True)
    run.add_argument("--methods", nargs="+", default=["AM1", "PM3", "PM7"])
    run.add_argument("--limit", type=int)
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--work-root", type=Path, default=Path("runs"))
    run.add_argument("--results-dir", type=Path, default=Path("results"))
    run.add_argument("--crest-executable", default=None)
    run.add_argument("--mopac-executable", default=None)
    run.add_argument("--xtb-executable", default=None)
    run.add_argument("--no-xtb", action="store_true", help="Skip xTB pre-optimization")
    run.add_argument("--threads", type=int, default=None)
    run.add_argument("--crest-method", default=None)
    run.add_argument("--crest-ewin", type=float, default=None)
    run.add_argument("--crest-rthr", type=float, default=None)
    run.add_argument(
        "--crest-quick-mode", choices=["quick", "squick", "mquick"], default=None
    )
    run.add_argument("--max-conformers-to-optimize", type=int, default=None)

    export = sub.add_parser("export-summary")
    export.add_argument("--input", type=Path, required=True)
    export.add_argument("--reactions", type=Path, required=True)
    export.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.verbose)

    if args.command == "validate-data":
        issues = validate_workbook(args.xlsx)
        if issues:
            print("INVALID")
            for issue in issues:
                print(f"- {issue}")
            return 1
        print(f"OK: {args.xlsx}")
        return 0

    if args.command == "run":
        config = _config_from_args(args)
        run_pipeline(
            args.xlsx,
            methods=[m.upper() for m in args.methods],
            config=config,
            limit=args.limit,
            dry_run=args.dry_run,
        )
        print(f"Wrote results to {config.results_dir}")
        return 0

    if args.command == "export-summary":
        formation = read_csv_rows(args.input)
        reactions = read_csv_rows(args.reactions)
        write_simple_xlsx(
            args.output,
            {
                "formacao_grimperium_mini": (
                    list(formation[0]) if formation else [],
                    [dict(row) for row in formation],
                ),
                "reacoes_grimperium_mini": (
                    list(reactions[0]) if reactions else [],
                    [dict(row) for row in reactions],
                ),
            },
        )
        print(f"Wrote {args.output}")
        return 0

    parser.error("unknown command")
    return 2


def _config_from_args(args: argparse.Namespace) -> MiniConfig:
    _defaults = MiniConfig()
    config = MiniConfig(
        work_root=args.work_root,
        results_dir=args.results_dir,
        crest_executable=args.crest_executable or _defaults.crest_executable,
        mopac_executable=args.mopac_executable or _defaults.mopac_executable,
        xtb_executable=args.xtb_executable or _defaults.xtb_executable,
        xtb_enabled=not args.no_xtb,
        threads=args.threads if args.threads is not None else _defaults.threads,
        crest_method=args.crest_method or _defaults.crest_method,
        crest_ewin=args.crest_ewin if args.crest_ewin is not None else _defaults.crest_ewin,
        crest_rthr=args.crest_rthr if args.crest_rthr is not None else _defaults.crest_rthr,
        crest_quick_mode=args.crest_quick_mode,
        max_conformers_to_optimize=args.max_conformers_to_optimize,
    )
    return config


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
