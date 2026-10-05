"""Conformer-selection validation (spec section 5): sample, run, summarize.

Commands::

    # 1. Draw the versioned sample (reads the dataset only).
    poetry run python scripts/conformer_validation sample \
        --dataset data/thermo_cbs_chon_v2.csv --size 200 --seed 20261004 \
        --out reports/conformer_validation/sample.csv

    # 2. Plan, then run CREST once per molecule and arms A/B/C on it.
    poetry run python scripts/conformer_validation run \
        --sample reports/conformer_validation/sample.csv \
        --work-dir runs/conformer_validation \
        --results reports/conformer_validation/results.jsonl \
        --workers 4 --crest-threads 4 --crest-timeout 7200 --dry-run
    #    (drop --dry-run to execute; --limit N for a pilot; --arms A,B to
    #    restrict arms; --ensemble-dir DIR to smoke-test from stored
    #    <DIR>/<mol_id>/crest_conformers.xyz without running CREST;
    #    --retry-errors to rerun pairs whose latest line is an error)

    # 3. Per-arm metrics and the pre-registered decision.
    poetry run python scripts/conformer_validation summarize \
        --results reports/conformer_validation/results.jsonl \
        --out reports/conformer_validation/

Arms: A = CREST energy Top-10; B = CONFPASS port Top-10; C = folding
filter (k=6, f=1.00, c_min=1, no r_max) + B. PM7 with
``require_minimum``. ``run`` is resumable: finished (mol_id, arm) pairs
in the results file are skipped and finished CREST runs are reused.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from conformer_validation.arms import ARM_NAMES, parse_arms  # noqa: E402
from conformer_validation.run import (  # noqa: E402
    MoleculeTask,
    build_tasks,
    done_pairs,
    load_references,
    plan_text,
    read_results,
    run_tasks,
)
from conformer_validation.sampling import (  # noqa: E402
    DEFAULT_SEED,
    DEFAULT_SIZE,
    draw,
    load_dataset,
    render_report,
    sha256_of,
)
from conformer_validation.summarize import (  # noqa: E402
    render_report as render_summary,
)
from conformer_validation.summarize import (  # noqa: E402
    summarize,
    write_outputs,
)

DEFAULT_DATASET = Path("data/thermo_cbs_chon_v2.csv")


def _sample(args: argparse.Namespace) -> int:
    frame = load_dataset(args.dataset)
    sample, report = draw(
        frame,
        size=args.size,
        seed=args.seed,
        dataset_sha256=sha256_of(args.dataset),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    sample.to_csv(args.out, index=False)
    print(render_report(report, sample))
    print(f"wrote {args.out}")
    return 0


def _run(args: argparse.Namespace) -> int:
    import pandas as pd

    sample = pd.read_csv(args.sample, dtype={"mol_id": str, "smiles": str})
    references = load_references(args.dataset, sample)
    done = done_pairs(read_results(args.results), retry_errors=args.retry_errors)
    template = MoleculeTask(
        mol_id="-",
        smiles="-",
        arms=(),
        reference_xyz=None,
        h298_cbs=None,
        work_dir=str(args.work_dir.resolve()),
        ensemble_dir=(
            None if args.ensemble_dir is None else str(args.ensemble_dir.resolve())
        ),
        crest_threads=args.crest_threads,
        crest_timeout=args.crest_timeout,
        crest_executable=args.crest_executable,
        xtb_executable=args.xtb_executable,
        mopac_executable=args.mopac_executable,
    )
    tasks = build_tasks(
        sample,
        references,
        arms=parse_arms(args.arms),
        done=done,
        limit=args.limit,
        template=template,
    )
    missing = sorted(set(sample["mol_id"]) - set(references))
    if missing:
        print(
            f"warning: {len(missing)} sampled molecules have no reference in "
            f"{args.dataset}; their delta and RMSD will be empty"
        )
    print(plan_text(tasks, workers=args.workers))
    if args.dry_run:
        return 0
    written = run_tasks(tasks, args.results, workers=args.workers)
    print(f"appended {written} lines to {args.results}")
    return 0


def _summarize(args: argparse.Namespace) -> int:
    summaries, decision = summarize(read_results(args.results))
    summary_path, report_path = write_outputs(summaries, decision, args.out)
    print(render_summary(summaries, decision))
    print(f"wrote {summary_path} and {report_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Argument parser for the three commands."""
    parser = argparse.ArgumentParser(
        prog="conformer_validation",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    commands = parser.add_subparsers(dest="command", required=True)

    sample = commands.add_parser("sample", help="draw the stratified sample")
    sample.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    sample.add_argument("--size", type=int, default=DEFAULT_SIZE)
    sample.add_argument("--seed", type=int, default=DEFAULT_SEED)
    sample.add_argument("--out", type=Path, required=True)
    sample.set_defaults(handler=_sample)

    run = commands.add_parser("run", help="run CREST and the arms")
    run.add_argument("--sample", type=Path, required=True)
    run.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    run.add_argument("--work-dir", type=Path, required=True)
    run.add_argument("--results", type=Path, required=True)
    run.add_argument("--workers", type=int, default=4)
    run.add_argument("--crest-threads", type=int, default=4)
    run.add_argument("--crest-timeout", type=float, default=7200.0)
    run.add_argument("--limit", type=int, default=None)
    run.add_argument("--arms", default=",".join(ARM_NAMES))
    run.add_argument("--ensemble-dir", type=Path, default=None)
    run.add_argument("--retry-errors", action="store_true")
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--crest-executable", default="crest")
    run.add_argument("--xtb-executable", default="xtb")
    run.add_argument("--mopac-executable", default="mopac")
    run.set_defaults(handler=_run)

    summary = commands.add_parser("summarize", help="per-arm metrics and decision")
    summary.add_argument("--results", type=Path, required=True)
    summary.add_argument("--out", type=Path, required=True)
    summary.set_defaults(handler=_summarize)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    args = build_parser().parse_args(argv)
    if getattr(args, "workers", 1) < 1:
        raise SystemExit("--workers must be >= 1")
    return int(args.handler(args))


if __name__ == "__main__":
    raise SystemExit(main())
