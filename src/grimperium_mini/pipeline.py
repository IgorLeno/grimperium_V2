"""End-to-end CREST + MOPAC pipeline orchestration."""

from __future__ import annotations

import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

from .config import MiniConfig
from .crest import build_crest_command, run_crest, use_single_conformer
from .io import load_pipeline_tasks, write_csv
from .metrics import absolute_deviation, relative_deviation_percent
from .models import ConformerResult, PipelineTask
from .mopac import run_mopac
from .rdkit_seed import generate_initial_xyz
from .reactions import calculate_reactions, load_reaction_rows, write_reactions_csv

FORMATION_COLUMNS = [
    "mol_id",
    "molecule",
    "smiles",
    "formula",
    "group",
    "method",
    "status",
    "crest_status",
    "mopac_status",
    "n_conformers_found",
    "n_conformers_optimized",
    "selected_conformer_rank",
    "selected_conformer_file",
    "hof_selected_kcalmol",
    "hof_selected_kJmol",
    "hf_exp_kJmol",
    "hf_tcc_am1_kJmol",
    "hf_tcc_pm3_kJmol",
    "hf_tcc_pm7_kJmol",
    "hf_joback_kJmol",
    "hf_constantinou_gani_kJmol",
    "ad_kJmol",
    "rd_percent",
    "work_dir",
    "error_message",
    "started_at",
    "finished_at",
    "total_time_s",
]

DETAIL_COLUMNS = [
    "mol_id",
    "molecule",
    "method",
    "conformer_rank",
    "conformer_file",
    "mopac_status",
    "hof_kcalmol",
    "hof_kJmol",
    "mopac_output_file",
    "error_message",
]


def run_pipeline(
    xlsx: Path,
    methods: list[str],
    config: MiniConfig,
    limit: int | None = None,
    dry_run: bool = False,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Run formation calculations and reaction post-processing."""
    tasks = load_pipeline_tasks(xlsx, methods=methods, limit=limit)
    formation_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    for task in tasks:
        row, details = process_task(task, config, dry_run=dry_run)
        formation_rows.append(row)
        detail_rows.extend(_detail_to_rows(details))

    formation_path = config.results_dir / "formacao_grimperium_mini.csv"
    detail_path = config.results_dir / "conformers_detail.csv"
    _backup_if_exists(formation_path)
    _backup_if_exists(detail_path)
    write_csv(formation_path, formation_rows, FORMATION_COLUMNS)
    write_csv(detail_path, detail_rows, DETAIL_COLUMNS)

    reactions = calculate_reactions(load_reaction_rows(xlsx), formation_rows)
    reactions_path = config.results_dir / "reacoes_grimperium_mini.csv"
    _backup_if_exists(reactions_path)
    write_reactions_csv(reactions_path, reactions)
    return formation_rows, detail_rows


def process_task(
    task: PipelineTask,
    config: MiniConfig,
    dry_run: bool = False,
) -> tuple[dict[str, object], list[ConformerResult]]:
    """Process one molecule + method task."""
    started = datetime.now(timezone.utc)
    start = time.time()
    work_dir = config.work_root / _safe_name(task.mol_id) / task.method
    conformer_details: list[ConformerResult] = []
    error = ""
    crest_status = "not_attempted"
    mopac_status = "not_attempted"
    selected: ConformerResult | None = None
    n_found = 0
    n_opt = 0

    if dry_run:
        work_dir.mkdir(parents=True, exist_ok=True)
        planned_input = work_dir / "input.xyz"
        crest_cmd = build_crest_command(planned_input, config)
        print("DRY-RUN CREST:", " ".join(crest_cmd))
        print(
            "DRY-RUN MOPAC:",
            " ".join(
                [
                    config.mopac_executable,
                    str(work_dir / f"conformer_0001_{task.method}.mop"),
                ]
            ),
        )
        crest_status = "dry_run"
        mopac_status = "dry_run"
    else:
        try:
            input_xyz = generate_initial_xyz(
                task.smiles, work_dir / "rdkit_initial.xyz", task.molecule
            )
            if task.run_crest:
                crest_result = run_crest(input_xyz, work_dir / "crest", config)
            else:
                crest_result = use_single_conformer(
                    input_xyz, work_dir / "crest_skipped"
                )
            crest_status = crest_result.status
            n_found = crest_result.n_conformers_found
            if crest_result.status not in {"success", "skipped"}:
                error = crest_result.error_message
            else:
                max_conf = (
                    task.max_conformers_to_process
                    or config.max_conformers_to_optimize
                    or len(crest_result.conformer_files)
                )
                for rank, conformer in enumerate(
                    crest_result.conformer_files[:max_conf], start=1
                ):
                    mop = run_mopac(
                        conformer,
                        work_dir / "mopac" / f"conformer_{rank:04d}",
                        task.method,
                        config,
                        charge=task.charge,
                        multiplicity=task.multiplicity,
                        stem=f"conformer_{rank:04d}",
                    )
                    detail = ConformerResult(
                        mol_id=task.mol_id,
                        molecule=task.molecule,
                        method=task.method,
                        conformer_rank=rank,
                        conformer_file=conformer,
                        mopac_status=mop.status,
                        hof_kcalmol=mop.hof_kcalmol,
                        hof_kJmol=mop.hof_kJmol,
                        mopac_output_file=mop.output_file,
                        error_message=mop.error_message,
                    )
                    conformer_details.append(detail)
                optimized = [d for d in conformer_details if d.hof_kcalmol is not None]
                n_opt = len(optimized)
                selected = select_lowest_hof(optimized)
                mopac_status = _aggregate_mopac_status(conformer_details)
                if selected is None and not error:
                    error = "No conformer with extractable heat of formation"
        except Exception as exc:
            error = str(exc)
            crest_status = "failed" if crest_status == "not_attempted" else crest_status

    finished = datetime.now(timezone.utc)
    status = "dry_run" if dry_run else ("success" if selected is not None else "failed")
    hof_kcal = selected.hof_kcalmol if selected else None
    hof_kj = selected.hof_kJmol if selected else None
    row: dict[str, object] = {
        "mol_id": task.mol_id,
        "molecule": task.molecule,
        "smiles": task.smiles,
        "formula": task.formula,
        "group": task.group,
        "method": task.method,
        "status": status,
        "crest_status": crest_status,
        "mopac_status": mopac_status,
        "n_conformers_found": n_found,
        "n_conformers_optimized": n_opt,
        "selected_conformer_rank": selected.conformer_rank if selected else "",
        "selected_conformer_file": str(selected.conformer_file) if selected else "",
        "hof_selected_kcalmol": hof_kcal if hof_kcal is not None else "",
        "hof_selected_kJmol": hof_kj if hof_kj is not None else "",
        "hf_exp_kJmol": task.hf_exp_kJmol if task.hf_exp_kJmol is not None else "",
        "hf_tcc_am1_kJmol": (
            task.hf_tcc_am1_kJmol if task.hf_tcc_am1_kJmol is not None else ""
        ),
        "hf_tcc_pm3_kJmol": (
            task.hf_tcc_pm3_kJmol if task.hf_tcc_pm3_kJmol is not None else ""
        ),
        "hf_tcc_pm7_kJmol": (
            task.hf_tcc_pm7_kJmol if task.hf_tcc_pm7_kJmol is not None else ""
        ),
        "hf_joback_kJmol": (
            task.hf_joback_kJmol if task.hf_joback_kJmol is not None else ""
        ),
        "hf_constantinou_gani_kJmol": (
            task.hf_constantinou_gani_kJmol
            if task.hf_constantinou_gani_kJmol is not None
            else ""
        ),
        "ad_kJmol": absolute_deviation(hof_kj, task.hf_exp_kJmol) or "",
        "rd_percent": relative_deviation_percent(hof_kj, task.hf_exp_kJmol) or "",
        "work_dir": str(work_dir),
        "error_message": error,
        "started_at": started.isoformat(),
        "finished_at": finished.isoformat(),
        "total_time_s": round(time.time() - start, 3),
    }
    return row, conformer_details


def select_lowest_hof(conformers: list[ConformerResult]) -> ConformerResult | None:
    """Select the conformer with the lowest heat of formation."""
    valid = [c for c in conformers if c.hof_kcalmol is not None]
    if not valid:
        return None
    return min(
        valid,
        key=lambda c: c.hof_kcalmol if c.hof_kcalmol is not None else float("inf"),
    )


def _aggregate_mopac_status(details: list[ConformerResult]) -> str:
    if not details:
        return "not_attempted"
    if any(
        d.mopac_status in {"success", "warning"} and d.hof_kcalmol is not None
        for d in details
    ):
        return "success"
    return "failed"


def _detail_to_rows(details: list[ConformerResult]) -> list[dict[str, object]]:
    return [
        {
            "mol_id": d.mol_id,
            "molecule": d.molecule,
            "method": d.method,
            "conformer_rank": d.conformer_rank,
            "conformer_file": str(d.conformer_file),
            "mopac_status": d.mopac_status,
            "hof_kcalmol": d.hof_kcalmol if d.hof_kcalmol is not None else "",
            "hof_kJmol": d.hof_kJmol if d.hof_kJmol is not None else "",
            "mopac_output_file": (
                str(d.mopac_output_file) if d.mopac_output_file else ""
            ),
            "error_message": d.error_message,
        }
        for d in details
    ]


def _backup_if_exists(path: Path) -> None:
    if not path.exists():
        return
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy(path, path.with_name(f"{path.stem}.{stamp}.bak{path.suffix}"))


def _safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in value)
