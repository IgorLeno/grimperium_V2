"""End-to-end CREST + MOPAC pipeline orchestration."""

from __future__ import annotations

import logging
import shutil
import time
from collections.abc import Callable
from datetime import datetime
from pathlib import Path
from typing import Any, cast

from .config import MiniConfig
from .crest import build_crest_command, run_crest, use_single_conformer
from .io import append_summary_row, load_completed_mol_ids, load_molecules
from .models import ConformerResult, SemiempiricalMethod
from .mopac import run_mopac
from .rdkit_seed import generate_initial_xyz
from .xtb import run_xtb_preopt, verify_xtb_runtime

logger = logging.getLogger(__name__)

SUMMARY_COLUMNS = [
    "mol_id",
    "molecule",
    "smiles",
    "n_conformers_founded",
    "hf_am1_kJmol",
    "hf_pm3_kJmol",
    "hf_pm7_kJmol",
    "status",
    "crest_status",
    "mopac_status",
    "hf_exp_kJmol",
    "total_time_s",
    "total_time_min",
]

_METHODS: list[SemiempiricalMethod] = ["AM1", "PM3", "PM7"]


def run_pipeline(
    xlsx: Path,
    config: MiniConfig,
    limit: int | None = None,
    dry_run: bool = False,
    on_progress: Callable[[str, Any], None] | None = None,
) -> list[dict[str, object]]:
    """Run formation calculations for each molecule (CREST 1× → MOPAC 3×)."""
    tasks = load_molecules(xlsx, limit=limit)
    summary_path = config.results_dir / "grimperium_mini_summary.xlsx"

    completed_ids: set[str] = set()
    if not dry_run:
        completed_ids = load_completed_mol_ids(summary_path)
        if completed_ids:
            logger.info(
                "Resume mode: skipping %d already-successful molecules: %s",
                len(completed_ids),
                sorted(completed_ids),
            )

    if config.xtb_enabled and not dry_run:
        ok, reason = verify_xtb_runtime(config.xtb_executable, threads=config.threads)
        if not ok:
            logger.warning(
                "xTB pre-optimization disabled: binary check failed (%s). "
                "If you are using the Fedora-packaged xtb 6.7.1-4 (fc43), this is a "
                "known Fortran format-string bug in optimizer.f90:852. "
                "Workarounds: install conda-forge xtb=6.7.1, pass --no-xtb, "
                "or set GRIMPERIUM_MINI_XTB_ENABLED=false.",
                reason,
            )
            config.xtb_enabled = False

    rows: list[dict[str, object]] = []
    for mol in tasks:
        mol_id = str(mol.get("mol_id") or "")

        if mol_id in completed_ids:
            logger.debug("Skipping %s (already success)", mol_id)
            skipped_row: dict[str, object] = {
                **mol,  # type: ignore[arg-type]
                "status": "skipped",
                "crest_status": "skipped",
                "mopac_status": "skipped",
                "n_conformers_founded": "",
                "hf_am1_kJmol": "",
                "hf_pm3_kJmol": "",
                "hf_pm7_kJmol": "",
                "total_time_s": 0,
                "total_time_min": 0,
            }
            if on_progress:
                on_progress("done", skipped_row)
            continue

        if on_progress:
            on_progress("start", mol)
        row = process_molecule(mol, config, dry_run=dry_run)
        append_summary_row(summary_path, row, SUMMARY_COLUMNS)
        rows.append(row)
        if on_progress:
            on_progress("done", row)

    return rows


def process_molecule(
    mol: dict[str, object],
    config: MiniConfig,
    dry_run: bool = False,
) -> dict[str, object]:
    """Process one molecule: RDKit → xTB → CREST → MOPAC×3."""
    mol_id = str(mol.get("mol_id") or "")
    smiles = str(mol.get("smiles") or "")
    molecule = str(mol.get("molecule") or "")
    charge = cast(int, mol.get("charge", 0))
    multiplicity = cast(int, mol.get("multiplicity", 1))
    run_crest_flag = cast(bool, mol.get("run_crest", True))

    work_dir = config.work_root / _safe_name(mol_id)
    start = time.time()

    if dry_run:
        work_dir.mkdir(parents=True, exist_ok=True)
        rdkit_xyz = work_dir / "rdkit_initial.xyz"
        xtb_xyz = (work_dir / "xtb" / f"{mol_id}_preopt.xyz").resolve()
        crest_input = xtb_xyz if config.xtb_enabled else rdkit_xyz.resolve()
        if config.xtb_enabled:
            print(
                "DRY-RUN XTB:",
                config.xtb_executable,
                str(rdkit_xyz.resolve()),
                "--opt --gfn2 --temp 300 --T",
                str(config.threads),
            )
        crest_cmd = build_crest_command(crest_input, config)
        print("DRY-RUN CREST:", " ".join(crest_cmd))
        for method in _METHODS:
            mop_file = work_dir / "mopac" / method.lower() / f"{mol_id}_{method}.mop"
            print("DRY-RUN MOPAC:", config.mopac_executable, str(mop_file.resolve()))
        return _dry_run_row(mol, start)

    try:
        input_xyz = generate_initial_xyz(
            smiles, work_dir / "rdkit_initial.xyz", molecule
        )

        if config.xtb_enabled:
            xtb_result = run_xtb_preopt(
                input_xyz,
                work_dir / "xtb",
                mol_id,
                config.xtb_executable,
                config.timeout_xtb_s,
                threads=config.threads,
            )
            if not xtb_result.success or xtb_result.output_xyz is None:
                return _failed_row(mol, "xtb_failed", xtb_result.error_message, start)
            crest_input_path: Path = xtb_result.output_xyz
        else:
            crest_input_path = input_xyz

        if run_crest_flag:
            crest_result = run_crest(crest_input_path, work_dir, config)
        else:
            crest_result = use_single_conformer(crest_input_path, work_dir)

        _CREST_RECOVERABLE = {"success", "skipped", "timeout_partial", "partial"}
        if crest_result.status not in _CREST_RECOVERABLE:
            return _failed_row(
                mol, crest_result.status, crest_result.error_message, start
            )

        best_conformer = crest_result.conformer_files[0]
        n_conformers = crest_result.n_conformers_found
        final_crest_status = crest_result.status

        hof: dict[str, float | None] = {}
        for method in _METHODS:
            mop_result = run_mopac(
                best_conformer,
                work_dir / "mopac" / method.lower(),
                method,
                config,
                charge=charge,
                multiplicity=multiplicity,
                stem=f"{mol_id}_{method}",
            )
            hof[method] = mop_result.hof_kJmol

        all_ok = all(v is not None for v in hof.values())
        any_ok = any(v is not None for v in hof.values())
        mopac_status = "success" if all_ok else ("partial" if any_ok else "failed")
        overall = "success" if all_ok else ("partial" if any_ok else "failed")
        elapsed = time.time() - start
        exp = mol.get("hf_exp_kJmol")

        return {
            "mol_id": mol_id,
            "molecule": molecule,
            "smiles": smiles,
            "n_conformers_founded": n_conformers,
            "hf_am1_kJmol": hof["AM1"] if hof["AM1"] is not None else "",
            "hf_pm3_kJmol": hof["PM3"] if hof["PM3"] is not None else "",
            "hf_pm7_kJmol": hof["PM7"] if hof["PM7"] is not None else "",
            "status": overall,
            "crest_status": final_crest_status,
            "mopac_status": mopac_status,
            "hf_exp_kJmol": exp if exp is not None else "",
            "total_time_s": round(elapsed, 1),
            "total_time_min": round(elapsed / 60, 2),
        }

    except Exception as exc:
        return _failed_row(mol, "failed", str(exc), start)


def _dry_run_row(mol: dict[str, object], start: float) -> dict[str, object]:
    elapsed = time.time() - start
    exp = mol.get("hf_exp_kJmol")
    return {
        "mol_id": str(mol.get("mol_id") or ""),
        "molecule": str(mol.get("molecule") or ""),
        "smiles": str(mol.get("smiles") or ""),
        "n_conformers_founded": 0,
        "hf_am1_kJmol": "",
        "hf_pm3_kJmol": "",
        "hf_pm7_kJmol": "",
        "status": "dry_run",
        "crest_status": "dry_run",
        "mopac_status": "dry_run",
        "hf_exp_kJmol": exp if exp is not None else "",
        "total_time_s": round(elapsed, 1),
        "total_time_min": round(elapsed / 60, 2),
    }


def _failed_row(
    mol: dict[str, object], crest_status: str, error: str, start: float
) -> dict[str, object]:
    elapsed = time.time() - start
    exp = mol.get("hf_exp_kJmol")
    return {
        "mol_id": str(mol.get("mol_id") or ""),
        "molecule": str(mol.get("molecule") or ""),
        "smiles": str(mol.get("smiles") or ""),
        "n_conformers_founded": 0,
        "hf_am1_kJmol": "",
        "hf_pm3_kJmol": "",
        "hf_pm7_kJmol": "",
        "status": "failed",
        "crest_status": crest_status,
        "mopac_status": "not_attempted",
        "hf_exp_kJmol": exp if exp is not None else "",
        "total_time_s": round(elapsed, 1),
        "total_time_min": round(elapsed / 60, 2),
    }


def select_lowest_hof(conformers: list[ConformerResult]) -> ConformerResult | None:
    """Select the conformer with the lowest heat of formation."""
    valid = [c for c in conformers if c.hof_kcalmol is not None]
    if not valid:
        return None
    return min(
        valid,
        key=lambda c: c.hof_kcalmol if c.hof_kcalmol is not None else float("inf"),
    )


def _backup_if_exists(path: Path) -> None:
    if not path.exists():
        return
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy(path, path.with_name(f"{path.stem}.{stamp}.bak{path.suffix}"))


def _safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in value)
