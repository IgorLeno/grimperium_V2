"""Multi-conformer MOPAC analysis pipeline.

Reuses existing CREST conformer files (conformer_NNNN.xyz) and runs MOPAC
on the best N conformers for AM1, PM3, and PM7.  Never re-runs CREST.
"""

from __future__ import annotations

import logging
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

from .config import MiniConfig
from .io import load_molecules, read_csv_rows, write_csv
from .metrics import absolute_deviation
from .models import SemiempiricalMethod
from .mopac import run_mopac

logger = logging.getLogger(__name__)

_METHODS: list[SemiempiricalMethod] = ["AM1", "PM3", "PM7"]


def _safe_name(value: str) -> str:
    """Convert a string to a filesystem-safe name (mirrors pipeline.py)."""
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in value)


def select_existing_conformers(run_dir: Path, max_n: int = 10) -> list[Path]:
    """Return up to *max_n* existing conformer XYZ files sorted by energy.

    CREST writes conformers in ascending energy order, so
    ``conformer_0001.xyz`` is always the lowest-energy structure.
    Lexicographic sort of the filenames preserves that ordering.

    Returns an empty list when the directory does not exist or contains
    no ``conformer_*.xyz`` files.
    """
    if not run_dir.is_dir():
        return []
    paths = sorted(run_dir.glob("conformer_*.xyz"))
    return paths[:max_n]


def _build_columns(max_conformers: int) -> list[str]:
    """Build the ordered column list for the wide CSV layout."""
    meta = [
        "mol_id",
        "molecule",
        "smiles",
        "hf_exp_kJmol",
        "n_conformers_used",
        "status",
        "mopac_status",
        "total_time_s",
        "total_time_min",
    ]
    method_cols: list[str] = []
    for method in _METHODS:
        m = method.lower()
        for n in range(1, max_conformers + 1):
            method_cols.append(f"hof_{m}_conf{n:02d}_kJmol")
        for n in range(1, max_conformers + 1):
            method_cols.append(f"absdev_{m}_conf{n:02d}_kJmol")
        method_cols.append(f"best_conf_index_{m}")
        method_cols.append(f"best_hof_{m}_kJmol")
    return meta + method_cols


def process_molecule_multi(
    mol: dict[str, object],
    config: MiniConfig,
    max_conformers: int = 10,
) -> dict[str, object]:
    """Process one molecule in multi-conformer mode.

    Reads up to *max_conformers* pre-computed CREST conformer files from
    ``config.work_root / mol_id`` and runs MOPAC (AM1 + PM3 + PM7) on each.
    CREST is **never** invoked.

    Result columns per method *M*:
    - ``hof_{M}_conf{N:02d}_kJmol``   — heat of formation (blank on MOPAC failure)
    - ``absdev_{M}_conf{N:02d}_kJmol``— |HoF − HoF_exp| (blank if no exp. data)
    - ``best_conf_index_{M}``          — 1-based index of closest conformer to exp.
    - ``best_hof_{M}_kJmol``           — HoF of that conformer
    """
    mol_id = str(mol.get("mol_id") or "")
    smiles = str(mol.get("smiles") or "")
    molecule = str(mol.get("molecule") or "")
    charge = cast(int, mol.get("charge", 0))
    multiplicity = cast(int, mol.get("multiplicity", 1))
    hf_exp = mol.get("hf_exp_kJmol")
    hf_exp_float: float | None = cast(float, hf_exp) if hf_exp is not None else None

    run_dir = config.work_root / _safe_name(mol_id)
    start = time.time()

    conformers = select_existing_conformers(run_dir, max_n=max_conformers)
    if not conformers:
        logger.warning("No conformer files found for %s in %s", mol_id, run_dir)
        return _no_conformers_row(mol, max_conformers, start)

    n_confs = len(conformers)
    logger.info(
        "Multi-conformer: processing %s with %d conformers (max %d)",
        mol_id,
        n_confs,
        max_conformers,
    )

    # hof_results[method] = list of kJmol values (None on MOPAC failure)
    hof_results: dict[str, list[float | None]] = {m: [] for m in _METHODS}
    any_mopac_ok = False
    all_mopac_ok = True

    for method in _METHODS:
        for idx, conf_path in enumerate(conformers, start=1):
            mop_dir = run_dir / "mopac_multi" / method.lower() / f"conf{idx:02d}"
            result = run_mopac(
                conf_path,
                mop_dir,
                method,
                config,
                charge=charge,
                multiplicity=multiplicity,
                stem=f"{mol_id}_conf{idx:02d}",
            )
            hof = result.hof_kJmol
            hof_results[method].append(hof)
            if hof is not None:
                any_mopac_ok = True
            else:
                all_mopac_ok = False
                logger.debug(
                    "MOPAC failed for %s %s conf %d: %s",
                    mol_id,
                    method,
                    idx,
                    result.error_message,
                )

    mopac_status = (
        "success" if all_mopac_ok else ("partial" if any_mopac_ok else "failed")
    )
    overall = "success" if any_mopac_ok else "failed"
    elapsed = time.time() - start

    row: dict[str, object] = {
        "mol_id": mol_id,
        "molecule": molecule,
        "smiles": smiles,
        "hf_exp_kJmol": hf_exp if hf_exp is not None else "",
        "n_conformers_used": n_confs,
        "status": overall,
        "mopac_status": mopac_status,
        "total_time_s": round(elapsed, 1),
        "total_time_min": round(elapsed / 60, 2),
    }

    for method in _METHODS:
        m = method.lower()
        hofs = hof_results[method]
        # Pad to max_conformers so every column is present
        padded: list[float | None] = hofs + [None] * (max_conformers - len(hofs))

        for n, hof in enumerate(padded, start=1):
            row[f"hof_{m}_conf{n:02d}_kJmol"] = hof if hof is not None else ""

        for n, hof in enumerate(padded, start=1):
            dev = absolute_deviation(hof, hf_exp_float)
            row[f"absdev_{m}_conf{n:02d}_kJmol"] = dev if dev is not None else ""

        # Best = conformer with minimum absolute deviation vs. experimental
        # Blank when experimental HoF is unavailable.
        best_idx: int | None = None
        best_hof: float | None = None
        if hf_exp_float is not None:
            best_dev: float | None = None
            for n, hof in enumerate(hofs, start=1):
                if hof is None:
                    continue
                dev_val = absolute_deviation(hof, hf_exp_float)
                if dev_val is not None and (best_dev is None or dev_val < best_dev):
                    best_dev = dev_val
                    best_idx = n
                    best_hof = hof

        row[f"best_conf_index_{m}"] = best_idx if best_idx is not None else ""
        row[f"best_hof_{m}_kJmol"] = best_hof if best_hof is not None else ""

    return row


def _no_conformers_row(
    mol: dict[str, object], max_conformers: int, start: float
) -> dict[str, object]:
    """Return a blank result row for a molecule whose run directory has no conformers."""
    elapsed = time.time() - start
    hf_exp = mol.get("hf_exp_kJmol")
    row: dict[str, object] = {
        "mol_id": str(mol.get("mol_id") or ""),
        "molecule": str(mol.get("molecule") or ""),
        "smiles": str(mol.get("smiles") or ""),
        "hf_exp_kJmol": hf_exp if hf_exp is not None else "",
        "n_conformers_used": 0,
        "status": "no_conformers",
        "mopac_status": "not_attempted",
        "total_time_s": round(elapsed, 1),
        "total_time_min": round(elapsed / 60, 2),
    }
    for method in _METHODS:
        m = method.lower()
        for n in range(1, max_conformers + 1):
            row[f"hof_{m}_conf{n:02d}_kJmol"] = ""
            row[f"absdev_{m}_conf{n:02d}_kJmol"] = ""
        row[f"best_conf_index_{m}"] = ""
        row[f"best_hof_{m}_kJmol"] = ""
    return row


def run_multi_conformer_pipeline(
    xlsx: Path,
    config: MiniConfig,
    limit: int | None = None,
    max_conformers: int = 10,
    on_progress: Callable[[str, Any], None] | None = None,
) -> list[dict[str, object]]:
    """Run multi-conformer MOPAC analysis for every molecule in *xlsx*.

    Reuses existing CREST conformer files.  Never runs CREST.
    Writes a wide CSV to ``config.results_dir/grimperium_mini_multiconf_summary.csv``.

    Resume mode: molecules already present with ``status == 'success'`` in the
    output CSV are skipped and their existing rows are preserved as-is.

    Args:
        xlsx: Path to the pipeline workbook.
        config: Runtime configuration (work_root and results_dir are used).
        limit: Optional cap on the number of molecules to process.
        max_conformers: Maximum number of conformers to consider per molecule.
        on_progress: Optional callback called with ``("start", mol_dict)`` and
            ``("done", row_dict)`` for each molecule.  Skipped molecules fire
            only a ``"done"`` event with ``status="skipped"``.

    Returns:
        List of result row dicts (one per molecule in *tasks*, including
        preserved rows for skipped molecules).
    """
    tasks = load_molecules(xlsx, limit=limit)
    columns = _build_columns(max_conformers)
    output_path = config.results_dir / "grimperium_mini_multiconf_summary.csv"

    # Resume: find already-successful rows from a prior run
    completed_ids: set[str] = set()
    existing_by_id: dict[str, dict[str, object]] = {}
    if output_path.exists():
        prior = read_csv_rows(output_path)
        for r in prior:
            existing_by_id[r["mol_id"]] = dict(r)
        completed_ids = {r["mol_id"] for r in prior if r.get("status") == "success"}
        if completed_ids:
            logger.info(
                "Resume mode: skipping %d already-successful molecules: %s",
                len(completed_ids),
                sorted(completed_ids),
            )

    rows: list[dict[str, object]] = []
    for mol in tasks:
        mol_id = str(mol.get("mol_id") or "")
        if mol_id in completed_ids:
            logger.debug("Skipping %s (already success)", mol_id)
            if on_progress:
                on_progress("done", {"mol_id": mol_id, "status": "skipped"})
            rows.append(existing_by_id[mol_id])
            continue
        if on_progress:
            on_progress("start", mol)
        row = process_molecule_multi(mol, config, max_conformers=max_conformers)
        rows.append(row)
        if on_progress:
            on_progress("done", row)

    write_csv(output_path, rows, columns)
    logger.info("Wrote multi-conformer summary to %s", output_path)
    return rows
