"""CREST execution and multi-XYZ parsing."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from .config import MiniConfig
from .models import CrestResult


def build_crest_command(input_xyz: Path, config: MiniConfig) -> list[str]:
    """Build the CREST command for an XYZ file."""
    method_flag = {
        "gfn2": "--gfn2",
        "gfnff": "--gfnff",
        "gfn2//gfnff": "--gfn2//gfnff",
    }[config.crest_method]
    cmd = [
        config.crest_executable,
        str(input_xyz.resolve()),
        method_flag,
        "--v3",
        "--ewin",
        str(config.crest_ewin),
        "--rthr",
        str(config.crest_rthr),
        "--T",
        str(config.threads),
    ]
    if config.crest_quick_mode:
        cmd.append(f"--{config.crest_quick_mode}")
    if config.crest_max_structures and config.crest_max_structures > 0:
        cmd.extend(["--mstruct", str(config.crest_max_structures)])
    return cmd


def split_crest_conformers(ensemble_file: Path, output_dir: Path) -> list[Path]:
    """Split a CREST multi-XYZ file into conformer_0001.xyz files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    lines = ensemble_file.read_text(encoding="utf-8").splitlines()
    conformers: list[list[str]] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n_atoms = int(lines[i].strip())
        except ValueError:
            i += 1
            continue
        block = lines[i : i + n_atoms + 2]
        if len(block) == n_atoms + 2:
            conformers.append(block)
            i += n_atoms + 2
        else:
            break

    paths: list[Path] = []
    for idx, block in enumerate(conformers, start=1):
        path = output_dir / f"conformer_{idx:04d}.xyz"
        path.write_text("\n".join(block) + "\n", encoding="utf-8")
        paths.append(path)
    return paths


def use_single_conformer(input_xyz: Path, work_dir: Path) -> CrestResult:
    """Use an existing XYZ as the only conformer."""
    work_dir.mkdir(parents=True, exist_ok=True)
    target = work_dir / "conformer_0001.xyz"
    if input_xyz.resolve() != target.resolve():
        shutil.copy(input_xyz, target)
    return CrestResult("skipped", [target], 1, work_dir)


def _find_ensemble(work_dir: Path) -> Path | None:
    """Return the best available CREST output file, or None if nothing exists."""
    for candidate in (
        work_dir / "crest_conformers.xyz",
        work_dir / "crest_best.xyz",
    ):
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    return None


def run_crest(input_xyz: Path, work_dir: Path, config: MiniConfig) -> CrestResult:
    """Run CREST; on timeout or early stop, attempt to recover partial results.

    If crest_conformers.xyz or crest_best.xyz already exists in work_dir (e.g. from
    an interrupted previous run), CREST is skipped and the existing ensemble is used.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    existing = _find_ensemble(work_dir)
    if existing is not None:
        conformers = split_crest_conformers(existing, work_dir)
        if conformers:
            return CrestResult("success", conformers, len(conformers), work_dir)

    input_copy = work_dir / "input.xyz"
    shutil.copy(input_xyz, input_copy)
    cmd = build_crest_command(input_copy, config)
    # Prevent OpenBLAS from spawning its own threads inside each xTB call,
    # which would compete with CREST's OpenMP thread pool and degrade performance.
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1"}

    timed_out = False
    non_zero_exit = False

    try:
        proc = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=config.timeout_crest_s,
            env=env,
        )
        if proc.returncode != 0:
            non_zero_exit = True
    except subprocess.TimeoutExpired:
        timed_out = True
    except OSError as exc:
        return CrestResult("failed", [], 0, work_dir, str(exc), cmd)

    ensemble = _find_ensemble(work_dir)

    if ensemble is None:
        if timed_out:
            reason = "CREST timeout — no partial results found"
        elif non_zero_exit:
            reason = f"CREST exited with {proc.returncode}: {proc.stderr[:500]}"
        else:
            reason = "CREST conformers not found"
        return CrestResult("failed", [], 0, work_dir, reason, cmd)

    conformers = split_crest_conformers(ensemble, work_dir)
    if not conformers:
        return CrestResult("failed", [], 0, work_dir, "No conformers parsed", cmd)

    if timed_out:
        status = "timeout_partial"
    elif non_zero_exit:
        status = "partial"
    else:
        status = "success"

    return CrestResult(status, conformers, len(conformers), work_dir, command=cmd)
