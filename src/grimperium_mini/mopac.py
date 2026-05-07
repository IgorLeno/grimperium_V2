"""MOPAC input generation, execution, and heat-of-formation parsing."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .config import MiniConfig
from .models import KCAL_TO_KJ, MopacResult, SemiempiricalMethod

_MULTIPLICITY = {
    2: "DOUBLET",
    3: "TRIPLET",
    4: "QUARTET",
    5: "QUINTET",
    6: "SEXTET",
}

_HOF_PATTERNS = [
    re.compile(
        r"FINAL\s+HEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+(?:\.\d+)?)\s*KCAL/MOL",
        re.IGNORECASE,
    ),
    re.compile(
        r"\bHEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+(?:\.\d+)?)\s*KCAL/MOL",
        re.IGNORECASE,
    ),
]


def build_mopac_keywords(
    method: SemiempiricalMethod,
    charge: int = 0,
    multiplicity: int = 1,
    extra: list[str] | None = None,
) -> list[str]:
    """Build method-independent MOPAC keywords."""
    keywords = [method, *(extra if extra is not None else ["PRECISE", "EF", "AUX"])]
    if charge != 0:
        keywords.append(f"CHARGE={charge}")
    if multiplicity != 1:
        if multiplicity not in _MULTIPLICITY:
            raise ValueError(f"Unsupported multiplicity: {multiplicity}")
        keywords.append(_MULTIPLICITY[multiplicity])
    return keywords


def create_mopac_input(
    xyz_file: Path,
    output_mop: Path,
    method: SemiempiricalMethod,
    charge: int = 0,
    multiplicity: int = 1,
    extra_keywords: list[str] | None = None,
) -> None:
    """Create a MOPAC .mop file from XYZ coordinates."""
    lines = xyz_file.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError("XYZ file must contain atom count and comment line")
    n_atoms = int(lines[0].strip())
    if len(lines) < n_atoms + 2:
        raise ValueError("XYZ file has fewer coordinate lines than declared")
    keywords = build_mopac_keywords(method, charge, multiplicity, extra_keywords)
    out = [" ".join(keywords), lines[1].strip(), ""]
    for line in lines[2 : 2 + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ coordinate line: {line}")
        element, x, y, z = parts[:4]
        out.append(f"  {element:2s}  {x:>14s}  1  {y:>14s}  1  {z:>14s}  1")
    output_mop.parent.mkdir(parents=True, exist_ok=True)
    output_mop.write_text("\n".join(out) + "\n", encoding="utf-8")


def parse_hof_kcalmol(content: str) -> float | None:
    """Extract heat of formation in kcal/mol from MOPAC text."""
    for pattern in _HOF_PATTERNS:
        match = pattern.search(content)
        if match:
            return float(match.group(1))
    return None


def parse_hof_file(path: Path) -> float | None:
    """Extract heat of formation from an output or archive file."""
    return parse_hof_kcalmol(path.read_text(encoding="utf-8", errors="replace"))


def run_mopac(
    xyz_file: Path,
    work_dir: Path,
    method: SemiempiricalMethod,
    config: MiniConfig,
    charge: int = 0,
    multiplicity: int = 1,
    stem: str = "conformer",
) -> MopacResult:
    """Run MOPAC and parse heat of formation."""
    work_dir.mkdir(parents=True, exist_ok=True)
    mop_file = work_dir / f"{stem}_{method}.mop"
    create_mopac_input(
        xyz_file,
        mop_file,
        method=method,
        charge=charge,
        multiplicity=multiplicity,
        extra_keywords=config.mopac_keywords_extra,
    )
    cmd = [config.mopac_executable, str(mop_file)]
    try:
        proc = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=config.timeout_mopac_s,
        )
    except subprocess.TimeoutExpired:
        return MopacResult(
            "failed", input_file=mop_file, error_message="MOPAC timeout", command=cmd
        )
    except OSError as exc:
        return MopacResult(
            "failed", input_file=mop_file, error_message=str(exc), command=cmd
        )

    output = mop_file.with_suffix(".out")
    if not output.exists():
        output = mop_file.with_suffix(".arc")
    if not output.exists():
        return MopacResult(
            "failed",
            input_file=mop_file,
            error_message="No MOPAC .out or .arc file found",
            command=cmd,
        )

    hof = parse_hof_file(output)
    if hof is None:
        return MopacResult(
            "failed",
            input_file=mop_file,
            output_file=output,
            error_message="Heat of formation not found",
            command=cmd,
        )
    status = "warning" if proc.returncode != 0 else "success"
    error = f"MOPAC exited with {proc.returncode}" if proc.returncode != 0 else ""
    return MopacResult(
        status,
        hof_kcalmol=hof,
        hof_kJmol=hof * KCAL_TO_KJ,
        output_file=output,
        input_file=mop_file,
        error_message=error,
        command=cmd,
    )
