"""xTB pre-optimization step for grimperium_mini."""

from __future__ import annotations

import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass
class XTBResult:
    status: str  # "success", "failed", "skipped"
    output_xyz: Path | None
    error_message: str = ""
    time_s: float = 0.0

    @property
    def success(self) -> bool:
        return self.status in {"success", "skipped"}


def run_xtb_preopt(
    input_xyz: Path,
    work_dir: Path,
    mol_id: str,
    xtb_executable: str,
    timeout_s: int,
    threads: int = 1,
) -> XTBResult:
    """Run GFN2-xTB geometry pre-optimization before CREST."""
    work_dir.mkdir(parents=True, exist_ok=True)
    output_xyz = work_dir / f"{mol_id}_preopt.xyz"

    cmd = [
        xtb_executable,
        str(input_xyz.resolve()),
        "--opt",
        "--gfn2",
        "--temp",
        "300",
        "--T",
        str(threads),
    ]

    start = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
        elapsed = time.perf_counter() - start

        xtbopt = work_dir / "xtbopt.xyz"
        if proc.returncode == 0 and xtbopt.exists():
            xtbopt.replace(output_xyz)
            return XTBResult("success", output_xyz, time_s=elapsed)

        return XTBResult(
            "failed",
            None,
            f"xTB rc={proc.returncode}: {proc.stderr[:300]}",
            elapsed,
        )

    except FileNotFoundError as exc:
        elapsed = time.perf_counter() - start
        return XTBResult("failed", None, f"xTB executable not found: {exc}", elapsed)

    except subprocess.TimeoutExpired:
        elapsed = time.perf_counter() - start
        return XTBResult("failed", None, f"xTB timed out after {timeout_s}s", elapsed)
