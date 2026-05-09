"""xTB pre-optimization step for grimperium_mini."""

from __future__ import annotations

import subprocess
import tempfile
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
        "--silent",
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


_H2O_XYZ = """\
3
water preflight
O  0.0000  0.0000  0.0000
H  0.9572  0.0000  0.0000
H -0.2399  0.9267  0.0000
"""

_PREFLIGHT_TIMEOUT_S = 60


def verify_xtb_runtime(xtb_executable: str, threads: int = 1) -> tuple[bool, str]:
    """Test the xTB binary with a tiny molecule before the main run.

    Returns (True, "") if xtb works, (False, reason) if it is broken or absent.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        xyz = tmp / "h2o.xyz"
        xyz.write_text(_H2O_XYZ)
        cmd = [
            xtb_executable,
            str(xyz),
            "--opt",
            "--gfn2",
            "--silent",
            "--T",
            str(threads),
        ]
        try:
            proc = subprocess.run(
                cmd,
                cwd=tmp,
                capture_output=True,
                text=True,
                timeout=_PREFLIGHT_TIMEOUT_S,
                check=False,
            )
        except FileNotFoundError:
            return False, f"xTB executable not found: {xtb_executable}"
        except subprocess.TimeoutExpired:
            return False, f"xTB preflight timed out after {_PREFLIGHT_TIMEOUT_S}s"

        xtbopt = tmp / "xtbopt.xyz"
        if proc.returncode == 0 and xtbopt.exists():
            return True, ""

        stderr_snippet = proc.stderr[:200].strip()
        return (
            False,
            (
                f"rc={proc.returncode}: {stderr_snippet}"
                if stderr_snippet
                else f"rc={proc.returncode}, xtbopt.xyz not produced"
            ),
        )
