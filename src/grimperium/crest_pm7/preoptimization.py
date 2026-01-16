"""xTB Pre-optimization for molecular structures.

This module provides xTBPreOptimizer for optional GFN2-xTB pre-optimization
before CREST conformer generation. Pre-optimization can improve convergence
and reduce overall computation time for complex molecules.

Usage:
    >>> from grimperium.crest_pm7.preoptimization import xTBPreOptimizer
    >>> preopt = xTBPreOptimizer(enabled=True, timeout_seconds=300)
    >>> result = preopt.preoptimize_structure("mol001", Path("input.xyz"), Path("work"))
    >>> if result.success:
    ...     print(f"Pre-optimized: {result.output_xyz}")
"""

import logging
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

LOG = logging.getLogger("grimperium.crest_pm7.preoptimization")

ESTIMATED_TIME_PER_MOLECULE_SECONDS = 15.0


@dataclass
class xTBPreOptResult:
    """Result from xTB pre-optimization.

    Attributes:
        success: Whether pre-optimization completed successfully
        output_xyz: Path to pre-optimized XYZ file (None if failed)
        error_message: Error description if failed
        time_seconds: Execution time in seconds
    """

    success: bool
    output_xyz: Path | None = None
    error_message: str = ""
    time_seconds: float = 0.0


@dataclass
class xTBPreOptimizer:
    """GFN2-xTB pre-optimizer for molecular structures.

    Performs optional geometry pre-optimization using xTB before
    CREST conformer generation. This can improve convergence for
    poorly initialized structures.

    Attributes:
        enabled: Whether pre-optimization is active (default: False)
        timeout_seconds: Maximum execution time per molecule

    Example:
        >>> preopt = xTBPreOptimizer(enabled=True)
        >>> result = preopt.preoptimize_structure("benzene", input_xyz, work_dir)
        >>> if result.success:
        ...     # Use pre-optimized structure
        ...     optimized_xyz = result.output_xyz
    """

    enabled: bool = False
    timeout_seconds: int = 300
    _xtb_command: list[str] = field(default_factory=lambda: ["xtb"])

    def preoptimize_structure(
        self,
        mol_id: str,
        xyz_path: Path,
        work_dir: Path,
    ) -> xTBPreOptResult:
        """Run GFN2-xTB geometry optimization on input structure.

        Args:
            mol_id: Molecule identifier for naming output files
            xyz_path: Path to input XYZ file
            work_dir: Working directory for xTB output

        Returns:
            xTBPreOptResult with success status and output path

        Behavior:
            - If disabled: Returns original xyz_path as output
            - If enabled: Runs `xtb input.xyz --opt --gfn2 --temp 300`
            - Output file: {mol_id}_preopt.xyz in work_dir
        """
        if not self.enabled:
            LOG.debug(f"xTB pre-optimization disabled, returning original: {xyz_path}")
            return xTBPreOptResult(
                success=True,
                output_xyz=xyz_path,
                error_message="",
                time_seconds=0.0,
            )

        if not xyz_path.exists():
            error_msg = f"Input XYZ file not found: {xyz_path}"
            LOG.error(error_msg)
            return xTBPreOptResult(
                success=False,
                output_xyz=None,
                error_message=error_msg,
                time_seconds=0.0,
            )

        work_dir.mkdir(parents=True, exist_ok=True)
        output_xyz = work_dir / f"{mol_id}_preopt.xyz"

        cmd = [
            *self._xtb_command,
            str(xyz_path),
            "--opt",
            "--gfn2",
            "--temp",
            "300",
        ]

        LOG.info(f"Running xTB pre-optimization for {mol_id}")
        LOG.debug(f"Command: {' '.join(cmd)}")

        start_time = time.perf_counter()
        try:
            result = subprocess.run(
                cmd,
                cwd=work_dir,
                capture_output=True,
                text=True,
                timeout=self.timeout_seconds,
                check=False,
            )
            elapsed = time.perf_counter() - start_time

            xtbopt_output = work_dir / "xtbopt.xyz"
            if result.returncode == 0 and xtbopt_output.exists():
                xtbopt_output.rename(output_xyz)
                LOG.info(f"xTB pre-optimization complete: {output_xyz} ({elapsed:.1f}s)")
                return xTBPreOptResult(
                    success=True,
                    output_xyz=output_xyz,
                    error_message="",
                    time_seconds=elapsed,
                )
            else:
                error_msg = f"xTB failed (rc={result.returncode}): {result.stderr[:200]}"
                LOG.warning(error_msg)
                return xTBPreOptResult(
                    success=False,
                    output_xyz=None,
                    error_message=error_msg,
                    time_seconds=elapsed,
                )

        except FileNotFoundError as e:
            elapsed = time.perf_counter() - start_time
            error_msg = f"xTB executable not found: {e}"
            LOG.error(error_msg)
            return xTBPreOptResult(
                success=False,
                output_xyz=None,
                error_message=error_msg,
                time_seconds=elapsed,
            )

        except subprocess.TimeoutExpired:
            elapsed = time.perf_counter() - start_time
            error_msg = f"xTB timed out after {self.timeout_seconds}s"
            LOG.warning(error_msg)
            return xTBPreOptResult(
                success=False,
                output_xyz=None,
                error_message=error_msg,
                time_seconds=elapsed,
            )

    def get_estimated_time_hours(self, n_molecules: int) -> float:
        """Estimate total pre-optimization time.

        Args:
            n_molecules: Number of molecules to process

        Returns:
            Estimated hours (0.0 if disabled)
        """
        if not self.enabled:
            return 0.0
        return (n_molecules * ESTIMATED_TIME_PER_MOLECULE_SECONDS) / 3600.0
