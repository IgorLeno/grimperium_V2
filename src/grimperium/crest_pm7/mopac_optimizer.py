"""MOPAC PM7 optimization wrapper.

Handles MOPAC execution and robust output parsing.
"""

import logging
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .config import HOFConfidence, MOPACStatus, PM7Config
from .energy_extractor import extract_hof

LOG = logging.getLogger("grimperium.crest_pm7.mopac_optimizer")

# Mapeamento de multiplicidade para keywords do MOPAC
_MULTIPLICITY_KEYWORDS = {
    2: "DOUBLET",
    3: "TRIPLET",
    4: "QUARTET",
    5: "QUINTET",
    6: "SEXTET",
}


@dataclass
class MOPACResult:
    """Result of MOPAC PM7 optimization.

    Attributes:
        status: Execution status
        hof: Heat of formation (kcal/mol)
        hof_confidence: Confidence level of HOF extraction
        hof_method: Method used for HOF extraction
        execution_time: Execution time in seconds
        timeout_used: Timeout value used
        output_file: Path to MOPAC output file
        error_message: Error message if failed
    """

    status: MOPACStatus = MOPACStatus.NOT_ATTEMPTED
    hof: Optional[float] = None
    hof_confidence: HOFConfidence = HOFConfidence.LOW
    hof_method: Optional[str] = None
    execution_time: float = 0.0
    timeout_used: float = 0.0
    output_file: Optional[Path] = None
    error_message: Optional[str] = None


def _create_mopac_input(
    xyz_file: Path,
    output_mop: Path,
    charge: int = 0,
    multiplicity: int = 1,
) -> bool:
    """Create MOPAC input file from XYZ.

    Args:
        xyz_file: Input XYZ file
        output_mop: Output .mop file
        charge: Molecular charge
        multiplicity: Spin multiplicity (1=singlet, 2=doublet, etc.)

    Returns:
        True if input file created successfully

    Raises:
        ValueError: If multiplicity is not supported
    """
    try:
        with open(xyz_file, encoding="utf-8") as f:
            lines = f.readlines()

        n_atoms = int(lines[0].strip())
        comment = lines[1].strip() if len(lines) > 1 else ""

        # Build MOPAC keywords
        keywords = ["PM7", "PRECISE"]
        if charge != 0:
            keywords.append(f"CHARGE={charge}")

        if multiplicity != 1:
            if multiplicity not in _MULTIPLICITY_KEYWORDS:
                raise ValueError(
                    f"Unsupported multiplicity {multiplicity}. "
                    f"Supported values: 1 (singlet), {', '.join(f'{k} ({v.lower()})' for k, v in sorted(_MULTIPLICITY_KEYWORDS.items()))}"
                )
            keywords.append(_MULTIPLICITY_KEYWORDS[multiplicity])

        with open(output_mop, "w", encoding="utf-8") as f:
            f.write(" ".join(keywords) + "\n")
            f.write(f"{comment}\n")
            f.write("\n")  # Blank line

            # Write coordinates
            for line in lines[2:2 + n_atoms]:
                parts = line.split()
                if len(parts) >= 4:
                    element = parts[0]
                    x, y, z = parts[1], parts[2], parts[3]
                    # MOPAC format: element x flag y flag z flag
                    f.write(f"  {element:2s}  {x:>12s}  1  {y:>12s}  1  {z:>12s}  1\n")

        return True

    except Exception as e:
        LOG.warning(f"Failed to create MOPAC input: {e}")
        return False


def _detect_scf_failure(output_content: str) -> bool:
    """Detect SCF convergence failure in MOPAC output.

    Args:
        output_content: MOPAC output file content

    Returns:
        True if SCF failure detected
    """
    scf_patterns = [
        r"SCF\s+FAILED",
        r"UNABLE\s+TO\s+ACHIEVE\s+SCF",
        r"NO\s+CONVERGENCE",
        r"SCF\s+NOT\s+CONVERGED",
    ]
    for pattern in scf_patterns:
        if re.search(pattern, output_content, re.IGNORECASE):
            return True
    return False


def _detect_geometry_error(output_content: str) -> bool:
    """Detect geometry optimization error in MOPAC output.

    Args:
        output_content: MOPAC output file content

    Returns:
        True if geometry error detected
    """
    geo_patterns = [
        r"GEOMETRY\s+OPTIMIZATION\s+FAILED",
        r"ATOMS\s+TOO\s+CLOSE",
        r"GRADIENT\s+TOO\s+LARGE",
        r"ABNORMAL\s+TERMINATION",
    ]
    for pattern in geo_patterns:
        if re.search(pattern, output_content, re.IGNORECASE):
            return True
    return False


def run_mopac(
    mol_id: str,
    xyz_file: Path,
    config: PM7Config,
    timeout: float,
    nheavy: Optional[int] = None,
    work_dir: Optional[Path] = None,
    conf_index: int = 0,
) -> MOPACResult:
    """Run MOPAC PM7 optimization on a conformer.

    Args:
        mol_id: Molecule identifier
        xyz_file: Path to input XYZ file
        config: Pipeline configuration
        timeout: Timeout in seconds
        nheavy: Number of heavy atoms (for HOF validation)
        work_dir: Working directory (default: config.temp_dir/mol_id)
        conf_index: Conformer index

    Returns:
        MOPACResult with status and extracted energy
    """
    result = MOPACResult()
    result.timeout_used = timeout

    # Set up working directory
    if work_dir is None:
        work_dir = config.temp_dir / mol_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # Create input file
    mop_file = work_dir / f"{mol_id}_conf{conf_index:03d}.mop"
    out_file = work_dir / f"{mol_id}_conf{conf_index:03d}.out"

    if not _create_mopac_input(xyz_file, mop_file):
        result.status = MOPACStatus.NOT_ATTEMPTED
        result.error_message = "Failed to create MOPAC input file"
        return result

    start_time = time.time()

    try:
        # Run MOPAC
        cmd = [config.mopac_executable, str(mop_file)]
        LOG.debug(f"Running MOPAC for {mol_id} conf{conf_index}: {' '.join(cmd)}")

        proc = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        result.execution_time = time.time() - start_time

        # Check returncode of subprocess
        had_nonzero_exit = False
        if proc.returncode != 0:
            had_nonzero_exit = True
            error_msg = (
                f"MOPAC returned non-zero exit code {proc.returncode}. "
                f"cmd: {' '.join(cmd)}, stdout: {proc.stdout[:200] if proc.stdout else 'N/A'}, "
                f"stderr: {proc.stderr[:200] if proc.stderr else 'N/A'}"
            )
            LOG.warning(f"MOPAC non-zero exit for {mol_id} conf{conf_index}: {error_msg}")
            # Continue to attempt HOF extraction instead of returning immediately

        # MOPAC creates output with same name but .out extension
        if not out_file.exists():
            # Try arc file
            arc_file = mop_file.with_suffix(".arc")
            if arc_file.exists():
                LOG.info(
                    f"Using fallback arc file for {mol_id} conf{conf_index}: "
                    f"expected {out_file}, using {arc_file}"
                )
                out_file = arc_file

        if not out_file.exists():
            result.status = MOPACStatus.NOT_ATTEMPTED
            result.error_message = "No MOPAC output file found"
            return result

        result.output_file = out_file

        # Read and analyze output
        with open(out_file, encoding="utf-8", errors="replace") as f:
            output_content = f.read()

        # Check for errors
        if _detect_scf_failure(output_content):
            result.status = MOPACStatus.SCF_FAILED
            result.error_message = "SCF convergence failure"
            LOG.warning(f"MOPAC SCF failed for {mol_id} conf{conf_index}")
            # Still try to extract HOF if available
            hof, method, confidence = extract_hof(output_content, nheavy)
            if hof is not None:
                result.hof = hof
                result.hof_method = method
                result.hof_confidence = confidence
            return result

        if _detect_geometry_error(output_content):
            result.status = MOPACStatus.GEOMETRY_ERROR
            result.error_message = "Geometry optimization error"
            LOG.warning(f"MOPAC geometry error for {mol_id} conf{conf_index}")
            # Try to extract HOF even with geometry error (consistent with SCF branch)
            hof, method, confidence = extract_hof(output_content, nheavy)
            if hof is not None:
                result.hof = hof
                result.hof_method = method
                result.hof_confidence = confidence
            return result

        # Extract HOF
        hof, method, confidence = extract_hof(output_content, nheavy)

        if hof is None:
            # No HOF extracted - check if we had a non-zero exit
            if had_nonzero_exit:
                result.status = MOPACStatus.ERROR
                result.error_message = error_msg
            elif method is not None:
                # Pattern matched but validation failed
                result.status = MOPACStatus.SUCCESS
                result.hof_method = method
                result.hof_confidence = confidence
                result.error_message = "HOF extraction succeeded but validation failed"
            else:
                result.status = MOPACStatus.SUCCESS
                result.error_message = "No HOF value found in output"
            return result

        # HOF extracted successfully
        if had_nonzero_exit:
            # Had non-zero exit but managed to extract HOF - still an error
            result.status = MOPACStatus.ERROR
            result.error_message = error_msg
        else:
            # Success with HOF
            result.status = MOPACStatus.SUCCESS
        
        result.hof = hof
        result.hof_method = method
        result.hof_confidence = confidence

        LOG.info(
            f"MOPAC completed for {mol_id} conf{conf_index}: "
            f"HOF={hof:.2f} kcal/mol ({method}, {confidence.value}) "
            f"in {result.execution_time:.1f}s"
        )

    except subprocess.TimeoutExpired:
        result.execution_time = time.time() - start_time
        result.status = MOPACStatus.TIMEOUT
        result.error_message = f"MOPAC timeout after {timeout:.0f}s"
        LOG.warning(f"MOPAC timeout for {mol_id} conf{conf_index}")

    except Exception as e:
        result.execution_time = time.time() - start_time
        result.status = MOPACStatus.ERROR
        result.error_message = f"MOPAC error: {e}"
        LOG.error(f"MOPAC error for {mol_id} conf{conf_index}: {e}")

    return result


def optimize_conformer(
    mol_id: str,
    xyz_file: Path,
    config: PM7Config,
    timeout: float,
    nheavy: Optional[int] = None,
    conf_index: int = 0,
) -> MOPACResult:
    """Wrapper for run_mopac.

    This function exists for API stability and potential future extension.
    Currently forwards all arguments to run_mopac.

    Args:
        mol_id: Molecule identifier
        xyz_file: Path to conformer XYZ file
        config: Pipeline configuration
        timeout: Timeout in seconds
        nheavy: Number of heavy atoms
        conf_index: Conformer index

    Returns:
        MOPACResult
    """
    return run_mopac(
        mol_id=mol_id,
        xyz_file=xyz_file,
        config=config,
        timeout=timeout,
        nheavy=nheavy,
        conf_index=conf_index,
    )
