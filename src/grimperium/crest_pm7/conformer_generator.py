"""Conformer generation using CREST.

Wraps CREST executable and handles XYZ to SDF conversion.
"""

import logging
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .config import CRESTStatus, PM7Config

LOG = logging.getLogger("grimperium.crest_pm7.conformer_generator")


@dataclass
class CRESTResult:
    """Result of CREST conformer generation.

    Attributes:
        status: Execution status
        conformers_found: Number of conformers found
        conformer_files: List of XYZ file paths
        sdf_files: List of SDF file paths (after conversion)
        execution_time: Execution time in seconds
        error_message: Error message if failed
        work_dir: Working directory used
    """

    status: CRESTStatus = CRESTStatus.NOT_ATTEMPTED
    conformers_found: int = 0
    conformer_files: list[Path] = field(default_factory=list)
    sdf_files: list[Path] = field(default_factory=list)
    execution_time: float = 0.0
    error_message: Optional[str] = None
    work_dir: Optional[Path] = None


def _parse_xyz_file(xyz_path: Path) -> tuple[list[tuple[float, float, float]], int]:
    """Parse XYZ file to extract coordinates.

    Args:
        xyz_path: Path to XYZ file

    Returns:
        Tuple of (coordinates list, atom count)

    Raises:
        ValueError: If file is malformed or cannot be parsed
    """
    coords = []
    with open(xyz_path, encoding="utf-8") as f:
        lines = f.readlines()

    # Verify file has minimum lines
    if len(lines) < 2:
        raise ValueError(
            f"XYZ file '{xyz_path}' has insufficient lines ({len(lines)}). "
            "Expected at least 2 lines (atom count + comment)."
        )

    # Parse atom count
    try:
        n_atoms = int(lines[0].strip())
    except ValueError as e:
        raise ValueError(
            f"XYZ file '{xyz_path}' has invalid atom count header: '{lines[0].strip()}'. "
            f"Expected an integer. Error: {e}"
        ) from e

    # Validate n_atoms is positive
    if n_atoms <= 0:
        raise ValueError(
            f"XYZ file '{xyz_path}' has invalid atom count: {n_atoms}. "
            "Expected a positive integer."
        )

    # Skip comment line (line 1) and parse coordinates
    for line_idx, line in enumerate(lines[2 : 2 + n_atoms], start=2):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(
                f"XYZ file '{xyz_path}' has malformed coordinate line at index {line_idx}: "
                f"'{line.strip()}'. Expected at least 4 parts (element, x, y, z)."
            )
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            coords.append((x, y, z))
        except ValueError as e:
            raise ValueError(
                f"XYZ file '{xyz_path}' has invalid coordinate values at line {line_idx}: "
                f"'{line.strip()}'. Error: {e}"
            ) from e

    # Verify coordinate count matches
    if len(coords) != n_atoms:
        raise ValueError(
            f"XYZ file '{xyz_path}' coordinate count mismatch: "
            f"expected {n_atoms} atoms but parsed {len(coords)} coordinates."
        )

    return coords, n_atoms


def _convert_xyz_to_sdf_obabel(
    xyz_file: Path,
    sdf_file: Path,
) -> bool:
    """Convert XYZ to SDF using Open Babel.

    Args:
        xyz_file: Input XYZ file
        sdf_file: Output SDF file

    Returns:
        True if conversion succeeded
    """
    try:
        result = subprocess.run(
            ["obabel", str(xyz_file), "-O", str(sdf_file), "-h"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode == 0 and sdf_file.exists():
            LOG.debug(f"Converted {xyz_file} to {sdf_file} via obabel")
            return True
        LOG.debug(f"obabel conversion failed: {result.stderr}")
        return False
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        LOG.debug(f"obabel conversion error: {e}")
        return False


def _convert_xyz_to_sdf_rdkit(
    xyz_file: Path,
    sdf_file: Path,
    smiles: str,
) -> bool:
    """Convert XYZ to SDF using RDKit fallback.

    Requires SMILES to provide connectivity information.

    Args:
        xyz_file: Input XYZ file
        sdf_file: Output SDF file
        smiles: SMILES string for connectivity

    Returns:
        True if conversion succeeded
    """
    writer = None
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D

        # Create RDKit mol from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            LOG.warning(f"RDKit cannot parse SMILES: {smiles}")
            return False

        mol = Chem.AddHs(mol)

        # Parse XYZ coordinates
        coords, n_atoms_xyz = _parse_xyz_file(xyz_file)

        # Verify atom count
        if mol.GetNumAtoms() != n_atoms_xyz:
            LOG.warning(
                f"Atom count mismatch: RDKit={mol.GetNumAtoms()}, XYZ={n_atoms_xyz}"
            )
            return False

        # Create conformer with XYZ coordinates
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, Point3D(x, y, z))
        mol.AddConformer(conf, assignId=True)

        # Export SDF with try/finally to ensure writer is closed
        writer = Chem.SDWriter(str(sdf_file))
        try:
            writer.write(mol)
        finally:
            if writer is not None:
                writer.close()

        LOG.debug(f"Converted {xyz_file} to {sdf_file} via RDKit")
        return True

    except ImportError:
        LOG.warning("RDKit not available for XYZ->SDF conversion")
        return False
    except ValueError as e:
        LOG.warning(f"XYZ parsing error: {e}")
        return False
    except Exception as e:
        LOG.warning(f"RDKit conversion error: {e}")
        return False
    finally:
        # Garantir que writer seja fechado mesmo em caso de exceção
        if writer is not None:
            try:
                writer.close()
            except Exception:
                pass


def convert_xyz_to_sdf(
    xyz_file: Path,
    sdf_file: Path,
    smiles: Optional[str] = None,
) -> bool:
    """Convert XYZ to SDF with Open Babel + RDKit fallback.

    Args:
        xyz_file: Input XYZ file
        sdf_file: Output SDF file
        smiles: SMILES string (required for RDKit fallback)

    Returns:
        True if conversion succeeded
    """
    # Try Open Babel first
    if _convert_xyz_to_sdf_obabel(xyz_file, sdf_file):
        return True

    # Fallback to RDKit (requires SMILES)
    if smiles:
        return _convert_xyz_to_sdf_rdkit(xyz_file, sdf_file, smiles)

    LOG.warning("XYZ->SDF conversion failed: no obabel and no SMILES for RDKit")
    return False


def _split_crest_conformers(
    ensemble_file: Path,
    output_dir: Path,
    mol_id: str,
) -> list[Path]:
    """Split CREST ensemble file into individual XYZ files.

    Args:
        ensemble_file: crest_conformers.xyz with all conformers
        output_dir: Directory for output files
        mol_id: Molecule identifier

    Returns:
        List of paths to individual conformer XYZ files
    """
    output_files = []

    with open(ensemble_file, encoding="utf-8") as f:
        content = f.read()

    # Split by blank lines between conformers
    conformers = []
    lines = content.strip().split("\n")

    i = 0
    while i < len(lines):
        try:
            n_atoms = int(lines[i].strip())
            # Collect n_atoms + 2 lines (count, comment, atoms)
            conformer_lines = lines[i : i + n_atoms + 2]
            if len(conformer_lines) == n_atoms + 2:
                conformers.append("\n".join(conformer_lines))
                i += n_atoms + 2
            else:
                # Log quando conformer é pulado por ter número errado de linhas
                LOG.warning(
                    f"Skipping malformed conformer block at index {i} in {ensemble_file}: "
                    f"expected {n_atoms + 2} lines but got {len(conformer_lines)}. "
                    f"n_atoms={n_atoms}, lines preview: {lines[i:i+3]}"
                )
                i += 1
        except ValueError as e:
            # Log quando conformer é pulado por erro de parsing
            line_preview = lines[i][:50] if i < len(lines) else "N/A"
            LOG.warning(
                f"Skipping malformed conformer block at index {i} in {ensemble_file}: "
                f"could not parse atom count. Line: '{line_preview}'. Error: {e}"
            )
            i += 1

    for idx, conf_content in enumerate(conformers):
        output_file = output_dir / f"{mol_id}_conf{idx:03d}.xyz"
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(conf_content + "\n")
        output_files.append(output_file)

    LOG.debug(f"Split ensemble into {len(output_files)} conformer files")
    return output_files


def run_crest(
    mol_id: str,
    input_xyz: Path,
    config: PM7Config,
    smiles: Optional[str] = None,
) -> CRESTResult:
    """Run CREST conformer generation.

    Args:
        mol_id: Molecule identifier
        input_xyz: Path to input XYZ file
        config: Pipeline configuration
        smiles: SMILES string (for XYZ->SDF conversion)

    Returns:
        CRESTResult with status and output files
    """
    result = CRESTResult()

    # Create work directory
    work_dir = config.temp_dir / mol_id
    work_dir.mkdir(parents=True, exist_ok=True)
    result.work_dir = work_dir

    # Copy input file
    input_copy = work_dir / "input.xyz"
    shutil.copy(input_xyz, input_copy)

    start_time = time.time()

    try:
        # Run CREST
        cmd = [
            config.crest_executable,
            str(input_copy),
            "--quick",
            f"--ewin",
            str(config.energy_window),
        ]

        LOG.info(f"Running CREST for {mol_id}: {' '.join(cmd)}")

        proc = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=config.crest_timeout,
        )

        result.execution_time = time.time() - start_time

        if proc.returncode != 0:
            result.status = CRESTStatus.FAILED
            result.error_message = (
                f"CREST returned {proc.returncode}: {proc.stderr[:500]}"
            )
            LOG.warning(f"CREST failed for {mol_id}: {result.error_message}")
            return result

        # Check for output file
        ensemble_file = work_dir / "crest_conformers.xyz"
        if not ensemble_file.exists():
            # Try alternative name
            ensemble_file = work_dir / "crest_best.xyz"

        if not ensemble_file.exists():
            result.status = CRESTStatus.FAILED
            result.error_message = "No conformer output file found"
            LOG.warning(f"CREST output not found for {mol_id}")
            return result

        # Split conformers
        conformer_files = _split_crest_conformers(ensemble_file, work_dir, mol_id)

        if not conformer_files:
            result.status = CRESTStatus.FAILED
            result.error_message = "No conformers extracted from output"
            return result

        result.conformers_found = len(conformer_files)
        result.conformer_files = conformer_files

        # Convert to SDF
        sdf_dir = work_dir / "sdf"
        sdf_dir.mkdir(exist_ok=True)

        for xyz_file in conformer_files:
            sdf_file = sdf_dir / xyz_file.with_suffix(".sdf").name
            if convert_xyz_to_sdf(xyz_file, sdf_file, smiles):
                result.sdf_files.append(sdf_file)

        result.status = CRESTStatus.SUCCESS
        LOG.info(
            f"CREST completed for {mol_id}: "
            f"{result.conformers_found} conformers in {result.execution_time:.1f}s"
        )

    except subprocess.TimeoutExpired:
        result.execution_time = time.time() - start_time
        result.status = CRESTStatus.FAILED
        result.error_message = f"CREST timeout after {config.crest_timeout}s"
        LOG.warning(f"CREST timeout for {mol_id}")

    except Exception as e:
        result.execution_time = time.time() - start_time
        result.status = CRESTStatus.FAILED
        result.error_message = f"CREST error: {e}"
        LOG.error(f"CREST error for {mol_id}: {e}")

    return result
