"""Molecule processing for CREST-PM7 Pipeline.

Defines result dataclasses and the main MoleculeProcessor class.
"""

import logging
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from .config import (
    CRESTStatus,
    HOFConfidence,
    MOPACStatus,
    PM7Config,
    QualityGrade,
    TimeoutConfidence,
)
from .conformer_generator import run_crest
from .conformer_selector import calculate_delta_e, get_num_conformers
from .mopac_optimizer import optimize_conformer
from .timeout_predictor import TimeoutPredictor

LOG = logging.getLogger("grimperium.crest_pm7.molecule_processor")


@dataclass
class ConformerData:
    """Data for a single conformer.

    Attributes:
        index: Conformer index
        mol_id: Parent molecule ID
        crest_status: CREST generation status
        crest_geometry_file: Path to XYZ file
        crest_error_message: CREST error if any
        mopac_status: MOPAC optimization status
        mopac_output_file: Path to MOPAC output
        mopac_execution_time: MOPAC execution time
        mopac_timeout_used: Timeout value used
        mopac_error_message: MOPAC error if any
        energy_hof: Heat of formation
        hof_confidence: HOF extraction confidence
        hof_extraction_method: Method used for HOF extraction
        hof_extraction_successful: Whether HOF was extracted
    """

    index: int
    mol_id: str

    # CREST
    crest_status: CRESTStatus = CRESTStatus.NOT_ATTEMPTED
    crest_geometry_file: Path | None = None
    crest_error_message: str | None = None

    # MOPAC
    mopac_status: MOPACStatus = MOPACStatus.NOT_ATTEMPTED
    mopac_output_file: Path | None = None
    mopac_execution_time: float = 0.0
    mopac_timeout_used: float | None = None
    mopac_error_message: str | None = None

    # Energy
    energy_hof: float | None = None
    hof_confidence: HOFConfidence = HOFConfidence.LOW
    hof_extraction_method: str | None = None
    hof_extraction_successful: bool = False

    @property
    def is_successful(self) -> bool:
        """Whether this conformer was successfully processed."""
        return (
            self.mopac_status == MOPACStatus.SUCCESS
            and self.hof_extraction_successful
            and self.energy_hof is not None
        )

    def to_dict(self) -> dict[str, Any]:
        """Convert to JSON-safe dictionary."""
        return {
            "index": self.index,
            "mol_id": self.mol_id,
            "crest_status": self.crest_status.value,
            "crest_geometry_file": (
                str(self.crest_geometry_file) if self.crest_geometry_file else None
            ),
            "crest_error_message": self.crest_error_message,
            "mopac_status": self.mopac_status.value,
            "mopac_output_file": (
                str(self.mopac_output_file) if self.mopac_output_file else None
            ),
            "mopac_execution_time": self.mopac_execution_time,
            "mopac_timeout_used": self.mopac_timeout_used,
            "mopac_error_message": self.mopac_error_message,
            "energy_hof": self.energy_hof,
            # hof_confidence is non-Optional with a default, always serializes its value
            "hof_confidence": self.hof_confidence.value,
            "hof_extraction_method": self.hof_extraction_method,
            "hof_extraction_successful": self.hof_extraction_successful,
            "is_successful": self.is_successful,
        }


@dataclass
class PM7Result:
    """Complete result for a molecule processing.

    Attributes:
        mol_id: Molecule identifier
        smiles: SMILES string
        timestamp: Processing timestamp (timezone-aware UTC)
        phase: Processing phase
        nheavy: Number of heavy atoms
        nrotbonds: Number of rotatable bonds
        tpsa: Topological polar surface area
        aromatic_rings: Number of aromatic rings
        has_heteroatoms: Whether molecule has heteroatoms
        crest_status: Overall CREST status
        crest_conformers_generated: Number of conformers from CREST
        crest_time: CREST execution time
        crest_error: CREST error message
        conformers: List of ConformerData
        num_conformers_selected: Number of conformers selected for MOPAC
        total_execution_time: Total processing time
        delta_e_12: Energy difference conf1-conf2
        delta_e_13: Energy difference conf1-conf3
        delta_e_15: Energy difference conf1-conf5
        timeout_predicted: Predicted timeout
        timeout_confidence: Timeout prediction confidence
        decisions: List of decision audit trail
        quality_grade: Quality grade
        issues: List of issues found
        success: Whether processing succeeded
        error_message: Error message if failed
    """

    # Identification
    mol_id: str
    smiles: str
    timestamp: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    phase: str = "A"

    # Metadata
    nheavy: int | None = None
    nrotbonds: int | None = None
    tpsa: float | None = None
    aromatic_rings: int | None = None
    has_heteroatoms: bool | None = None

    # CREST
    crest_status: CRESTStatus = CRESTStatus.NOT_ATTEMPTED
    crest_conformers_generated: int = 0
    crest_time: float | None = None
    crest_error: str | None = None

    # MOPAC
    conformers: list[ConformerData] = field(default_factory=list)
    num_conformers_selected: int | None = None
    total_execution_time: float | None = None

    # Energy differences
    delta_e_12: float | None = None
    delta_e_13: float | None = None
    delta_e_15: float | None = None

    # Timeout
    timeout_predicted: float | None = None
    timeout_confidence: TimeoutConfidence | None = None

    # Decisions (audit)
    decisions: list[str] = field(default_factory=list)

    # Quality
    quality_grade: QualityGrade = QualityGrade.FAILED
    issues: list[str] = field(default_factory=list)

    # Final
    success: bool = False
    error_message: str | None = None

    @property
    def most_stable_hof(self) -> float | None:
        """Get HOF of most stable conformer."""
        successful = [c for c in self.conformers if c.is_successful]
        if not successful:
            return None
        # is_successful guarantees energy_hof is not None
        energies: list[float] = [c.energy_hof for c in successful]  # type: ignore[misc]
        return min(energies)

    @property
    def successful_conformers(self) -> list[ConformerData]:
        """Get list of successfully processed conformers."""
        return [c for c in self.conformers if c.is_successful]

    def to_dict(self) -> dict[str, Any]:
        """Convert to JSON-safe dictionary."""
        return {
            "mol_id": self.mol_id,
            "smiles": self.smiles,
            "timestamp": self.timestamp.isoformat(),
            "phase": self.phase,
            "nheavy": self.nheavy,
            "nrotbonds": self.nrotbonds,
            "tpsa": self.tpsa,
            "aromatic_rings": self.aromatic_rings,
            "has_heteroatoms": self.has_heteroatoms,
            "crest_status": self.crest_status.value,
            "crest_conformers_generated": self.crest_conformers_generated,
            "crest_time": self.crest_time,
            "crest_error": self.crest_error,
            "conformers": [c.to_dict() for c in self.conformers],
            "num_conformers_selected": self.num_conformers_selected,
            "total_execution_time": self.total_execution_time,
            "delta_e_12": self.delta_e_12,
            "delta_e_13": self.delta_e_13,
            "delta_e_15": self.delta_e_15,
            "timeout_predicted": self.timeout_predicted,
            "timeout_confidence": (
                self.timeout_confidence.value if self.timeout_confidence else None
            ),
            "decisions": self.decisions,
            "quality_grade": self.quality_grade.value,
            "issues": self.issues,
            "success": self.success,
            "error_message": self.error_message,
            "most_stable_hof": self.most_stable_hof,
        }


def compute_molecular_descriptors(smiles: str) -> dict[str, Any]:
    """Compute molecular descriptors from SMILES.

    Args:
        smiles: SMILES string

    Returns:
        Dictionary with nheavy, nrotbonds, tpsa, aromatic_rings, has_heteroatoms
    """
    result: dict[str, Any] = {
        "nheavy": None,
        "nrotbonds": None,
        "tpsa": None,
        "aromatic_rings": None,
        "has_heteroatoms": None,
    }

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return result

        result["nheavy"] = Descriptors.HeavyAtomCount(mol)
        result["nrotbonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        result["tpsa"] = Descriptors.TPSA(mol)
        result["aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)

        # Check for heteroatoms (N, O, S, P, etc.)
        hetero_atoms = {"N", "O", "S", "P", "F", "Cl", "Br", "I"}
        has_hetero: bool = any(
            atom.GetSymbol() in hetero_atoms for atom in mol.GetAtoms()
        )
        result["has_heteroatoms"] = has_hetero

    except Exception as e:
        LOG.warning(f"Failed to compute descriptors for {smiles}: {e}")

    return result


def _collect_issues(result: PM7Result) -> list[str]:
    """Collect quality issues from a PM7Result.

    Pure helper function that analyzes result and returns list of issues.

    Args:
        result: PM7Result to analyze

    Returns:
        List of issue strings
    """
    if not result.success:
        return []  # No detailed issues for failed results

    successful = result.successful_conformers
    if not successful:
        return []

    issues = []

    # Check HOF confidence
    high_conf = sum(1 for c in successful if c.hof_confidence == HOFConfidence.HIGH)
    if high_conf == 0:
        issues.append("no_high_confidence_hof")

    # Check conformer coverage
    if (
        result.num_conformers_selected
        and len(successful) < result.num_conformers_selected
    ):
        issues.append("incomplete_conformer_coverage")

    # Check timeout
    if result.timeout_confidence == TimeoutConfidence.LOW:
        issues.append("low_timeout_confidence")

    return issues


def _grade_from_issues(
    issues: list[str], success: bool, has_conformers: bool
) -> QualityGrade:
    """Determine quality grade based on issues.

    Pure helper function that returns grade based solely on input parameters.

    Args:
        issues: List of issue strings
        success: Whether processing succeeded
        has_conformers: Whether any conformers were successful

    Returns:
        QualityGrade
    """
    if not success:
        return QualityGrade.FAILED

    if not has_conformers:
        return QualityGrade.FAILED

    if len(issues) == 0:
        return QualityGrade.A
    elif len(issues) == 1:
        return QualityGrade.B
    else:
        return QualityGrade.C


class MoleculeProcessor:
    """Processes a single molecule through CREST + MOPAC pipeline."""

    def __init__(
        self,
        config: PM7Config,
        timeout_predictor: TimeoutPredictor | None = None,
    ) -> None:
        """Initialize processor.

        Args:
            config: Pipeline configuration
            timeout_predictor: Optional timeout predictor
        """
        self.config = config
        self.timeout_predictor = timeout_predictor or TimeoutPredictor(
            recalibrate_interval=config.timeout_predictor_recalibrate_interval
        )

    def _assign_quality_grade(self, result: PM7Result) -> QualityGrade:
        """Assign quality grade based on results.

        Uses pure helper functions for issue collection and grading.

        Args:
            result: PM7Result to grade

        Returns:
            QualityGrade
        """
        # Collect issues using pure helper
        issues = _collect_issues(result)

        # Assign issues to result (caller may also do this)
        result.issues = issues

        # Determine grade using pure helper
        return _grade_from_issues(
            issues,
            success=result.success,
            has_conformers=len(result.successful_conformers) > 0,
        )

    def process(
        self,
        mol_id: str,
        smiles: str,
        input_xyz: Path | None = None,
    ) -> PM7Result:
        """Process a single molecule.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string
            input_xyz: Optional input XYZ (if None, generate from SMILES)

        Returns:
            PM7Result with all processing results
        """
        start_time = time.time()

        # Get phase value for string representation
        phase_value = (
            self.config.phase.value
            if hasattr(self.config.phase, "value")
            else str(self.config.phase)
        )

        result = PM7Result(
            mol_id=mol_id,
            smiles=smiles,
            phase=phase_value,
        )

        LOG.info(f"Processing molecule: {mol_id} ({smiles})")

        # Compute descriptors
        descriptors = compute_molecular_descriptors(smiles)
        result.nheavy = descriptors["nheavy"]
        result.nrotbonds = descriptors["nrotbonds"]
        result.tpsa = descriptors["tpsa"]
        result.aromatic_rings = descriptors["aromatic_rings"]
        result.has_heteroatoms = descriptors["has_heteroatoms"]

        if result.nheavy is None:
            result.error_message = "Failed to parse SMILES"
            result.total_execution_time = time.time() - start_time
            return result

        result.decisions.append(f"nheavy={result.nheavy}, nrotbonds={result.nrotbonds}")

        # Generate initial XYZ if not provided
        if input_xyz is None:
            input_xyz = self._generate_xyz_from_smiles(mol_id, smiles)
            if input_xyz is None:
                result.error_message = "Failed to generate initial 3D coordinates"
                result.total_execution_time = time.time() - start_time
                return result

        # Determine number of conformers
        num_conformers, reason = get_num_conformers(
            result.nrotbonds or 0,
            self.config,
        )
        result.num_conformers_selected = num_conformers
        result.decisions.append(f"conformers={num_conformers}: {reason}")

        # Predict timeout
        timeout, confidence = self.timeout_predictor.predict(
            result.nheavy,
            num_conformers,
        )
        result.timeout_predicted = timeout
        result.timeout_confidence = confidence
        result.decisions.append(f"timeout={timeout:.0f}s ({confidence.value})")

        # Run CREST
        crest_result = run_crest(
            mol_id=mol_id,
            input_xyz=input_xyz,
            config=self.config,
            smiles=smiles,
        )

        result.crest_status = crest_result.status
        result.crest_conformers_generated = crest_result.conformers_found
        result.crest_time = crest_result.execution_time
        result.crest_error = crest_result.error_message

        if crest_result.status != CRESTStatus.SUCCESS:
            result.error_message = f"CREST failed: {crest_result.error_message}"
            result.num_conformers_selected = 0  # No conformers when CREST fails
            result.total_execution_time = time.time() - start_time
            return result

        # Select conformers
        conformers_to_process = crest_result.conformer_files[:num_conformers]
        # Update num_conformers_selected to reflect actual count (may be less than requested)
        result.num_conformers_selected = len(conformers_to_process)
        result.decisions.append(
            f"processing {len(conformers_to_process)} of {crest_result.conformers_found} conformers"
        )

        # Run MOPAC on each conformer with dynamic timeout redistribution
        remaining_timeout = timeout
        remaining_conformers = len(conformers_to_process)

        for idx, xyz_file in enumerate(conformers_to_process):
            # Calculate timeout for this conformer based on remaining time
            per_conformer_timeout = remaining_timeout / max(1, remaining_conformers)

            conf_data = ConformerData(index=idx, mol_id=mol_id)
            conf_data.crest_status = CRESTStatus.SUCCESS
            conf_data.crest_geometry_file = xyz_file

            mopac_result = optimize_conformer(
                mol_id=mol_id,
                xyz_file=xyz_file,
                config=self.config,
                timeout=per_conformer_timeout,
                nheavy=result.nheavy,
                conf_index=idx,
            )

            conf_data.mopac_status = mopac_result.status
            conf_data.mopac_output_file = mopac_result.output_file
            conf_data.mopac_execution_time = mopac_result.execution_time
            conf_data.mopac_timeout_used = mopac_result.timeout_used
            conf_data.mopac_error_message = mopac_result.error_message

            if mopac_result.hof is not None:
                conf_data.energy_hof = mopac_result.hof
                conf_data.hof_confidence = mopac_result.hof_confidence
                conf_data.hof_extraction_method = mopac_result.hof_method
                conf_data.hof_extraction_successful = True

                # Record for timeout predictor
                self.timeout_predictor.add_observation(
                    result.nheavy,
                    mopac_result.execution_time,
                )

            result.conformers.append(conf_data)

            # Update remaining time for next conformers
            remaining_timeout = max(0, remaining_timeout - mopac_result.execution_time)
            remaining_conformers -= 1

        # Calculate energy differences - is_successful guarantees energy_hof is not None
        successful_energies: list[float] = [
            c.energy_hof
            for c in result.conformers
            if c.is_successful and c.energy_hof is not None
        ]
        if successful_energies:
            deltas = calculate_delta_e(successful_energies)
            result.delta_e_12 = deltas["delta_e_12"]
            result.delta_e_13 = deltas["delta_e_13"]
            result.delta_e_15 = deltas["delta_e_15"]

        # Determine success
        result.success = len(result.successful_conformers) > 0

        # Assign quality grade
        result.quality_grade = self._assign_quality_grade(result)

        result.total_execution_time = time.time() - start_time

        LOG.info(
            f"Completed {mol_id}: success={result.success}, "
            f"grade={result.quality_grade.value}, "
            f"HOF={result.most_stable_hof}, "
            f"time={result.total_execution_time:.1f}s"
        )

        return result

    def _generate_xyz_from_smiles(
        self,
        mol_id: str,
        smiles: str,
    ) -> Path | None:
        """Generate initial XYZ coordinates from SMILES.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string

        Returns:
            Path to XYZ file or None if failed
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            if embed_result < 0:
                # Try with random seed
                embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
                if embed_result < 0:
                    return None

            # Optimize geometry and check result
            optim_result = AllChem.MMFFOptimizeMolecule(mol)
            if optim_result == -1:
                # -1 means missing force field parameters
                mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else smiles[:50]
                LOG.warning(
                    f"MMFF optimization failed for '{mol_name}' (return code -1): "
                    "missing force field parameters. Geometry may be unreliable."
                )
            elif optim_result != 0:
                mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else smiles[:50]
                LOG.warning(
                    f"MMFF optimization returned non-zero code {optim_result} for '{mol_name}'. "
                    "Geometry may not be fully optimized."
                )

            # Write XYZ
            work_dir = self.config.temp_dir / mol_id
            work_dir.mkdir(parents=True, exist_ok=True)
            xyz_path = work_dir / f"{mol_id}_input.xyz"

            conf = mol.GetConformer()
            with open(xyz_path, "w", encoding="utf-8") as f:
                f.write(f"{mol.GetNumAtoms()}\n")
                f.write(f"{mol_id} generated from SMILES\n")
                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

            return xyz_path

        except Exception as e:
            LOG.warning(f"Failed to generate XYZ from SMILES: {e}")
            return None
