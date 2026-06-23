"""Semiempirical Method A runner for standard enthalpy of formation."""

from __future__ import annotations

import os
from collections.abc import Callable
from datetime import datetime, timezone
from importlib import import_module
from pathlib import Path
from typing import Any, Protocol, cast

from grimperium.calculation.contracts.enums import (
    ArtifactType,
    OverallStatus,
    PropertyRole,
    StageExecutionStatus,
)
from grimperium.calculation.contracts.models import (
    CalculationArtifact,
    CalculationMethodReference,
    ConformerResult,
    HamiltonianResult,
    MoleculeCalculationResult,
    MoleculeData,
    PropertyEstimate,
    RunMetadata,
    StageExecutionRecord,
)
from grimperium.calculation.contracts.quantity import Quantity
from grimperium.crest_pm7.config import CRESTStatus, MOPACStatus, PM7Config
from grimperium.crest_pm7.mopac_descriptors import extract_mopac_descriptors
from grimperium.crest_pm7.mopac_optimizer import MOPACResult, optimize_conformer
from grimperium.crest_pm7.preoptimization import xTBPreOptimizer

PROPERTY_ID = "standard_enthalpy_of_formation"
METHOD_ID = "semiempirical_am1_pm3_pm7"
METHOD_VERSION = "0.1.0"
HAMILTONIANS = ("AM1", "PM3", "PM7")


class XtbResultLike(Protocol):
    """Protocol for injected xTB pre-optimization results."""

    success: bool
    output_xyz: Path | None
    error_message: str
    time_seconds: float


class CrestResultLike(Protocol):
    """Protocol for CREST conformer generation results."""

    status: CRESTStatus
    conformers_found: int
    conformer_files: list[Path]
    execution_time: float
    error_message: str | None


GeometryGenerator = Callable[[str, Path, str], Path]
XtbPreoptimizer = Callable[[Path, Path, str], XtbResultLike]
CrestRunner = Callable[[str, Path, PM7Config, str | None], CrestResultLike]
MopacRunner = Callable[..., MOPACResult]
DescriptorExtractor = Callable[[Path], dict[str, Any]]


def _generate_initial_xyz(smiles: str, output_xyz: Path, name: str) -> Path:
    """Generate an initial 3D XYZ geometry using RDKit."""
    try:
        chem: Any = import_module("rdkit.Chem")
        all_chem: Any = import_module("rdkit.Chem.AllChem")
    except ImportError as exc:  # pragma: no cover - depends on optional environment
        raise RuntimeError("RDKit is required to generate initial geometries") from exc

    mol = chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse SMILES: {smiles}")
    mol = chem.AddHs(mol)
    status = all_chem.EmbedMolecule(mol, randomSeed=0xF00D)
    if status != 0:
        status = all_chem.EmbedMolecule(
            mol,
            useRandomCoords=True,
            randomSeed=0xF00D,
        )
    if status != 0:
        raise ValueError(f"RDKit failed to embed SMILES: {smiles}")
    all_chem.UFFOptimizeMolecule(mol, maxIters=500)

    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), name]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {pos.x:.8f} {pos.y:.8f} {pos.z:.8f}")

    output_xyz.parent.mkdir(parents=True, exist_ok=True)
    output_xyz.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_xyz


def _run_crest_default(
    mol_id: str,
    input_xyz: Path,
    config: PM7Config,
    smiles: str | None,
) -> CrestResultLike:
    module: Any = import_module("grimperium.crest_pm7.conformer_generator")
    return cast(CrestResultLike, module.run_crest(mol_id, input_xyz, config, smiles))


def _stage_record(
    *,
    stage_id: str,
    program: str,
    role: str,
    status: StageExecutionStatus,
    requested: bool,
    execution_time_s: float | None,
    settings: dict[str, Any],
    artifact_ids: list[str] | None = None,
    error_message: str | None = None,
    started_at: datetime | None = None,
    completed_at: datetime | None = None,
) -> StageExecutionRecord:
    return StageExecutionRecord(
        stage_id=stage_id,
        program=program,
        role=role,
        status=status,
        requested=requested,
        started_at=started_at,
        completed_at=completed_at,
        execution_time_s=execution_time_s,
        program_version=None,
        settings=settings,
        artifact_ids=artifact_ids or [],
        error_message=error_message,
    )


def _mopac_status_to_overall(status: MOPACStatus, hof: float | None) -> OverallStatus:
    if status is MOPACStatus.SUCCESS and hof is not None:
        return OverallStatus.SUCCESS
    return OverallStatus.FAILED


class SemiempiricalFormationEnthalpyRunner:
    """Run Method A: AM1/PM3/PM7 HOF on the lowest-CREST-energy conformer."""

    def __init__(
        self,
        *,
        config: PM7Config | None = None,
        work_root: Path | None = None,
        xtb_enabled: bool = True,
        geometry_generator: GeometryGenerator = _generate_initial_xyz,
        xtb_preoptimizer: XtbPreoptimizer | None = None,
        crest_runner: CrestRunner = _run_crest_default,
        mopac_runner: MopacRunner = optimize_conformer,
        descriptor_extractor: DescriptorExtractor = extract_mopac_descriptors,
    ) -> None:
        self.config = config or PM7Config()
        self.work_root = work_root or self.config.temp_dir
        self.xtb_enabled = xtb_enabled
        self.geometry_generator: GeometryGenerator = geometry_generator
        self.crest_runner: CrestRunner = crest_runner
        self.mopac_runner: MopacRunner = mopac_runner
        self.descriptor_extractor: DescriptorExtractor = descriptor_extractor
        if xtb_preoptimizer is None:
            optimizer = xTBPreOptimizer(enabled=xtb_enabled)
            self.xtb_preoptimizer: XtbPreoptimizer = (
                lambda input_xyz, work_dir, mol_id: optimizer.preoptimize_structure(
                    mol_id,
                    input_xyz,
                    work_dir,
                )
            )
        else:
            self.xtb_preoptimizer = xtb_preoptimizer

    def calculate_single_smiles(
        self,
        smiles: str,
        *,
        molecule_id: str,
        name: str | None = None,
        charge: int = 0,
        multiplicity: int = 1,
    ) -> MoleculeCalculationResult:
        """Calculate AM1, PM3, and PM7 formation enthalpy estimates."""
        started_at = datetime.now(timezone.utc)
        run_dir = self.work_root / molecule_id
        run_dir.mkdir(parents=True, exist_ok=True)
        molecule_name = name or molecule_id

        stage_executions: list[StageExecutionRecord] = []
        artifacts: list[CalculationArtifact] = []
        initial_xyz = self.geometry_generator(
            smiles,
            run_dir / f"{molecule_id}_initial.xyz",
            molecule_name,
        )
        stage_executions.append(
            _stage_record(
                stage_id="rdkit_initial_geometry",
                program="rdkit",
                role="initial_geometry",
                status=StageExecutionStatus.SUCCESS,
                requested=True,
                execution_time_s=None,
                settings={},
            )
        )

        crest_input_xyz = initial_xyz
        if self.xtb_enabled:
            xtb_started_at = datetime.now(timezone.utc)
            xtb_result = self.xtb_preoptimizer(
                initial_xyz, run_dir / "xtb", molecule_id
            )
            xtb_completed_at = datetime.now(timezone.utc)
            xtb_status = (
                StageExecutionStatus.SUCCESS
                if xtb_result.success and xtb_result.output_xyz is not None
                else StageExecutionStatus.FAILED
            )
            stage_executions.append(
                _stage_record(
                    stage_id="xtb_preoptimization",
                    program="xtb",
                    role="preoptimization",
                    status=xtb_status,
                    requested=True,
                    execution_time_s=xtb_result.time_seconds,
                    settings={"enabled": True},
                    error_message=xtb_result.error_message or None,
                    started_at=xtb_started_at,
                    completed_at=xtb_completed_at,
                )
            )
            if xtb_status is StageExecutionStatus.SUCCESS and xtb_result.output_xyz:
                crest_input_xyz = xtb_result.output_xyz
        else:
            stage_executions.append(
                _stage_record(
                    stage_id="xtb_preoptimization",
                    program="xtb",
                    role="preoptimization",
                    status=StageExecutionStatus.NOT_REQUESTED,
                    requested=False,
                    execution_time_s=None,
                    settings={"enabled": False},
                )
            )

        crest_result = self.crest_runner(
            molecule_id,
            crest_input_xyz,
            self.config,
            smiles,
        )
        if (
            crest_result.status is not CRESTStatus.SUCCESS
            or not crest_result.conformer_files
        ):
            return self._failed_result(
                smiles=smiles,
                name=name,
                charge=charge,
                multiplicity=multiplicity,
                started_at=started_at,
                stage_executions=[
                    *stage_executions,
                    _stage_record(
                        stage_id="crest_conformer_search",
                        program="crest",
                        role="conformer_search",
                        status=StageExecutionStatus.FAILED,
                        requested=True,
                        execution_time_s=crest_result.execution_time,
                        settings={"selection_strategy": "lowest_crest_energy"},
                        error_message=crest_result.error_message,
                    ),
                ],
                error_message=crest_result.error_message or "CREST failed",
            )

        selected_xyz = crest_result.conformer_files[0]
        stage_executions.append(
            _stage_record(
                stage_id="crest_conformer_search",
                program="crest",
                role="conformer_search",
                status=StageExecutionStatus.SUCCESS,
                requested=True,
                execution_time_s=crest_result.execution_time,
                settings={
                    "selection_strategy": "lowest_crest_energy",
                    "max_selected": 1,
                    "conformers_found": crest_result.conformers_found,
                },
            )
        )

        hamiltonian_results: list[HamiltonianResult] = []
        estimates: list[PropertyEstimate] = []
        for hamiltonian in HAMILTONIANS:
            mopac_started_at = datetime.now(timezone.utc)
            mopac_result = self.mopac_runner(
                mol_id=molecule_id,
                xyz_file=selected_xyz,
                config=self.config,
                timeout=self.config.mopac_timeout_base,
                conf_index=0,
                hamiltonian=hamiltonian,
                charge=charge,
                multiplicity=multiplicity,
            )
            mopac_completed_at = datetime.now(timezone.utc)
            artifact_ids = self._artifacts_for_mopac_result(
                artifacts,
                mopac_result,
                hamiltonian,
            )
            descriptors = (
                self._numeric_descriptors(mopac_result.output_file)
                if mopac_result.output_file is not None
                else {}
            )
            overall = _mopac_status_to_overall(mopac_result.status, mopac_result.hof)
            energy = (
                Quantity(value=mopac_result.hof, unit="kcal/mol")
                if mopac_result.hof is not None
                else None
            )
            hamiltonian_results.append(
                HamiltonianResult(
                    hamiltonian=hamiltonian,
                    status=overall,
                    energy_hof=energy,
                    electronic_descriptors=descriptors,
                    conformer_index=0,
                    artifact_ids=artifact_ids,
                    error_message=mopac_result.error_message,
                )
            )
            stage_executions.append(
                _stage_record(
                    stage_id=f"mopac_{hamiltonian.lower()}",
                    program="mopac",
                    role="semiempirical_hof",
                    status=(
                        StageExecutionStatus.SUCCESS
                        if overall is OverallStatus.SUCCESS
                        else StageExecutionStatus.FAILED
                    ),
                    requested=True,
                    execution_time_s=mopac_result.execution_time,
                    settings={"hamiltonian": hamiltonian},
                    artifact_ids=artifact_ids,
                    error_message=mopac_result.error_message,
                    started_at=mopac_started_at,
                    completed_at=mopac_completed_at,
                )
            )
            if energy is not None and overall is OverallStatus.SUCCESS:
                estimates.append(
                    PropertyEstimate(
                        estimate_id=f"{molecule_id}_{hamiltonian.lower()}_hof",
                        property_id=PROPERTY_ID,
                        role=PropertyRole.FINAL,
                        method_id=METHOD_ID,
                        method_version=METHOD_VERSION,
                        hamiltonian=hamiltonian,
                        value=energy,
                        value_kcal_mol=energy.value,
                        value_kj_mol=None,
                        conformer_source_id=0,
                        uncertainty=None,
                        model_path=None,
                    )
                )

        overall_status = (
            OverallStatus.SUCCESS
            if len(estimates) == len(HAMILTONIANS)
            else OverallStatus.PARTIAL if estimates else OverallStatus.FAILED
        )
        completed_at = datetime.now(timezone.utc)
        return MoleculeCalculationResult(
            molecule=MoleculeData(
                smiles=smiles,
                name=name,
                charge=charge,
                multiplicity=multiplicity,
            ),
            run=RunMetadata(
                run_id=molecule_id,
                execution_phase=str(self.config.phase.value),
                method_ref=CalculationMethodReference(
                    method_id=METHOD_ID,
                    method_version=METHOD_VERSION,
                    property_id=PROPERTY_ID,
                ),
                started_at=started_at,
                completed_at=completed_at,
                grimperium_version=None,
            ),
            overall_status=overall_status,
            conformers=[
                ConformerResult(
                    conformer_index=0,
                    crest_energy_hartree=None,
                    hamiltonian_results=hamiltonian_results,
                )
            ],
            molecular_descriptors=None,
            estimates=estimates,
            artifacts=artifacts,
            stage_executions=stage_executions,
            error_message=None if estimates else "No successful MOPAC estimates",
        )

    def _failed_result(
        self,
        *,
        smiles: str,
        name: str | None,
        charge: int,
        multiplicity: int,
        started_at: datetime,
        stage_executions: list[StageExecutionRecord],
        error_message: str,
    ) -> MoleculeCalculationResult:
        return MoleculeCalculationResult(
            molecule=MoleculeData(
                smiles=smiles,
                name=name,
                charge=charge,
                multiplicity=multiplicity,
            ),
            run=RunMetadata(
                run_id="",
                execution_phase=str(self.config.phase.value),
                method_ref=CalculationMethodReference(
                    method_id=METHOD_ID,
                    method_version=METHOD_VERSION,
                    property_id=PROPERTY_ID,
                ),
                started_at=started_at,
                completed_at=datetime.now(timezone.utc),
                grimperium_version=None,
            ),
            overall_status=OverallStatus.FAILED,
            conformers=[],
            molecular_descriptors=None,
            estimates=[],
            artifacts=[],
            stage_executions=stage_executions,
            error_message=error_message,
        )

    def _artifacts_for_mopac_result(
        self,
        artifacts: list[CalculationArtifact],
        mopac_result: MOPACResult,
        hamiltonian: str,
    ) -> list[str]:
        if mopac_result.output_file is None:
            return []
        artifact_id = f"mopac_{hamiltonian.lower()}_out"
        artifacts.append(
            CalculationArtifact(
                artifact_id=artifact_id,
                artifact_type=ArtifactType.MOPAC_OUT,
                relative_path=self._relative_artifact_path(mopac_result.output_file),
                hamiltonian=hamiltonian,
                conformer_index=0,
            )
        )
        return [artifact_id]

    def _relative_artifact_path(self, path: Path) -> str:
        try:
            return path.resolve().relative_to(self.work_root.resolve()).as_posix()
        except ValueError:
            return os.path.relpath(path, self.work_root).replace(os.sep, "/")

    def _numeric_descriptors(self, output_file: Path) -> dict[str, float]:
        descriptors = self.descriptor_extractor(output_file)
        numeric: dict[str, float] = {}
        for key, value in descriptors.items():
            if isinstance(value, int | float):
                numeric[key] = float(value)
        return numeric
