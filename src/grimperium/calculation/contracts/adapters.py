"""Adapters from legacy result models to canonical calculation contracts."""

from __future__ import annotations

from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Protocol
from uuid import uuid4

from .enums import ArtifactType, OverallStatus, PropertyRole, StageExecutionStatus
from .models import (
    CalculationArtifact,
    CalculationMethodReference,
    ConformerResult,
    HamiltonianResult,
    MolecularDescriptors,
    MoleculeCalculationResult,
    MoleculeData,
    PropertyEstimate,
    RunMetadata,
    StageExecutionRecord,
)
from .quantity import Quantity

PROPERTY_ID = "standard_enthalpy_of_formation"
LEGACY_METHOD_ID = "pm7_delta_legacy"
LEGACY_METHOD_VERSION = "0.0.0"
HAMILTONIAN = "PM7"


class ConformerData(Protocol):
    """Structural subset consumed from legacy ConformerData."""

    index: int
    mol_id: str
    crest_rank: int | None
    crest_geometry_file: Path | None
    mopac_output_file: Path | None
    mopac_error_message: str | None
    energy_hof: float | None
    is_successful: bool


class PM7Result(Protocol):
    """Structural subset consumed from legacy PM7Result."""

    mol_id: str
    smiles: str
    timestamp: datetime
    phase: str
    nrotbonds: int | None
    tpsa: float | None
    rdkit_descriptors: dict[str, float]
    conformers: Sequence[ConformerData]
    k_selected_pm7: int | None
    total_execution_time: float | None
    success: bool
    error_message: str | None

    @property
    def successful_conformers(self) -> Sequence[ConformerData]:
        """Return successfully processed conformers (read-only)."""

    def get_selected_conformer(self) -> ConformerData | None:
        """Return the legacy PM7-selected conformer."""


_RDKIT_DESCRIPTOR_MAP = {
    "rdkit_nrotbonds": "nrotbonds",
    "rdkit_tpsa": "tpsa",
    "rdkit_num_rings": "num_rings",
    "rdkit_fsp3": "fsp3",
    "rdkit_mol_weight": "mol_weight",
    "rdkit_hbond_donors": "hbond_donors",
    "rdkit_hbond_acceptors": "hbond_acceptors",
    "rdkit_nC": "nC",
    "rdkit_nH": "nH",
    "rdkit_nO": "nO",
    "rdkit_nN": "nN",
    "rdkit_bonds_single": "bonds_single",
    "rdkit_bonds_double": "bonds_double",
    "rdkit_bonds_triple": "bonds_triple",
    "rdkit_bonds_aromatic": "bonds_aromatic",
}

_ELECTRONIC_DESCRIPTOR_KEYS = [
    "mopac_dipole_debye",
    "mopac_ionization_potential_ev",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
    "mopac_cosmo_area_a2",
    "mopac_cosmo_volume_a3",
    "mopac_gradient_norm",
    "mopac_num_scf_cycles",
    "mopac_point_group",
    "mopac_time_s",
]


def _relative_path(path: Path) -> str:
    return path.as_posix().lstrip("/")


def _resolved_conformer_index(conformer: ConformerData) -> int:
    return conformer.crest_rank if conformer.crest_rank is not None else conformer.index


def _artifact_id(prefix: str, conformer: ConformerData) -> str:
    conformer_key = _resolved_conformer_index(conformer)
    return f"{prefix}-{conformer.mol_id}-{conformer_key}"


def _artifact_from_path(
    *,
    artifact_id: str,
    artifact_type: ArtifactType,
    path: Path | None,
    conformer: ConformerData,
    hamiltonian: str | None,
) -> CalculationArtifact | None:
    if path is None:
        return None
    return CalculationArtifact(
        artifact_id=artifact_id,
        artifact_type=artifact_type,
        relative_path=_relative_path(path),
        hamiltonian=hamiltonian,
        conformer_index=_resolved_conformer_index(conformer),
    )


def _overall_status(result: PM7Result) -> OverallStatus:
    if result.success:
        if result.error_message:
            return OverallStatus.PARTIAL
        return OverallStatus.SUCCESS
    if result.successful_conformers:
        return OverallStatus.PARTIAL
    return OverallStatus.FAILED


def _hamiltonian_status(conformer: ConformerData) -> OverallStatus:
    if conformer.is_successful:
        return OverallStatus.SUCCESS
    return OverallStatus.FAILED


def _molecular_descriptors(result: PM7Result) -> MolecularDescriptors | None:
    values: dict[str, float | None] = dict.fromkeys(_RDKIT_DESCRIPTOR_MAP.values())
    for legacy_key, target_key in _RDKIT_DESCRIPTOR_MAP.items():
        raw_value = result.rdkit_descriptors.get(legacy_key)
        values[target_key] = None if raw_value is None else float(raw_value)

    if result.nrotbonds is not None and values["nrotbonds"] is None:
        values["nrotbonds"] = float(result.nrotbonds)
    if result.tpsa is not None and values["tpsa"] is None:
        values["tpsa"] = float(result.tpsa)

    if all(value is None for value in values.values()):
        return None
    return MolecularDescriptors(**values)


def _electronic_descriptors(conformer: ConformerData) -> dict[str, float]:
    descriptors: dict[str, float] = {}
    for key in _ELECTRONIC_DESCRIPTOR_KEYS:
        value = getattr(conformer, key, None)
        if isinstance(value, int | float):
            descriptors[key] = float(value)
    return descriptors


def _stage_records(
    result: PM7Result, artifact_ids: list[str]
) -> list[StageExecutionRecord]:
    stage_status = (
        StageExecutionStatus.SUCCESS
        if result.successful_conformers
        else StageExecutionStatus.FAILED
    )
    return [
        StageExecutionRecord(
            stage_id="legacy_crest_pm7",
            program="grimperium.crest_pm7",
            role="baseline",
            status=stage_status,
            requested=True,
            started_at=result.timestamp,
            completed_at=result.timestamp,
            execution_time_s=result.total_execution_time,
            program_version=None,
            settings={"hamiltonian": HAMILTONIAN},
            artifact_ids=artifact_ids,
            error_message=result.error_message,
        )
    ]


def _conformer_result(
    conformer: ConformerData,
    artifact_ids: list[str],
) -> ConformerResult:
    energy = (
        Quantity(value=conformer.energy_hof, unit="kcal/mol")
        if conformer.energy_hof is not None
        else None
    )
    return ConformerResult(
        conformer_index=_resolved_conformer_index(conformer),
        crest_energy_hartree=None,
        hamiltonian_results=[
            HamiltonianResult(
                hamiltonian=HAMILTONIAN,
                status=_hamiltonian_status(conformer),
                energy_hof=energy,
                electronic_descriptors=_electronic_descriptors(conformer),
                conformer_index=_resolved_conformer_index(conformer),
                artifact_ids=artifact_ids,
                error_message=conformer.mopac_error_message,
            )
        ],
    )


def _artifacts_for_conformer(conformer: ConformerData) -> list[CalculationArtifact]:
    artifacts = []
    crest_artifact = _artifact_from_path(
        artifact_id=_artifact_id("crest-xyz", conformer),
        artifact_type=ArtifactType.CREST_XYZ,
        path=conformer.crest_geometry_file,
        conformer=conformer,
        hamiltonian=None,
    )
    if crest_artifact is not None:
        artifacts.append(crest_artifact)

    mopac_artifact = _artifact_from_path(
        artifact_id=_artifact_id("mopac-out", conformer),
        artifact_type=ArtifactType.MOPAC_OUT,
        path=conformer.mopac_output_file,
        conformer=conformer,
        hamiltonian=HAMILTONIAN,
    )
    if mopac_artifact is not None:
        artifacts.append(mopac_artifact)
    return artifacts


def pm7result_to_canonical(
    result: PM7Result,
    *,
    method_id: str = LEGACY_METHOD_ID,
    method_version: str = LEGACY_METHOD_VERSION,
    property_role: PropertyRole = PropertyRole.BASELINE,
) -> MoleculeCalculationResult:
    """Convert a legacy PM7Result to a canonical calculation result.

    The optional ``method_id``/``method_version``/``property_role`` parameters let
    higher-level runners (e.g. the PM7 + Delta Learning runner or the batch
    executor) attribute the PM7 estimate to the true method that produced it. The
    defaults preserve the original legacy-baseline behavior for existing callers.
    """
    run_id = str(uuid4())
    artifacts: list[CalculationArtifact] = []
    conformers: list[ConformerResult] = []

    for conformer in result.conformers:
        conformer_artifacts = _artifacts_for_conformer(conformer)
        artifacts.extend(conformer_artifacts)
        conformers.append(
            _conformer_result(
                conformer,
                [artifact.artifact_id for artifact in conformer_artifacts],
            )
        )

    estimates: list[PropertyEstimate] = []
    selected = result.get_selected_conformer()
    if selected is not None and selected.energy_hof is not None:
        estimates.append(
            PropertyEstimate(
                estimate_id=str(uuid4()),
                property_id=PROPERTY_ID,
                role=property_role,
                method_id=method_id,
                method_version=method_version,
                hamiltonian=HAMILTONIAN,
                value=Quantity(value=selected.energy_hof, unit="kcal/mol"),
                value_kcal_mol=selected.energy_hof,
                value_kj_mol=None,
                conformer_source_id=(
                    result.k_selected_pm7
                    if result.k_selected_pm7 is not None
                    else selected.crest_rank
                ),
                uncertainty=None,
                model_path=None,
            )
        )

    method_ref = CalculationMethodReference(
        method_id=method_id,
        method_version=method_version,
        property_id=PROPERTY_ID,
    )
    timestamp = result.timestamp
    if timestamp.tzinfo is None:
        timestamp = timestamp.replace(tzinfo=timezone.utc)

    artifact_ids = [artifact.artifact_id for artifact in artifacts]
    return MoleculeCalculationResult(
        molecule=MoleculeData(
            smiles=result.smiles,
            name=result.mol_id,
        ),
        run=RunMetadata(
            run_id=run_id,
            execution_phase=result.phase or "A",
            method_ref=method_ref,
            started_at=timestamp,
            completed_at=timestamp,
            grimperium_version=None,
        ),
        overall_status=_overall_status(result),
        conformers=conformers,
        molecular_descriptors=_molecular_descriptors(result),
        estimates=estimates,
        artifacts=artifacts,
        stage_executions=_stage_records(result, artifact_ids),
        error_message=result.error_message,
    )
