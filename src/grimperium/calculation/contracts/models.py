"""Dataclass models for the canonical calculation contract."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import PurePosixPath, PureWindowsPath
from typing import Any

from .enums import ArtifactType, OverallStatus, PropertyRole, StageExecutionStatus
from .quantity import Quantity


@dataclass(kw_only=True)
class StageExecutionRecord:
    """Execution metadata for one calculation stage."""

    stage_id: str
    program: str
    role: str
    status: StageExecutionStatus
    requested: bool
    started_at: datetime | None
    completed_at: datetime | None
    execution_time_s: float | None
    program_version: str | None
    settings: dict[str, Any]
    artifact_ids: list[str]
    error_message: str | None


@dataclass(kw_only=True)
class CalculationArtifact:
    """A file produced by a calculation run."""

    artifact_id: str
    artifact_type: ArtifactType
    relative_path: str
    hamiltonian: str | None
    conformer_index: int | None

    def __post_init__(self) -> None:
        if PurePosixPath(self.relative_path).is_absolute() or PureWindowsPath(
            self.relative_path
        ).is_absolute():
            raise ValueError("CalculationArtifact.relative_path must not be absolute")


@dataclass(kw_only=True)
class MoleculeData:
    """Molecule identity and charge state."""

    smiles: str
    name: str | None
    charge: int = 0
    multiplicity: int = 1


@dataclass(kw_only=True)
class MolecularDescriptors:
    """Legacy RDKit descriptors at molecule level."""

    nrotbonds: float | None = None
    tpsa: float | None = None
    num_rings: float | None = None
    fsp3: float | None = None
    mol_weight: float | None = None
    hbond_donors: float | None = None
    hbond_acceptors: float | None = None
    nC: float | None = None
    nH: float | None = None
    nO: float | None = None
    nN: float | None = None
    bonds_single: float | None = None
    bonds_double: float | None = None
    bonds_triple: float | None = None
    bonds_aromatic: float | None = None


@dataclass(kw_only=True)
class HamiltonianResult:
    """Result for one Hamiltonian on one conformer."""

    hamiltonian: str
    status: OverallStatus
    energy_hof: Quantity | None
    electronic_descriptors: dict[str, float]
    conformer_index: int | None
    artifact_ids: list[str]
    error_message: str | None


@dataclass(kw_only=True)
class ConformerResult:
    """Canonical conformer result."""

    conformer_index: int
    crest_energy_hartree: float | None
    hamiltonian_results: list[HamiltonianResult]


@dataclass(kw_only=True)
class PropertyEstimate:
    """Canonical property estimate."""

    estimate_id: str
    property_id: str
    role: PropertyRole
    method_id: str
    method_version: str
    hamiltonian: str | None
    value: Quantity
    value_kcal_mol: float | None
    value_kj_mol: float | None
    conformer_source_id: int | None
    uncertainty: float | None
    model_path: str | None


@dataclass(kw_only=True)
class CalculationMethodReference:
    """Reference to the calculation method definition."""

    method_id: str
    method_version: str
    property_id: str


@dataclass(kw_only=True)
class RunMetadata:
    """Run metadata for a molecule calculation."""

    run_id: str
    execution_phase: str
    method_ref: CalculationMethodReference
    started_at: datetime | None
    completed_at: datetime | None
    grimperium_version: str | None


@dataclass(kw_only=True)
class MoleculeCalculationResult:
    """Canonical calculation result for one molecule."""

    molecule: MoleculeData
    run: RunMetadata
    overall_status: OverallStatus
    conformers: list[ConformerResult]
    molecular_descriptors: MolecularDescriptors | None
    estimates: list[PropertyEstimate]
    artifacts: list[CalculationArtifact]
    stage_executions: list[StageExecutionRecord]
    error_message: str | None = None
    schema_version: int = field(default=1)
