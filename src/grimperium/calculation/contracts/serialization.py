"""Manual serialization for canonical calculation contracts."""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime
from typing import Any

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

SUPPORTED_SCHEMA_VERSION = 1


def _datetime_to_str(value: datetime | None) -> str | None:
    return value.isoformat() if value is not None else None


def _datetime_from_str(value: str | None) -> datetime | None:
    return datetime.fromisoformat(value) if value is not None else None


def _float_or_none(value: Any) -> float | None:
    return None if value is None else float(value)


def _int_or_none(value: Any) -> int | None:
    return None if value is None else int(value)


def _bool_from_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().casefold()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no"}:
            return False
        msg = f"Cannot convert string to bool: {value!r}"
        raise ValueError(msg)
    return bool(value)


def _stage_execution_to_dict(record: StageExecutionRecord) -> dict[str, Any]:
    return {
        "stage_id": record.stage_id,
        "program": record.program,
        "role": record.role,
        "status": record.status.value,
        "requested": record.requested,
        "started_at": _datetime_to_str(record.started_at),
        "completed_at": _datetime_to_str(record.completed_at),
        "execution_time_s": record.execution_time_s,
        "program_version": record.program_version,
        "settings": deepcopy(record.settings),
        "artifact_ids": list(record.artifact_ids),
        "error_message": record.error_message,
    }


def _stage_execution_from_dict(data: dict[str, Any]) -> StageExecutionRecord:
    return StageExecutionRecord(
        stage_id=str(data["stage_id"]),
        program=str(data["program"]),
        role=str(data["role"]),
        status=StageExecutionStatus(data["status"]),
        requested=_bool_from_value(data["requested"]),
        started_at=_datetime_from_str(data["started_at"]),
        completed_at=_datetime_from_str(data["completed_at"]),
        execution_time_s=_float_or_none(data["execution_time_s"]),
        program_version=data["program_version"],
        settings=dict(data["settings"]),
        artifact_ids=list(data["artifact_ids"]),
        error_message=data["error_message"],
    )


def _artifact_to_dict(artifact: CalculationArtifact) -> dict[str, Any]:
    return {
        "artifact_id": artifact.artifact_id,
        "artifact_type": artifact.artifact_type.value,
        "relative_path": artifact.relative_path,
        "hamiltonian": artifact.hamiltonian,
        "conformer_index": artifact.conformer_index,
    }


def _artifact_from_dict(data: dict[str, Any]) -> CalculationArtifact:
    return CalculationArtifact(
        artifact_id=str(data["artifact_id"]),
        artifact_type=ArtifactType(data["artifact_type"]),
        relative_path=str(data["relative_path"]),
        hamiltonian=data["hamiltonian"],
        conformer_index=_int_or_none(data["conformer_index"]),
    )


def _molecular_descriptors_to_dict(
    descriptors: MolecularDescriptors | None,
) -> dict[str, Any] | None:
    if descriptors is None:
        return None
    return {
        "nrotbonds": descriptors.nrotbonds,
        "tpsa": descriptors.tpsa,
        "num_rings": descriptors.num_rings,
        "fsp3": descriptors.fsp3,
        "mol_weight": descriptors.mol_weight,
        "hbond_donors": descriptors.hbond_donors,
        "hbond_acceptors": descriptors.hbond_acceptors,
        "nC": descriptors.nC,
        "nH": descriptors.nH,
        "nO": descriptors.nO,
        "nN": descriptors.nN,
        "bonds_single": descriptors.bonds_single,
        "bonds_double": descriptors.bonds_double,
        "bonds_triple": descriptors.bonds_triple,
        "bonds_aromatic": descriptors.bonds_aromatic,
    }


def _molecular_descriptors_from_dict(
    data: dict[str, Any] | None,
) -> MolecularDescriptors | None:
    if data is None:
        return None
    return MolecularDescriptors(
        nrotbonds=_float_or_none(data["nrotbonds"]),
        tpsa=_float_or_none(data["tpsa"]),
        num_rings=_float_or_none(data["num_rings"]),
        fsp3=_float_or_none(data["fsp3"]),
        mol_weight=_float_or_none(data["mol_weight"]),
        hbond_donors=_float_or_none(data["hbond_donors"]),
        hbond_acceptors=_float_or_none(data["hbond_acceptors"]),
        nC=_float_or_none(data["nC"]),
        nH=_float_or_none(data["nH"]),
        nO=_float_or_none(data["nO"]),
        nN=_float_or_none(data["nN"]),
        bonds_single=_float_or_none(data["bonds_single"]),
        bonds_double=_float_or_none(data["bonds_double"]),
        bonds_triple=_float_or_none(data["bonds_triple"]),
        bonds_aromatic=_float_or_none(data["bonds_aromatic"]),
    )


def _hamiltonian_to_dict(result: HamiltonianResult) -> dict[str, Any]:
    return {
        "hamiltonian": result.hamiltonian,
        "status": result.status.value,
        "energy_hof": result.energy_hof.to_dict() if result.energy_hof else None,
        "electronic_descriptors": dict(result.electronic_descriptors),
        "conformer_index": result.conformer_index,
        "artifact_ids": list(result.artifact_ids),
        "error_message": result.error_message,
    }


def _hamiltonian_from_dict(data: dict[str, Any]) -> HamiltonianResult:
    return HamiltonianResult(
        hamiltonian=str(data["hamiltonian"]),
        status=OverallStatus(data["status"]),
        energy_hof=(
            Quantity.from_dict(data["energy_hof"])
            if data["energy_hof"] is not None
            else None
        ),
        electronic_descriptors={
            str(key): float(value)
            for key, value in dict(data["electronic_descriptors"]).items()
        },
        conformer_index=_int_or_none(data["conformer_index"]),
        artifact_ids=list(data["artifact_ids"]),
        error_message=data["error_message"],
    )


def _conformer_to_dict(result: ConformerResult) -> dict[str, Any]:
    return {
        "conformer_index": result.conformer_index,
        "crest_energy_hartree": result.crest_energy_hartree,
        "hamiltonian_results": [
            _hamiltonian_to_dict(hamiltonian)
            for hamiltonian in result.hamiltonian_results
        ],
    }


def _conformer_from_dict(data: dict[str, Any]) -> ConformerResult:
    return ConformerResult(
        conformer_index=int(data["conformer_index"]),
        crest_energy_hartree=_float_or_none(data["crest_energy_hartree"]),
        hamiltonian_results=[
            _hamiltonian_from_dict(item) for item in data["hamiltonian_results"]
        ],
    )


def _estimate_to_dict(estimate: PropertyEstimate) -> dict[str, Any]:
    return {
        "estimate_id": estimate.estimate_id,
        "property_id": estimate.property_id,
        "role": estimate.role.value,
        "method_id": estimate.method_id,
        "method_version": estimate.method_version,
        "hamiltonian": estimate.hamiltonian,
        "value": estimate.value.to_dict(),
        "value_kcal_mol": estimate.value_kcal_mol,
        "value_kj_mol": estimate.value_kj_mol,
        "conformer_source_id": estimate.conformer_source_id,
        "uncertainty": estimate.uncertainty,
        "model_path": estimate.model_path,
    }


def _estimate_from_dict(data: dict[str, Any]) -> PropertyEstimate:
    return PropertyEstimate(
        estimate_id=str(data["estimate_id"]),
        property_id=str(data["property_id"]),
        role=PropertyRole(data["role"]),
        method_id=str(data["method_id"]),
        method_version=str(data["method_version"]),
        hamiltonian=data["hamiltonian"],
        value=Quantity.from_dict(data["value"]),
        value_kcal_mol=_float_or_none(data["value_kcal_mol"]),
        value_kj_mol=_float_or_none(data["value_kj_mol"]),
        conformer_source_id=_int_or_none(data["conformer_source_id"]),
        uncertainty=_float_or_none(data["uncertainty"]),
        model_path=data["model_path"],
    )


def _method_ref_to_dict(method_ref: CalculationMethodReference) -> dict[str, Any]:
    return {
        "method_id": method_ref.method_id,
        "method_version": method_ref.method_version,
        "property_id": method_ref.property_id,
    }


def _method_ref_from_dict(data: dict[str, Any]) -> CalculationMethodReference:
    return CalculationMethodReference(
        method_id=str(data["method_id"]),
        method_version=str(data["method_version"]),
        property_id=str(data["property_id"]),
    )


def _run_to_dict(run: RunMetadata) -> dict[str, Any]:
    return {
        "run_id": run.run_id,
        "execution_phase": run.execution_phase,
        "method_ref": _method_ref_to_dict(run.method_ref),
        "started_at": _datetime_to_str(run.started_at),
        "completed_at": _datetime_to_str(run.completed_at),
        "grimperium_version": run.grimperium_version,
    }


def _run_from_dict(data: dict[str, Any]) -> RunMetadata:
    return RunMetadata(
        run_id=str(data["run_id"]),
        execution_phase=str(data["execution_phase"]),
        method_ref=_method_ref_from_dict(data["method_ref"]),
        started_at=_datetime_from_str(data["started_at"]),
        completed_at=_datetime_from_str(data["completed_at"]),
        grimperium_version=data["grimperium_version"],
    )


def to_dict(result: MoleculeCalculationResult) -> dict[str, Any]:
    """Serialize a canonical result to a JSON-safe dictionary."""
    return {
        "schema_version": result.schema_version,
        "molecule": {
            "smiles": result.molecule.smiles,
            "name": result.molecule.name,
            "charge": result.molecule.charge,
            "multiplicity": result.molecule.multiplicity,
        },
        "run": _run_to_dict(result.run),
        "overall_status": result.overall_status.value,
        "conformers": [_conformer_to_dict(item) for item in result.conformers],
        "molecular_descriptors": _molecular_descriptors_to_dict(
            result.molecular_descriptors
        ),
        "estimates": [_estimate_to_dict(item) for item in result.estimates],
        "artifacts": [_artifact_to_dict(item) for item in result.artifacts],
        "stage_executions": [
            _stage_execution_to_dict(item) for item in result.stage_executions
        ],
        "error_message": result.error_message,
    }


def from_dict(data: dict[str, Any]) -> MoleculeCalculationResult:
    """Deserialize a canonical result from a dictionary."""
    schema_version = int(data["schema_version"])
    if schema_version > SUPPORTED_SCHEMA_VERSION:
        msg = f"Unsupported schema_version {schema_version}; max supported is 1"
        raise ValueError(msg)

    molecule = data["molecule"]
    return MoleculeCalculationResult(
        schema_version=schema_version,
        molecule=MoleculeData(
            smiles=str(molecule["smiles"]),
            name=molecule["name"],
            charge=int(molecule["charge"]),
            multiplicity=int(molecule["multiplicity"]),
        ),
        run=_run_from_dict(data["run"]),
        overall_status=OverallStatus(data["overall_status"]),
        conformers=[_conformer_from_dict(item) for item in data["conformers"]],
        molecular_descriptors=_molecular_descriptors_from_dict(
            data["molecular_descriptors"]
        ),
        estimates=[_estimate_from_dict(item) for item in data["estimates"]],
        artifacts=[_artifact_from_dict(item) for item in data["artifacts"]],
        stage_executions=[
            _stage_execution_from_dict(item) for item in data["stage_executions"]
        ],
        error_message=data["error_message"],
    )
