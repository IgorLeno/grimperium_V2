from datetime import datetime, timezone

import pytest

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
    MolecularDescriptors,
    MoleculeCalculationResult,
    MoleculeData,
    PropertyEstimate,
    RunMetadata,
    StageExecutionRecord,
)
from grimperium.calculation.contracts.quantity import Quantity
from grimperium.calculation.contracts.serialization import from_dict, to_dict


def _sample_result() -> MoleculeCalculationResult:
    started_at = datetime(2026, 6, 22, 10, 30, 1, 123456, tzinfo=timezone.utc)
    completed_at = datetime(2026, 6, 22, 10, 31, 2, 654321, tzinfo=timezone.utc)
    method_ref = CalculationMethodReference(
        method_id="pm7_delta_legacy",
        method_version="0.0.0",
        property_id="standard_enthalpy_of_formation",
    )
    return MoleculeCalculationResult(
        molecule=MoleculeData(smiles="CCO", name="ethanol"),
        run=RunMetadata(
            run_id="run-001",
            execution_phase="B",
            method_ref=method_ref,
            started_at=started_at,
            completed_at=completed_at,
            grimperium_version="0.2.0",
        ),
        overall_status=OverallStatus.SUCCESS,
        conformers=[
            ConformerResult(
                conformer_index=2,
                crest_energy_hartree=-154.12345678901234,
                hamiltonian_results=[
                    HamiltonianResult(
                        hamiltonian="PM7",
                        status=OverallStatus.SUCCESS,
                        energy_hof=Quantity(value=-56.7890123456789, unit="kcal/mol"),
                        electronic_descriptors={"mopac_homo_ev": -10.123456789},
                        conformer_index=2,
                        artifact_ids=["artifact-001"],
                        error_message=None,
                    )
                ],
            )
        ],
        molecular_descriptors=MolecularDescriptors(
            nrotbonds=1.0,
            tpsa=20.23,
            mol_weight=46.069,
        ),
        estimates=[
            PropertyEstimate(
                estimate_id="estimate-001",
                property_id="standard_enthalpy_of_formation",
                role=PropertyRole.BASELINE,
                method_id="pm7_delta_legacy",
                method_version="0.0.0",
                hamiltonian="PM7",
                value=Quantity(value=-56.7890123456789, unit="kcal/mol"),
                value_kcal_mol=-56.7890123456789,
                value_kj_mol=None,
                conformer_source_id=2,
                uncertainty=None,
                model_path=None,
            )
        ],
        artifacts=[
            CalculationArtifact(
                artifact_id="artifact-001",
                artifact_type=ArtifactType.MOPAC_OUT,
                relative_path="molecule/conf_2/pm7.out",
                hamiltonian="PM7",
                conformer_index=2,
            )
        ],
        stage_executions=[
            StageExecutionRecord(
                stage_id="mopac_pm7",
                program="mopac",
                role="baseline",
                status=StageExecutionStatus.SUCCESS,
                requested=True,
                started_at=started_at,
                completed_at=completed_at,
                execution_time_s=61.530865,
                program_version="MOPAC2016",
                settings={"hamiltonian": "PM7"},
                artifact_ids=["artifact-001"],
                error_message=None,
            )
        ],
    )


def test_quantity_converts_kj_mol_to_kcal_mol() -> None:
    assert Quantity(value=4.184, unit="kJ/mol").to_kcal_mol() == pytest.approx(1.0)
    assert Quantity(value=1.25, unit="kcal/mol").to_kcal_mol() == pytest.approx(1.25)


def test_quantity_rejects_unknown_units() -> None:
    with pytest.raises(ValueError, match="Unsupported energy unit"):
        Quantity(value=1.0, unit="hartree").to_kcal_mol()


def test_result_serialization_uses_strings_and_iso_datetimes() -> None:
    payload = to_dict(_sample_result())

    assert payload["schema_version"] == 1
    assert payload["overall_status"] == "success"
    assert payload["estimates"][0]["role"] == "baseline"
    assert payload["artifacts"][0]["artifact_type"] == "mopac_out"
    assert payload["stage_executions"][0]["status"] == "success"
    assert payload["run"]["started_at"] == "2026-06-22T10:30:01.123456+00:00"


def test_result_round_trip_preserves_numeric_precision() -> None:
    original = _sample_result()
    restored = from_dict(to_dict(original))

    assert restored.conformers[0].crest_energy_hartree == -154.12345678901234
    assert restored.estimates[0].value.value == -56.7890123456789
    assert restored.stage_executions[0].execution_time_s == 61.530865
    assert to_dict(restored) == to_dict(original)


def test_from_dict_rejects_future_schema_version() -> None:
    payload = to_dict(_sample_result())
    payload["schema_version"] = 2

    with pytest.raises(ValueError, match="Unsupported schema_version"):
        from_dict(payload)
