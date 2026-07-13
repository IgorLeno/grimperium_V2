from pathlib import Path

from grimperium.calculation.contracts.adapters import pm7result_to_canonical
from grimperium.calculation.contracts.enums import OverallStatus, PropertyRole
from grimperium.crest_pm7.config import (
    CRESTStatus,
    HOFConfidence,
    MOPACStatus,
)
from grimperium.crest_pm7.molecule_processor import ConformerData, PM7Result


def _pm7_result() -> PM7Result:
    failed = ConformerData(index=0, mol_id="mol-001", crest_rank=1)
    failed.crest_status = CRESTStatus.SUCCESS
    failed.crest_geometry_file = Path("work/mol-001/conf_1.xyz")
    failed.mopac_status = MOPACStatus.ERROR
    failed.mopac_error_message = "SCF failed"

    selected = ConformerData(index=1, mol_id="mol-001", crest_rank=2)
    selected.crest_status = CRESTStatus.SUCCESS
    selected.crest_geometry_file = Path("work/mol-001/conf_2.xyz")
    selected.mopac_status = MOPACStatus.SUCCESS
    selected.mopac_output_file = Path("work/mol-001/conf_2.out")
    selected.mopac_execution_time = 12.5
    selected.energy_hof = -42.25
    selected.hof_confidence = HOFConfidence.HIGH
    selected.hof_extraction_method = "FINAL_HEAT_OF_FORMATION"
    selected.hof_extraction_successful = True

    higher_hof = ConformerData(index=2, mol_id="mol-001", crest_rank=3)
    higher_hof.crest_status = CRESTStatus.SUCCESS
    higher_hof.crest_geometry_file = Path("work/mol-001/conf_3.xyz")
    higher_hof.mopac_status = MOPACStatus.SUCCESS
    higher_hof.mopac_output_file = Path("work/mol-001/conf_3.out")
    higher_hof.mopac_execution_time = 11.0
    higher_hof.energy_hof = -40.0
    higher_hof.hof_confidence = HOFConfidence.HIGH
    higher_hof.hof_extraction_method = "FINAL_HEAT_OF_FORMATION"
    higher_hof.hof_extraction_successful = True

    return PM7Result(
        mol_id="mol-001",
        smiles="CCO",
        phase="B",
        nheavy=3,
        nrotbonds=1,
        tpsa=20.23,
        aromatic_rings=0,
        has_heteroatoms=True,
        rdkit_descriptors={
            "rdkit_nrotbonds": 1.0,
            "rdkit_tpsa": 20.23,
            "rdkit_num_rings": 0.0,
            "rdkit_fsp3": 1.0,
            "rdkit_mol_weight": 46.069,
            "rdkit_hbond_donors": 1.0,
            "rdkit_hbond_acceptors": 1.0,
            "rdkit_nC": 2.0,
            "rdkit_nH": 6.0,
            "rdkit_nO": 1.0,
            "rdkit_nN": 0.0,
            "rdkit_bonds_single": 8.0,
            "rdkit_bonds_double": 0.0,
            "rdkit_bonds_triple": 0.0,
            "rdkit_bonds_aromatic": 0.0,
            "not_legacy_schema": 999.0,
        },
        crest_status=CRESTStatus.SUCCESS,
        crest_conformers_generated=3,
        crest_time=7.5,
        conformers=[failed, selected, higher_hof],
        num_conformers_selected=3,
        k_selected_pm7=2,
        total_execution_time=31.0,
        success=True,
    )


def test_pm7_adapter_maps_selected_conformer_and_preserves_hof() -> None:
    canonical = pm7result_to_canonical(_pm7_result())

    assert canonical.run.execution_phase == "B"
    assert canonical.run.method_ref.method_id == "crest_pm7"
    assert canonical.run.method_ref.method_version == "0.0.0"
    assert canonical.overall_status is OverallStatus.SUCCESS
    assert len(canonical.estimates) == 1
    estimate = canonical.estimates[0]
    assert estimate.role is PropertyRole.BASELINE
    assert estimate.hamiltonian == "PM7"
    assert estimate.value.value == -42.25
    assert estimate.value.unit == "kcal/mol"
    assert estimate.value_kcal_mol == -42.25
    assert estimate.value_kj_mol is None
    assert estimate.conformer_source_id == 2
    assert estimate.model_path is None


def test_pm7_adapter_maps_only_legacy_rdkit_descriptor_fields() -> None:
    canonical = pm7result_to_canonical(_pm7_result())

    assert canonical.molecular_descriptors is not None
    descriptors = canonical.molecular_descriptors
    assert descriptors.nrotbonds == 1.0
    assert descriptors.tpsa == 20.23
    assert descriptors.mol_weight == 46.069
    assert not hasattr(descriptors, "not_legacy_schema")


def test_pm7_adapter_maps_electronic_descriptors_from_selected_conformer() -> None:
    canonical = pm7result_to_canonical(_pm7_result())
    hamiltonian = canonical.conformers[1].hamiltonian_results[0]

    assert hamiltonian.hamiltonian == "PM7"
    assert hamiltonian.conformer_index == 2
    assert hamiltonian.energy_hof is not None
    assert hamiltonian.energy_hof.value == -42.25
    assert hamiltonian.electronic_descriptors == {}


def test_pm7_adapter_keeps_artifact_paths_relative() -> None:
    canonical = pm7result_to_canonical(_pm7_result())

    assert canonical.artifacts
    for artifact in canonical.artifacts:
        assert not Path(artifact.relative_path).is_absolute()


def test_pm7_adapter_uses_authoritative_run_id_when_provided() -> None:
    canonical = pm7result_to_canonical(
        _pm7_result(),
        run_id="run_authoritative",
    )

    assert canonical.run.run_id == "run_authoritative"
