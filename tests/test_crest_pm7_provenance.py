from pathlib import Path

from grimperium.calculation.contracts.adapters import (
    canonical_pm7_method_id,
    pm7result_to_canonical,
)
from grimperium.calculation.methods.registry import get_calculation_method
from grimperium.crest_pm7.batch.models import Batch
from grimperium.crest_pm7.config import CRESTStatus, MOPACStatus
from grimperium.crest_pm7.molecule_processor import ConformerData, PM7Result


def _pm7_result() -> PM7Result:
    conformer = ConformerData(index=0, mol_id="mol_001", crest_rank=1)
    conformer.crest_status = CRESTStatus.SUCCESS
    conformer.crest_geometry_file = Path("work/mol_001/conf_1.xyz")
    conformer.mopac_status = MOPACStatus.SUCCESS
    conformer.mopac_output_file = Path("work/mol_001/conf_1.out")
    conformer.energy_hof = -55.0
    conformer.hof_extraction_successful = True
    return PM7Result(
        mol_id="mol_001",
        smiles="CCO",
        phase="A",
        conformers=[conformer],
        k_selected_pm7=1,
        success=True,
    )


def test_crest_pm7_method_definition_is_registered() -> None:
    method = get_calculation_method(
        "crest_pm7", property_id="standard_enthalpy_of_formation"
    )

    assert method.method_id == "crest_pm7"
    assert method.version == "1.0"
    assert method.model_requirement.model_required is False
    assert method.compatibility.baseline_hamiltonian == "PM7"


def test_batch_defaults_to_crest_pm7_method() -> None:
    batch = Batch(
        batch_id="batch_0001",
        crest_timeout_minutes=30,
        mopac_timeout_minutes=60,
    )

    assert batch.method_id == "crest_pm7"


def test_pm7_delta_method_ids_map_to_crest_pm7_in_canonical_reads() -> None:
    canonical = pm7result_to_canonical(
        _pm7_result(), method_id=canonical_pm7_method_id("pm7_delta_learning")
    )

    assert canonical_pm7_method_id("pm7_delta_learning") == "crest_pm7"
    assert canonical_pm7_method_id("pm7_delta_legacy") == "crest_pm7"
    assert canonical.run.method_ref.method_id == "crest_pm7"
    assert canonical.estimates[0].method_id == "crest_pm7"
