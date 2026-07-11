from pathlib import Path
from unittest.mock import patch

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager
from grimperium.crest_pm7.batch.output_contracts import BatchOutputLayout
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import CRESTStatus, MOPACStatus, PM7Config
from grimperium.crest_pm7.csv_enhancements import CSVManagerExtensions
from grimperium.crest_pm7.molecule_processor import ConformerData, PM7Result


class _FakeProcessor:
    def __init__(self) -> None:
        self.calls: list[dict[str, object]] = []

    def update_timeouts(self, **_kwargs: object) -> None:
        pass

    def process_with_fixed_timeout(
        self,
        mol_id: str,
        smiles: str,
        input_xyz: Path | None = None,
        progress_callback: object = None,
        charge: int = 0,
        multiplicity: int = 1,
    ) -> PM7Result:
        self.calls.append(
            {
                "mol_id": mol_id,
                "charge": charge,
                "multiplicity": multiplicity,
            }
        )
        if mol_id == "mol_fail":
            return PM7Result(
                mol_id=mol_id,
                smiles=smiles,
                conformers=[],
                success=False,
                error_message="MOPAC timeout",
            )
        conformer = ConformerData(index=0, mol_id=mol_id, crest_rank=1)
        conformer.crest_status = CRESTStatus.SUCCESS
        conformer.mopac_status = MOPACStatus.SUCCESS
        conformer.energy_hof = -55.0
        conformer.hof_extraction_successful = True
        return PM7Result(
            mol_id=mol_id,
            smiles=smiles,
            phase="A",
            nheavy=3,
            rdkit_descriptors={},
            crest_status=CRESTStatus.SUCCESS,
            crest_conformers_generated=1,
            conformers=[conformer],
            num_conformers_selected=1,
            k_selected_pm7=1,
            total_execution_time=1.0,
            success=True,
        )


def test_batch_e2e_fake_pipeline_dual_write_canonical_and_recovery(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        [
            {
                "mol_id": "mol_ok",
                "smiles": "CCO",
                "nheavy": 3,
                "charge": -1,
                "multiplicity": 2,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            },
            {
                "mol_id": "mol_fail",
                "smiles": "CCC",
                "nheavy": 3,
                "charge": 0,
                "multiplicity": 1,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            },
            {
                "mol_id": "mol_stuck",
                "smiles": "C",
                "nheavy": 1,
                "charge": 0,
                "multiplicity": 1,
                "status": MoleculeStatus.RUNNING.value,
                "reruns": 0,
            },
        ]
    ).to_csv(csv_path, index=False)
    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()
    state_manager = BatchStateManager(tmp_path / "batch_state.csv", PM7Config())
    state_manager.reconcile_molecules(csv_manager.state_seed_rows())
    recovered = state_manager.reset_stuck_running_molecules()
    for mol_id in recovered:
        csv_manager.apply_operational_status(mol_id, MoleculeStatus.PENDING.value)

    batch = csv_manager.select_batch(
        batch_id="batch_0001",
        batch_size=2,
        crest_timeout_minutes=30,
        mopac_timeout_minutes=60,
    )
    state_manager.mark_selected_from_batch(batch)
    fake_processor = _FakeProcessor()
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=type(
            "DetailStub",
            (),
            {
                "pm7result_to_detail": lambda self, **kwargs: {},
                "save_detail": lambda self, detail: None,
            },
        )(),
        pm7_config=PM7Config(temp_dir=tmp_path),
        processor_adapter=fake_processor,
        output_layout=BatchOutputLayout(tmp_path / "out"),
        write_xlsx=False,
    )

    with patch.object(
        CSVManagerExtensions, "update_molecule_with_mopac_results", return_value=True
    ):
        result = manager.execute_batch(batch)

    state = pd.read_csv(tmp_path / "batch_state.csv", keep_default_na=False).set_index(
        "mol_id"
    )
    scientific = pd.read_csv(csv_path).set_index("mol_id")
    canonical = pd.read_csv(tmp_path / "out" / "calculation_results.csv")
    assert state.loc["mol_ok", "status"] == MoleculeStatus.OK.value
    assert scientific.loc["mol_ok", "status"] == MoleculeStatus.OK.value
    assert state.loc["mol_fail", "status"] == MoleculeStatus.RERUN.value
    assert scientific.loc["mol_fail", "status"] == MoleculeStatus.RERUN.value
    assert result.success_count == 1
    assert result.rerun_count == 1
    assert set(canonical["molecule_name"]) == {"mol_ok"}
    mol_ok_call = next(
        call for call in fake_processor.calls if call["mol_id"] == "mol_ok"
    )
    assert mol_ok_call["charge"] == -1
    assert mol_ok_call["multiplicity"] == 2
