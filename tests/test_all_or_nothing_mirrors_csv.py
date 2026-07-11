from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import BatchFailurePolicy, MoleculeStatus
from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.crest_pm7.config import CRESTStatus, PM7Config, QualityGrade
from grimperium.crest_pm7.csv_enhancements import CSVManagerExtensions


class _FakeProcessor:
    def update_timeouts(self, **_kwargs: object) -> None:
        pass

    def process_with_fixed_timeout(
        self, mol_id: str, smiles: str, progress_callback: object = None
    ) -> MagicMock:
        result = MagicMock()
        result.mol_id = mol_id
        result.smiles = smiles
        result.crest_status = CRESTStatus.SUCCESS
        result.crest_conformers_generated = 1
        result.crest_time = 1.0
        result.conformers = []
        result.num_conformers_selected = 0
        result.k_selected_pm7 = None
        result.total_execution_time = 1.0
        result.rdkit_descriptors = {}
        result.get_selected_conformer.return_value = None
        if mol_id == "mol_ok":
            result.success = True
            result.most_stable_hof = -55.0
            result.quality_grade = QualityGrade.A
            result.error_message = None
        else:
            result.success = False
            result.most_stable_hof = None
            result.quality_grade = QualityGrade.FAILED
            result.error_message = "MOPAC failed"
        return result


def test_all_or_nothing_reset_mirrors_pending_to_scientific_csv(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        [
            {
                "mol_id": "mol_ok",
                "smiles": "CCO",
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
                "H298_pm7": "",
            },
            {
                "mol_id": "mol_fail",
                "smiles": "CCC",
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
                "H298_pm7": "",
            },
        ]
    ).to_csv(csv_path, index=False)
    csv_manager = BatchCSVManager(csv_path)
    csv_manager.load_csv()
    state_manager = BatchStateManager(tmp_path / "batch_state.csv", PM7Config())
    state_manager.reconcile_molecules(csv_manager.state_seed_rows())
    batch = csv_manager.select_batch(
        batch_id="batch_0001",
        batch_size=2,
        crest_timeout_minutes=30,
        mopac_timeout_minutes=60,
        failure_policy=BatchFailurePolicy.ALL_OR_NOTHING,
    )
    state_manager.mark_selected_from_batch(batch)
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=MagicMock(),
        pm7_config=PM7Config(temp_dir=tmp_path),
        processor_adapter=_FakeProcessor(),
    )

    with patch.object(
        CSVManagerExtensions, "update_molecule_with_mopac_results", return_value=True
    ):
        result = manager.execute_batch(batch)

    scientific = pd.read_csv(csv_path).set_index("mol_id")
    state = pd.read_csv(tmp_path / "batch_state.csv", keep_default_na=False).set_index(
        "mol_id"
    )
    assert result.invalidated is True
    assert scientific.loc["mol_ok", "status"] == MoleculeStatus.PENDING.value
    assert scientific.loc["mol_fail", "status"] == MoleculeStatus.PENDING.value
    assert float(scientific.loc["mol_ok", "H298_pm7"]) == -55.0
    assert state.loc["mol_ok", "status"] == MoleculeStatus.PENDING.value
    assert state.loc["mol_fail", "status"] == MoleculeStatus.PENDING.value
