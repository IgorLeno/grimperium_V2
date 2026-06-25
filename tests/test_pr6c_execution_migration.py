from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from grimperium.crest_pm7.batch.enums import BatchFailurePolicy, MoleculeStatus
from grimperium.crest_pm7.batch.execution_manager import (
    BatchExecutionManager,
    create_execution_manager,
)
from grimperium.crest_pm7.batch.models import Batch, BatchMolecule, BatchResult
from grimperium.crest_pm7.csv_enhancements import CSVManagerExtensions


def _pm7_config(tmp_path: Path) -> MagicMock:
    config = MagicMock()
    config.temp_dir = tmp_path
    return config


def _success_pm7_result() -> MagicMock:
    pm7_result = MagicMock()
    pm7_result.success = True
    pm7_result.most_stable_hof = -55.0
    pm7_result.quality_grade.value = "A"
    pm7_result.error_message = None
    pm7_result.get_selected_conformer.return_value = None
    pm7_result.k_selected_pm7 = 1
    return pm7_result


def _failed_pm7_result() -> MagicMock:
    pm7_result = MagicMock()
    pm7_result.success = False
    pm7_result.most_stable_hof = None
    pm7_result.quality_grade.value = "F"
    pm7_result.error_message = "MOPAC failed"
    return pm7_result


def _result() -> BatchResult:
    return BatchResult(batch_id="batch_0001", total_count=1)


def _batch(failure_policy: BatchFailurePolicy = BatchFailurePolicy.PARTIAL_OK) -> Batch:
    return Batch(
        batch_id="batch_0001",
        molecules=[
            BatchMolecule(
                mol_id="mol_001",
                smiles="CCO",
                batch_order=1,
                nheavy=3,
            )
        ],
        crest_timeout_minutes=30,
        mopac_timeout_minutes=60,
        failure_policy=failure_policy,
        method_id="pm7_delta_learning",
        method_version="1.0.0",
        method_snapshot={"property_id": "standard_enthalpy_of_formation"},
    )


def test_execution_manager_accepts_required_state_manager(tmp_path: Path) -> None:
    manager = BatchExecutionManager(
        csv_manager=MagicMock(),
        state_manager=MagicMock(),
        detail_manager=MagicMock(),
        pm7_config=_pm7_config(tmp_path),
        processor_adapter=MagicMock(),
    )

    assert manager.state_manager is not None


def test_process_molecule_marks_running_in_state_manager_not_csv(
    tmp_path: Path,
) -> None:
    csv_manager = MagicMock()
    state_manager = MagicMock()
    processor = MagicMock()
    processor.process_with_fixed_timeout.side_effect = RuntimeError("boom")
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=MagicMock(),
        pm7_config=_pm7_config(tmp_path),
        processor_adapter=processor,
    )

    manager._process_molecule(
        mol_id="mol_001",
        smiles="CCO",
        batch_id="batch_0001",
        batch_order=1,
        crest_timeout=30.0,
        mopac_timeout=60.0,
        result=_result(),
        hof_values=[],
        method_id="pm7_delta_learning",
        method_version="1.0.0",
        method_snapshot={},
    )

    running_call = state_manager.update_molecule_status.call_args_list[0]
    assert running_call.args == ("mol_001", MoleculeStatus.RUNNING.value)
    assert (
        running_call.kwargs["extra_fields"]
        | {
            "method_id": "pm7_delta_learning",
            "method_version": "1.0.0",
            "method_definition_snapshot": "{}",
        }
        == running_call.kwargs["extra_fields"]
    )
    csv_manager.mark_running.assert_not_called()


def test_process_molecule_success_updates_state_and_scientific_csv(
    tmp_path: Path,
) -> None:
    csv_manager = MagicMock()
    csv_manager.pm7result_to_csv_update.return_value = {"H298_pm7": -55.0}
    csv_manager.get_reference_hof.return_value = -50.0
    state_manager = MagicMock()
    processor = MagicMock()
    processor.process_with_fixed_timeout.return_value = _success_pm7_result()
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=MagicMock(),
        pm7_config=_pm7_config(tmp_path),
        processor_adapter=processor,
    )
    manager._batch_settings = {}
    manager._batch_logger = MagicMock()

    with patch.object(
        CSVManagerExtensions,
        "update_molecule_with_mopac_results",
        return_value=True,
    ):
        manager._process_molecule(
            mol_id="mol_001",
            smiles="CCO",
            batch_id="batch_0001",
            batch_order=1,
            crest_timeout=30.0,
            mopac_timeout=60.0,
            result=_result(),
            hof_values=[],
            method_id="pm7_delta_learning",
            method_version="1.0.0",
            method_snapshot={},
        )

    state_manager.update_molecule_status.assert_any_call(
        "mol_001",
        MoleculeStatus.OK.value,
    )
    csv_manager.pm7result_to_csv_update.assert_called_once()
    csv_manager.mark_success.assert_called_once_with("mol_001", {"H298_pm7": -55.0})


def test_process_molecule_failure_updates_state_without_scientific_csv(
    tmp_path: Path,
) -> None:
    csv_manager = MagicMock()
    state_manager = MagicMock()
    state_manager.record_rerun_or_skip.return_value = MoleculeStatus.RERUN.value
    processor = MagicMock()
    processor.process_with_fixed_timeout.return_value = _failed_pm7_result()
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=MagicMock(),
        pm7_config=_pm7_config(tmp_path),
        processor_adapter=processor,
    )

    result = _result()
    manager._process_molecule(
        mol_id="mol_001",
        smiles="CCO",
        batch_id="batch_0001",
        batch_order=1,
        crest_timeout=30.0,
        mopac_timeout=60.0,
        result=result,
        hof_values=[],
        method_id="pm7_delta_learning",
        method_version="1.0.0",
        method_snapshot={},
    )

    state_manager.record_rerun_or_skip.assert_called_once_with(
        "mol_001", "MOPAC failed"
    )
    csv_manager.pm7result_to_csv_update.assert_not_called()
    csv_manager.mark_success.assert_not_called()
    assert result.rerun_count == 1


def test_execute_batch_all_or_nothing_resets_state_manager_not_csv(
    tmp_path: Path,
) -> None:
    csv_manager = MagicMock()
    state_manager = MagicMock()
    state_manager.record_rerun_or_skip.return_value = MoleculeStatus.RERUN.value
    processor = MagicMock()
    processor.process_with_fixed_timeout.return_value = _failed_pm7_result()
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=MagicMock(),
        pm7_config=_pm7_config(tmp_path),
        processor_adapter=processor,
    )

    manager.execute_batch(_batch(BatchFailurePolicy.ALL_OR_NOTHING))

    state_manager.reset_all_or_nothing.assert_called_once_with("batch_0001")
    csv_manager.reset_batch.assert_not_called()


def test_execute_batch_records_batch_method_metadata_when_molecule_starts(
    tmp_path: Path,
) -> None:
    csv_manager = MagicMock()
    state_manager = MagicMock()
    state_manager.record_rerun_or_skip.return_value = MoleculeStatus.RERUN.value
    processor = MagicMock()
    processor.process_with_fixed_timeout.return_value = _failed_pm7_result()
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=MagicMock(),
        pm7_config=_pm7_config(tmp_path),
        processor_adapter=processor,
    )

    manager.execute_batch(_batch())

    running_call = state_manager.update_molecule_status.call_args_list[0]
    assert running_call.args == ("mol_001", MoleculeStatus.RUNNING.value)
    assert (
        running_call.kwargs["extra_fields"]
        | {
            "method_id": "pm7_delta_learning",
            "method_version": "1.0.0",
            "method_definition_snapshot": (
                '{"property_id": "standard_enthalpy_of_formation"}'
            ),
        }
        == running_call.kwargs["extra_fields"]
    )


def test_create_execution_manager_wires_state_manager_from_state_csv_path(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    state_csv_path = tmp_path / "batch_state.csv"
    detail_dir = tmp_path / "details"

    manager = create_execution_manager(
        csv_path=csv_path,
        state_csv_path=state_csv_path,
        detail_dir=detail_dir,
        pm7_config=_pm7_config(tmp_path),
    )

    assert manager.csv_manager.csv_path == csv_path
    assert manager.state_manager.state_csv_path == state_csv_path
