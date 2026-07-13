from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pandas as pd
from rich.console import Console

from grimperium.calculation.methods import get_calculation_method
from grimperium.cli.views.batch_view import BatchView
from grimperium.crest_pm7.batch.enums import MoleculeStatus


def test_batch_view_wires_unique_canonical_output_layout(tmp_path: Path) -> None:
    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        [
            {
                "mol_id": "mol_001",
                "smiles": "CCO",
                "nheavy": 3,
                "status": MoleculeStatus.PENDING.value,
                "reruns": 0,
            }
        ]
    ).to_csv(csv_path, index=False)

    controller = SimpleNamespace(
        console=Console(record=True),
        settings_manager=None,
    )
    view = BatchView(controller)  # type: ignore[arg-type]
    view.csv_path = csv_path
    view.detail_dir = tmp_path / "details"

    exec_manager, batch = view._prepare_batch()

    assert batch.batch_id
    assert exec_manager._output_layout is not None
    assert exec_manager._output_layout.output_dir.parent == tmp_path / "batch_outputs"
    assert batch.batch_id in exec_manager._output_layout.output_dir.name
    assert exec_manager.state_manager.state_csv_path == tmp_path / "batch_state.csv"


def test_batch_view_method_fields_always_crest_pm7_with_delta_session() -> None:
    delta = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )
    controller = MagicMock()
    controller.current_method_definition = delta
    controller.console = Console(record=True)
    view = BatchView(controller)

    method_id, method_version, snapshot = view._method_run_fields()

    assert method_id == "crest_pm7"
    assert method_version == "1.0.0"
    assert snapshot == {"method_id": "crest_pm7"}
