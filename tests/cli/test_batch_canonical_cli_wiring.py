from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from rich.console import Console

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
