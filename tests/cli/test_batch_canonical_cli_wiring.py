from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pandas as pd
import pytest
from rich.console import Console

from grimperium.calculation.methods import get_calculation_method
from grimperium.cli.session import SessionContext
from grimperium.cli.views.batch_view import BatchView
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.runs.models import RunStatus
from grimperium.runs.persistence import MANIFEST_FILENAME
from grimperium.runs.service import RunService


def _pending_csv(path: Path) -> None:
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
    ).to_csv(path, index=False)


def _ok_csv(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "mol_id": "mol_001",
                "smiles": "CCO",
                "nheavy": 3,
                "status": MoleculeStatus.OK.value,
                "reruns": 0,
            }
        ]
    ).to_csv(path, index=False)


def _view(tmp_path: Path, csv_path: Path, run_service: RunService) -> BatchView:
    controller = SimpleNamespace(
        console=Console(record=True),
        settings_manager=None,
        session=SessionContext(),
        run_service=run_service,
    )
    view = BatchView(controller)  # type: ignore[arg-type]
    view.csv_path = csv_path
    view.detail_dir = tmp_path / "details"
    return view


def test_batch_view_wires_output_under_runs_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    csv_path = tmp_path / "thermo_pm7.csv"
    _pending_csv(csv_path)
    service = RunService(runs_root)
    view = _view(tmp_path, csv_path, service)

    exec_manager, batch, manifest = view._prepare_batch()

    assert batch.batch_id
    assert exec_manager._output_layout is not None
    assert exec_manager._output_layout.output_dir == runs_root / manifest.run_id
    assert exec_manager.canonical_run_id == manifest.run_id
    assert not (tmp_path / "batch_outputs").exists()
    assert exec_manager.state_manager.state_csv_path == tmp_path / "batch_state.csv"


def test_batch_view_manifest_and_canonical_share_run_id(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    csv_path = tmp_path / "thermo_pm7.csv"
    _pending_csv(csv_path)
    service = RunService(runs_root)
    view = _view(tmp_path, csv_path, service)

    exec_manager, _batch, manifest = view._prepare_batch()

    assert exec_manager.canonical_run_id == manifest.run_id
    assert (runs_root / manifest.run_id / MANIFEST_FILENAME).exists()
    assert exec_manager._output_layout.output_dir.name == manifest.run_id


def test_batch_view_serialized_paths_relative_and_reloadable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    csv_path = tmp_path / "thermo_pm7.csv"
    _pending_csv(csv_path)
    service = RunService(runs_root)
    view = _view(tmp_path, csv_path, service)

    exec_manager, _batch, manifest = view._prepare_batch()
    layout = exec_manager._output_layout
    assert layout is not None
    layout.output_dir.mkdir(parents=True, exist_ok=True)
    layout.calculation_results_csv.write_text(
        "mol_id,H298_pm7\nmol_001,-10\n", encoding="utf-8"
    )
    layout.batch_state_csv.write_text("mol_id,status\nmol_001,OK\n", encoding="utf-8")

    service.start_run(manifest.run_id)
    attached = view._attach_existing_outputs(manifest, layout)
    service.complete_run(manifest.run_id, success_count=1, failure_count=0)

    raw = json.loads(
        (runs_root / manifest.run_id / MANIFEST_FILENAME).read_text(encoding="utf-8")
    )
    serialized = raw["output_paths"]["calculation_results_csv"]
    assert not Path(serialized).is_absolute()
    assert serialized == f"{manifest.run_id}/calculation_results.csv"

    reloaded = RunService(runs_root).get_run(manifest.run_id)
    assert reloaded.output_paths["calculation_results_csv"].exists()
    assert attached.output_paths["calculation_results_csv"].exists()


def test_batch_view_empty_batch_cancels_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    csv_path = tmp_path / "thermo_pm7.csv"
    _ok_csv(csv_path)
    service = RunService(runs_root)
    view = _view(tmp_path, csv_path, service)

    view._run_batch_with_tracker()

    assert view.controller.session.run.run_id is not None
    cancelled = service.get_run(view.controller.session.run.run_id)
    assert cancelled.status is RunStatus.CANCELLED
    assert cancelled.error == "No molecules available for processing"
    assert not (tmp_path / "batch_outputs").exists()


def test_batch_view_no_orphan_scientific_dirs_outside_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runs_root = tmp_path / "runs"
    monkeypatch.setenv("GRIMPERIUM_RUNS_DIR", str(runs_root))
    csv_path = tmp_path / "thermo_pm7.csv"
    _pending_csv(csv_path)
    service = RunService(runs_root)
    view = _view(tmp_path, csv_path, service)

    exec_manager, _batch, manifest = view._prepare_batch()

    assert exec_manager._output_layout is not None
    assert exec_manager._output_layout.output_dir == service.run_dir(manifest.run_id)
    assert not (tmp_path / "batch_outputs").exists()
    assert list(tmp_path.glob("**/calculation_results.csv")) == []


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
