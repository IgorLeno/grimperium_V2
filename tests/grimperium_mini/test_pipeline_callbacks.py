"""Tests for on_progress callback in run_pipeline."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from grimperium_mini.config import MiniConfig
from grimperium_mini.pipeline import run_pipeline

DATA = Path("data/grimperium_mini_pipeline_tcc.xlsx")

_DRY_ROW: dict[str, object] = {
    "mol_id": "m0",
    "molecule": "methane",
    "smiles": "C",
    "n_conformers_founded": 0,
    "hf_am1_kJmol": "",
    "hf_pm3_kJmol": "",
    "hf_pm7_kJmol": "",
    "status": "dry_run",
    "crest_status": "dry_run",
    "mopac_status": "dry_run",
    "hf_exp_kJmol": "",
    "total_time_s": 0.0,
    "total_time_min": 0.0,
}


@pytest.mark.skipif(not DATA.exists(), reason="xlsx fixture not available")
def test_callback_called_per_task(tmp_path: Path) -> None:
    config = MiniConfig(work_root=tmp_path / "runs", results_dir=tmp_path / "results")
    mock_cb = MagicMock()
    run_pipeline(DATA, config=config, limit=2, dry_run=True, on_progress=mock_cb)
    done_events = [c for c in mock_cb.call_args_list if c.args[0] == "done"]
    assert len(done_events) >= 1


@pytest.mark.skipif(not DATA.exists(), reason="xlsx fixture not available")
def test_callback_none_no_error(tmp_path: Path) -> None:
    config = MiniConfig(work_root=tmp_path / "runs", results_dir=tmp_path / "results")
    run_pipeline(DATA, config=config, limit=1, dry_run=True, on_progress=None)


def test_run_pipeline_accepts_callback_parameter() -> None:
    """Ensure the signature includes on_progress without runtime error."""
    import inspect

    from grimperium_mini.pipeline import run_pipeline as rp

    sig = inspect.signature(rp)
    assert "on_progress" in sig.parameters


def test_callback_invoked_with_mock_tasks() -> None:
    mols: list[dict[str, object]] = [
        {
            "mol_id": f"m{i}",
            "molecule": "methane",
            "smiles": "C",
            "formula": "CH4",
            "group": "alkane",
            "charge": 0,
            "multiplicity": 1,
            "run_crest": True,
            "hf_exp_kJmol": None,
        }
        for i in range(3)
    ]
    config = MiniConfig()
    mock_cb = MagicMock()

    with (
        patch("grimperium_mini.pipeline.load_molecules", return_value=mols),
        patch("grimperium_mini.pipeline.process_molecule", return_value=_DRY_ROW),
        patch("grimperium_mini.pipeline.append_summary_row"),
    ):
        run_pipeline(
            Path("dummy.xlsx"),
            config=config,
            dry_run=True,
            on_progress=mock_cb,
        )

    done_calls = [c for c in mock_cb.call_args_list if c.args[0] == "done"]
    assert len(done_calls) == 3
