"""Tests for on_progress callback in run_pipeline."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pytest

from grimperium_mini.config import MiniConfig
from grimperium_mini.pipeline import run_pipeline


DATA = Path("data/grimperium_mini_pipeline_tcc.xlsx")


@pytest.mark.skipif(not DATA.exists(), reason="xlsx fixture not available")
def test_callback_called_per_task() -> None:
    config = MiniConfig()
    mock_cb = MagicMock()
    run_pipeline(DATA, methods=["AM1"], config=config, limit=2, dry_run=True, on_progress=mock_cb)
    done_calls = [c for c in mock_cb.call_args_list if c == call("done", c.args[1])]
    # at least one "done" call per task processed
    done_events = [c for c in mock_cb.call_args_list if c.args[0] == "done"]
    assert len(done_events) >= 1


@pytest.mark.skipif(not DATA.exists(), reason="xlsx fixture not available")
def test_callback_none_no_error() -> None:
    config = MiniConfig()
    run_pipeline(DATA, methods=["AM1"], config=config, limit=1, dry_run=True, on_progress=None)


def test_run_pipeline_accepts_callback_parameter() -> None:
    """Ensure the signature includes on_progress without runtime error on mock data."""
    import inspect
    from grimperium_mini.pipeline import run_pipeline as rp

    sig = inspect.signature(rp)
    assert "on_progress" in sig.parameters


def test_callback_invoked_with_mock_tasks() -> None:
    from grimperium_mini.models import PipelineTask

    tasks = [
        PipelineTask(
            mol_id=f"m{i}",
            molecule="methane",
            smiles="C",
            formula="CH4",
            group="alkane",
            method="PM7",
        )
        for i in range(3)
    ]
    config = MiniConfig()
    mock_cb = MagicMock()

    with (
        patch("grimperium_mini.pipeline.load_pipeline_tasks", return_value=tasks),
        patch("grimperium_mini.pipeline.process_task") as mock_process,
        patch("grimperium_mini.pipeline.write_csv"),
        patch("grimperium_mini.pipeline._backup_if_exists"),
        patch("grimperium_mini.pipeline.load_reaction_rows", return_value=[]),
        patch("grimperium_mini.pipeline.calculate_reactions", return_value=[]),
        patch("grimperium_mini.pipeline.write_reactions_csv"),
    ):
        mock_process.return_value = ({"mol_id": "m0", "status": "dry_run", "method": "PM7"}, [])
        run_pipeline(
            Path("dummy.xlsx"),
            methods=["PM7"],
            config=config,
            dry_run=True,
            on_progress=mock_cb,
        )

    done_calls = [c for c in mock_cb.call_args_list if c.args[0] == "done"]
    assert len(done_calls) == 3
