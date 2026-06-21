from pathlib import Path

import pytest

from grimperium_mini.config import MiniConfig
from grimperium_mini.io import load_molecules
from grimperium_mini.models import ConformerResult
from grimperium_mini.pipeline import (
    process_molecule,
    run_pipeline,
    select_lowest_hof,
)

DATA = Path("data/grimperium_mini_pipeline_tcc.xlsx")


def test_select_lowest_hof() -> None:
    conformers = [
        ConformerResult(
            "m1", "mol", "PM7", 1, Path("c1.xyz"), "success", -1.0, -4.184, None
        ),
        ConformerResult(
            "m1", "mol", "PM7", 2, Path("c2.xyz"), "success", -3.0, -12.552, None
        ),
    ]
    result = select_lowest_hof(conformers)
    assert result is not None
    assert result.conformer_rank == 2


@pytest.mark.skipif(not DATA.exists(), reason="xlsx fixture not available")
def test_dry_run_does_not_call_external_subprocess(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(*args: object, **kwargs: object) -> None:
        raise AssertionError("subprocess should not be called in dry-run")

    monkeypatch.setattr("grimperium_mini.crest.subprocess.run", fail)
    monkeypatch.setattr("grimperium_mini.mopac.subprocess.run", fail)
    mols = load_molecules(DATA, limit=1)
    row = process_molecule(
        mols[0],
        MiniConfig(work_root=tmp_path / "runs", results_dir=tmp_path / "results"),
        dry_run=True,
    )
    assert row["status"] == "dry_run"


@pytest.mark.skipif(not DATA.exists(), reason="xlsx fixture not available")
def test_run_pipeline_dry_run_writes_summary_xlsx(tmp_path: Path) -> None:
    config = MiniConfig(work_root=tmp_path / "runs", results_dir=tmp_path / "results")
    rows = run_pipeline(DATA, config, limit=2, dry_run=True)
    assert len(rows) == 2
    assert (tmp_path / "results" / "grimperium_mini_summary.xlsx").exists()
