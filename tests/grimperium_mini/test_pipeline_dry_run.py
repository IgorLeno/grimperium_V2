from pathlib import Path

from grimperium_mini.config import MiniConfig
from grimperium_mini.io import load_pipeline_tasks
from grimperium_mini.models import ConformerResult
from grimperium_mini.pipeline import (
    process_task,
    run_pipeline,
    select_lowest_hof,
)

DATA = Path("data/grimperium_mini_pipeline_tcc.xlsx")


def test_select_lowest_hof():
    conformers = [
        ConformerResult(
            "m1", "mol", "PM7", 1, Path("c1.xyz"), "success", -1.0, -4.184, None
        ),
        ConformerResult(
            "m1", "mol", "PM7", 2, Path("c2.xyz"), "success", -3.0, -12.552, None
        ),
    ]
    assert select_lowest_hof(conformers).conformer_rank == 2


def test_dry_run_does_not_call_external_subprocess(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("subprocess should not be called in dry-run")

    monkeypatch.setattr("grimperium_mini.crest.subprocess.run", fail)
    monkeypatch.setattr("grimperium_mini.mopac.subprocess.run", fail)
    task = load_pipeline_tasks(DATA, methods=["PM7"], limit=1)[0]
    row, details = process_task(
        task,
        MiniConfig(work_root=tmp_path / "runs", results_dir=tmp_path / "results"),
        dry_run=True,
    )
    assert row["status"] == "dry_run"
    assert details == []


def test_run_pipeline_dry_run_writes_expected_csvs(tmp_path):
    config = MiniConfig(work_root=tmp_path / "runs", results_dir=tmp_path / "results")
    rows, _ = run_pipeline(DATA, ["PM7"], config, limit=2, dry_run=True)
    assert len(rows) == 2
    assert (tmp_path / "results" / "formacao_grimperium_mini.csv").exists()
    assert (tmp_path / "results" / "reacoes_grimperium_mini.csv").exists()
