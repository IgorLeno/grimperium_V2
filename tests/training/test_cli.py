from __future__ import annotations

from pathlib import Path

import pytest

from semi_imperium.training.__main__ import main
from semi_imperium.training.molecule_queue import TrainingQueue
from tests.training.conftest import CbsFactory, row


def _args(tmp_path: Path, cbs: Path, *rest: str) -> list[str]:
    return ["--csv", str(cbs), "--store", str(tmp_path / "store"), *rest]


def test_status_counts_against_the_queue(
    tmp_path: Path, make_cbs: CbsFactory, capsys: pytest.CaptureFixture[str]
) -> None:
    cbs = make_cbs([row("a"), row("b", smiles="smiles-a", conformer_rank_by_h298=2)])
    assert main([*_args(tmp_path, cbs), "status"]) == 0
    out = capsys.readouterr().out
    assert "na fila (pending): 1 de 1" in out  # one per SMILES
    assert "ocioso" in out


def test_status_after_a_result_counts_ok(
    tmp_path: Path, make_cbs: CbsFactory, capsys: pytest.CaptureFixture[str]
) -> None:
    cbs = make_cbs([row("a"), row("b")])
    queue = TrainingQueue(tmp_path / "store" / "training.sqlite3", cbs)
    queue.ensure_seeded()
    queue.record_success(
        "a", hof_kcal_mol=-1.0, conformers_evaluated=1, elapsed_s=1.0, batch_id="x"
    )
    assert main([*_args(tmp_path, cbs), "status"]) == 0
    assert "concluídas (ok): 1 de 2 (50.0%)" in capsys.readouterr().out


def test_status_with_a_missing_cbs_exits_2(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main([*_args(tmp_path, tmp_path / "absent.csv"), "status"]) == 2
    assert "CBS CSV not found" in capsys.readouterr().err


def test_run_with_a_missing_cbs_exits_2(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main([*_args(tmp_path, tmp_path / "absent.csv"), "run"]) == 2
    assert "CBS CSV not found" in capsys.readouterr().err


def test_there_is_no_migrate_command(tmp_path: Path, make_cbs: CbsFactory) -> None:
    with pytest.raises(SystemExit):
        main([*_args(tmp_path, make_cbs([row("a")])), "migrate"])


def test_stop_without_a_job_exits_1(tmp_path: Path, make_cbs: CbsFactory) -> None:
    assert main([*_args(tmp_path, make_cbs([row("a")])), "stop"]) == 1
