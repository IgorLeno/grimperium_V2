from __future__ import annotations

from pathlib import Path

import pytest

import semi_imperium.training.__main__ as cli
from semi_imperium.training.__main__ import main
from semi_imperium.training.molecule_queue import TrainingQueue
from semi_imperium.training.runner import RunSettings, RunSummary
from tests.training.conftest import CbsFactory, row, set_ok


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
    set_ok(tmp_path / "store" / "training.sqlite3", "a", hof=-1.0)
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


def _bin_dir(folder: Path) -> Path:
    folder.mkdir(parents=True)
    for name in ("crest", "xtb", "mopac"):
        (folder / name).write_text("")
    return folder


def _captured_settings(monkeypatch: pytest.MonkeyPatch) -> list[RunSettings]:
    seen: list[RunSettings] = []

    def fake_run(settings: RunSettings) -> RunSummary:
        seen.append(settings)
        return RunSummary(attempted=0, ok=0, failed=0, stopped=False)

    monkeypatch.setattr(cli, "run_training", fake_run)
    return seen


def test_run_uses_the_programs_of_the_given_bin_dir(
    tmp_path: Path, make_cbs: CbsFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen = _captured_settings(monkeypatch)
    bins = _bin_dir(tmp_path / "env" / "bin")
    args = [*_args(tmp_path, make_cbs([row("a")])), "run", "--bin-dir", str(bins)]
    assert main(args) == 0
    (settings,) = seen
    assert (
        settings.crest_executable,
        settings.xtb_executable,
        settings.mopac_executable,
    ) == (str(bins / "crest"), str(bins / "xtb"), str(bins / "mopac"))


def test_run_with_a_bin_dir_missing_a_program_exits_2(
    tmp_path: Path, make_cbs: CbsFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen = _captured_settings(monkeypatch)
    bins = _bin_dir(tmp_path / "env" / "bin")
    (bins / "mopac").unlink()
    args = [*_args(tmp_path, make_cbs([row("a")])), "run", "--bin-dir", str(bins)]
    assert main(args) == 2
    assert seen == []


def test_run_defaults_to_the_fixed_environment_when_it_exists(
    tmp_path: Path, make_cbs: CbsFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen = _captured_settings(monkeypatch)
    root = tmp_path / "mamba"
    bins = _bin_dir(root / "envs" / "grimperium-quimica" / "bin")
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(root))
    assert main([*_args(tmp_path, make_cbs([row("a")])), "run"]) == 0
    assert seen[0].crest_executable == str(bins / "crest")


def test_run_without_the_fixed_environment_falls_back_to_path_and_warns(
    tmp_path: Path,
    make_cbs: CbsFactory,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    seen = _captured_settings(monkeypatch)
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "no-mamba"))
    assert main([*_args(tmp_path, make_cbs([row("a")])), "run"]) == 0
    assert seen[0].crest_executable == "crest"
    assert "grimperium-quimica" in capsys.readouterr().err
