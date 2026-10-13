from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from semi_imperium.training.__main__ import main
from tests.training.conftest import CsvFactory, row


def _args(tmp_path: Path, csv: Path, *rest: str) -> list[str]:
    return ["--csv", str(csv), "--store", str(tmp_path / "store"), *rest]


def test_status_prints_progress(
    tmp_path: Path, make_csv: CsvFactory, capsys: pytest.CaptureFixture[str]
) -> None:
    csv = make_csv([row("a", "OK", protocol="1"), row("b", "Pending", protocol="")])
    assert main([*_args(tmp_path, csv), "status"]) == 0
    out = capsys.readouterr().out
    assert "protocolo 1" in out and "1 de 2" in out
    assert "parado" in out or "ocioso" in out


def test_migrate_is_a_dry_run_unless_applied(
    tmp_path: Path, make_csv: CsvFactory, capsys: pytest.CaptureFixture[str]
) -> None:
    csv = make_csv([row("a", "OK"), row("b", "Selected")])
    assert main([*_args(tmp_path, csv), "migrate"]) == 0
    assert "ensaio" in capsys.readouterr().out
    assert "protocol" not in pd.read_csv(csv, nrows=0).columns or (
        pd.read_csv(csv, dtype=str, keep_default_na=False)["protocol"].eq("").all()
    )
    assert main([*_args(tmp_path, csv), "migrate", "--apply"]) == 0
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["a", "protocol"] == "1" and saved.loc["b", "status"] == "Pending"


def test_run_refuses_an_unmigrated_csv_with_exit_2(
    tmp_path: Path, make_csv: CsvFactory, capsys: pytest.CaptureFixture[str]
) -> None:
    csv = make_csv([row("a", "OK"), row("b", "Pending")])
    assert main([*_args(tmp_path, csv), "run", "--limit", "1"]) == 2
    assert "migrate --apply" in capsys.readouterr().err


def test_stop_without_a_job_exits_1(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = make_csv([row("a", "Pending", protocol="")])
    assert main([*_args(tmp_path, csv), "stop"]) == 1
