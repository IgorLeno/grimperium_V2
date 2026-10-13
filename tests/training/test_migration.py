from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pytest

from semi_imperium.training.migration import migrate
from tests.training.conftest import CsvFactory, row


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(make_csv: CsvFactory) -> Path:
    return make_csv(
        [
            row("a", "OK", H298_pm7=-3.83),
            row("b", "OK", H298_pm7=-5.5),
            row("c", "Selected"),
            row("d", "Pending"),
            row("e", "Skip"),
        ]
    )


def test_dry_run_reports_and_writes_nothing(make_csv: CsvFactory) -> None:
    path = _rows(make_csv)
    before = _digest(path)
    report = migrate(path)
    assert (report.rows, report.legacy_marked, report.selected_reset) == (5, 2, 1)
    assert report.applied is False and report.backup_path is None
    assert _digest(path) == before
    assert not list(path.parent.glob("*.pre-protocol-*"))


def test_apply_marks_legacy_rows_resets_selected_and_backs_up(
    make_csv: CsvFactory,
) -> None:
    path = _rows(make_csv)
    original = pd.read_csv(path, dtype=str, keep_default_na=False)
    report = migrate(path, apply=True)
    assert report.applied is True
    assert report.backup_path is not None and report.backup_path.is_file()
    after = pd.read_csv(path, dtype=str, keep_default_na=False)
    by_id = after.set_index("mol_id")
    assert by_id.loc["a", "protocol"] == "1" and by_id.loc["b", "protocol"] == "1"
    assert by_id.loc["c", "status"] == "Pending"
    assert by_id.loc["d", "protocol"] == "" and by_id.loc["e", "status"] == "Skip"
    # The preliminary values stay put. The project's atomic writer prints floats
    # with six decimals, so compare numerically rather than textually.
    before = original.set_index("mol_id")
    for mol_id in ("a", "b"):
        assert float(by_id.loc[mol_id, "H298_pm7"]) == pytest.approx(
            float(before.loc[mol_id, "H298_pm7"]), abs=5e-7
        )


def test_apply_is_idempotent(make_csv: CsvFactory) -> None:
    path = _rows(make_csv)
    migrate(path, apply=True)
    once = _digest(path)
    report = migrate(path, apply=True)
    assert (report.legacy_marked, report.selected_reset) == (0, 0)
    assert report.already_marked == 2
    assert _digest(path) == once


def test_rows_already_in_protocol_2_keep_their_mark(make_csv: CsvFactory) -> None:
    path = make_csv([row("a", "OK", protocol="2"), row("b", "OK")])
    migrate(path, apply=True)
    by_id = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("mol_id")
    assert by_id.loc["a", "protocol"] == "2"
    assert by_id.loc["b", "protocol"] == "1"


def test_missing_file_is_an_error_not_an_empty_csv(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        migrate(tmp_path / "absent.csv")
    assert not (tmp_path / "absent.csv").exists()
