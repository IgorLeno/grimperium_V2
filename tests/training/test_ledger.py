from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from semi_imperium.training.ledger import (
    NotMigratedError,
    TrainingLedger,
    read_progress,
)
from tests.training.conftest import CsvFactory, row


def _migrated(make_csv: CsvFactory) -> Path:
    return make_csv(
        [
            row("big", "Pending", nheavy=9, protocol=""),
            row("small", "Pending", nheavy=2, protocol=""),
            row("mid", "Rerun", nheavy=5, protocol="", reruns=1),
            row("orphan", "Running", nheavy=3, protocol=""),
            row("sel", "Selected", nheavy=4, protocol=""),
            row("old", "OK", nheavy=1, protocol="1", H298_pm7=-3.83),
            row("done", "OK", nheavy=1, protocol="2", H298_pm7=-1.0),
            row("skip", "Skip", nheavy=1, protocol=""),
        ]
    )


def test_progress_counts_each_bucket(make_csv: CsvFactory) -> None:
    progress = read_progress(_migrated(make_csv))
    assert progress.total == 8
    assert (progress.protocol2, progress.protocol1, progress.ok_unmarked) == (1, 1, 0)
    assert (progress.pending, progress.selected, progress.rerun) == (2, 1, 1)
    assert (progress.running, progress.skipped) == (1, 1)
    assert progress.migrated is True


def test_progress_reports_unmigrated_csv(make_csv: CsvFactory) -> None:
    path = make_csv([row("a", "OK"), row("b", "Pending")])
    progress = read_progress(path)
    assert progress.migrated is False
    assert progress.ok_unmarked == 1


def test_ledger_refuses_an_unmigrated_csv(make_csv: CsvFactory) -> None:
    ledger = TrainingLedger(make_csv([row("a", "OK")]))
    with pytest.raises(NotMigratedError):
        ledger.load()


def test_queue_is_smallest_molecules_first_and_skips_finished_rows(
    make_csv: CsvFactory,
) -> None:
    ledger = TrainingLedger(_migrated(make_csv))
    ledger.load()
    assert [m.mol_id for m in ledger.queue()] == [
        "small",
        "orphan",
        "sel",
        "mid",
        "big",
    ]
    assert [m.mol_id for m in ledger.queue(limit=2)] == ["small", "orphan"]


def test_queue_casts_charge_and_multiplicity_to_int(make_csv: CsvFactory) -> None:
    path = make_csv([row("rad", "Pending", multiplicity=2.0, charge=-1.0, protocol="")])
    ledger = TrainingLedger(path)
    ledger.load()
    (molecule,) = ledger.queue()
    assert (molecule.charge, molecule.multiplicity) == (-1, 2)
    assert isinstance(molecule.charge, int) and isinstance(molecule.multiplicity, int)


def test_record_success_writes_value_delta_and_protocol(make_csv: CsvFactory) -> None:
    path = make_csv([row("a", "Pending", H298_cbs=-17.8, protocol="")])
    ledger = TrainingLedger(path)
    ledger.load()
    ledger.mark_running("a")
    ledger.record_success(
        "a",
        hof_kcal_mol=-16.5,
        conformers_evaluated=3,
        elapsed_s=120.0,
        batch_id="p2-x",
    )
    saved = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["a", "status"] == "OK"
    assert saved.loc["a", "protocol"] == "2"
    assert float(saved.loc["a", "H298_pm7"]) == pytest.approx(-16.5)
    assert float(saved.loc["a", "target_delta_kcalmol"]) == pytest.approx(-1.3)
    assert float(saved.loc["a", "total_time"]) == pytest.approx(2.0)
    assert saved.loc["a", "crest_status"] == "SUCCESS"
    assert saved.loc["a", "mopac_status"] == "OK"
    assert saved.loc["a", "batch_id"] == "p2-x"


def test_failures_become_rerun_until_the_limit_then_skip(make_csv: CsvFactory) -> None:
    path = make_csv([row("a", "Pending", protocol="")])
    ledger = TrainingLedger(path, max_reruns=2)
    ledger.load()
    assert ledger.record_failure("a", "saddle point") == "Rerun"
    assert ledger.record_failure("a", "saddle point") == "Skip"
    saved = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["a", "status"] == "Skip"
    assert saved.loc["a", "reruns"] in ("2", "2.0", "2.000000")
    assert saved.loc["a", "protocol"] == ""


def test_legacy_values_are_never_touched_by_a_run(make_csv: CsvFactory) -> None:
    path = make_csv(
        [
            row("old", "OK", protocol="1", H298_pm7=-3.83),
            row("new", "Pending", protocol=""),
        ]
    )
    ledger = TrainingLedger(path)
    ledger.load()
    ledger.mark_running("new")
    ledger.record_success(
        "new", hof_kcal_mol=-1.0, conformers_evaluated=1, elapsed_s=6.0, batch_id="p2"
    )
    saved = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("mol_id")
    assert (saved.loc["old", "status"], saved.loc["old", "protocol"]) == ("OK", "1")
    assert float(saved.loc["old", "H298_pm7"]) == pytest.approx(-3.83)
