from __future__ import annotations

import sqlite3
from contextlib import closing
from pathlib import Path

import pytest

from semi_imperium.training.db import database_path
from semi_imperium.training.molecule_queue import (
    CbsReadError,
    TrainingQueue,
    read_activity,
    read_progress,
    select_reference_rows,
)
from tests.training.conftest import CbsFactory, row


def _cbs(make_cbs: CbsFactory) -> Path:
    return make_cbs(
        [
            row("big", nheavy=9),
            row("small", nheavy=2),
            row("mid", nheavy=5),
            row("tiny", nheavy=1),
        ]
    )


def _queue(tmp_path: Path, cbs: Path, **kwargs: int) -> TrainingQueue:
    return TrainingQueue(database_path(tmp_path / "store"), cbs, **kwargs)


def _status(tmp_path: Path, mol_id: str) -> str:
    with closing(sqlite3.connect(database_path(tmp_path / "store"))) as connection:
        found = connection.execute(
            "SELECT status FROM training_queue WHERE mol_id = ?", (mol_id,)
        ).fetchone()
    return str(found[0])


# ----------------------------------------------------------------- selection


def test_one_row_per_smiles_prefers_the_reference_conformer(
    make_cbs: CbsFactory,
) -> None:
    cbs = make_cbs(
        [
            row("c2", smiles="CC", conformer_rank_by_h298=2, H298_cbs=-10.0),
            row("c1", smiles="CC", conformer_rank_by_h298=1, H298_cbs=-11.0),
            row("solo", smiles="C"),
        ]
    )
    selected = select_reference_rows(cbs)
    assert sorted(selected["mol_id"]) == ["c1", "solo"]


def test_selection_breaks_rank_ties_by_h298_then_mol_id(make_cbs: CbsFactory) -> None:
    cbs = make_cbs(
        [
            row("b", smiles="CC", conformer_rank_by_h298=1, H298_cbs=-10.0),
            row("a", smiles="CC", conformer_rank_by_h298=1, H298_cbs=-10.0),
            row("c", smiles="CC", conformer_rank_by_h298=1, H298_cbs=-12.0),
        ]
    )
    assert list(select_reference_rows(cbs)["mol_id"]) == ["c"]
    cbs = make_cbs(
        [
            row("b", smiles="CC", conformer_rank_by_h298=1, H298_cbs=-10.0),
            row("a", smiles="CC", conformer_rank_by_h298=1, H298_cbs=-10.0),
        ]
    )
    assert list(select_reference_rows(cbs)["mol_id"]) == ["a"]


def test_missing_cbs_file_and_columns_are_explicit(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    with pytest.raises(FileNotFoundError):
        select_reference_rows(tmp_path / "absent.csv")
    with pytest.raises(CbsReadError, match="conformer_rank_by_h298"):
        select_reference_rows(make_cbs([{"mol_id": "a", "smiles": "C"}]))


# --------------------------------------------------------------------- queue


def test_seeding_is_idempotent_and_never_resets_state(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _queue(tmp_path, _cbs(make_cbs))
    assert queue.ensure_seeded() == 4
    queue.mark_running("tiny")
    assert queue.ensure_seeded() == 0
    assert _status(tmp_path, "tiny") == "running"


def test_queue_is_smallest_molecules_first_and_skips_finished_rows(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _queue(tmp_path, _cbs(make_cbs))
    queue.ensure_seeded()
    queue.mark_running("small")
    queue.record_success(
        "tiny", hof_kcal_mol=-1.0, conformers_evaluated=1, elapsed_s=6.0, batch_id="b"
    )
    queue.record_failure("big", "boom")
    ids = [m.mol_id for m in queue.queue()]
    # ok is gone; running (orphan) and rerun stay; smallest first
    assert ids == ["small", "mid", "big"]
    assert [m.mol_id for m in queue.queue(limit=2)] == ["small", "mid"]


def test_queued_molecule_carries_identity_and_cbs_value(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = make_cbs([row("r", smiles="[CH3]", multiplicity=2, charge=0, H298_cbs=34.8)])
    queue = _queue(tmp_path, cbs)
    queue.ensure_seeded()
    (molecule,) = queue.queue()
    assert (molecule.smiles, molecule.multiplicity, molecule.charge) == ("[CH3]", 2, 0)
    assert molecule.h298_cbs == pytest.approx(34.8)


def test_success_stores_value_delta_and_provenance(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = make_cbs([row("a", H298_cbs=-20.0)])
    queue = _queue(tmp_path, cbs)
    queue.ensure_seeded()
    queue.record_success(
        "a", hof_kcal_mol=-15.0, conformers_evaluated=3, elapsed_s=90.0, batch_id="b1"
    )
    with closing(sqlite3.connect(database_path(tmp_path / "store"))) as connection:
        stored = connection.execute(
            "SELECT status, h298_pm7, delta_kcal_mol, conformers_evaluated,"
            " elapsed_s, batch_id FROM training_queue WHERE mol_id = 'a'"
        ).fetchone()
    assert stored == ("ok", -15.0, -5.0, 3, 90.0, "b1")  # delta = CBS - PM7


def test_failures_rerun_until_the_limit_then_skip(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _queue(tmp_path, make_cbs([row("a")]), max_reruns=2)
    queue.ensure_seeded()
    assert queue.record_failure("a", "x") == "rerun"
    assert queue.record_failure("a", "x") == "skip"
    assert queue.queue() == []


# ------------------------------------------------------------------ progress


def test_progress_before_any_run_is_derived_from_the_cbs_without_writing(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    db = database_path(tmp_path / "store")
    progress = read_progress(db, _cbs(make_cbs))
    assert (progress.total, progress.pending, progress.ok) == (4, 4, 0)
    assert not db.exists()


def test_progress_counts_each_bucket_from_the_database(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = make_cbs([row(m) for m in ("a", "b", "c", "d", "e")])
    queue = _queue(tmp_path, cbs, max_reruns=1)
    queue.ensure_seeded()
    queue.record_success(
        "a", hof_kcal_mol=-1.0, conformers_evaluated=1, elapsed_s=1.0, batch_id="b"
    )
    queue.mark_running("b")
    queue.record_failure("c", "x")  # max_reruns=1 -> skip
    progress = read_progress(database_path(tmp_path / "store"), cbs)
    assert (progress.total, progress.ok, progress.running) == (5, 1, 1)
    assert (progress.pending, progress.rerun, progress.skipped) == (2, 0, 1)


# ------------------------------------------------------------------ activity


def test_activity_without_a_database_is_empty(tmp_path: Path) -> None:
    activity = read_activity(database_path(tmp_path / "store"))
    assert (activity.in_flight, activity.recent, activity.mean_elapsed_s) == (
        (),
        (),
        None,
    )
    assert not database_path(tmp_path / "store").exists()


def test_activity_mean_only_counts_molecules_since_the_run_started(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _queue(tmp_path, _cbs(make_cbs))
    queue.ensure_seeded()
    queue.record_success(
        "tiny",
        hof_kcal_mol=-1.0,
        conformers_evaluated=1,
        elapsed_s=100.0,
        batch_id="old",
    )
    db = database_path(tmp_path / "store")
    with closing(sqlite3.connect(db)) as connection, connection:
        connection.execute(
            "UPDATE training_queue SET updated_at = '2000-01-01T00:00:00+00:00'"
            " WHERE mol_id = 'tiny'"
        )
    since = "2020-01-01T00:00:00+00:00"
    queue.record_success(
        "small",
        hof_kcal_mol=-2.0,
        conformers_evaluated=2,
        elapsed_s=10.0,
        batch_id="new",
    )
    queue.mark_running("mid")
    activity = read_activity(db, since=since, recent=5)
    assert activity.mean_elapsed_s == 10.0
    assert [m.mol_id for m in activity.recent] == ["small", "tiny"]  # newest first
    (running,) = activity.in_flight
    assert (running.mol_id, running.nheavy) == ("mid", 5)
    assert running.started_at is not None
