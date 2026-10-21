from __future__ import annotations

import importlib
import sqlite3
from contextlib import closing
from pathlib import Path

import pytest

from semi_imperium.domain.identity import MolecularIdentity
from semi_imperium.training.db import database_path
from semi_imperium.training.molecule_queue import (
    CbsReadError,
    TrainingQueue,
    _identity_key,
    read_activity,
    read_progress,
    select_reference_rows,
)
from tests.training.conftest import (
    CbsFactory,
    row,
    set_ok,
    set_running,
    set_status,
)


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


#: Two CBS spellings of buta-1,3-dien-2-yl (cbs_52946/cbs_52947): different
#: SMILES, one ``MolecularIdentity``, CBS values 8.2 kcal/mol apart.
BUTADIENYL = (
    row("cbs_52946", smiles="C=C[C]=C", multiplicity=2, nheavy=4, H298_cbs=78.63363),
    row("cbs_52947", smiles="C=[C]C=C", multiplicity=2, nheavy=4, H298_cbs=86.79913),
)


def test_spellings_of_one_identity_keep_only_the_lowest_h298(
    make_cbs: CbsFactory,
) -> None:
    cbs = make_cbs([BUTADIENYL[1], BUTADIENYL[0], row("solo", smiles="C")])
    assert sorted(select_reference_rows(cbs)["mol_id"]) == ["cbs_52946", "solo"]


@pytest.mark.parametrize(
    ("smiles", "charge", "multiplicity"),
    [
        ("C=C[C]=C", 0, 2),
        ("[c]1ccc2c(c1)C=CC=C2", 0, 2),
        ("C[NH3+]", 1, 1),
        ("CCO", 0, 1),
    ],
)
def test_the_queue_key_is_the_molecule_id_the_worker_uses(
    smiles: str, charge: int, multiplicity: int
) -> None:
    chem = importlib.import_module("rdkit.Chem")
    expected = MolecularIdentity.from_smiles(
        smiles, charge=charge, multiplicity=multiplicity
    ).molecule_id
    assert _identity_key(chem, smiles, charge, multiplicity) == expected


def test_smiles_rdkit_cannot_parse_stay_one_row_each(make_cbs: CbsFactory) -> None:
    cbs = make_cbs([row("x", smiles="not-a-smiles"), row("y", smiles="also-not")])
    assert sorted(select_reference_rows(cbs)["mol_id"]) == ["x", "y"]


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
    (claimed,) = queue.claim("central", 1)
    assert queue.ensure_seeded() == 0
    assert _status(tmp_path, claimed.mol_id) == "running"


def _insert_legacy(tmp_path: Path, mol_id: str, status: str, reruns: int) -> None:
    """A row seeded before the identity dedup (one row per raw SMILES)."""
    source = next(r for r in BUTADIENYL if r["mol_id"] == mol_id)
    with closing(sqlite3.connect(database_path(tmp_path / "store"))) as connection:
        connection.execute(
            "INSERT INTO training_queue (mol_id, smiles, charge, multiplicity,"
            " nheavy, h298_cbs, status, reruns) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (
                mol_id,
                source["smiles"],
                source["charge"],
                source["multiplicity"],
                source["nheavy"],
                source["H298_cbs"],
                status,
                reruns,
            ),
        )
        connection.commit()


def _failures(tmp_path: Path) -> list[tuple[str, str]]:
    with closing(sqlite3.connect(database_path(tmp_path / "store"))) as connection:
        return [
            (str(m), str(r))
            for m, r in connection.execute(
                "SELECT mol_id, reason FROM training_failure ORDER BY id"
            )
        ]


def test_seeding_skips_an_open_duplicate_seeded_before_the_dedup(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _queue(tmp_path, make_cbs(list(BUTADIENYL)))
    # The production state: both spellings in ``rerun`` after colliding.
    queue.ensure_seeded()
    _insert_legacy(tmp_path, "cbs_52947", "rerun", 1)
    set_status(database_path(tmp_path / "store"), "cbs_52946", "rerun", reruns=1)

    assert queue.ensure_seeded() == 0
    assert _status(tmp_path, "cbs_52947") == "skip"
    assert _status(tmp_path, "cbs_52946") == "rerun"
    ((mol_id, reason),) = _failures(tmp_path)
    assert mol_id == "cbs_52947"
    assert "cbs_52946" in reason
    # A second seeding finds nothing left to skip and logs nothing more.
    queue.ensure_seeded()
    assert len(_failures(tmp_path)) == 1


@pytest.mark.parametrize("status", ["ok", "running", "skip"])
def test_seeding_never_touches_a_duplicate_that_is_not_open(
    tmp_path: Path, make_cbs: CbsFactory, status: str
) -> None:
    queue = _queue(tmp_path, make_cbs(list(BUTADIENYL)))
    queue.ensure_seeded()
    _insert_legacy(tmp_path, "cbs_52947", status, 0)
    queue.ensure_seeded()
    assert _status(tmp_path, "cbs_52947") == status
    assert _failures(tmp_path) == []


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
    _queue(tmp_path, cbs).ensure_seeded()
    db = database_path(tmp_path / "store")
    set_ok(db, "a", hof=-1.0)
    set_running(db, "b")
    set_status(db, "c", "skip", reruns=1)
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
    _queue(tmp_path, _cbs(make_cbs)).ensure_seeded()
    db = database_path(tmp_path / "store")
    set_ok(db, "tiny", hof=-1.0, elapsed_s=100.0, batch_id="old")
    with closing(sqlite3.connect(db)) as connection, connection:
        connection.execute(
            "UPDATE training_queue SET updated_at = '2000-01-01T00:00:00+00:00'"
            " WHERE mol_id = 'tiny'"
        )
    since = "2020-01-01T00:00:00+00:00"
    set_ok(db, "small", hof=-2.0, conformers=2, elapsed_s=10.0, batch_id="new")
    set_running(db, "mid", cell_id="cell-x")
    activity = read_activity(db, since=since, recent=5)
    assert activity.mean_elapsed_s == 10.0
    assert [m.mol_id for m in activity.recent] == ["small", "tiny"]  # newest first
    (running,) = activity.in_flight
    assert (running.mol_id, running.nheavy, running.cell_id) == ("mid", 5, "cell-x")
    assert running.started_at is not None
