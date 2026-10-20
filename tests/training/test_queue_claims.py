"""Claims, leases and attempt-based acceptance of the training queue."""

from __future__ import annotations

import sqlite3
import threading
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

import pytest

from semi_imperium.training.db import database_path
from semi_imperium.training.geometry_store import GeometryRecord, GeometryStore
from semi_imperium.training.molecule_queue import TrainingQueue
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


def _db(tmp_path: Path) -> Path:
    return database_path(tmp_path / "store")


def _seeded(tmp_path: Path, cbs: Path, **kwargs: int) -> TrainingQueue:
    queue = TrainingQueue(_db(tmp_path), cbs, **kwargs)
    queue.ensure_seeded()
    return queue


def _row(tmp_path: Path, mol_id: str) -> sqlite3.Row:
    with closing(sqlite3.connect(_db(tmp_path))) as connection:
        connection.row_factory = sqlite3.Row
        found = connection.execute(
            "SELECT * FROM training_queue WHERE mol_id = ?", (mol_id,)
        ).fetchone()
    assert found is not None
    return found


def _attempts(tmp_path: Path) -> list[sqlite3.Row]:
    with closing(sqlite3.connect(_db(tmp_path))) as connection:
        connection.row_factory = sqlite3.Row
        return connection.execute(
            "SELECT * FROM training_attempt ORDER BY claimed_at, mol_id"
        ).fetchall()


def _set_status(tmp_path: Path, mol_id: str, status: str, reruns: int = 0) -> None:
    with closing(sqlite3.connect(_db(tmp_path))) as connection:
        connection.execute(
            "UPDATE training_queue SET status = ?, reruns = ? WHERE mol_id = ?",
            (status, reruns, mol_id),
        )
        connection.commit()


# --------------------------------------------------------------------- claim


def test_claim_hands_out_the_smallest_molecules_and_marks_them(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, _cbs(make_cbs))
    before = datetime.now(timezone.utc)
    claimed = queue.claim("central", 2, lease_s=300.0)
    assert [m.mol_id for m in claimed] == ["tiny", "small"]
    assert len({m.attempt_id for m in claimed}) == 2
    for molecule in claimed:
        stored = _row(tmp_path, molecule.mol_id)
        assert stored["status"] == "running"
        assert stored["cell_id"] == "central"
        assert stored["attempt_id"] == molecule.attempt_id
        assert datetime.fromisoformat(stored["lease_until"]) > before
    assert _row(tmp_path, "mid")["status"] == "pending"


def test_claim_carries_identity_and_cbs_value(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = make_cbs([row("r", smiles="[CH3]", multiplicity=2, charge=0, H298_cbs=34.8)])
    (molecule,) = _seeded(tmp_path, cbs).claim("central", 5)
    assert (molecule.smiles, molecule.multiplicity, molecule.charge) == ("[CH3]", 2, 0)
    assert molecule.h298_cbs == pytest.approx(34.8)


def test_claim_logs_one_attempt_per_molecule(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    claimed = _seeded(tmp_path, _cbs(make_cbs)).claim("cell-a", 3)
    attempts = _attempts(tmp_path)
    assert {(a["attempt_id"], a["mol_id"], a["cell_id"]) for a in attempts} == {
        (m.attempt_id, m.mol_id, "cell-a") for m in claimed
    }
    assert all(a["outcome"] is None and a["claimed_at"] for a in attempts)


def test_claim_takes_pending_and_rerun_but_never_running_ok_or_skip(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, _cbs(make_cbs))
    _set_status(tmp_path, "tiny", "ok")
    _set_status(tmp_path, "small", "skip", reruns=3)
    _set_status(tmp_path, "mid", "running")
    _set_status(tmp_path, "big", "rerun", reruns=1)
    assert [m.mol_id for m in queue.claim("central", 10)] == ["big"]
    assert queue.claim("central", 10) == []


def test_claim_of_zero_is_empty_and_negative_is_refused(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, _cbs(make_cbs))
    assert queue.claim("central", 0) == []
    with pytest.raises(ValueError):
        queue.claim("central", -1)
    with pytest.raises(ValueError):
        queue.claim("", 1)


def test_concurrent_claimers_never_share_a_molecule(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    cbs = make_cbs([row(f"m{i:02d}", nheavy=i % 7 + 1) for i in range(40)])
    _seeded(tmp_path, cbs)
    taken: list[str] = []
    lock = threading.Lock()

    def claimer(cell: str) -> None:
        queue = TrainingQueue(_db(tmp_path), cbs)  # own connections, like processes
        while batch := queue.claim(cell, 1):
            with lock:
                taken.extend(m.mol_id for m in batch)

    threads = [threading.Thread(target=claimer, args=(f"c{i}",)) for i in range(6)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert sorted(taken) == sorted(f"m{i:02d}" for i in range(40))


# ------------------------------------------------------------ lease, release


def test_an_expired_lease_never_returns_the_molecule_by_itself(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a")]))
    queue.claim("cell-a", 1, lease_s=-60.0)  # already lapsed
    assert queue.claim("cell-b", 1) == []
    assert _row(tmp_path, "a")["cell_id"] == "cell-a"


def test_renew_extends_only_the_cells_own_running_attempts(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, _cbs(make_cbs))
    (mine,) = queue.claim("cell-a", 1, lease_s=10.0)
    (theirs,) = queue.claim("cell-b", 1, lease_s=10.0)
    before = _row(tmp_path, theirs.mol_id)["lease_until"]
    started = _row(tmp_path, mine.mol_id)["updated_at"]

    renewed = queue.renew(
        "cell-a", [mine.attempt_id, theirs.attempt_id, "unknown"], lease_s=600.0
    )

    assert renewed == 1
    assert _row(tmp_path, mine.mol_id)["lease_until"] > before
    assert _row(tmp_path, theirs.mol_id)["lease_until"] == before
    # the dashboard reads updated_at as "started at": renewal must not move it
    assert _row(tmp_path, mine.mol_id)["updated_at"] == started


def test_release_returns_a_cells_molecules_to_the_queue(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, _cbs(make_cbs))
    _set_status(tmp_path, "small", "rerun", reruns=1)
    mine = queue.claim("cell-a", 2)  # tiny (pending) and small (rerun)
    (theirs,) = queue.claim("cell-b", 1)

    released = queue.release("cell-a")

    assert sorted(released) == ["small", "tiny"]
    tiny, small = _row(tmp_path, "tiny"), _row(tmp_path, "small")
    assert (tiny["status"], tiny["reruns"]) == ("pending", 0)
    assert (small["status"], small["reruns"]) == ("rerun", 1)  # no rerun spent
    for stored in (tiny, small):
        assert (stored["cell_id"], stored["attempt_id"], stored["lease_until"]) == (
            None,
            None,
            None,
        )
    assert _row(tmp_path, theirs.mol_id)["cell_id"] == "cell-b"
    by_id = {a["attempt_id"]: a for a in _attempts(tmp_path)}
    for molecule in mine:
        assert by_id[molecule.attempt_id]["outcome"] == "released"
        assert by_id[molecule.attempt_id]["released_at"]
    assert by_id[theirs.attempt_id]["outcome"] is None
    # released molecules are claimable again, by anyone
    assert {m.mol_id for m in queue.claim("cell-c", 5)} == {"tiny", "small", "big"}


def test_release_can_take_back_legacy_rows_with_no_cell(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, _cbs(make_cbs))
    _set_status(tmp_path, "mid", "running")  # orphan from before the migration
    queue.claim("cell-b", 1)
    assert queue.release("central") == []
    assert queue.release("central", include_unowned=True) == ["mid"]
    assert _row(tmp_path, "mid")["status"] == "pending"
    assert _row(tmp_path, "tiny")["cell_id"] == "cell-b"


# ---------------------------------------------------------------- acceptance


def _geometry(mol_id: str, hof: float = -15.0) -> GeometryRecord:
    return GeometryRecord(
        mol_id=mol_id,
        smiles=f"smiles-{mol_id}",
        hamiltonian="PM7",
        protocol="2",
        signature="sig",
        calculation_id=f"calc-{mol_id}-{hof}",
        conformer_index=1,
        elements=("C",),
        coordinates=((0.0, 0.0, 0.0),),
        hof_kcal_mol=hof,
        lowest_frequency_cm1=100.0,
        minimum_confirmed=True,
    )


def _succeed(queue: TrainingQueue, attempt_id: str, mol_id: str, hof: float) -> str:
    return queue.accept_success(
        attempt_id,
        geometry=_geometry(mol_id, hof),
        conformers_evaluated=3,
        elapsed_s=90.0,
        batch_id="b1",
    )


def _attempt(tmp_path: Path, attempt_id: str) -> sqlite3.Row:
    (found,) = [a for a in _attempts(tmp_path) if a["attempt_id"] == attempt_id]
    return found


def test_a_valid_result_is_accepted_with_its_geometry(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a", H298_cbs=-20.0)]))
    (claimed,) = queue.claim("cell-a", 1)

    assert _succeed(queue, claimed.attempt_id, "a", -15.0) == "accepted"

    stored = _row(tmp_path, "a")
    assert (stored["status"], stored["h298_pm7"], stored["delta_kcal_mol"]) == (
        "ok",
        -15.0,
        -5.0,  # delta = CBS - PM7
    )
    assert (stored["conformers_evaluated"], stored["elapsed_s"]) == (3, 90.0)
    assert (stored["batch_id"], stored["cell_id"], stored["attempt_id"]) == (
        "b1",
        "cell-a",
        claimed.attempt_id,
    )
    assert stored["lease_until"] is None
    geometry = GeometryStore(_db(tmp_path)).get("a")
    assert geometry is not None and geometry.hof_kcal_mol == -15.0
    assert _attempt(tmp_path, claimed.attempt_id)["outcome"] == "accepted"


def test_the_same_attempt_reported_twice_changes_nothing(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a")]))
    (claimed,) = queue.claim("cell-a", 1)
    _succeed(queue, claimed.attempt_id, "a", -15.0)

    assert _succeed(queue, claimed.attempt_id, "a", -99.0) == "duplicate"
    assert queue.accept_failure(claimed.attempt_id, "late") == "duplicate"

    assert _row(tmp_path, "a")["h298_pm7"] == -15.0
    geometry = GeometryStore(_db(tmp_path)).get("a")
    assert geometry is not None and geometry.hof_kcal_mol == -15.0


def test_a_result_after_release_is_accepted_while_the_molecule_is_open(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a")]))
    (claimed,) = queue.claim("cell-a", 1)
    queue.release("cell-a")

    assert _succeed(queue, claimed.attempt_id, "a", -15.0) == "accepted"
    assert _row(tmp_path, "a")["status"] == "ok"
    assert queue.claim("cell-b", 1) == []


def test_the_first_valid_result_wins_and_the_late_one_is_logged(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a")]))
    (first,) = queue.claim("cell-a", 1)
    queue.release("cell-a")
    (second,) = queue.claim("cell-b", 1)

    assert _succeed(queue, second.attempt_id, "a", -15.0) == "accepted"
    assert _succeed(queue, first.attempt_id, "a", -16.0) == "discarded"

    stored = _row(tmp_path, "a")
    assert (stored["h298_pm7"], stored["cell_id"]) == (-15.0, "cell-b")
    late = _attempt(tmp_path, first.attempt_id)
    assert late["outcome"] == "discarded" and late["finished_at"]
    assert second.attempt_id in late["note"]


def test_a_released_attempt_can_beat_the_cell_that_took_over(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a")]))
    (first,) = queue.claim("cell-a", 1)
    queue.release("cell-a")
    (second,) = queue.claim("cell-b", 1)

    assert _succeed(queue, first.attempt_id, "a", -16.0) == "accepted"
    assert queue.accept_failure(second.attempt_id, "crest timeout") == "discarded"
    assert _succeed(queue, second.attempt_id, "a", -15.0) == "duplicate"

    stored = _row(tmp_path, "a")
    assert (stored["status"], stored["h298_pm7"], stored["reruns"]) == (
        "ok",
        -16.0,
        0,
    )


def test_a_failure_of_the_current_attempt_spends_a_rerun_then_skips(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a")]), max_reruns=2)
    (claimed,) = queue.claim("cell-a", 1)
    assert queue.accept_failure(claimed.attempt_id, "saddle") == "rerun"
    stored = _row(tmp_path, "a")
    assert (stored["status"], stored["reruns"], stored["cell_id"]) == (
        "rerun",
        1,
        None,
    )
    failed = _attempt(tmp_path, claimed.attempt_id)
    assert (failed["outcome"], failed["note"]) == ("failed", "saddle")

    (again,) = queue.claim("cell-a", 1)
    assert queue.accept_failure(again.attempt_id, "saddle") == "skip"
    assert queue.claim("cell-a", 1) == []


def test_a_failure_of_a_released_attempt_spends_nothing(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a")]))
    (claimed,) = queue.claim("cell-a", 1)
    queue.release("cell-a")

    assert queue.accept_failure(claimed.attempt_id, "power cut") == "discarded"

    stored = _row(tmp_path, "a")
    assert (stored["status"], stored["reruns"]) == ("pending", 0)
    assert "power cut" in _attempt(tmp_path, claimed.attempt_id)["note"]


def test_unknown_attempts_and_inconsistent_results_are_refused(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    queue = _seeded(tmp_path, make_cbs([row("a"), row("b")]))
    (claimed,) = queue.claim("cell-a", 1)
    with pytest.raises(KeyError):
        _succeed(queue, "nope", "a", -15.0)
    with pytest.raises(KeyError):
        queue.accept_failure("nope", "x")
    with pytest.raises(ValueError, match="mol_id"):
        _succeed(queue, claimed.attempt_id, "b", -15.0)
    with pytest.raises(ValueError, match="non-finite"):
        _succeed(queue, claimed.attempt_id, "a", float("nan"))
    assert _row(tmp_path, "a")["status"] == "running"
