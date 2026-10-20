"""Claims, leases and attempt-based acceptance of the training queue."""

from __future__ import annotations

import sqlite3
import threading
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

import pytest

from semi_imperium.training.db import database_path
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
