"""Admission test: a new cell's five reference HOFs against the central's."""

from __future__ import annotations

from pathlib import Path

import pytest

from semi_imperium.training.cells import (
    ADMISSION_BLOCKED,
    ADMISSION_PASSED,
    ADMISSION_PENDING,
    ADMISSION_PREFIX,
    ADMISSION_REFERENCES,
    AdmissionTest,
    AdmissionUnavailable,
    CellRegistry,
)
from semi_imperium.training.db import database_path
from semi_imperium.training.molecule_queue import TrainingQueue
from tests.training.conftest import (
    CbsFactory,
    accept_references,
    reference_rows,
    row,
)

HOF = -48.9388


@pytest.fixture
def db(tmp_path: Path, make_cbs: CbsFactory) -> Path:
    path = database_path(tmp_path / "store")
    TrainingQueue(path, make_cbs([*reference_rows(), row("a")])).ensure_seeded()
    accept_references(path, hof=HOF)
    return path


def _cell(db: Path) -> str:
    cell, _token = CellRegistry(db).create("Lab")
    return cell.cell_id


def _all(test: AdmissionTest, cell_id: str) -> dict[str, str]:
    """Claim every reference; ``mol_id -> attempt_id``."""
    return {m.mol_id: m.attempt_id for m in test.claim(cell_id, 10)}


def test_claim_hands_the_five_references_once(db: Path) -> None:
    test, cell_id = AdmissionTest(db), _cell(db)
    first = test.claim(cell_id, 2)
    rest = test.claim(cell_id, 10)
    assert [m.mol_id for m in first + rest] == list(ADMISSION_REFERENCES)
    assert all(m.attempt_id.startswith(ADMISSION_PREFIX) for m in first + rest)
    assert (first[1].smiles, first[1].multiplicity) == ("[H][C]=C", 2)
    assert test.claim(cell_id, 10) == []
    assert test.held(cell_id) == {m.attempt_id for m in first + rest}
    assert test.claim(_cell(db), 1)[0].mol_id == ADMISSION_REFERENCES[0]


def test_five_within_tolerance_admit_the_cell(db: Path) -> None:
    test, cell_id = AdmissionTest(db), _cell(db)
    attempts = _all(test, cell_id)
    for n, attempt_id in enumerate(attempts.values()):
        # 0.1 kcal/mol is inside, by either side.
        hof = HOF + (0.1 if n % 2 else -0.1)
        assert test.record(attempt_id, hof_kcal_mol=hof, failure=None) == "accepted"
        expected = ADMISSION_PASSED if n == 4 else ADMISSION_PENDING
        assert CellRegistry(db).get(cell_id).admission == expected
    assert test.held(cell_id) == set()
    assert test.note(cell_id, ADMISSION_PASSED) is None


def test_progress_while_pending(db: Path) -> None:
    test, cell_id = AdmissionTest(db), _cell(db)
    assert test.note(cell_id, ADMISSION_PENDING) is None
    attempts = _all(test, cell_id)
    test.record(attempts["cbs_09698"], hof_kcal_mol=HOF, failure=None)
    assert "1 de 5" in str(test.note(cell_id, ADMISSION_PENDING))


def test_one_outside_tolerance_blocks_with_the_reason(db: Path) -> None:
    test, cell_id = AdmissionTest(db), _cell(db)
    attempts = _all(test, cell_id)
    test.record(attempts["cbs_09698"], hof_kcal_mol=HOF + 0.1001, failure=None)
    assert CellRegistry(db).get(cell_id).admission == ADMISSION_BLOCKED
    note = str(test.note(cell_id, ADMISSION_BLOCKED))
    assert note.startswith("cbs_09698 (CO): HOF -48.83870 contra -48.93880")
    assert "0.1 kcal/mol" in note
    # The test is over: what is still running no longer counts.
    assert test.held(cell_id) == set()
    late = test.record(attempts["cbs_00187"], hof_kcal_mol=HOF, failure=None)
    assert late == "discarded"


def test_a_failed_reference_blocks_with_the_failure(db: Path) -> None:
    test, cell_id = AdmissionTest(db), _cell(db)
    attempts = _all(test, cell_id)
    test.record(attempts["cbs_24671"], hof_kcal_mol=None, failure="crest timeout")
    assert CellRegistry(db).get(cell_id).admission == ADMISSION_BLOCKED
    assert test.note(cell_id, ADMISSION_BLOCKED) == (
        "cbs_24671 (C[N+](=O)[O-]) falhou: crest timeout"
    )


def test_a_repeated_result_is_a_duplicate(db: Path) -> None:
    test, cell_id = AdmissionTest(db), _cell(db)
    attempt_id = _all(test, cell_id)["cbs_09698"]
    test.record(attempt_id, hof_kcal_mol=HOF, failure=None)
    assert test.record(attempt_id, hof_kcal_mol=HOF, failure=None) == "duplicate"


def test_a_restart_hands_the_open_references_again(db: Path) -> None:
    test, cell_id = AdmissionTest(db), _cell(db)
    attempts = _all(test, cell_id)
    test.record(attempts["cbs_09698"], hof_kcal_mol=HOF, failure=None)
    test.release(cell_id)
    assert test.held(cell_id) == set()
    again = _all(test, cell_id)
    assert list(again) == list(ADMISSION_REFERENCES[1:])
    # The released attempt was computed: it still counts, once.
    old = attempts["cbs_00187"]
    assert test.record(old, hof_kcal_mol=HOF, failure=None) == "accepted"
    new = again["cbs_00187"]
    assert test.record(new, hof_kcal_mol=HOF, failure=None) == "duplicate"


def test_back_to_pending_starts_the_test_over(db: Path) -> None:
    test, registry, cell_id = AdmissionTest(db), CellRegistry(db), _cell(db)
    attempts = _all(test, cell_id)
    test.record(attempts["cbs_09698"], hof_kcal_mol=HOF + 5, failure=None)
    registry.set_admission(cell_id, ADMISSION_PENDING)
    assert test.note(cell_id, ADMISSION_PENDING) is None
    assert test.record(attempts["cbs_00187"], hof_kcal_mol=HOF, failure=None) == (
        "discarded"
    )
    assert list(_all(test, cell_id)) == list(ADMISSION_REFERENCES)


def test_a_manual_decision_ends_the_open_attempts(db: Path) -> None:
    test, registry, cell_id = AdmissionTest(db), CellRegistry(db), _cell(db)
    attempts = _all(test, cell_id)
    registry.set_admission(cell_id, ADMISSION_PASSED)
    assert test.held(cell_id) == set()
    assert test.record(attempts["cbs_09698"], hof_kcal_mol=HOF, failure=None) == (
        "discarded"
    )


def test_without_the_central_values_there_is_no_test(
    tmp_path: Path, make_cbs: CbsFactory
) -> None:
    path = database_path(tmp_path / "store")
    TrainingQueue(path, make_cbs(reference_rows())).ensure_seeded()
    with pytest.raises(AdmissionUnavailable, match="cbs_09698"):
        AdmissionTest(path).claim(_cell(path), 1)


def test_unknown_attempt(db: Path) -> None:
    with pytest.raises(KeyError):
        AdmissionTest(db).owner(f"{ADMISSION_PREFIX}nope")
    with pytest.raises(KeyError):
        AdmissionTest(db).record("adm-nope", hof_kcal_mol=HOF, failure=None)
