"""Cell registry: who may compute for the central, and how each cell stands."""

from __future__ import annotations

import hashlib
from contextlib import closing
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from semi_imperium.training.cells import (
    ADMISSION_BLOCKED,
    ADMISSION_PASSED,
    ADMISSION_PENDING,
    Cell,
    CellRegistry,
    cell_state,
)
from semi_imperium.training.db import connect

NOW = datetime(2026, 10, 20, 12, 0, tzinfo=timezone.utc)
VERSIONS = {"crest": "3.0.2", "crest_sha256": "abc", "xtb": "6.7.1"}


def _registry(tmp_path: Path) -> CellRegistry:
    return CellRegistry(tmp_path / "training.sqlite3")


def _register(registry: CellRegistry, cell_id: str, refusal: str | None = None) -> Cell:
    return registry.record_registration(
        cell_id,
        host="lab-pc",
        cores=16,
        workers=4,
        crest_threads=4,
        program_versions=VERSIONS,
        protocol_signature="sig",
        refusal=refusal,
    )


def test_create_returns_the_token_once_and_stores_only_its_hash(
    tmp_path: Path,
) -> None:
    registry = _registry(tmp_path)
    cell, token = registry.create("PC do Prof. Silva")
    assert len(token) >= 40
    assert cell.name == "PC do Prof. Silva"
    assert cell.admission == ADMISSION_PENDING
    assert cell.registered_at is None and cell.last_seen_at is None
    # WAL mode: recent pages live in the -wal file, so read every database file.
    for path in tmp_path.glob("training.sqlite3*"):
        assert token.encode() not in path.read_bytes()
    with closing(connect(tmp_path / "training.sqlite3")) as connection:
        stored = connection.execute(
            "SELECT token_hash FROM cells WHERE cell_id = ?", (cell.cell_id,)
        ).fetchone()[0]
    assert stored == hashlib.sha256(token.encode()).hexdigest()


def test_each_cell_gets_its_own_id_and_token_never_central(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    first, token_a = registry.create("a")
    second, token_b = registry.create("a")
    assert first.cell_id != second.cell_id and token_a != token_b
    assert "central" not in (first.cell_id, second.cell_id)


def test_authenticate_maps_a_token_to_its_cell(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    cell, token = registry.create("a")
    found = registry.authenticate(token)
    assert found is not None and found.cell_id == cell.cell_id
    assert registry.authenticate("not-a-token") is None
    assert registry.authenticate("") is None


def test_revoked_token_no_longer_authenticates(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    cell, token = registry.create("a")
    revoked = registry.revoke(cell.cell_id)
    assert revoked.revoked_at is not None
    assert registry.authenticate(token) is None
    with closing(connect(tmp_path / "training.sqlite3")) as connection:
        assert (
            connection.execute(
                "SELECT token_hash FROM cells WHERE cell_id = ?", (cell.cell_id,)
            ).fetchone()[0]
            is None
        )


def test_registration_records_the_machine_and_clears_an_old_refusal(
    tmp_path: Path,
) -> None:
    registry = _registry(tmp_path)
    cell, _ = registry.create("a")
    refused = _register(registry, cell.cell_id, refusal="crest: 3.0.1 ≠ 3.0.2")
    assert refused.refusal == "crest: 3.0.1 ≠ 3.0.2"
    assert refused.registered_at is None and refused.last_seen_at is not None
    accepted = _register(registry, cell.cell_id)
    assert accepted.refusal is None and accepted.registered_at is not None
    assert accepted.host == "lab-pc" and accepted.cores == 16
    assert accepted.workers == 4 and accepted.crest_threads == 4
    assert accepted.program_versions == VERSIONS
    assert accepted.protocol_signature == "sig"


def test_rename_touch_admission_and_listing(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    cell, _ = registry.create("a")
    other, _ = registry.create("b")
    assert registry.rename(cell.cell_id, "Lab 2").name == "Lab 2"
    registry.touch(other.cell_id)
    assert registry.get(other.cell_id).last_seen_at is not None
    assert (
        registry.set_admission(cell.cell_id, ADMISSION_PASSED).admission
        == ADMISSION_PASSED
    )
    assert [c.cell_id for c in registry.list()] == [cell.cell_id, other.cell_id]


def test_unknown_cell_and_bad_values_are_explicit(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    cell, _ = registry.create("a")
    for call in (
        lambda: registry.get("nope"),
        lambda: registry.rename("nope", "x"),
        lambda: registry.revoke("nope"),
        lambda: registry.touch("nope"),
    ):
        with pytest.raises(KeyError):
            call()
    with pytest.raises(ValueError):
        registry.create("  ")
    with pytest.raises(ValueError):
        registry.rename(cell.cell_id, "")
    with pytest.raises(ValueError):
        registry.set_admission(cell.cell_id, "maybe")


def _cell(**changes: object) -> Cell:
    base: dict[str, object] = {
        "cell_id": "cell-1",
        "name": "a",
        "host": "h",
        "cores": 8,
        "workers": 2,
        "crest_threads": 4,
        "program_versions": {},
        "protocol_signature": "sig",
        "admission": ADMISSION_PASSED,
        "created_at": NOW.isoformat(),
        "registered_at": NOW.isoformat(),
        "last_seen_at": NOW.isoformat(),
        "revoked_at": None,
        "refusal": None,
    }
    base.update(changes)
    return Cell(**base)  # type: ignore[arg-type]


def _seen(seconds_ago: float) -> str:
    return (NOW - timedelta(seconds=seconds_ago)).isoformat()


@pytest.mark.parametrize(
    ("seconds_ago", "expected"),
    [(0, "online"), (89, "online"), (90, "late"), (300, "late"), (301, "silent")],
)
def test_state_follows_the_central_clock(seconds_ago: float, expected: str) -> None:
    assert cell_state(_cell(last_seen_at=_seen(seconds_ago)), NOW) == expected


def test_state_precedence() -> None:
    assert cell_state(_cell(registered_at=None, last_seen_at=None), NOW) == (
        "never_seen"
    )
    assert cell_state(_cell(refusal="x"), NOW) == "version_mismatch"
    assert cell_state(_cell(admission=ADMISSION_BLOCKED), NOW) == "blocked"
    assert cell_state(_cell(revoked_at=NOW.isoformat(), refusal="x"), NOW) == "revoked"
