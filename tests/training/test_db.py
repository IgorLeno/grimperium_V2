"""Schema migration of the training database (it already holds production data)."""

from __future__ import annotations

import sqlite3
from contextlib import closing
from pathlib import Path

from semi_imperium.training.db import connect

# ``training_queue`` as released before the cell columns existed.
_OLD_QUEUE = """
CREATE TABLE training_queue (
    mol_id TEXT PRIMARY KEY,
    smiles TEXT NOT NULL,
    charge INTEGER NOT NULL,
    multiplicity INTEGER NOT NULL,
    nheavy INTEGER NOT NULL,
    h298_cbs REAL NOT NULL,
    status TEXT NOT NULL DEFAULT 'pending',
    reruns INTEGER NOT NULL DEFAULT 0,
    h298_pm7 REAL,
    delta_kcal_mol REAL,
    conformers_evaluated INTEGER,
    elapsed_s REAL,
    batch_id TEXT,
    updated_at TEXT
);
INSERT INTO training_queue (mol_id, smiles, charge, multiplicity, nheavy,
    h298_cbs, status, reruns, h298_pm7, delta_kcal_mol, updated_at)
VALUES ('done', 'C', 0, 1, 1, -17.8, 'ok', 0, -12.0, -5.8, '2026-10-20T10:00:00'),
       ('todo', 'CC', 0, 1, 2, -20.0, 'rerun', 1, NULL, NULL, NULL);
"""

_NEW_COLUMNS = {"cell_id", "attempt_id", "lease_until"}


def _old_database(tmp_path: Path) -> Path:
    path = tmp_path / "training.sqlite3"
    with closing(sqlite3.connect(path)) as connection:
        connection.executescript(_OLD_QUEUE)
    return path


def _columns(connection: sqlite3.Connection, table: str) -> set[str]:
    return {str(r[1]) for r in connection.execute(f"PRAGMA table_info({table})")}


def test_connect_adds_the_cell_columns_and_keeps_every_row(tmp_path: Path) -> None:
    path = _old_database(tmp_path)
    with closing(connect(path)) as connection:
        assert _columns(connection, "training_queue") >= _NEW_COLUMNS
        rows = connection.execute(
            "SELECT mol_id, status, reruns, h298_pm7, cell_id, attempt_id, lease_until"
            " FROM training_queue ORDER BY mol_id"
        ).fetchall()
    assert [tuple(r) for r in rows] == [
        ("done", "ok", 0, -12.0, None, None, None),
        ("todo", "rerun", 1, None, None, None, None),
    ]


def test_migration_is_idempotent(tmp_path: Path) -> None:
    path = _old_database(tmp_path)
    for _ in range(3):
        with closing(connect(path)) as connection:
            columns = [
                str(r[1])
                for r in connection.execute("PRAGMA table_info(training_queue)")
            ]
    assert len(columns) == len(set(columns))
    assert set(columns) >= _NEW_COLUMNS


def test_a_fresh_database_has_the_cell_columns_and_the_attempt_log(
    tmp_path: Path,
) -> None:
    with closing(connect(tmp_path / "fresh.sqlite3")) as connection:
        assert _columns(connection, "training_queue") >= _NEW_COLUMNS
        assert _columns(connection, "training_attempt") >= {
            "attempt_id",
            "mol_id",
            "cell_id",
            "claimed_at",
            "released_at",
            "finished_at",
            "outcome",
            "note",
        }
        indexes = {
            str(r[1]) for r in connection.execute("PRAGMA index_list(training_queue)")
        }
    assert "training_queue_cell" in indexes


def test_connect_adds_the_cells_table_to_an_old_database(tmp_path: Path) -> None:
    path = _old_database(tmp_path)
    for _ in range(2):  # idempotent
        with closing(connect(path)) as connection:
            assert _columns(connection, "cells") >= {
                "cell_id",
                "name",
                "token_hash",
                "admission",
                "host",
                "cores",
                "workers",
                "crest_threads",
                "program_versions",
                "protocol_signature",
                "refusal",
                "created_at",
                "registered_at",
                "last_seen_at",
                "revoked_at",
            }
            count = connection.execute("SELECT COUNT(*) FROM training_queue")
            assert count.fetchone()[0] == 2
