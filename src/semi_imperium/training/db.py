"""SQLite plumbing shared by the training stores (one file, WAL mode).

Several processes touch this file at once (the coordinator, the dashboard API
and the CLI), so every call opens a short-lived connection instead of sharing
one, and the schema is created idempotently on connect. Columns added after a
database already held data are migrated there too (see ``_migrate``).
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

DB_FILENAME = "training.sqlite3"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS final_geometry (
    mol_id TEXT PRIMARY KEY,
    smiles TEXT NOT NULL,
    hamiltonian TEXT NOT NULL,
    protocol TEXT NOT NULL,
    signature TEXT NOT NULL,
    calculation_id TEXT NOT NULL,
    conformer_index INTEGER,
    n_atoms INTEGER NOT NULL,
    elements TEXT NOT NULL,
    coordinates TEXT NOT NULL,
    hof_kcal_mol REAL,
    lowest_frequency_cm1 REAL,
    minimum_confirmed INTEGER NOT NULL,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS training_queue (
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
CREATE INDEX IF NOT EXISTS training_queue_status ON training_queue (status);
CREATE TABLE IF NOT EXISTS training_job (
    id INTEGER PRIMARY KEY CHECK (id = 1),
    state TEXT NOT NULL,
    pid INTEGER,
    started_at TEXT,
    heartbeat_at TEXT,
    finished_at TEXT,
    workers INTEGER NOT NULL DEFAULT 0,
    limit_count INTEGER,
    done INTEGER NOT NULL DEFAULT 0,
    failed INTEGER NOT NULL DEFAULT 0,
    message TEXT NOT NULL DEFAULT '',
    stop_requested INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE IF NOT EXISTS training_attempt (
    attempt_id TEXT PRIMARY KEY,
    mol_id TEXT NOT NULL,
    cell_id TEXT NOT NULL,
    claimed_at TEXT NOT NULL,
    released_at TEXT,
    finished_at TEXT,
    outcome TEXT,
    note TEXT
);
CREATE INDEX IF NOT EXISTS training_attempt_mol ON training_attempt (mol_id);
CREATE TABLE IF NOT EXISTS training_failure (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mol_id TEXT NOT NULL,
    reason TEXT NOT NULL,
    at TEXT NOT NULL
);
"""

#: ``training_queue`` columns that arrived after the table held real data; a
#: ``CREATE TABLE IF NOT EXISTS`` cannot add them, so ``_migrate`` does.
_QUEUE_COLUMNS = (
    ("cell_id", "TEXT"),
    ("attempt_id", "TEXT"),
    ("lease_until", "TEXT"),
)


def database_path(store_root: Path) -> Path:
    """Location of the training database inside a Semi-Imperium store."""
    return store_root / DB_FILENAME


def connect(path: Path) -> sqlite3.Connection:
    """Open ``path`` (creating it and its schema); autocommit, WAL, 30 s wait."""
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path, timeout=30.0, isolation_level=None)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA busy_timeout=30000")
    connection.executescript(_SCHEMA)
    _migrate(connection)
    return connection


def _migrate(connection: sqlite3.Connection) -> None:
    """Add the missing ``training_queue`` columns, keeping every row.

    The check is repeated under the write lock: the dashboard and the
    coordinator may connect at the same moment, and the second ``ALTER``
    would otherwise fail with a duplicate column.
    """
    if not _missing_columns(connection):
        return
    connection.execute("BEGIN IMMEDIATE")
    try:
        for name, kind in _missing_columns(connection):
            connection.execute(f"ALTER TABLE training_queue ADD COLUMN {name} {kind}")
        connection.execute(
            "CREATE INDEX IF NOT EXISTS training_queue_cell"
            " ON training_queue (cell_id)"
        )
        connection.execute("COMMIT")
    except BaseException:
        connection.execute("ROLLBACK")
        raise


def _missing_columns(connection: sqlite3.Connection) -> list[tuple[str, str]]:
    present = {
        str(r["name"]) for r in connection.execute("PRAGMA table_info(training_queue)")
    }
    return [(name, kind) for name, kind in _QUEUE_COLUMNS if name not in present]
