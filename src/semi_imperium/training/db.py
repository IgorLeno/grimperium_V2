"""SQLite plumbing shared by the training stores (one file, WAL mode).

Several processes touch this file at once (the coordinator, the dashboard API
and the CLI), so every call opens a short-lived connection instead of sharing
one, and the schema is created idempotently on connect.
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
CREATE TABLE IF NOT EXISTS training_failure (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mol_id TEXT NOT NULL,
    reason TEXT NOT NULL,
    at TEXT NOT NULL
);
"""


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
    return connection
