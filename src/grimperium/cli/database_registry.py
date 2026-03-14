"""
Database registry for GRIMPERIUM CLI.

Manages database definitions via a JSON registry file,
replacing hardcoded mock data with dynamic, enriched entries.
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class DatabaseInfo:
    """Represents a molecular database with metadata."""

    name: str
    alias: str
    csv_path: str  # Relative to data_dir
    description: str
    pipeline: str  # "reference" | "crest_pm7" | "pm7_only"
    created_at: str  # ISO date
    max_conformers: int = 10
    reference_db: str | None = None
    status: str = "ready"  # "ready" | "in_development" | "empty"
    # Computed at load time (not persisted):
    molecules: int = 0
    last_updated: date | None = None
    properties: list[str] = field(default_factory=list)


class DatabaseRegistry:
    """Manages database entries via a JSON registry file.

    Auto-creates ``databases_registry.json`` with default entries
    if the file does not exist. Enriches each entry at load time
    by counting OK molecules from the corresponding CSV.
    """

    REGISTRY_FILENAME = "databases_registry.json"

    def __init__(self, data_dir: Path) -> None:
        self.data_dir = data_dir
        self._registry_path = data_dir / self.REGISTRY_FILENAME
        self._cache: list[DatabaseInfo] | None = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load(self) -> list[DatabaseInfo]:
        """Read JSON registry, enrich entries, and return them.

        Creates the registry with defaults if the file is missing.
        Uses a cache so repeated calls within the same session are fast.
        """
        if self._cache is not None:
            return self._cache

        if not self._registry_path.exists():
            self._write_defaults()

        raw = self._read_json()
        entries = [self._dict_to_info(d) for d in raw]
        self._cache = [self._enrich_entry(e) for e in entries]
        return self._cache

    def save(self, entries: list[DatabaseInfo]) -> None:
        """Persist a list of DatabaseInfo back to the JSON file."""
        data = [self._info_to_dict(e) for e in entries]
        self._registry_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._registry_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        self._cache = None  # Invalidate cache

    def reload(self) -> list[DatabaseInfo]:
        """Force re-read from disk (invalidates cache)."""
        self._cache = None
        return self.load()

    def add_entry(self, entry: DatabaseInfo) -> None:
        """Add a new entry and persist. (For wizard — Session 2.)"""
        entries = self.load()
        entries.append(entry)
        self.save(entries)

    def get_by_alias(self, alias: str) -> DatabaseInfo | None:
        """Look up a database by its short alias (case-insensitive)."""
        alias_lower = alias.lower()
        for entry in self.load():
            if entry.alias.lower() == alias_lower:
                return entry
        return None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _write_defaults(self) -> None:
        """Create the registry file with default entries."""
        self._registry_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._registry_path, "w", encoding="utf-8") as f:
            json.dump(self._create_defaults(), f, indent=2, ensure_ascii=False)

    @staticmethod
    def _create_defaults() -> list[dict[str, Any]]:
        """Default entries: CBS, PM7, NIST."""
        return [
            {
                "name": "CBS Reference (CHON-only)",
                "alias": "CBS",
                "csv_path": "thermo_cbs_chon.csv",
                "description": (
                    "CHON-filtered CBS reference energies (C, H, O, N only). "
                    "Recommended for Phase A."
                ),
                "pipeline": "reference",
                "created_at": "2026-01-15",
                "max_conformers": 10,
                "reference_db": None,
                "status": "ready",
                "properties": [
                    "H298_cbs",
                    "H298_b3",
                    "smiles",
                    "charge",
                    "multiplicity",
                ],
            },
            {
                "name": "CREST PM7",
                "alias": "PM7",
                "csv_path": "thermo_pm7.csv",
                "description": "CREST conformer search with PM7 optimization",
                "pipeline": "crest_pm7",
                "created_at": "2026-01-01",
                "max_conformers": 10,
                "reference_db": "CBS",
                "status": "in_development",
                "properties": [
                    "H298_pm7",
                    "conformers",
                    "smiles",
                    "quality_grade",
                ],
            },
            {
                "name": "NIST Experimental",
                "alias": "NIST",
                "csv_path": "",
                "description": "NIST experimental thermochemistry data",
                "pipeline": "reference",
                "created_at": "2026-01-01",
                "max_conformers": 10,
                "reference_db": None,
                "status": "in_development",
                "properties": ["H298_exp", "uncertainty", "reference"],
            },
        ]

    def _dict_to_info(self, d: dict[str, Any]) -> DatabaseInfo:
        """Convert a raw JSON dict into a DatabaseInfo."""
        return DatabaseInfo(
            name=d["name"],
            alias=d["alias"],
            csv_path=d.get("csv_path", ""),
            description=d.get("description", ""),
            pipeline=d.get("pipeline", "reference"),
            created_at=d.get("created_at", ""),
            max_conformers=d.get("max_conformers", 10),
            reference_db=d.get("reference_db"),
            status=d.get("status", "in_development"),
            properties=d.get("properties", []),
        )

    @staticmethod
    def _info_to_dict(info: DatabaseInfo) -> dict[str, Any]:
        """Convert a DatabaseInfo to a JSON-serializable dict.

        Only persists fields that belong in the registry (not computed ones).
        """
        return {
            "name": info.name,
            "alias": info.alias,
            "csv_path": info.csv_path,
            "description": info.description,
            "pipeline": info.pipeline,
            "created_at": info.created_at,
            "max_conformers": info.max_conformers,
            "reference_db": info.reference_db,
            "status": info.status,
            "properties": info.properties,
        }

    def _enrich_entry(self, entry: DatabaseInfo) -> DatabaseInfo:
        """Enrich a DatabaseInfo with computed fields from the CSV.

        For entries with a CSV path, counts molecules and determines
        last_updated from file modification time. For crest_pm7 pipeline,
        counts only status=="OK" rows.
        """
        if not entry.csv_path:
            return entry

        csv_full = self.data_dir / entry.csv_path
        if not csv_full.exists():
            return entry

        try:
            df = pd.read_csv(csv_full, low_memory=False)

            if entry.pipeline == "crest_pm7" and "status" in df.columns:
                entry.molecules = int((df["status"] == "OK").sum())
            else:
                entry.molecules = len(df)

            # Determine status based on molecule count
            if entry.molecules > 0 and entry.status == "in_development":
                entry.status = "ready"

            # Last updated from file mtime
            try:
                mtime = os.path.getmtime(csv_full)
                entry.last_updated = date.fromtimestamp(mtime)
            except OSError:
                pass

        except (pd.errors.EmptyDataError, OSError, ValueError) as exc:
            logger.debug("Failed to enrich %s: %s", entry.alias, exc)

        return entry

    def _read_json(self) -> list[dict[str, Any]]:
        """Read the JSON registry file."""
        try:
            with open(self._registry_path, encoding="utf-8") as f:
                data = json.load(f)
            if isinstance(data, list):
                return data
            logger.warning("Registry is not a list, recreating defaults")
            self._write_defaults()
            return self._create_defaults()
        except (json.JSONDecodeError, OSError) as exc:
            logger.warning("Failed to read registry: %s — recreating", exc)
            self._write_defaults()
            return self._create_defaults()
