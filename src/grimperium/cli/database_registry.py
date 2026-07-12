"""Dynamic database registry for the GRIMPERIUM CLI."""

from __future__ import annotations

import json
import logging
import os
import shutil
import uuid
from dataclasses import dataclass, field
from datetime import date
from importlib import resources
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

SCHEMA_VERSION = 2
OFFICIAL_DATABASE_PACKAGE = "grimperium.resources.databases"
USER_ID_PREFIX = "user."
CAPABILITIES = frozenset(
    {
        "readable",
        "writable",
        "batch_input",
        "analysis_input",
        "model_training",
        "reference_values",
    }
)


@dataclass
class DatabaseInfo:
    """Catalogued database plus computed availability."""

    schema_version: int
    database_id: str
    origin: str
    manifest_version: str
    role: str
    capabilities: frozenset[str]
    path: Path
    metadata: dict[str, Any] = field(default_factory=dict)
    status: str = "missing"
    molecules: int = 0
    last_updated: date | None = None
    properties: list[str] = field(default_factory=list)

    @property
    def name(self) -> str:
        return str(self.metadata.get("name", self.database_id))

    @property
    def alias(self) -> str:
        return str(self.metadata.get("alias", self.database_id))

    @property
    def description(self) -> str:
        return str(self.metadata.get("description", ""))

    @property
    def pipeline(self) -> str:
        return self.role

    @property
    def csv_path(self) -> str:
        return "" if self.path == Path("") else str(self.path)

    @property
    def created_at(self) -> str:
        return str(self.metadata.get("created_at", ""))

    @property
    def reference_db(self) -> str | None:
        value = self.metadata.get("reference_db")
        return str(value) if value is not None else None

    @property
    def max_conformers(self) -> int:
        value = self.metadata.get("max_conformers", 10)
        if isinstance(value, int):
            return value
        if isinstance(value, float):
            return int(value)
        if isinstance(value, str) and value.isdecimal():
            return int(value)
        return 10

    def has_capability(self, capability: str) -> bool:
        """Return whether this database supports a catalog action."""
        return capability in self.capabilities


@dataclass(frozen=True)
class RegistryScanResult:
    """CSV discovery result from the data directory."""

    known: list[Path]
    missing: list[DatabaseInfo]
    unregistered: list[Path]


@dataclass
class _RegistryPayload:
    schema_version: int
    entries: list[dict[str, Any]]
    overrides: dict[str, dict[str, Any]]


class DatabaseRegistry:
    """Loads official manifests and a user overlay catalog."""

    REGISTRY_FILENAME = "databases_registry.json"

    def __init__(self, data_dir: Path, config_dir: Path | None = None) -> None:
        self.data_dir = data_dir
        self.config_dir = config_dir or self._default_config_dir()
        self.overlay_path = self.config_dir / self.REGISTRY_FILENAME
        self._legacy_registry_path = data_dir / self.REGISTRY_FILENAME
        self._cache: list[DatabaseInfo] | None = None

    def load(self) -> list[DatabaseInfo]:
        """Load official manifests plus overlay entries."""
        if self._cache is not None:
            return self._cache

        self._migrate_legacy_registry()
        official_entries = self._load_official_entries()
        payload = self._read_overlay()
        merged = self._apply_overrides(official_entries, payload.overrides)
        merged.extend(payload.entries)
        entries = [self._dict_to_info(raw) for raw in merged]
        self._cache = [self._enrich_entry(entry) for entry in entries]
        return self._cache

    def reload(self) -> list[DatabaseInfo]:
        """Force a re-read from disk."""
        self._cache = None
        return self.load()

    def save(self, entries: list[DatabaseInfo]) -> None:
        """Persist user entries and official overrides to the overlay."""
        payload = _RegistryPayload(
            schema_version=SCHEMA_VERSION, entries=[], overrides={}
        )
        for entry in entries:
            raw = self._info_to_dict(entry)
            if entry.origin == "official":
                payload.overrides[entry.database_id] = self._official_override(raw)
            else:
                payload.entries.append(raw)
        self._write_overlay(payload)
        self._cache = None

    def add_entry(self, entry: DatabaseInfo) -> None:
        """Add a user database entry to the overlay."""
        if entry.origin == "official":
            raise ValueError("Official database manifests cannot be user entries")
        payload = self._read_overlay()
        payload.entries.append(self._info_to_dict(entry))
        self._write_overlay(payload)
        self._cache = None

    def add_user_database(
        self,
        *,
        path: Path,
        name: str,
        alias: str,
        description: str,
        role: str,
        capabilities: set[str] | frozenset[str],
    ) -> DatabaseInfo:
        """Create and persist a user database with a stable user UUID."""
        self._validate_capabilities(capabilities)
        database_id = f"{USER_ID_PREFIX}{uuid.uuid4()}"
        entry = DatabaseInfo(
            schema_version=SCHEMA_VERSION,
            database_id=database_id,
            origin="user",
            manifest_version="1.0.0",
            role=role,
            capabilities=frozenset(capabilities),
            path=path,
            metadata={
                "name": name,
                "alias": alias,
                "description": description,
            },
        )
        self.add_entry(entry)
        return self.get_by_id(database_id) or entry

    def update_entry(
        self,
        database_id: str,
        *,
        path: Path | None = None,
        metadata: dict[str, Any] | None = None,
        role: str | None = None,
        capabilities: set[str] | frozenset[str] | None = None,
    ) -> None:
        """Update a user entry or store limited official overrides."""
        current = self.get_by_id(database_id)
        if current is None:
            raise ValueError(f"Unknown database: {database_id}")
        if capabilities is not None:
            self._validate_capabilities(capabilities)

        payload = self._read_overlay()
        if current.origin == "official":
            override = payload.overrides.get(database_id, {})
            if path is not None:
                override["path"] = self._serialize_path(path)
            if metadata:
                override_metadata = dict(override.get("metadata", {}))
                override_metadata.update(metadata)
                override["metadata"] = override_metadata
            payload.overrides[database_id] = override
        else:
            payload.entries = [
                self._updated_user_entry(
                    raw, database_id, path, metadata, role, capabilities
                )
                for raw in payload.entries
            ]
        self._write_overlay(payload)
        self._cache = None

    def remove_from_catalog(self, database_id: str) -> None:
        """Remove an overlay entry or official override without deleting files."""
        payload = self._read_overlay()
        payload.entries = [
            raw for raw in payload.entries if raw.get("database_id") != database_id
        ]
        payload.overrides.pop(database_id, None)
        self._write_overlay(payload)
        self._cache = None

    def get_by_alias(self, alias: str) -> DatabaseInfo | None:
        """Look up a database by short alias, case-insensitively."""
        alias_lower = alias.lower()
        for entry in self.load():
            if entry.alias.lower() == alias_lower:
                return entry
        return None

    def get_by_id(self, database_id: str) -> DatabaseInfo | None:
        """Look up a database by stable ID."""
        for entry in self.load():
            if entry.database_id == database_id:
                return entry
        return None

    def scan_data_dir(self) -> RegistryScanResult:
        """Detect known, missing and unregistered CSV files."""
        entries = self.load()
        known_paths = {
            entry.path.resolve() for entry in entries if self._path_is_set(entry.path)
        }
        known = [
            entry.path
            for entry in entries
            if self._path_is_set(entry.path) and entry.path.exists()
        ]
        missing = [
            entry
            for entry in entries
            if self._path_is_set(entry.path) and not entry.path.exists()
        ]
        csv_files = (
            sorted(self.data_dir.glob("*.csv")) if self.data_dir.exists() else []
        )
        unregistered = [path for path in csv_files if path.resolve() not in known_paths]
        return RegistryScanResult(
            known=known, missing=missing, unregistered=unregistered
        )

    @staticmethod
    def _default_config_dir() -> Path:
        env_dir = os.environ.get("GRIMPERIUM_CONFIG_DIR")
        if env_dir:
            return Path(env_dir)
        return Path.home() / ".grimperium"

    def _load_official_entries(self) -> list[dict[str, Any]]:
        package_files = resources.files(OFFICIAL_DATABASE_PACKAGE)
        entries: list[dict[str, Any]] = []
        for resource in sorted(package_files.iterdir(), key=lambda item: item.name):
            if resource.name.endswith(".json"):
                raw = json.loads(resource.read_text(encoding="utf-8"))
                self._validate_manifest(raw, source=resource.name)
                entries.append(raw)
        return entries

    def _read_overlay(self) -> _RegistryPayload:
        if not self.overlay_path.exists():
            return _RegistryPayload(
                schema_version=SCHEMA_VERSION,
                entries=[],
                overrides={},
            )
        try:
            data = json.loads(self.overlay_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Invalid database overlay JSON: {self.overlay_path}"
            ) from exc

        if isinstance(data, list):
            return _RegistryPayload(
                schema_version=SCHEMA_VERSION,
                entries=[self._legacy_to_v2(raw) for raw in data],
                overrides={},
            )
        if not isinstance(data, dict):
            raise ValueError("Database overlay must be a mapping")
        entries = data.get("entries", [])
        overrides = data.get("overrides", {})
        if not isinstance(entries, list) or not isinstance(overrides, dict):
            raise ValueError("Database overlay has invalid entries/overrides shape")
        return _RegistryPayload(
            schema_version=int(data.get("schema_version", SCHEMA_VERSION)),
            entries=[self._normalize_entry(raw) for raw in entries],
            overrides={
                str(key): value
                for key, value in overrides.items()
                if isinstance(value, dict)
            },
        )

    def _write_overlay(self, payload: _RegistryPayload) -> None:
        self.overlay_path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "schema_version": payload.schema_version,
            "entries": payload.entries,
            "overrides": payload.overrides,
        }
        self.overlay_path.write_text(
            json.dumps(data, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

    def _migrate_legacy_registry(self) -> None:
        if not self._legacy_registry_path.exists():
            return
        raw_entries = self._read_legacy_entries(self._legacy_registry_path)
        payload = self._read_overlay()
        for raw in raw_entries:
            mapped_id = self._legacy_alias_to_official_id(str(raw.get("alias", "")))
            if mapped_id is None:
                payload.entries.append(self._legacy_to_v2(raw))
                continue
            override = payload.overrides.get(mapped_id, {})
            path = raw.get("csv_path")
            if isinstance(path, str) and path:
                override["path"] = path
            metadata = dict(override.get("metadata", {}))
            for key in ("name", "alias", "description", "properties", "created_at"):
                if key in raw:
                    metadata[key] = raw[key]
            override["metadata"] = metadata
            payload.overrides[mapped_id] = override
        self._write_overlay(payload)
        backup_path = self._legacy_registry_path.with_suffix(
            self._legacy_registry_path.suffix + ".bak"
        )
        if not backup_path.exists():
            shutil.copy2(self._legacy_registry_path, backup_path)
        self._legacy_registry_path.unlink()

    @staticmethod
    def _read_legacy_entries(path: Path) -> list[dict[str, Any]]:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return []
        if isinstance(data, list):
            return [raw for raw in data if isinstance(raw, dict)]
        return []

    @staticmethod
    def _legacy_alias_to_official_id(alias: str) -> str | None:
        return {
            "cbs": "official.cbs_chon",
            "pm7": "official.crest_pm7",
            "nist": "official.nist_experimental",
        }.get(alias.lower())

    def _apply_overrides(
        self,
        official_entries: list[dict[str, Any]],
        overrides: dict[str, dict[str, Any]],
    ) -> list[dict[str, Any]]:
        merged: list[dict[str, Any]] = []
        for raw in official_entries:
            database_id = str(raw["database_id"])
            override = overrides.get(database_id, {})
            merged.append(self._merge_entry(raw, override))
        return merged

    @staticmethod
    def _merge_entry(raw: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
        merged = dict(raw)
        if "path" in override:
            merged["path"] = override["path"]
        metadata = dict(raw.get("metadata", {}))
        override_metadata = override.get("metadata", {})
        if isinstance(override_metadata, dict):
            metadata.update(override_metadata)
        merged["metadata"] = metadata
        return merged

    def _dict_to_info(self, raw: dict[str, Any]) -> DatabaseInfo:
        normalized = self._normalize_entry(raw)
        self._validate_manifest(normalized, source=str(normalized.get("database_id")))
        metadata = dict(normalized.get("metadata", {}))
        properties = metadata.get("properties", [])
        return DatabaseInfo(
            schema_version=int(normalized["schema_version"]),
            database_id=str(normalized["database_id"]),
            origin=str(normalized["origin"]),
            manifest_version=str(normalized["manifest_version"]),
            role=str(normalized["role"]),
            capabilities=frozenset(normalized["capabilities"]),
            path=self._resolve_path(str(normalized.get("path", ""))),
            metadata=metadata,
            properties=(
                [str(item) for item in properties]
                if isinstance(properties, list)
                else []
            ),
        )

    def _info_to_dict(self, info: DatabaseInfo) -> dict[str, Any]:
        return {
            "schema_version": info.schema_version,
            "database_id": info.database_id,
            "origin": info.origin,
            "manifest_version": info.manifest_version,
            "role": info.role,
            "capabilities": sorted(info.capabilities),
            "path": self._serialize_path(info.path),
            "metadata": dict(info.metadata),
        }

    @staticmethod
    def _official_override(raw: dict[str, Any]) -> dict[str, Any]:
        return {
            "path": raw.get("path", ""),
            "metadata": raw.get("metadata", {}),
        }

    def _updated_user_entry(
        self,
        raw: dict[str, Any],
        database_id: str,
        path: Path | None,
        metadata: dict[str, Any] | None,
        role: str | None,
        capabilities: set[str] | frozenset[str] | None,
    ) -> dict[str, Any]:
        if raw.get("database_id") != database_id:
            return raw
        updated = dict(raw)
        if path is not None:
            updated["path"] = self._serialize_path(path)
        if metadata:
            updated_metadata = dict(updated.get("metadata", {}))
            updated_metadata.update(metadata)
            updated["metadata"] = updated_metadata
        if role is not None:
            updated["role"] = role
        if capabilities is not None:
            updated["capabilities"] = sorted(capabilities)
        return updated

    def _normalize_entry(self, raw: dict[str, Any]) -> dict[str, Any]:
        if "database_id" not in raw:
            return self._legacy_to_v2(raw)
        normalized = dict(raw)
        normalized.setdefault("schema_version", SCHEMA_VERSION)
        normalized.setdefault("origin", "user")
        normalized.setdefault("manifest_version", "1.0.0")
        normalized.setdefault("metadata", {})
        normalized.setdefault("capabilities", ["readable"])
        normalized.setdefault("path", "")
        return normalized

    @staticmethod
    def _legacy_to_v2(raw: dict[str, Any]) -> dict[str, Any]:
        alias = str(raw.get("alias") or raw.get("name") or "database")
        role = str(raw.get("pipeline", "reference"))
        capabilities = {"readable"}
        if role == "crest_pm7":
            capabilities.update({"writable", "batch_input", "analysis_input"})
        if role == "reference":
            capabilities.update({"analysis_input", "reference_values"})
        return {
            "schema_version": SCHEMA_VERSION,
            "database_id": f"{USER_ID_PREFIX}{uuid.uuid4()}",
            "origin": "user",
            "manifest_version": "1.0.0",
            "role": role,
            "capabilities": sorted(capabilities),
            "path": raw.get("csv_path", ""),
            "metadata": {
                "name": raw.get("name", alias),
                "alias": alias,
                "description": raw.get("description", ""),
                "properties": raw.get("properties", []),
                "created_at": raw.get("created_at", ""),
                "reference_db": raw.get("reference_db"),
                "max_conformers": raw.get("max_conformers", 10),
            },
        }

    def _resolve_path(self, raw_path: str) -> Path:
        if not raw_path:
            return Path("")
        path = Path(raw_path)
        if path.is_absolute():
            return path
        return self.data_dir / path

    def _serialize_path(self, path: Path) -> str:
        if not self._path_is_set(path):
            return ""
        try:
            return str(path.relative_to(self.data_dir))
        except ValueError:
            return str(path)

    def _enrich_entry(self, entry: DatabaseInfo) -> DatabaseInfo:
        if not self._path_is_set(entry.path):
            entry.status = "missing"
            return entry
        if not entry.path.exists():
            entry.status = "missing"
            return entry
        try:
            header = pd.read_csv(entry.path, nrows=0)
        except (pd.errors.EmptyDataError, OSError, UnicodeError) as exc:
            logger.debug("Failed to read database header %s: %s", entry.path, exc)
            entry.status = "unreadable"
            return entry
        if not self._has_valid_csv_schema(header):
            entry.status = "invalid_schema"
            return entry
        try:
            df = pd.read_csv(entry.path, low_memory=False)
        except (pd.errors.EmptyDataError, OSError, UnicodeError, ValueError) as exc:
            logger.debug("Failed to enrich database %s: %s", entry.path, exc)
            entry.status = "unreadable"
            return entry
        if entry.role == "crest_pm7" and "status" in df.columns:
            entry.molecules = int((df["status"] == "OK").sum())
        else:
            entry.molecules = len(df)
        try:
            entry.last_updated = date.fromtimestamp(entry.path.stat().st_mtime)
        except OSError:
            pass
        entry.status = "available"
        return entry

    @staticmethod
    def _has_valid_csv_schema(df: pd.DataFrame) -> bool:
        columns = set(df.columns)
        return bool(columns) and ("smiles" in columns or "mol_id" in columns)

    @staticmethod
    def _path_is_set(path: Path) -> bool:
        return path != Path("")

    @staticmethod
    def _validate_capabilities(capabilities: set[str] | frozenset[str]) -> None:
        unknown = set(capabilities) - CAPABILITIES
        if unknown:
            raise ValueError(f"Unknown database capabilities: {sorted(unknown)}")

    @staticmethod
    def _validate_manifest(raw: dict[str, Any], *, source: str) -> None:
        required = {
            "schema_version",
            "database_id",
            "origin",
            "manifest_version",
            "role",
            "capabilities",
            "path",
            "metadata",
        }
        missing = required - set(raw)
        if missing:
            raise ValueError(
                f"Database manifest {source} missing fields: {sorted(missing)}"
            )
        if raw["schema_version"] != SCHEMA_VERSION:
            raise ValueError(
                f"Database manifest {source} has unsupported schema_version"
            )
        capabilities = raw["capabilities"]
        if not isinstance(capabilities, list) or not all(
            isinstance(item, str) for item in capabilities
        ):
            raise ValueError(f"Database manifest {source} capabilities must be strings")
        unknown = set(capabilities) - CAPABILITIES
        if unknown:
            raise ValueError(
                f"Database manifest {source} has unknown capabilities: {sorted(unknown)}"
            )
