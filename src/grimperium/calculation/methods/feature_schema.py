"""Feature schema catalog for calculation method compatibility."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from hashlib import sha256

from grimperium.ml.features import FEATURE_COLUMNS

FEATURE_SCHEMA_ID = "grimperium_dhf_v1"
FEATURE_SCHEMA_HASH = "623d64b75be5f6e68930e432745659abcf0bb5e23633586b58a5a1a9a95d9c56"


@dataclass(frozen=True)
class FeatureSchemaDefinition:
    """Ordered feature schema definition used for model compatibility checks."""

    schema_id: str
    columns: tuple[str, ...]
    schema_hash: str


def calculate_feature_schema_hash(columns: Sequence[str]) -> str:
    """Calculate the stable hash for an ordered feature column list."""
    return sha256("\n".join(columns).encode("utf-8")).hexdigest()


_FEATURE_SCHEMAS: dict[str, FeatureSchemaDefinition] = {
    FEATURE_SCHEMA_ID: FeatureSchemaDefinition(
        schema_id=FEATURE_SCHEMA_ID,
        columns=tuple(FEATURE_COLUMNS),
        schema_hash=FEATURE_SCHEMA_HASH,
    )
}


def get_feature_schema(schema_id: str) -> FeatureSchemaDefinition:
    """Return a known feature schema definition by ID."""
    try:
        return _FEATURE_SCHEMAS[schema_id]
    except KeyError as exc:
        raise ValueError(f"Unknown feature schema: {schema_id}") from exc


def validate_feature_schema(
    schema_id: str,
    columns: Sequence[str],
    schema_hash: str | None = None,
) -> None:
    """Validate feature columns and optional hash against the schema catalog."""
    schema = get_feature_schema(schema_id)
    provided_columns = tuple(columns)
    if provided_columns != schema.columns:
        raise ValueError(
            f"Feature schema {schema_id} column order does not match catalog"
        )

    actual_hash = calculate_feature_schema_hash(provided_columns)
    if actual_hash != schema.schema_hash:
        raise ValueError(f"Feature schema {schema_id} hash does not match catalog")

    if schema_hash is not None and schema_hash != schema.schema_hash:
        raise ValueError(f"Feature schema {schema_id} hash mismatch: {schema_hash}")
