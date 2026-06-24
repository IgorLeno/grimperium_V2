"""Model bundle compatibility checks for calculation methods."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from grimperium.calculation.methods import (
    CalculationMethodDefinition,
    validate_feature_schema,
)
from grimperium.ml.persistence import load_model, load_model_metadata

_REQUIRED_METADATA_KEYS = (
    "property_id",
    "baseline_hamiltonian",
    "feature_schema_id",
    "feature_schema_hash",
    "feature_columns",
)


class ModelCompatibilityError(ValueError):
    """Raised when a model bundle cannot satisfy a method requirement."""


def _has_complete_metadata(metadata: dict[str, Any]) -> bool:
    return all(metadata.get(key) is not None for key in _REQUIRED_METADATA_KEYS)


def _validate_metadata(
    metadata: dict[str, Any],
    method: CalculationMethodDefinition,
) -> None:
    compatibility = method.compatibility
    if metadata.get("property_id") != compatibility.property:
        raise ModelCompatibilityError(
            "Model property_id does not match method compatibility"
        )

    if metadata.get("baseline_hamiltonian") != compatibility.baseline_hamiltonian:
        raise ModelCompatibilityError(
            "Model baseline_hamiltonian does not match method compatibility"
        )

    if compatibility.feature_schema is None:
        return

    if metadata.get("feature_schema_id") != compatibility.feature_schema:
        raise ModelCompatibilityError(
            "Model feature_schema_id does not match method compatibility"
        )

    feature_columns = metadata.get("feature_columns")
    if not isinstance(feature_columns, list):
        raise ModelCompatibilityError("Model feature_columns metadata is missing")

    try:
        validate_feature_schema(
            compatibility.feature_schema,
            feature_columns,
            metadata.get("feature_schema_hash"),
        )
    except ValueError as exc:
        raise ModelCompatibilityError(str(exc)) from exc


def _validate_legacy_columns(
    model_path: Path,
    method: CalculationMethodDefinition,
) -> None:
    if method.compatibility.feature_schema is None:
        return

    try:
        _learner, feature_pipeline = load_model(model_path)
    except Exception as exc:
        raise ModelCompatibilityError(f"Could not inspect legacy model: {exc}") from exc

    feature_columns = getattr(feature_pipeline, "feature_cols", None)
    if not isinstance(feature_columns, list):
        raise ModelCompatibilityError("Legacy model feature columns are unavailable")

    try:
        validate_feature_schema(method.compatibility.feature_schema, feature_columns)
    except ValueError as exc:
        raise ModelCompatibilityError(
            f"legacy model feature columns do not match method: {exc}"
        ) from exc


def validate_model_for_method(
    model_path: Path,
    method: CalculationMethodDefinition,
) -> None:
    """Validate a model bundle against a calculation method definition."""
    if not method.model_requirement.model_required:
        return

    if not model_path.exists():
        raise ModelCompatibilityError(f"Model path does not exist: {model_path}")

    try:
        metadata = load_model_metadata(model_path)
    except Exception as exc:
        raise ModelCompatibilityError(f"Could not load model metadata: {exc}") from exc

    if _has_complete_metadata(metadata):
        _validate_metadata(metadata, method)
        return

    _validate_legacy_columns(model_path, method)
