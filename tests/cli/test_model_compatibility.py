"""Tests for model bundle metadata and method compatibility validation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import joblib
import pytest

from grimperium.calculation.methods import get_calculation_method, get_feature_schema
from grimperium.ml.features import FeaturePipeline
from grimperium.ml.persistence import load_model_metadata, save_model


def test_save_model_adds_default_compatibility_metadata(tmp_path: Path) -> None:
    model_path = tmp_path / "model.joblib"
    schema = get_feature_schema("grimperium_dhf_v1")

    save_model(
        {
            "learner": object(),
            "pipeline": object(),
            "metrics": {"train": {}, "test": {}},
        },
        model_path,
    )

    metadata = load_model_metadata(model_path)

    assert metadata["property_id"] == "standard_enthalpy_of_formation"
    assert metadata["baseline_hamiltonian"] == "PM7"
    assert metadata["feature_schema_id"] == schema.schema_id
    assert metadata["feature_schema_hash"] == schema.schema_hash
    assert tuple(metadata["feature_columns"]) == schema.columns


def _write_model_bundle(path: Path, **overrides: Any) -> None:
    schema = get_feature_schema("grimperium_dhf_v1")
    payload: dict[str, Any] = {
        "learner": object(),
        "pipeline": FeaturePipeline(list(schema.columns)),
        "metrics": {"train": {}, "test": {}},
        "version": "1.0.0",
        "trained_at": "2026-06-24T00:00:00+00:00",
        "property_id": "standard_enthalpy_of_formation",
        "baseline_hamiltonian": "PM7",
        "feature_schema_id": schema.schema_id,
        "feature_schema_hash": schema.schema_hash,
        "feature_columns": list(schema.columns),
    }
    payload.update(overrides)
    joblib.dump(payload, path)


def test_validate_model_for_method_rejects_baseline_mismatch(
    tmp_path: Path,
) -> None:
    from grimperium.cli.model_compatibility import (
        ModelCompatibilityError,
        validate_model_for_method,
    )

    model_path = tmp_path / "model.joblib"
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )
    _write_model_bundle(model_path, baseline_hamiltonian="AM1")

    with pytest.raises(ModelCompatibilityError, match="baseline_hamiltonian"):
        validate_model_for_method(model_path, method)


def test_validate_model_for_method_accepts_legacy_bundle_with_matching_columns(
    tmp_path: Path,
) -> None:
    from grimperium.cli.model_compatibility import validate_model_for_method

    model_path = tmp_path / "legacy.joblib"
    schema = get_feature_schema("grimperium_dhf_v1")
    joblib.dump(
        {
            "learner": object(),
            "pipeline": FeaturePipeline(list(schema.columns)),
            "metrics": {"train": {}, "test": {}},
            "version": "legacy",
            "trained_at": "unknown",
        },
        model_path,
    )
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )

    validate_model_for_method(model_path, method)


def test_validate_model_for_method_rejects_legacy_bundle_column_mismatch(
    tmp_path: Path,
) -> None:
    from grimperium.cli.model_compatibility import (
        ModelCompatibilityError,
        validate_model_for_method,
    )

    model_path = tmp_path / "legacy.joblib"
    schema = get_feature_schema("grimperium_dhf_v1")
    bad_columns = list(schema.columns)
    bad_columns[-1] = "wrong_baseline"
    joblib.dump(
        {
            "learner": object(),
            "pipeline": FeaturePipeline(bad_columns),
            "metrics": {"train": {}, "test": {}},
            "version": "legacy",
            "trained_at": "unknown",
        },
        model_path,
    )
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )

    with pytest.raises(ModelCompatibilityError, match="legacy model feature columns"):
        validate_model_for_method(model_path, method)
