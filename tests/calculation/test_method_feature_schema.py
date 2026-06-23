"""Tests for calculation method feature schema catalog."""

from __future__ import annotations

import pytest

from grimperium.calculation.methods.feature_schema import (
    FEATURE_SCHEMA_HASH,
    FEATURE_SCHEMA_ID,
    calculate_feature_schema_hash,
    get_feature_schema,
    validate_feature_schema,
)
from grimperium.ml.features import FEATURE_COLUMNS


def test_grimperium_dhf_v1_matches_feature_columns() -> None:
    schema = get_feature_schema(FEATURE_SCHEMA_ID)

    assert schema.schema_id == "grimperium_dhf_v1"
    assert schema.columns == tuple(FEATURE_COLUMNS)
    assert schema.schema_hash == FEATURE_SCHEMA_HASH
    assert schema.schema_hash == calculate_feature_schema_hash(FEATURE_COLUMNS)


def test_validate_feature_schema_rejects_unknown_schema() -> None:
    with pytest.raises(ValueError, match="Unknown feature schema"):
        validate_feature_schema("unknown_schema", FEATURE_COLUMNS)


def test_validate_feature_schema_rejects_column_order_mismatch() -> None:
    shuffled = list(FEATURE_COLUMNS)
    shuffled[0], shuffled[1] = shuffled[1], shuffled[0]

    with pytest.raises(ValueError, match="column order"):
        validate_feature_schema(FEATURE_SCHEMA_ID, shuffled)


def test_validate_feature_schema_rejects_hash_mismatch() -> None:
    with pytest.raises(ValueError, match="hash"):
        validate_feature_schema(FEATURE_SCHEMA_ID, FEATURE_COLUMNS, "not-the-hash")
