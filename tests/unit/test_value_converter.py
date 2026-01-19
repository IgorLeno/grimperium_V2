"""
Tests for MoleculeValueConverter.

Includes tests for Ajuste #1 (to_string allow_empty fix).
"""

from datetime import datetime
from enum import Enum

import pytest

from grimperium.core.value_converter import MoleculeValueConverter


class TestToFloat:
    """Tests for to_float conversion."""

    def test_valid_float(self):
        """Valid float returns (value, None)."""
        value, error = MoleculeValueConverter.to_float(5.3)
        assert value == 5.3
        assert error is None

    def test_valid_float_string(self):
        """Float as string is converted."""
        value, error = MoleculeValueConverter.to_float("5.3")
        assert value == 5.3
        assert error is None

    def test_valid_zero(self):
        """Zero is a valid float (not falsy!)."""
        value, error = MoleculeValueConverter.to_float(0.0)
        assert value == 0.0  # NOT None
        assert error is None

    def test_valid_zero_string(self):
        """Zero as string is valid."""
        value, error = MoleculeValueConverter.to_float("0.0")
        assert value == 0.0
        assert error is None

    def test_empty_string_returns_empty_error(self):
        """Empty string returns (None, 'empty')."""
        value, error = MoleculeValueConverter.to_float("")
        assert value is None
        assert error == "empty"

    def test_none_returns_empty_error(self):
        """None returns (None, 'empty')."""
        value, error = MoleculeValueConverter.to_float(None)
        assert value is None
        assert error == "empty"

    def test_nan_marker_returns_invalid_error(self):
        """'nan' string returns (None, 'invalid')."""
        value, error = MoleculeValueConverter.to_float("nan")
        assert value is None
        assert error == "invalid"

    def test_na_markers(self):
        """Various NA markers are invalid."""
        markers = ["nan", "none", "n/a", "na", "-", "?", "null"]
        for marker in markers:
            value, error = MoleculeValueConverter.to_float(marker)
            assert value is None, f"'{marker}' should return None"
            assert error == "invalid", f"'{marker}' should be invalid"

    def test_invalid_string(self):
        """Non-numeric string is invalid."""
        value, error = MoleculeValueConverter.to_float("not-a-number")
        assert value is None
        assert error == "invalid"

    def test_negative_float(self):
        """Negative float is valid."""
        value, error = MoleculeValueConverter.to_float(-123.456)
        assert value == -123.456
        assert error is None


class TestToInt:
    """Tests for to_int conversion."""

    def test_valid_int(self):
        """Valid int returns (value, None)."""
        value, error = MoleculeValueConverter.to_int(42)
        assert value == 42
        assert error is None

    def test_float_truncated(self):
        """Float is truncated to int."""
        value, error = MoleculeValueConverter.to_int(42.7)
        assert value == 42
        assert error is None

    def test_valid_zero(self):
        """Zero is valid."""
        value, error = MoleculeValueConverter.to_int(0)
        assert value == 0
        assert error is None

    def test_empty_returns_empty_error(self):
        """Empty returns (None, 'empty')."""
        value, error = MoleculeValueConverter.to_int("")
        assert value is None
        assert error == "empty"


class TestToStringAjuste1:
    """Tests for to_string - specifically for Ajuste #1 fix."""

    def test_valid_string(self):
        """Valid string returns (value, None)."""
        value, error = MoleculeValueConverter.to_string("hello")
        assert value == "hello"
        assert error is None

    def test_allow_empty_true_accepts_empty_string(self):
        """AJUSTE #1: allow_empty=True accepts empty string as valid."""
        value, error = MoleculeValueConverter.to_string("", allow_empty=True)
        assert value == ""  # ← Empty string, NOT None
        assert error is None  # ← No error

    def test_allow_empty_false_rejects_empty_string(self):
        """AJUSTE #1: allow_empty=False rejects empty string."""
        value, error = MoleculeValueConverter.to_string("", allow_empty=False)
        assert value is None
        assert error == "empty"

    def test_none_with_allow_empty_true(self):
        """AJUSTE #1: None with allow_empty=True becomes empty string."""
        value, error = MoleculeValueConverter.to_string(None, allow_empty=True)
        assert value == ""
        assert error is None

    def test_none_with_allow_empty_false(self):
        """AJUSTE #1: None with allow_empty=False is invalid."""
        value, error = MoleculeValueConverter.to_string(None, allow_empty=False)
        assert value is None
        assert error == "empty"

    def test_nan_marker_with_allow_empty_true(self):
        """AJUSTE #1: NA marker with allow_empty=True becomes empty."""
        value, error = MoleculeValueConverter.to_string("nan", allow_empty=True)
        assert value == ""  # ← Treat marker as empty
        assert error is None

    def test_nan_marker_with_allow_empty_false(self):
        """AJUSTE #1: NA marker with allow_empty=False is invalid."""
        value, error = MoleculeValueConverter.to_string("nan", allow_empty=False)
        assert value is None
        assert error == "invalid"

    def test_whitespace_stripped(self):
        """Whitespace is stripped from strings."""
        value, error = MoleculeValueConverter.to_string("  hello  ")
        assert value == "hello"
        assert error is None

    def test_default_allow_empty_is_false(self):
        """Default allow_empty is False."""
        value, error = MoleculeValueConverter.to_string("")
        assert value is None
        assert error == "empty"


class TestToEnum:
    """Tests for to_enum conversion."""

    class TestStatus(str, Enum):
        PENDING = "pending"
        COMPLETE = "complete"

    def test_valid_enum(self):
        """Valid enum value is converted."""
        value, error = MoleculeValueConverter.to_enum(
            "pending", self.TestStatus
        )
        assert value == self.TestStatus.PENDING
        assert error is None

    def test_case_insensitive(self):
        """Enum matching is case-insensitive."""
        value, error = MoleculeValueConverter.to_enum(
            "PENDING", self.TestStatus
        )
        assert value == self.TestStatus.PENDING
        assert error is None

    def test_invalid_value(self):
        """Invalid enum value returns error."""
        value, error = MoleculeValueConverter.to_enum(
            "invalid", self.TestStatus
        )
        assert value is None
        assert error == "invalid"

    def test_empty_returns_empty_error(self):
        """Empty value returns empty error."""
        value, error = MoleculeValueConverter.to_enum("", self.TestStatus)
        assert value is None
        assert error == "empty"


class TestToDatetime:
    """Tests for to_datetime conversion."""

    def test_valid_iso8601(self):
        """Valid ISO 8601 is parsed."""
        value, error = MoleculeValueConverter.to_datetime("2026-01-19T23:55:00")
        assert value == datetime(2026, 1, 19, 23, 55, 0)
        assert error is None

    def test_date_only(self):
        """Date-only ISO 8601 is parsed."""
        value, error = MoleculeValueConverter.to_datetime("2026-01-19")
        assert value == datetime(2026, 1, 19, 0, 0, 0)
        assert error is None

    def test_empty_returns_empty_error(self):
        """Empty returns empty error."""
        value, error = MoleculeValueConverter.to_datetime("")
        assert value is None
        assert error == "empty"

    def test_invalid_format(self):
        """Invalid format returns invalid error."""
        value, error = MoleculeValueConverter.to_datetime("not-a-date")
        assert value is None
        assert error == "invalid"


class TestToBool:
    """Tests for to_bool conversion."""

    def test_true_values(self):
        """Truthy values return True."""
        for val in ["true", "True", "TRUE", "yes", "1", "on", True, 1]:
            value, error = MoleculeValueConverter.to_bool(val)
            assert value is True, f"'{val}' should be True"
            assert error is None

    def test_false_values(self):
        """Falsy values return False."""
        for val in ["false", "False", "FALSE", "no", "0", "off", False, 0]:
            value, error = MoleculeValueConverter.to_bool(val)
            assert value is False, f"'{val}' should be False"
            assert error is None

    def test_empty_returns_empty_error(self):
        """Empty returns empty error."""
        value, error = MoleculeValueConverter.to_bool("")
        assert value is None
        assert error == "empty"

    def test_invalid_value(self):
        """Invalid value returns invalid error."""
        value, error = MoleculeValueConverter.to_bool("maybe")
        assert value is None
        assert error == "invalid"
