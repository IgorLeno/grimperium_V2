"""
CSV value conversion with clarity on type handling.

Distinguish:
- "field absent"   (column doesn't exist)
- "value empty"    (empty cell or None)
- "value invalid"  (corrupted data like "nan")
- "value zero"     (valid zero, not falsy!)

**AJUSTE v2.2:** Fixed to_string() logic
- allow_empty=True  → empty string is VALID ("", None)
- allow_empty=False → empty string is INVALID (None, "empty")
"""

from __future__ import annotations

from datetime import datetime
from enum import Enum
from typing import Any, Optional, Tuple, Type, TypeVar

import logging

import pandas as pd

logger = logging.getLogger(__name__)

T = TypeVar("T", bound=Enum)


class ConversionError(Exception):
    """Value conversion error."""

    pass


class MoleculeValueConverter:
    """
    Convert CSV values to Python types with clear semantics.

    Every conversion returns (value, error_reason):
    - value: Converted value (or None if failed)
    - error_reason: "empty" | "invalid" | None

    This allows calling code to decide: skip row, log warning, or fail.

    Examples:
        >>> MoleculeValueConverter.to_float("5.3")
        (5.3, None)

        >>> MoleculeValueConverter.to_float("")
        (None, "empty")

        >>> MoleculeValueConverter.to_float("nan")
        (None, "invalid")

        >>> MoleculeValueConverter.to_float("0.0")
        (0.0, None)  # Valid zero, NOT None!
    """

    # Special markers that mean "no data"
    NA_MARKERS: frozenset[str] = frozenset(
        {"nan", "none", "n/a", "na", "-", "?", "null"}
    )

    @staticmethod
    def to_float(
        value: Any,
        field_name: str = "unknown",
        allow_zero: bool = True,
    ) -> Tuple[Optional[float], Optional[str]]:
        """
        Convert to float with clear error semantics.

        Args:
            value: Value to convert
            field_name: Field name (for logging)
            allow_zero: If False, log debug when value is 0.0 (but still valid!)

        Returns:
            (float_value, error_reason) where:
            - (5.3, None) - valid float
            - (0.0, None) - valid zero (even if allow_zero=False)
            - (None, "empty") - missing value
            - (None, "invalid") - corrupted value

        Examples:
            >>> MoleculeValueConverter.to_float(5.3)
            (5.3, None)

            >>> MoleculeValueConverter.to_float("")
            (None, "empty")

            >>> MoleculeValueConverter.to_float("nan")
            (None, "invalid")

            >>> MoleculeValueConverter.to_float(0.0)
            (0.0, None)  # Zero is valid!
        """
        # Handle None and empty
        if value is None or value == "" or pd.isna(value):
            return None, "empty"

        # Handle string NA markers
        if isinstance(value, str):
            value_stripped = value.strip().lower()
            if value_stripped in MoleculeValueConverter.NA_MARKERS:
                return None, "invalid"

        # Try conversion
        try:
            result = float(value)

            # Log zero values (noteworthy but valid)
            if result == 0.0 and not allow_zero:
                logger.debug(f"Field '{field_name}' = 0.0 (valid but noteworthy)")

            return result, None

        except (ValueError, TypeError):
            return None, "invalid"

    @staticmethod
    def to_int(
        value: Any,
        field_name: str = "unknown",
    ) -> Tuple[Optional[int], Optional[str]]:
        """
        Convert to int with clear error semantics.

        Args:
            value: Value to convert
            field_name: Field name (for logging)

        Returns:
            (int_value, error_reason)

        Examples:
            >>> MoleculeValueConverter.to_int(42)
            (42, None)

            >>> MoleculeValueConverter.to_int("42")
            (42, None)

            >>> MoleculeValueConverter.to_int(42.7)
            (42, None)  # Truncates

            >>> MoleculeValueConverter.to_int("")
            (None, "empty")
        """
        # Try to convert to float first, then to int
        result, error = MoleculeValueConverter.to_float(
            value, field_name, allow_zero=True
        )

        if result is not None:
            return int(result), error

        return None, error

    @staticmethod
    def to_string(
        value: Any,
        field_name: str = "unknown",
        allow_empty: bool = False,
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Convert to string with clear error semantics.

        **AJUSTE v2.2:** Fixed allow_empty logic

        Args:
            value: Value to convert
            field_name: Field name (for logging)
            allow_empty:
                - True: empty string is VALID ("", None)
                - False: empty string is INVALID (None, "empty")

        Returns:
            (str_value, error_reason)

        Examples:
            >>> MoleculeValueConverter.to_string("hello")
            ("hello", None)

            >>> MoleculeValueConverter.to_string("", allow_empty=True)
            ("", None)  # ← Valid empty string

            >>> MoleculeValueConverter.to_string("", allow_empty=False)
            (None, "empty")  # ← Invalid, rejected

            >>> MoleculeValueConverter.to_string(None, allow_empty=True)
            ("", None)  # ← Convert None to empty

            >>> MoleculeValueConverter.to_string("nan", allow_empty=True)
            ("", None)  # ← Treat marker as empty

            >>> MoleculeValueConverter.to_string("nan", allow_empty=False)
            (None, "invalid")  # ← Reject marker
        """
        # Handle None
        if value is None or pd.isna(value):
            if allow_empty:
                return "", None  # ← Convert None to empty string
            else:
                return None, "empty"

        # Handle empty string (AJUSTE: Fixed logic)
        if value == "":
            if allow_empty:
                return "", None  # ← VALID: keep empty string
            else:
                return None, "empty"  # ← INVALID: reject empty

        # Convert to string and strip
        str_val = str(value).strip()

        # Check for NA markers
        if str_val.lower() in MoleculeValueConverter.NA_MARKERS:
            if allow_empty:
                return "", None  # ← Treat marker as empty (valid)
            else:
                return None, "invalid"  # ← Reject marker

        # Valid string
        return str_val, None

    @staticmethod
    def to_enum(
        value: Any,
        enum_class: Type[T],
        field_name: str = "unknown",
    ) -> Tuple[Optional[T], Optional[str]]:
        """
        Convert to enum with validation.

        Args:
            value: Value to convert
            enum_class: Enum class (e.g., MoleculeStatus)
            field_name: Field name (for logging)

        Returns:
            (enum_value, error_reason)

        Examples:
            >>> from enum import Enum
            >>> class Status(Enum):
            ...     PENDING = "pending"
            ...     COMPLETE = "complete"

            >>> MoleculeValueConverter.to_enum("pending", Status)
            (Status.PENDING, None)

            >>> MoleculeValueConverter.to_enum("PENDING", Status)
            (Status.PENDING, None)  # Case-insensitive

            >>> MoleculeValueConverter.to_enum("invalid", Status)
            (None, "invalid")
        """
        if value is None or value == "":
            return None, "empty"

        str_val = str(value).strip().lower()

        # Try to match enum value
        for member in enum_class:
            if member.value.lower() == str_val:
                return member, None

        # No match
        valid_values = [m.value for m in enum_class]
        logger.warning(f"Field '{field_name}': '{value}' not in {valid_values}")
        return None, "invalid"

    @staticmethod
    def to_datetime(
        value: Any,
        field_name: str = "unknown",
    ) -> Tuple[Optional[datetime], Optional[str]]:
        """
        Convert to datetime (ISO 8601) with validation.

        Args:
            value: Value to convert
            field_name: Field name (for logging)

        Returns:
            (datetime_value, error_reason)

        Examples:
            >>> MoleculeValueConverter.to_datetime("2026-01-19T23:55:00")
            (datetime(2026, 1, 19, 23, 55), None)

            >>> MoleculeValueConverter.to_datetime("2026-01-19")
            (datetime(2026, 1, 19, 0, 0), None)

            >>> MoleculeValueConverter.to_datetime("")
            (None, "empty")

            >>> MoleculeValueConverter.to_datetime("not-a-date")
            (None, "invalid")
        """
        if value is None or value == "":
            return None, "empty"

        if pd.isna(value):
            return None, "empty"

        str_val = str(value).strip()

        # Try ISO 8601 parsing
        try:
            dt = datetime.fromisoformat(str_val)
            return dt, None
        except ValueError:
            logger.warning(
                f"Field '{field_name}': '{value}' not ISO 8601. "
                f"Expected: 2026-01-19T23:55:00 or with timezone"
            )
            return None, "invalid"

    @staticmethod
    def to_bool(
        value: Any,
        field_name: str = "unknown",
    ) -> Tuple[Optional[bool], Optional[str]]:
        """
        Convert to boolean with validation.

        Args:
            value: Value to convert
            field_name: Field name (for logging)

        Returns:
            (bool_value, error_reason)

        Truthy: "true", "yes", "1", "on", True, 1
        Falsy: "false", "no", "0", "off", False, 0
        """
        if value is None or value == "":
            return None, "empty"

        if pd.isna(value):
            return None, "empty"

        # Direct bool
        if isinstance(value, bool):
            return value, None

        # Direct int (0 or 1)
        if isinstance(value, int):
            if value == 0:
                return False, None
            elif value == 1:
                return True, None
            else:
                return None, "invalid"

        # String conversion
        str_val = str(value).strip().lower()

        truthy = {"true", "yes", "1", "on"}
        falsy = {"false", "no", "0", "off"}

        if str_val in truthy:
            return True, None
        elif str_val in falsy:
            return False, None
        else:
            logger.warning(
                f"Field '{field_name}': '{value}' not a valid boolean. "
                f"Expected: true/false, yes/no, 1/0, on/off"
            )
            return None, "invalid"
