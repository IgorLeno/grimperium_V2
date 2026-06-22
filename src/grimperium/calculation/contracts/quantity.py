"""Quantity value object for canonical calculation contracts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

KJ_PER_KCAL = 4.184


@dataclass(kw_only=True)
class Quantity:
    """Numeric value with unit metadata."""

    value: float
    unit: str

    def to_kcal_mol(self) -> float:
        """Return this quantity in kcal/mol."""
        if self.unit == "kcal/mol":
            return self.value
        if self.unit == "kJ/mol":
            return self.value / KJ_PER_KCAL
        raise ValueError(f"Unsupported energy unit: {self.unit}")

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-safe dictionary."""
        return {"value": self.value, "unit": self.unit}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Quantity:
        """Deserialize from dictionary."""
        return cls(value=float(data["value"]), unit=str(data["unit"]))
