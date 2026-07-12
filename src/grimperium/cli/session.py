"""Typed CLI session context for Grimperium."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from grimperium.calculation.methods import CalculationMethodDefinition


class ModelState(str, Enum):
    """Compatibility state for the session-selected model."""

    NOT_SELECTED = "not_selected"
    NOT_REQUIRED = "not_required"
    READY = "ready"
    MISSING = "missing"
    INCOMPATIBLE = "incompatible"


@dataclass
class ModelRef:
    """Reference to a model selected for the current CLI session."""

    name: str | None = None
    path: Path | None = None
    state: ModelState = ModelState.NOT_SELECTED


@dataclass(frozen=True)
class DatasetRef:
    """Reference to a catalogued scientific dataset."""

    database_id: str
    alias: str
    name: str
    path: Path
    role: str
    capabilities: frozenset[str] = field(default_factory=frozenset)


@dataclass(frozen=True)
class AnalysisSourceRef:
    """Explicit ad-hoc source for analysis views outside the data catalog."""

    source_type: str
    name: str
    path: Path


@dataclass(frozen=True)
class RunRef:
    """Placeholder for a future persisted run reference."""

    run_id: str | None = None


@dataclass
class MethodExecutionOverrides:
    """Future method execution override schema.

    The schema is intentionally not wired into calculation execution yet.
    """

    n_conformers: int | None = None
    crest_timeout_minutes: float | None = None
    mopac_timeout_minutes: float | None = None
    xtb_timeout_seconds: float | None = None
    batch_size: int | None = None


@dataclass
class SessionContext:
    """Typed scientific context shared across CLI views."""

    property_id: str | None = None
    method_definition: CalculationMethodDefinition | None = None
    overrides: MethodExecutionOverrides = field(
        default_factory=MethodExecutionOverrides
    )
    dataset: DatasetRef | None = None
    model: ModelRef = field(default_factory=ModelRef)
    run: RunRef = field(default_factory=RunRef)
    analysis_source: AnalysisSourceRef | None = None

    @property
    def analysis_path(self) -> Path | None:
        """Return the active analysis CSV path from the session source of truth."""
        if self.analysis_source is not None:
            return self.analysis_source.path
        if self.dataset is not None:
            return self.dataset.path
        return None
