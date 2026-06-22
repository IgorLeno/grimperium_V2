"""Enums for canonical calculation contracts."""

from enum import Enum


class StageExecutionStatus(str, Enum):
    """Execution status for an individual calculation stage."""

    NOT_REQUESTED = "not_requested"
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    SKIPPED = "skipped"
    UNAVAILABLE = "unavailable"
    FAILED = "failed"
    TIMEOUT = "timeout"


class OverallStatus(str, Enum):
    """Overall scientific result status."""

    SUCCESS = "success"
    PARTIAL = "partial"
    FAILED = "failed"


class PropertyRole(str, Enum):
    """Role of a property estimate in a calculation method."""

    BASELINE = "baseline"
    CORRECTION = "correction"
    FINAL = "final"


class ArtifactType(str, Enum):
    """Known calculation artifact types."""

    CREST_XYZ = "crest_xyz"
    MOPAC_INPUT = "mopac_input"
    MOPAC_OUT = "mopac_out"
    OTHER = "other"
