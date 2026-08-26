"""Persisted records: runs, individual calculations and their provenance.

The records here are frozen dataclasses. State changes produce a new
record through :meth:`CalculationRecord.transition_to`, which validates
the transition and the state/verification coherence before returning it,
so an inconsistent record cannot be constructed by accident.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import PurePosixPath, PureWindowsPath
from typing import Any

from semi_imperium.domain.configuration import (
    CalculationSignature,
    EffectiveConfiguration,
)
from semi_imperium.domain.enums import (
    ALLOWED_CALCULATION_TRANSITIONS,
    ALLOWED_RUN_TRANSITIONS,
    COHERENT_VERIFICATION_OUTCOMES,
    CalculationState,
    RunState,
    VerificationOutcome,
)
from semi_imperium.domain.hashing import stable_digest
from semi_imperium.domain.identity import MolecularIdentity

RECORD_SCHEMA_VERSION = 1

CALCULATION_ID_PREFIX = "calc"
CALCULATION_ID_LENGTH = 16


def utc_now() -> datetime:
    """Return the current time as a timezone-aware UTC timestamp."""
    return datetime.now(timezone.utc)


def _require_aware(value: datetime, label: str) -> datetime:
    """Reject naive datetimes, which cannot be compared across machines."""
    if value.tzinfo is None or value.tzinfo.utcoffset(value) is None:
        raise ValueError(f"{label} must be timezone-aware, got {value!r}")
    return value


def _optional_datetime(value: Any, label: str) -> datetime | None:
    """Parse an optional ISO-8601 timestamp, enforcing timezone awareness."""
    if value is None:
        return None
    return _require_aware(datetime.fromisoformat(str(value)), label)


def _iso(value: datetime | None) -> str | None:
    """Render an optional datetime as ISO-8601."""
    return None if value is None else value.isoformat()


@dataclass(frozen=True)
class Timestamps:
    """Creation, start and completion instants of one unit of work.

    ``started_at`` and ``completed_at`` stay ``None`` until the event they
    describe has actually happened. That is the one place where ``None``
    is meaningful rather than ambiguous, because the accompanying state
    already says whether the work has started or finished.
    """

    created_at: datetime = field(default_factory=utc_now)
    started_at: datetime | None = None
    completed_at: datetime | None = None

    def __post_init__(self) -> None:
        _require_aware(self.created_at, "Timestamps.created_at")
        if self.started_at is not None:
            _require_aware(self.started_at, "Timestamps.started_at")
        if self.completed_at is not None:
            _require_aware(self.completed_at, "Timestamps.completed_at")
        if self.completed_at is not None and self.started_at is None:
            raise ValueError(
                "Timestamps.completed_at requires started_at to be set as well"
            )
        if (
            self.started_at is not None
            and self.completed_at is not None
            and self.completed_at < self.started_at
        ):
            raise ValueError("Timestamps.completed_at must not precede started_at")

    @property
    def duration_seconds(self) -> float | None:
        """Wall-clock duration, or ``None`` while the work is unfinished."""
        if self.started_at is None or self.completed_at is None:
            return None
        return (self.completed_at - self.started_at).total_seconds()

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "created_at": self.created_at.isoformat(),
            "started_at": _iso(self.started_at),
            "completed_at": _iso(self.completed_at),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> Timestamps:
        """Deserialize from JSON-compatible primitives."""
        created = datetime.fromisoformat(str(payload["created_at"]))
        return cls(
            created_at=_require_aware(created, "Timestamps.created_at"),
            started_at=_optional_datetime(
                payload.get("started_at"), "Timestamps.started_at"
            ),
            completed_at=_optional_datetime(
                payload.get("completed_at"), "Timestamps.completed_at"
            ),
        )


@dataclass(frozen=True)
class ScientificProvenance:
    """Who produced a number, with which software, and under which method.

    This is what makes a stored value defensible months later. It records
    the acting software rather than any interpretation of the result.
    """

    method_id: str
    method_version: str
    property_id: str
    semi_imperium_version: str
    grimperium_version: str
    program_versions: dict[str, str] = field(default_factory=dict)
    """Version string per external program, e.g. ``{"crest": "3.0.1"}``."""
    hostname: str | None = None
    source: str = "semi_imperium"
    """Which component produced the record. Never blank."""
    notes: str | None = None

    def __post_init__(self) -> None:
        if not self.source.strip():
            raise ValueError("ScientificProvenance.source must not be empty")

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "method_id": self.method_id,
            "method_version": self.method_version,
            "property_id": self.property_id,
            "semi_imperium_version": self.semi_imperium_version,
            "grimperium_version": self.grimperium_version,
            "program_versions": dict(self.program_versions),
            "hostname": self.hostname,
            "source": self.source,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ScientificProvenance:
        """Deserialize from JSON-compatible primitives."""
        versions = payload.get("program_versions") or {}
        return cls(
            method_id=str(payload["method_id"]),
            method_version=str(payload["method_version"]),
            property_id=str(payload["property_id"]),
            semi_imperium_version=str(payload["semi_imperium_version"]),
            grimperium_version=str(payload["grimperium_version"]),
            program_versions={str(k): str(v) for k, v in versions.items()},
            hostname=_optional_text(payload.get("hostname")),
            source=str(payload["source"]),
            notes=_optional_text(payload.get("notes")),
        )


@dataclass(frozen=True)
class CalculationResultData:
    """The scientific payload of a finished calculation.

    Kept separate from the record so that lifecycle bookkeeping and the
    numbers themselves never get tangled: a ``FAILED`` record simply has
    no result attached.
    """

    energy_hof_kcal_mol: float | None = None
    conformer_index: int | None = None
    conformers_evaluated: int = 0
    lowest_frequency_cm1: float | None = None
    artifact_paths: tuple[str, ...] = ()
    """Paths relative to the store root; absolute paths are rejected."""

    def __post_init__(self) -> None:
        if self.conformers_evaluated < 0:
            raise ValueError(
                "CalculationResultData.conformers_evaluated must be >= 0, "
                f"got {self.conformers_evaluated}"
            )
        for path in self.artifact_paths:
            if PurePosixPath(path).is_absolute() or PureWindowsPath(path).is_absolute():
                raise ValueError(
                    "CalculationResultData.artifact_paths must be relative to the "
                    f"store root, got {path!r}"
                )

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "energy_hof_kcal_mol": self.energy_hof_kcal_mol,
            "conformer_index": self.conformer_index,
            "conformers_evaluated": self.conformers_evaluated,
            "lowest_frequency_cm1": self.lowest_frequency_cm1,
            "artifact_paths": list(self.artifact_paths),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> CalculationResultData:
        """Deserialize from JSON-compatible primitives."""
        return cls(
            energy_hof_kcal_mol=_optional_float(payload.get("energy_hof_kcal_mol")),
            conformer_index=_optional_int(payload.get("conformer_index")),
            conformers_evaluated=int(payload.get("conformers_evaluated", 0)),
            lowest_frequency_cm1=_optional_float(payload.get("lowest_frequency_cm1")),
            artifact_paths=tuple(
                str(item) for item in payload.get("artifact_paths") or ()
            ),
        )


@dataclass(frozen=True)
class CalculationRecord:
    """One molecule computed once under one effective configuration.

    Identity has three layers:

    * ``molecule.molecule_id`` — which species was computed;
    * ``signature.digest`` — under which effective configuration;
    * ``calculation_id`` — this individual attempt, inside ``run_id``.

    The first two form :attr:`reuse_key`, the only thing reuse consults.
    """

    run_id: str
    molecule: MolecularIdentity
    signature: CalculationSignature
    provenance: ScientificProvenance
    state: CalculationState = CalculationState.PENDING
    verification: VerificationOutcome = VerificationOutcome.NOT_REQUESTED
    timestamps: Timestamps = field(default_factory=Timestamps)
    result: CalculationResultData | None = None
    error_message: str | None = None
    calculation_id: str = ""
    """Left blank on construction; derived deterministically in ``__post_init__``."""
    schema_version: int = RECORD_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.run_id.strip():
            raise ValueError("CalculationRecord.run_id must not be empty")
        expected_id = build_calculation_id(
            run_id=self.run_id,
            molecule_id=self.molecule.molecule_id,
            signature_digest=self.signature.digest,
        )
        if not self.calculation_id:
            object.__setattr__(self, "calculation_id", expected_id)
        elif self.calculation_id != expected_id:
            raise ValueError(
                "CalculationRecord.calculation_id does not match its run, "
                f"molecule and signature: {self.calculation_id!r} != {expected_id!r}"
            )
        _validate_coherence(self.state, self.verification)
        if self.state.has_result and self.result is None:
            raise ValueError(
                f"CalculationRecord in state {self.state.value!r} must carry a result"
            )
        if self.state is CalculationState.FAILED and not self.error_message:
            raise ValueError(
                "CalculationRecord in state 'failed' must carry an error_message"
            )

    @property
    def reuse_key(self) -> str:
        """Molecular identity plus configuration signature, as one key."""
        return f"{self.molecule.molecule_id}/{self.signature.digest}"

    @property
    def is_terminal(self) -> bool:
        """Whether this record can still change state."""
        return self.state.is_terminal

    def transition_to(
        self,
        state: CalculationState,
        *,
        verification: VerificationOutcome | None = None,
        result: CalculationResultData | None = None,
        error_message: str | None = None,
        at: datetime | None = None,
    ) -> CalculationRecord:
        """Return a new record in ``state``, validating the transition.

        Raises:
            ValueError: If the transition is not allowed, or if the new
                state and verification outcome contradict each other.
        """
        allowed = ALLOWED_CALCULATION_TRANSITIONS[self.state]
        if state not in allowed:
            allowed_names = ", ".join(sorted(item.value for item in allowed)) or "none"
            raise ValueError(
                f"Illegal calculation transition {self.state.value!r} -> "
                f"{state.value!r}; allowed: {allowed_names}"
            )

        moment = _require_aware(at or utc_now(), "CalculationRecord transition time")
        new_verification = (
            verification if verification is not None else self.verification
        )
        _validate_coherence(state, new_verification)

        timestamps = self.timestamps
        if state is CalculationState.RUNNING:
            timestamps = replace(timestamps, started_at=moment)
        elif state.is_terminal:
            timestamps = replace(
                timestamps,
                started_at=timestamps.started_at or moment,
                completed_at=moment,
            )

        return replace(
            self,
            state=state,
            verification=new_verification,
            result=result if result is not None else self.result,
            error_message=(
                error_message if error_message is not None else self.error_message
            ),
            timestamps=timestamps,
        )

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "schema_version": self.schema_version,
            "calculation_id": self.calculation_id,
            "run_id": self.run_id,
            "reuse_key": self.reuse_key,
            "molecule": self.molecule.to_dict(),
            "signature": self.signature.to_dict(),
            "state": self.state.value,
            "verification": self.verification.value,
            "timestamps": self.timestamps.to_dict(),
            "provenance": self.provenance.to_dict(),
            "result": None if self.result is None else self.result.to_dict(),
            "error_message": self.error_message,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> CalculationRecord:
        """Deserialize from JSON-compatible primitives.

        ``calculation_id`` and ``reuse_key`` are recomputed from the
        identity fields, so a tampered file fails loudly instead of
        quietly pointing reuse at the wrong result.
        """
        result_payload = payload.get("result")
        return cls(
            run_id=str(payload["run_id"]),
            molecule=MolecularIdentity.from_dict(payload["molecule"]),
            signature=CalculationSignature.from_dict(payload["signature"]),
            provenance=ScientificProvenance.from_dict(payload["provenance"]),
            state=CalculationState(str(payload["state"])),
            verification=VerificationOutcome(str(payload["verification"])),
            timestamps=Timestamps.from_dict(payload["timestamps"]),
            result=(
                None
                if result_payload is None
                else CalculationResultData.from_dict(result_payload)
            ),
            error_message=_optional_text(payload.get("error_message")),
            calculation_id=str(payload.get("calculation_id") or ""),
            schema_version=int(payload.get("schema_version", RECORD_SCHEMA_VERSION)),
        )


@dataclass(frozen=True)
class RunRecord:
    """One execution of Semi-Imperium over a set of molecules.

    The run owns the effective configuration and its signature; the
    individual calculations reference the same signature so that a run's
    results are never a mix of incomparable settings.
    """

    run_id: str
    configuration: EffectiveConfiguration
    provenance: ScientificProvenance
    state: RunState = RunState.PENDING
    timestamps: Timestamps = field(default_factory=Timestamps)
    label: str | None = None
    molecule_ids: tuple[str, ...] = ()
    error_message: str | None = None
    schema_version: int = RECORD_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.run_id.strip():
            raise ValueError("RunRecord.run_id must not be empty")
        if self.state is RunState.FAILED and not self.error_message:
            raise ValueError("RunRecord in state 'failed' must carry an error_message")

    @property
    def signature(self) -> CalculationSignature:
        """Signature of this run's effective configuration."""
        return self.configuration.signature()

    def transition_to(
        self,
        state: RunState,
        *,
        error_message: str | None = None,
        at: datetime | None = None,
    ) -> RunRecord:
        """Return a new run record in ``state``, validating the transition."""
        allowed = ALLOWED_RUN_TRANSITIONS[self.state]
        if state not in allowed:
            allowed_names = ", ".join(sorted(item.value for item in allowed)) or "none"
            raise ValueError(
                f"Illegal run transition {self.state.value!r} -> {state.value!r}; "
                f"allowed: {allowed_names}"
            )

        moment = _require_aware(at or utc_now(), "RunRecord transition time")
        timestamps = self.timestamps
        if state is RunState.RUNNING:
            timestamps = replace(timestamps, started_at=moment)
        elif state.is_terminal:
            timestamps = replace(
                timestamps,
                started_at=timestamps.started_at or moment,
                completed_at=moment,
            )

        return replace(
            self,
            state=state,
            error_message=(
                error_message if error_message is not None else self.error_message
            ),
            timestamps=timestamps,
        )

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "schema_version": self.schema_version,
            "run_id": self.run_id,
            "state": self.state.value,
            "label": self.label,
            "timestamps": self.timestamps.to_dict(),
            "configuration": self.configuration.to_dict(),
            "signature": self.signature.to_dict(),
            "provenance": self.provenance.to_dict(),
            "molecule_ids": list(self.molecule_ids),
            "error_message": self.error_message,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> RunRecord:
        """Deserialize from JSON-compatible primitives.

        The stored signature is verified against a freshly computed one so
        that an edited configuration can never keep an old digest.
        """
        configuration = EffectiveConfiguration.from_dict(payload["configuration"])
        record = cls(
            run_id=str(payload["run_id"]),
            configuration=configuration,
            provenance=ScientificProvenance.from_dict(payload["provenance"]),
            state=RunState(str(payload["state"])),
            timestamps=Timestamps.from_dict(payload["timestamps"]),
            label=_optional_text(payload.get("label")),
            molecule_ids=tuple(str(item) for item in payload.get("molecule_ids") or ()),
            error_message=_optional_text(payload.get("error_message")),
            schema_version=int(payload.get("schema_version", RECORD_SCHEMA_VERSION)),
        )
        stored_signature = payload.get("signature")
        if stored_signature is not None:
            stored_digest = str(stored_signature["digest"])
            if stored_digest != record.signature.digest:
                raise ValueError(
                    "Stored run signature does not match its configuration: "
                    f"{stored_digest!r} != {record.signature.digest!r}"
                )
        return record


def build_calculation_id(
    *, run_id: str, molecule_id: str, signature_digest: str
) -> str:
    """Return the deterministic id of one calculation inside one run.

    Determinism is what guarantees a run never holds two records for the
    same molecule under the same configuration.
    """
    digest = stable_digest(
        {
            "run_id": run_id,
            "molecule_id": molecule_id,
            "signature_digest": signature_digest,
        }
    )
    return f"{CALCULATION_ID_PREFIX}-{digest[:CALCULATION_ID_LENGTH]}"


def build_reuse_key(
    molecule: MolecularIdentity, signature: CalculationSignature
) -> str:
    """Return the reuse key for a molecule under a configuration signature."""
    return f"{molecule.molecule_id}/{signature.digest}"


def _validate_coherence(
    state: CalculationState, verification: VerificationOutcome
) -> None:
    """Reject state/verification pairs that would misreport the science."""
    coherent = COHERENT_VERIFICATION_OUTCOMES[state]
    if verification not in coherent:
        coherent_names = ", ".join(sorted(item.value for item in coherent))
        raise ValueError(
            f"Verification outcome {verification.value!r} is incoherent with "
            f"calculation state {state.value!r}; expected one of: {coherent_names}"
        )


def _optional_text(value: Any) -> str | None:
    """Return ``value`` as text, mapping absent values to ``None``."""
    return None if value is None else str(value)


def _optional_float(value: Any) -> float | None:
    """Return ``value`` as a float, mapping absent values to ``None``."""
    return None if value is None else float(value)


def _optional_int(value: Any) -> int | None:
    """Return ``value`` as an int, mapping absent values to ``None``."""
    return None if value is None else int(value)


__all__ = [
    "CALCULATION_ID_PREFIX",
    "RECORD_SCHEMA_VERSION",
    "CalculationRecord",
    "CalculationResultData",
    "RunRecord",
    "ScientificProvenance",
    "Timestamps",
    "build_calculation_id",
    "build_reuse_key",
    "utc_now",
]
