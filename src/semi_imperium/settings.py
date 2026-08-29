"""Focused CREST, MOPAC and runtime defaults for Semi-Imperium.

Everything here is a *default*: a starting point for the next run and
nothing more. The moment a run is prepared, the resolved
:class:`~semi_imperium.domain.configuration.EffectiveConfiguration` is
copied into that run's record, so editing a default later can never
retroactively re-describe a number that has already been computed.

The split between the two halves is deliberate:

* CREST, conformer selection, MOPAC and verification change *what* is
  computed, so they end up inside the calculation signature;
* :class:`RuntimeSettings` (executables, threads, timeouts, directories)
  changes only *how long* it takes and *where* the files land, so it is
  kept out of the signature and out of reuse decisions entirely.
"""

from __future__ import annotations

import shutil
import tempfile
from dataclasses import dataclass, field, replace
from pathlib import Path

from semi_imperium.domain.configuration import (
    ConformerSearchSettings,
    ConformerSelectionSettings,
    EffectiveConfiguration,
    SemiempiricalSettings,
    VerificationSettings,
)

#: Method identity of the focused Semi-Imperium pipeline.
METHOD_ID = "semi_imperium_conformer_mopac"
METHOD_VERSION = "1.0"

#: The single property this application computes.
PROPERTY_ID = "standard_enthalpy_of_formation"


@dataclass(frozen=True)
class ToolReadiness:
    """Whether one external program can actually be launched right now."""

    name: str
    executable: str
    available: bool
    detail: str

    @property
    def label(self) -> str:
        """Short human-readable state, never a blank cell."""
        return "available" if self.available else "not found"


@dataclass(frozen=True)
class RuntimeSettings:
    """Execution-only knobs: they never enter a calculation signature."""

    crest_executable: str = "crest"
    mopac_executable: str = "mopac"
    crest_threads: int = 4
    crest_timeout_seconds: float = 300.0
    mopac_timeout_seconds: float = 120.0
    work_dir: Path = field(
        default_factory=lambda: Path(tempfile.gettempdir()) / "semi_imperium"
    )
    store_root: Path = field(
        default_factory=lambda: Path.cwd() / "data" / "semi_imperium"
    )

    def __post_init__(self) -> None:
        if self.crest_threads < 1:
            raise ValueError(
                f"RuntimeSettings.crest_threads must be >= 1, got {self.crest_threads}"
            )
        for label in ("crest_timeout_seconds", "mopac_timeout_seconds"):
            value = float(getattr(self, label))
            if value <= 0:
                raise ValueError(f"RuntimeSettings.{label} must be > 0, got {value}")

    def readiness(self) -> tuple[ToolReadiness, ...]:
        """Report, per external program, whether it is on PATH right now.

        This is a live probe rather than stored state: a tool that was
        installed after the last run must not keep showing as missing.
        """
        return (
            _probe("CREST", self.crest_executable),
            _probe("MOPAC", self.mopac_executable),
        )

    @property
    def is_ready(self) -> bool:
        """Whether every external program needed by a full run was found."""
        return all(tool.available for tool in self.readiness())


@dataclass(frozen=True)
class SemiImperiumSettings:
    """The complete set of defaults behind the next prepared run."""

    conformer_search: ConformerSearchSettings = field(
        default_factory=ConformerSearchSettings
    )
    conformer_selection: ConformerSelectionSettings = field(
        default_factory=ConformerSelectionSettings
    )
    semiempirical: SemiempiricalSettings = field(default_factory=SemiempiricalSettings)
    verification: VerificationSettings = field(default_factory=VerificationSettings)
    runtime: RuntimeSettings = field(default_factory=RuntimeSettings)
    method_id: str = METHOD_ID
    method_version: str = METHOD_VERSION
    property_id: str = PROPERTY_ID

    def configuration_for(
        self,
        hamiltonian: str,
        *,
        crest_enabled: bool,
    ) -> EffectiveConfiguration:
        """Resolve these defaults into one run's frozen effective configuration.

        Args:
            hamiltonian: The semiempirical Hamiltonian this configuration
                computes. AM1, PM3 and PM7 are independent requests, so each
                gets its own configuration and therefore its own signature.
            crest_enabled: Whether the conformer search runs for this
                molecule. Disabling it is a configuration choice, not a
                different method, so it stays inside the same signature.
        """
        return EffectiveConfiguration(
            method_id=self.method_id,
            method_version=self.method_version,
            property_id=self.property_id,
            conformer_search=replace(self.conformer_search, enabled=crest_enabled),
            conformer_selection=self.conformer_selection,
            semiempirical=replace(self.semiempirical, hamiltonian=hamiltonian),
            verification=self.verification,
        )


def _probe(name: str, executable: str) -> ToolReadiness:
    """Locate one executable, accepting either a PATH name or a full path."""
    candidate = Path(executable)
    if candidate.is_absolute():
        available = candidate.is_file()
        detail = str(candidate)
    else:
        resolved = shutil.which(executable)
        available = resolved is not None
        detail = resolved or f"{executable!r} is not on PATH"
    return ToolReadiness(
        name=name,
        executable=executable,
        available=available,
        detail=detail,
    )


__all__ = [
    "METHOD_ID",
    "METHOD_VERSION",
    "PROPERTY_ID",
    "RuntimeSettings",
    "SemiImperiumSettings",
    "ToolReadiness",
]
