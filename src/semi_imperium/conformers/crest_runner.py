"""CREST executed as a local subprocess, one directory per molecule.

:class:`SubprocessCrestRunner` is the concrete :class:`CrestRunner` the
rest of the conformer stage was written against. It owns the only code
that spawns CREST (and the optional xTB pre-optimization), and it keeps
three guarantees that the science downstream depends on:

* the CREST input is the RDKit embedding of the request's SMILES with
  explicit hydrogens, so the ensemble keeps the atom order
  :class:`~semi_imperium.conformers.topology.SmilesTopology` derives;
* only ``crest_conformers.xyz`` is accepted as the ensemble. A run that
  does not write it is an error, never a silent fallback to
  ``crest_best.xyz`` (one structure is not an ensemble);
* a finished run is recorded next to its ensemble with the settings it
  was produced under, so asking again for the same molecule and the same
  settings reuses it instead of running CREST twice, and asking with
  different settings is refused instead of mixing ensembles.

RDKit is only needed to embed the input structure, so the default
initial-structure backend is imported lazily; tests inject both that
backend and the process runner and never spawn anything.
"""

from __future__ import annotations

import json
import re
import subprocess
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

from semi_imperium.conformers.backends import (
    ConformerBackendError,
    ConformerRequest,
    InitialStructureBackend,
)
from semi_imperium.conformers.crest import CrestRun
from semi_imperium.conformers.ensemble import (
    UNKNOWN_PROGRAM_VERSION,
    ConformerGeometry,
)
from semi_imperium.domain.configuration import ConformerSearchSettings
from semi_imperium.domain.hashing import stable_digest

#: Files the runner reads or writes inside each molecule directory.
ENSEMBLE_FILENAME = "crest_conformers.xyz"
RUN_RECORD_FILENAME = "crest_run.json"
INPUT_FILENAME = "input.xyz"
PREOPT_FILENAME = "preopt.xyz"
CREST_STDOUT_FILENAME = "crest.out"
CREST_STDERR_FILENAME = "crest.err"

#: Error codes raised by this module.
CREST_UNAVAILABLE = "crest_unavailable"
CREST_TIMEOUT = "crest_timeout"
CREST_MISSING_ENSEMBLE = "crest_missing_ensemble"
CREST_CACHE_SETTINGS_MISMATCH = "crest_cache_settings_mismatch"
CREST_CACHE_UNRECORDED = "crest_cache_unrecorded"
CREST_PREOPT_FAILED = "crest_preopt_failed"

#: Bumped when the meaning of a stored run record changes.
RUN_RECORD_VERSION = 1

_METHOD_FLAGS = {
    "gfn2": "--gfn2",
    "gfnff": "--gfnff",
    "gfn2//gfnff": "--gfn2//gfnff",
}
_QUICK_FLAGS = {
    "off": None,
    "quick": "--quick",
    "squick": "--squick",
    "mquick": "--mquick",
}
_PREOPTIMIZERS = ("none", "xtb")
_VERSION_PATTERN = re.compile(r"Version\s+(\d+(?:\.\d+)+)")
_SAFE_MOLECULE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
_STDERR_TAIL_CHARS = 2000


class ProcessRunner(Protocol):
    """Subprocess boundary, injectable for deterministic adapter tests."""

    def run(
        self, argv: list[str], *, cwd: Path, timeout: float
    ) -> subprocess.CompletedProcess[str]:
        """Execute one command and return its completed process."""
        ...


class SubprocessProcessRunner:
    """Default process runner: plain :func:`subprocess.run`, no shell."""

    def run(
        self, argv: list[str], *, cwd: Path, timeout: float
    ) -> subprocess.CompletedProcess[str]:
        # The executables are explicit local configuration and shell=False
        # keeps molecule data out of command interpretation.
        return subprocess.run(  # noqa: S603
            argv,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )


def _default_initial_structure() -> InitialStructureBackend:
    """Import the RDKit route only when a runner is actually built."""
    from semi_imperium.conformers.initial_structure import RDKitInitialStructure

    return RDKitInitialStructure()


@dataclass(frozen=True)
class SubprocessCrestRunner:
    """Runs CREST once per molecule and settings, reusing finished runs.

    Each molecule gets ``work_root / molecule_id``. Threads and timeouts
    are execution detail: they are recorded in the run record but do not
    take part in deciding whether a stored ensemble can be reused.
    """

    work_root: Path
    executable: str = "crest"
    xtb_executable: str = "xtb"
    threads: int = 4
    timeout_seconds: float = 7200.0
    xtb_timeout_seconds: float = 300.0
    initial_structure: InitialStructureBackend = field(
        default_factory=_default_initial_structure
    )
    process_runner: ProcessRunner = field(default_factory=SubprocessProcessRunner)

    def __post_init__(self) -> None:
        if self.threads < 1:
            raise ValueError(
                f"SubprocessCrestRunner.threads must be >= 1, got {self.threads}"
            )
        if self.timeout_seconds <= 0:
            raise ValueError(
                "SubprocessCrestRunner.timeout_seconds must be > 0, "
                f"got {self.timeout_seconds}"
            )
        if self.xtb_timeout_seconds <= 0:
            raise ValueError(
                "SubprocessCrestRunner.xtb_timeout_seconds must be > 0, "
                f"got {self.xtb_timeout_seconds}"
            )

    def molecule_dir(self, request: ConformerRequest) -> Path:
        """Return the directory that holds ``request``'s CREST run."""
        if not _SAFE_MOLECULE_ID.match(request.molecule_id):
            raise ValueError(
                f"Molecule id {request.molecule_id!r} cannot name a directory; "
                "use letters, digits, '_', '-' or '.'"
            )
        return self.work_root / request.molecule_id

    def run_record(self, request: ConformerRequest) -> dict[str, Any] | None:
        """Return the stored record of a finished run, if there is one."""
        path = self.molecule_dir(request) / RUN_RECORD_FILENAME
        if not path.exists():
            return None
        payload: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
        return payload

    def run(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> CrestRun:
        """Return the CREST ensemble for ``request``, running CREST if needed.

        Raises:
            ConformerBackendError: If a stored ensemble cannot be reused,
                CREST or xTB is missing, times out or fails to write the
                ensemble, or the pre-optimization fails.
            ValueError: If the settings name a method, quick mode or
                pre-optimizer this runner does not know.
        """
        mol_dir = self.molecule_dir(request)
        cached = self._cached(request, settings, mol_dir)
        if cached is not None:
            return cached

        mol_dir.mkdir(parents=True, exist_ok=True)
        input_name = self._prepare_input(request, settings, mol_dir)
        command = self.crest_command(request, settings, input_name=input_name)
        started = time.perf_counter()
        completed = self._execute(
            command,
            cwd=mol_dir,
            timeout=self.timeout_seconds,
            program="CREST",
            molecule_id=request.molecule_id,
        )
        elapsed = time.perf_counter() - started
        (mol_dir / CREST_STDOUT_FILENAME).write_text(
            completed.stdout or "", encoding="utf-8"
        )
        (mol_dir / CREST_STDERR_FILENAME).write_text(
            completed.stderr or "", encoding="utf-8"
        )
        version = parse_crest_version(completed.stdout or "")
        if completed.returncode != 0:
            return CrestRun(
                ensemble_xyz="",
                program_version=version,
                command=tuple(command),
                exit_code=completed.returncode,
                stderr=(completed.stderr or "")[-_STDERR_TAIL_CHARS:],
            )

        ensemble_path = mol_dir / ENSEMBLE_FILENAME
        if not ensemble_path.exists():
            raise ConformerBackendError(
                f"CREST finished for {request.molecule_id!r} without writing "
                f"{ENSEMBLE_FILENAME}; a single best structure is not accepted "
                "as an ensemble",
                code=CREST_MISSING_ENSEMBLE,
            )
        self._write_record(
            mol_dir,
            settings=settings,
            command=command,
            program_version=version,
            elapsed_seconds=elapsed,
        )
        return CrestRun(
            ensemble_xyz=ensemble_path.read_text(encoding="utf-8"),
            program_version=version,
            command=tuple(command),
        )

    def crest_command(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
        *,
        input_name: str,
    ) -> list[str]:
        """Build the CREST argument vector for ``settings``.

        Raises:
            ValueError: If the method or quick mode is unknown.
        """
        method = _METHOD_FLAGS.get(settings.method)
        if method is None:
            known = ", ".join(sorted(_METHOD_FLAGS))
            raise ValueError(
                f"Unknown CREST method {settings.method!r}; expected one of: {known}"
            )
        if settings.quick_mode not in _QUICK_FLAGS:
            known = ", ".join(sorted(_QUICK_FLAGS))
            raise ValueError(
                f"Unknown CREST quick mode {settings.quick_mode!r}; "
                f"expected one of: {known}"
            )
        command = [self.executable, input_name, method]
        if settings.use_v3:
            command.append("--v3")
        if settings.nci:
            command.append("--nci")
        quick = _QUICK_FLAGS[settings.quick_mode]
        if quick is not None:
            command.append(quick)
        command.extend(
            [
                "--ewin",
                str(settings.energy_window_kcal_mol),
                "--rthr",
                str(settings.rmsd_threshold),
                "--opt",
                str(settings.opt_level),
                "--chrg",
                str(request.charge),
                "--uhf",
                str(request.multiplicity - 1),
                "--T",
                str(self.threads),
            ]
        )
        return command

    def _cached(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
        mol_dir: Path,
    ) -> CrestRun | None:
        """Return the stored run when it was produced under ``settings``."""
        ensemble_path = mol_dir / ENSEMBLE_FILENAME
        record = self.run_record(request)
        if record is None:
            if ensemble_path.exists():
                raise ConformerBackendError(
                    f"{ensemble_path} exists without a {RUN_RECORD_FILENAME}; its "
                    "settings are unknown, so it is neither reused nor overwritten",
                    code=CREST_CACHE_UNRECORDED,
                )
            return None
        if record.get("settings_digest") != _settings_digest(settings):
            raise ConformerBackendError(
                f"The stored CREST run in {mol_dir} was produced under different "
                "search settings; use another work directory instead of mixing "
                "ensembles",
                code=CREST_CACHE_SETTINGS_MISMATCH,
            )
        if not ensemble_path.exists():
            raise ConformerBackendError(
                f"{mol_dir / RUN_RECORD_FILENAME} records a finished run but "
                f"{ENSEMBLE_FILENAME} is missing",
                code=CREST_MISSING_ENSEMBLE,
            )
        return CrestRun(
            ensemble_xyz=ensemble_path.read_text(encoding="utf-8"),
            program_version=str(
                record.get("program_version") or UNKNOWN_PROGRAM_VERSION
            ),
            command=tuple(str(part) for part in record.get("command", ())),
        )

    def _prepare_input(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
        mol_dir: Path,
    ) -> str:
        """Write the RDKit input (optionally xTB-relaxed); return its name."""
        if settings.preoptimizer not in _PREOPTIMIZERS:
            known = ", ".join(_PREOPTIMIZERS)
            raise ValueError(
                f"Unknown pre-optimizer {settings.preoptimizer!r}; "
                f"expected one of: {known}"
            )
        initial = self.initial_structure.build(request, settings)
        geometry = initial.conformers[0].geometry
        (mol_dir / INPUT_FILENAME).write_text(
            geometry_to_xyz(geometry, comment=request.molecule_id),
            encoding="utf-8",
        )
        if settings.preoptimizer == "none":
            return INPUT_FILENAME

        # xTB keeps the atom order of its input, so the relaxed structure
        # still matches the SMILES topology.
        command = [
            self.xtb_executable,
            INPUT_FILENAME,
            "--opt",
            "--gfn2",
            "--chrg",
            str(request.charge),
            "--uhf",
            str(request.multiplicity - 1),
        ]
        completed = self._execute(
            command,
            cwd=mol_dir,
            timeout=self.xtb_timeout_seconds,
            program="xTB",
            molecule_id=request.molecule_id,
        )
        optimized = mol_dir / "xtbopt.xyz"
        if completed.returncode != 0 or not optimized.exists():
            detail = (completed.stderr or "").strip()[-_STDERR_TAIL_CHARS:]
            raise ConformerBackendError(
                f"xTB pre-optimization failed for {request.molecule_id!r} "
                f"(exit code {completed.returncode}): {detail or 'no stderr'}",
                code=CREST_PREOPT_FAILED,
            )
        optimized.replace(mol_dir / PREOPT_FILENAME)
        return PREOPT_FILENAME

    def _execute(
        self,
        command: list[str],
        *,
        cwd: Path,
        timeout: float,
        program: str,
        molecule_id: str,
    ) -> subprocess.CompletedProcess[str]:
        """Run one external program, translating launch failures to codes."""
        try:
            return self.process_runner.run(command, cwd=cwd, timeout=timeout)
        except FileNotFoundError as exc:
            raise ConformerBackendError(
                f"{program} executable {command[0]!r} was not found: {exc}",
                code=CREST_UNAVAILABLE,
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise ConformerBackendError(
                f"{program} timed out after {timeout:g} s for {molecule_id!r}",
                code=CREST_TIMEOUT,
            ) from exc

    def _write_record(
        self,
        mol_dir: Path,
        *,
        settings: ConformerSearchSettings,
        command: Sequence[str],
        program_version: str,
        elapsed_seconds: float,
    ) -> None:
        """Persist the finished run atomically, after the ensemble exists."""
        record = {
            "record_version": RUN_RECORD_VERSION,
            "settings": settings.to_dict(),
            "settings_digest": _settings_digest(settings),
            "command": list(command),
            "program_version": program_version,
            "threads": self.threads,
            "timeout_seconds": self.timeout_seconds,
            "elapsed_seconds": elapsed_seconds,
        }
        target = mol_dir / RUN_RECORD_FILENAME
        partial = target.with_suffix(".json.partial")
        partial.write_text(json.dumps(record, indent=2), encoding="utf-8")
        partial.replace(target)


def parse_crest_version(stdout: str) -> str:
    """Read the version CREST prints in its banner, or say it is unknown."""
    match = _VERSION_PATTERN.search(stdout)
    return match.group(1) if match else UNKNOWN_PROGRAM_VERSION


def geometry_to_xyz(geometry: ConformerGeometry, *, comment: str = "") -> str:
    """Render ``geometry`` as one XYZ block, preserving its atom order."""
    lines = [str(geometry.atom_count), comment.replace("\n", " ")]
    for element, (x, y, z) in zip(geometry.elements, geometry.coordinates, strict=True):
        lines.append(f"{element:<2} {x:16.10f} {y:16.10f} {z:16.10f}")
    return "\n".join(lines) + "\n"


def _settings_digest(settings: ConformerSearchSettings) -> str:
    """Digest of the settings that decide which ensemble CREST samples."""
    return stable_digest(settings.to_dict())


__all__ = [
    "CREST_CACHE_SETTINGS_MISMATCH",
    "CREST_CACHE_UNRECORDED",
    "CREST_MISSING_ENSEMBLE",
    "CREST_PREOPT_FAILED",
    "CREST_TIMEOUT",
    "CREST_UNAVAILABLE",
    "ENSEMBLE_FILENAME",
    "RUN_RECORD_FILENAME",
    "ProcessRunner",
    "SubprocessCrestRunner",
    "SubprocessProcessRunner",
    "geometry_to_xyz",
    "parse_crest_version",
]
