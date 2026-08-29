"""Concrete MOPAC executable adapter for optimization and ``FORCE``.

The minimum workflow deliberately depends on a small protocol.  This module is
the production implementation of that protocol: it writes actual MOPAC data
sets, invokes the configured executable, reads the optimized Cartesian
geometry and heat of formation, and extracts normal-coordinate vectors needed
for bounded saddle recovery.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from grimperium.crest_pm7.config import PM7Config
from grimperium.crest_pm7.energy_extractor import extract_hof
from semi_imperium.conformers import ConformerGeometry
from semi_imperium.mopac.models import (
    SUPPORTED_HAMILTONIANS,
    DisplacementLineage,
    ForceRun,
    OptimizationRun,
)

_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][-+]?\d+)?"
_CARTESIAN_ROW = re.compile(
    rf"^\s*\d+\s+(?P<element>[A-Za-z]{{1,3}})\s+"
    rf"(?P<x>{_NUMBER})\s+(?P<y>{_NUMBER})\s+(?P<z>{_NUMBER})(?:\s|$)"
)
_NORMAL_SECTION = re.compile(
    r"\bNORMAL COORDINATE ANALYSIS\b(?P<body>.*?)"
    r"(?=\bMASS-WEIGHTED COORDINATE ANALYSIS\b|"
    r"\bDESCRIPTION OF VIBRATIONS\b|\Z)",
    re.IGNORECASE | re.DOTALL,
)
_ROOT_HEADER = re.compile(r"^\s*Root\s+No\.?(?P<inline>.*)$", re.IGNORECASE)
_INTEGER_LINE = re.compile(r"^\s*\d+(?:\s+\d+)*\s*$")
_OPTIMIZATION_FAILURES = (
    re.compile(r"GEOMETRY\s+OPTIMIZATION\s+FAILED", re.IGNORECASE),
    re.compile(r"EXCESS\s+NUMBER\s+OF\s+OPTIMIZATION\s+CYCLES", re.IGNORECASE),
    re.compile(r"UNABLE\s+TO\s+ACHIEVE\s+SELF[- ]CONSISTENCE", re.IGNORECASE),
    re.compile(r"SCF\s+(?:FAILED|NOT\s+CONVERGED)", re.IGNORECASE),
    re.compile(r"ABNORMAL\s+TERMINATION", re.IGNORECASE),
)
_MULTIPLICITY_KEYWORDS = {
    2: "DOUBLET",
    3: "TRIPLET",
    4: "QUARTET",
    5: "QUINTET",
    6: "SEXTET",
}


class CommandRunner(Protocol):
    """Subprocess boundary, injectable for deterministic adapter tests."""

    def run(
        self, argv: list[str], *, cwd: Path, timeout: float
    ) -> subprocess.CompletedProcess[str]:
        """Execute one MOPAC data set."""
        ...


class SubprocessCommandRunner:
    """Default command runner used by :class:`MopacExecutableBackend`."""

    def run(
        self, argv: list[str], *, cwd: Path, timeout: float
    ) -> subprocess.CompletedProcess[str]:
        # The executable is explicit local configuration and shell=False keeps
        # molecule/attempt data out of command interpretation.
        return subprocess.run(  # noqa: S603
            argv,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )


@dataclass(frozen=True)
class MopacExecutableSettings:
    """Execution details shared by every independent Hamiltonian run."""

    executable: str
    work_dir: Path
    timeout_seconds: float
    calculation_id: str
    scf_convergence: float = 1.0e-4
    charge: int = 0
    multiplicity: int = 1
    extra_optimization_keywords: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.executable.strip():
            raise ValueError("MopacExecutableSettings.executable must not be empty")
        if self.timeout_seconds <= 0:
            raise ValueError("MopacExecutableSettings.timeout_seconds must be > 0")
        _require_safe_component(self.calculation_id, label="calculation_id")
        if self.scf_convergence <= 0:
            raise ValueError("MopacExecutableSettings.scf_convergence must be > 0")
        if self.multiplicity < 1 or self.multiplicity > 6:
            raise ValueError("MopacExecutableSettings.multiplicity must be 1 through 6")

    @classmethod
    def from_pm7_config(
        cls,
        config: PM7Config,
        *,
        calculation_id: str,
        work_dir: Path | None = None,
        charge: int = 0,
        multiplicity: int = 1,
    ) -> MopacExecutableSettings:
        """Use Grimperium's installed MOPAC executable and execution limits."""
        return cls(
            executable=config.mopac_executable,
            work_dir=work_dir if work_dir is not None else config.temp_dir / "minima",
            timeout_seconds=config.mopac_timeout_base * config.mopac_timeout_margin,
            calculation_id=calculation_id,
            scf_convergence=config.mopac_scf_threshold,
            charge=charge,
            multiplicity=multiplicity,
        )


class MopacExecutableBackend:
    """Run real MOPAC optimization and FORCE jobs in auditable attempt folders."""

    def __init__(
        self,
        settings: MopacExecutableSettings,
        *,
        runner: CommandRunner | None = None,
    ) -> None:
        self.settings = settings
        self.runner = runner or SubprocessCommandRunner()

    @classmethod
    def from_pm7_config(
        cls,
        config: PM7Config,
        *,
        calculation_id: str,
        work_dir: Path | None = None,
        charge: int = 0,
        multiplicity: int = 1,
        runner: CommandRunner | None = None,
    ) -> MopacExecutableBackend:
        """Build the concrete backend from the active Grimperium config."""
        settings = MopacExecutableSettings.from_pm7_config(
            config,
            calculation_id=calculation_id,
            work_dir=work_dir,
            charge=charge,
            multiplicity=multiplicity,
        )
        return cls(settings, runner=runner)

    def optimize(
        self,
        *,
        hamiltonian: str,
        geometry: ConformerGeometry,
        source_conformer_index: int,
        attempt_id: str,
        displacement: DisplacementLineage | None,
    ) -> OptimizationRun:
        """Optimize one candidate and return only a genuinely usable result."""
        del source_conformer_index, displacement
        method = _normalize_hamiltonian(hamiltonian)
        attempt_dir = self._attempt_dir(method, attempt_id)
        input_path = attempt_dir / "optimization.mop"
        output_path = input_path.with_suffix(".out")
        keywords = (
            method,
            "EF",
            "PRECISE",
            f"SCFCRT={self.settings.scf_convergence:.1e}",
            *self.settings.extra_optimization_keywords,
        )
        _write_mopac_input(
            input_path,
            geometry,
            keywords=keywords,
            title=f"{attempt_id} geometry optimization",
            optimize=True,
            charge=self.settings.charge,
            multiplicity=self.settings.multiplicity,
        )
        execution_error = self._execute(input_path, output_path)
        if execution_error is not None:
            return OptimizationRun(
                converged=False,
                output_path=str(output_path) if output_path.exists() else None,
                error_message=execution_error,
            )

        output = output_path.read_text(encoding="utf-8", errors="replace")
        failure = next(
            (
                pattern.pattern
                for pattern in _OPTIMIZATION_FAILURES
                if pattern.search(output)
            ),
            None,
        )
        if failure is not None:
            return OptimizationRun(
                converged=False,
                output_path=str(output_path),
                error_message=f"MOPAC optimization diagnostic matched {failure!r}",
            )
        if "MOPAC DONE" not in output.upper():
            return OptimizationRun(
                converged=False,
                output_path=str(output_path),
                error_message="MOPAC optimization output has no completion marker",
            )
        optimized = parse_last_cartesian_geometry(
            output, expected_elements=geometry.elements
        )
        if optimized is None:
            return OptimizationRun(
                converged=False,
                output_path=str(output_path),
                error_message="MOPAC output has no complete final Cartesian geometry",
            )
        heat, _, _ = extract_hof(output, nheavy=None)
        if heat is None:
            return OptimizationRun(
                converged=False,
                output_path=str(output_path),
                error_message="MOPAC output has no parseable final heat of formation",
            )
        return OptimizationRun(
            converged=True,
            geometry=optimized,
            heat_of_formation_kcal_mol=heat,
            output_path=str(output_path),
        )

    def verify_force(
        self,
        *,
        hamiltonian: str,
        optimization: OptimizationRun,
        attempt_id: str,
    ) -> ForceRun:
        """Run FORCE on the exact optimized geometry and retain its mode vectors."""
        method = _normalize_hamiltonian(hamiltonian)
        geometry = optimization.geometry
        if not optimization.converged or geometry is None:
            return ForceRun(
                output="",
                execution_error="FORCE requires a converged optimized geometry",
            )
        attempt_dir = self._attempt_dir(method, attempt_id)
        input_path = attempt_dir / "verification-force.mop"
        output_path = input_path.with_suffix(".out")
        _write_mopac_input(
            input_path,
            geometry,
            keywords=(method, "FORCE", "NOREOR"),
            title=f"{attempt_id} minimum verification",
            optimize=False,
            charge=self.settings.charge,
            multiplicity=self.settings.multiplicity,
        )
        execution_error = self._execute(input_path, output_path)
        if not output_path.exists():
            return ForceRun(output="", execution_error=execution_error)
        output = output_path.read_text(encoding="utf-8", errors="replace")
        return ForceRun(
            output=output,
            output_path=str(output_path),
            normal_modes=parse_normal_coordinate_vectors(
                output, atom_count=geometry.atom_count
            ),
            execution_error=execution_error,
        )

    def _attempt_dir(self, hamiltonian: str, attempt_id: str) -> Path:
        _require_safe_component(attempt_id, label="attempt_id")
        path = (
            self.settings.work_dir
            / self.settings.calculation_id
            / hamiltonian.lower()
            / attempt_id
        )
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _execute(self, input_path: Path, output_path: Path) -> str | None:
        try:
            completed = self.runner.run(
                [self.settings.executable, str(input_path)],
                cwd=input_path.parent,
                timeout=self.settings.timeout_seconds,
            )
        except subprocess.TimeoutExpired:
            return f"MOPAC timed out after {self.settings.timeout_seconds:g} seconds"
        except OSError as exc:
            return f"MOPAC could not be executed: {exc}"
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "").strip()[:240]
            suffix = f": {detail}" if detail else ""
            return f"MOPAC exited with code {completed.returncode}{suffix}"
        if not output_path.exists():
            return f"MOPAC produced no output file at {output_path}"
        return None


def parse_last_cartesian_geometry(
    output: str, *, expected_elements: tuple[str, ...]
) -> ConformerGeometry | None:
    """Return the last complete MOPAC Cartesian-coordinate table."""
    tables: list[ConformerGeometry] = []
    lines = output.splitlines()
    for index, line in enumerate(lines):
        if "CARTESIAN COORDINATES" not in line.upper():
            continue
        elements: list[str] = []
        coordinates: list[tuple[float, float, float]] = []
        started = False
        for candidate in lines[index + 1 :]:
            match = _CARTESIAN_ROW.match(candidate)
            if match is None:
                if started:
                    break
                continue
            started = True
            elements.append(match.group("element"))
            coordinates.append(
                (
                    _as_float(match.group("x")),
                    _as_float(match.group("y")),
                    _as_float(match.group("z")),
                )
            )
        if tuple(elements) == expected_elements:
            tables.append(
                ConformerGeometry(
                    elements=tuple(elements), coordinates=tuple(coordinates)
                )
            )
    return tables[-1] if tables else None


def parse_normal_coordinate_vectors(
    output: str, *, atom_count: int
) -> dict[int, tuple[tuple[float, float, float], ...]]:
    """Parse MOPAC's projected normal-coordinate tables by one-based root."""
    if atom_count < 1:
        raise ValueError("atom_count must be >= 1")
    section = _NORMAL_SECTION.search(output)
    if section is None:
        return {}
    lines = section.group("body").splitlines()
    vectors: dict[int, tuple[tuple[float, float, float], ...]] = {}
    cursor = 0
    coordinate_count = 3 * atom_count
    while cursor < len(lines):
        header = _ROOT_HEADER.match(lines[cursor])
        if header is None:
            cursor += 1
            continue
        roots = _parse_roots(header.group("inline"))
        cursor += 1
        if not roots:
            while cursor < len(lines) and not lines[cursor].strip():
                cursor += 1
            if cursor < len(lines) and _INTEGER_LINE.match(lines[cursor]):
                roots = tuple(int(value) for value in lines[cursor].split())
                cursor += 1
        if not roots:
            continue
        rows: list[tuple[float, ...]] = []
        expected_row = 1
        scan = cursor
        while scan < len(lines) and expected_row <= coordinate_count:
            if _ROOT_HEADER.match(lines[scan]):
                break
            row = _coordinate_row(lines[scan], roots=len(roots))
            scan += 1
            if row is None or int(row[0]) != expected_row:
                continue
            rows.append(row[1:])
            expected_row += 1
        if len(rows) == coordinate_count:
            for column, root in enumerate(roots):
                flat = tuple(row[column] for row in rows)
                vectors[root] = tuple(
                    (flat[position], flat[position + 1], flat[position + 2])
                    for position in range(0, coordinate_count, 3)
                )
            cursor = scan
        else:
            cursor += 1
    return vectors


def _write_mopac_input(
    path: Path,
    geometry: ConformerGeometry,
    *,
    keywords: tuple[str, ...],
    title: str,
    optimize: bool,
    charge: int,
    multiplicity: int,
) -> None:
    words = list(keywords)
    if charge:
        words.append(f"CHARGE={charge}")
    if multiplicity != 1:
        words.append(_MULTIPLICITY_KEYWORDS[multiplicity])
    flag = 1 if optimize else 0
    lines = [" ".join(words), title, ""]
    for element, (x, y, z) in zip(geometry.elements, geometry.coordinates):
        lines.append(
            f"  {element:2s} {x:16.10f} {flag} {y:16.10f} {flag} " f"{z:16.10f} {flag}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _normalize_hamiltonian(hamiltonian: str) -> str:
    normalized = hamiltonian.strip().upper()
    if normalized not in SUPPORTED_HAMILTONIANS:
        raise ValueError(f"Unsupported Hamiltonian {hamiltonian!r}")
    return normalized


def _require_safe_component(value: str, *, label: str) -> None:
    if (
        not value
        or value in {".", ".."}
        or Path(value).name != value
        or not re.fullmatch(r"[A-Za-z0-9_.-]+", value)
    ):
        raise ValueError(f"{label} must be one safe path component")


def _parse_roots(text: str) -> tuple[int, ...]:
    stripped = text.strip()
    if not stripped or _INTEGER_LINE.match(stripped) is None:
        return ()
    return tuple(int(value) for value in stripped.split())


def _coordinate_row(line: str, *, roots: int) -> tuple[float, ...] | None:
    tokens = line.split()
    if len(tokens) != roots + 1 or not tokens[0].isdigit():
        return None
    try:
        return tuple(_as_float(value) for value in tokens)
    except ValueError:
        return None


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


__all__ = [
    "CommandRunner",
    "MopacExecutableBackend",
    "MopacExecutableSettings",
    "SubprocessCommandRunner",
    "parse_last_cartesian_geometry",
    "parse_normal_coordinate_vectors",
]
