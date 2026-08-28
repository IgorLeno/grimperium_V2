"""CREST adapter: raw search output in, typed ensemble out.

The process handling stays in the injected :class:`CrestRunner`, which
is the only thing that ever touches an executable. Everything in this
module is pure parsing and bookkeeping, so the CREST contract can be
exercised end to end from recorded output without CREST installed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from semi_imperium.conformers.backends import (
    ConformerBackendError,
    ConformerRequest,
)
from semi_imperium.conformers.ensemble import (
    UNKNOWN_PROGRAM_VERSION,
    Conformer,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerSearchProvenance,
)
from semi_imperium.domain.configuration import ConformerSearchSettings
from semi_imperium.domain.enums import ConformerSource

#: CREST reports total energies in Hartree on the XYZ comment line.
HARTREE_TO_KCAL_MOL = 627.5094740631

_ENERGY_FACTORS = {
    "hartree": HARTREE_TO_KCAL_MOL,
    "kcal/mol": 1.0,
}


@dataclass(frozen=True)
class CrestRun:
    """Raw output of one CREST execution, as the runner observed it."""

    ensemble_xyz: str
    program_version: str = UNKNOWN_PROGRAM_VERSION
    command: tuple[str, ...] = ()
    exit_code: int = 0
    stderr: str = ""


class CrestRunner(Protocol):
    """Executes CREST for one molecule and reports what it produced."""

    def run(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> CrestRun:
        """Return the raw CREST output for ``request``."""
        ...


@dataclass(frozen=True)
class CrestConformerSearch:
    """Turns one CREST execution into a provenance-carrying ensemble."""

    runner: CrestRunner
    energy_unit: str = "hartree"

    def search(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        """Run the search and parse its ensemble.

        Raises:
            ConformerBackendError: If the search is disabled, the runner
                reported a non-zero exit code, or the ensemble could not
                be parsed.
        """
        if not settings.enabled:
            raise ConformerBackendError(
                "The CREST search is disabled for this configuration; use the "
                "initial-3D route instead of calling the search backend",
                code="crest_disabled",
            )
        run = self.runner.run(request, settings)
        if run.exit_code != 0:
            detail = run.stderr.strip() or "no stderr reported"
            raise ConformerBackendError(
                f"CREST exited with code {run.exit_code} for "
                f"{request.molecule_id!r}: {detail}",
                code="crest_failed",
            )
        conformers = parse_crest_ensemble(
            run.ensemble_xyz,
            energy_unit=self.energy_unit,
        )
        provenance = ConformerSearchProvenance(
            source=ConformerSource.CREST,
            program=settings.program,
            program_version=run.program_version or UNKNOWN_PROGRAM_VERSION,
            settings=settings,
            run_id=request.run_id,
            command=run.command,
        )
        return ConformerEnsemble(conformers=conformers, provenance=provenance)


def parse_crest_ensemble(
    text: str,
    *,
    energy_unit: str = "hartree",
) -> tuple[Conformer, ...]:
    """Parse a multi-structure CREST XYZ ensemble.

    The whole ensemble is returned in the order CREST wrote it, with no
    conformer dropped: a malformed block raises instead of being skipped,
    because a silently shorter ensemble would change the science without
    telling anyone.

    Args:
        text: Contents of ``crest_conformers.xyz`` or an equivalent file.
        energy_unit: Unit of the energy on each comment line.

    Raises:
        ConformerBackendError: If the ensemble is empty or malformed.
        ValueError: If ``energy_unit`` is unknown.
    """
    factor = _ENERGY_FACTORS.get(energy_unit)
    if factor is None:
        known = ", ".join(sorted(_ENERGY_FACTORS))
        raise ValueError(
            f"Unknown CREST energy unit {energy_unit!r}; expected one of: {known}"
        )

    lines = text.strip().splitlines()
    conformers: list[Conformer] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        conformer, cursor = _parse_block(lines, cursor, len(conformers), factor)
        conformers.append(conformer)

    if not conformers:
        raise ConformerBackendError(
            "CREST produced no conformer in its ensemble output",
            code="crest_empty_ensemble",
        )
    return tuple(conformers)


def _parse_block(
    lines: list[str],
    cursor: int,
    index: int,
    factor: float,
) -> tuple[Conformer, int]:
    """Parse one XYZ block, returning the conformer and the next cursor."""
    header = lines[cursor].strip()
    try:
        atom_count = int(header)
    except ValueError as exc:
        raise ConformerBackendError(
            f"CREST ensemble line {cursor + 1} should hold an atom count, "
            f"got {header!r}",
            code="crest_parse_failed",
        ) from exc
    if atom_count < 1:
        raise ConformerBackendError(
            f"CREST ensemble line {cursor + 1} declares {atom_count} atoms",
            code="crest_parse_failed",
        )

    end = cursor + 2 + atom_count
    if end > len(lines):
        raise ConformerBackendError(
            f"CREST ensemble block starting at line {cursor + 1} is truncated: "
            f"{atom_count} atoms were declared",
            code="crest_parse_failed",
        )

    elements: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for offset, line in enumerate(lines[cursor + 2 : end]):
        parts = line.split()
        if len(parts) < 4:
            raise ConformerBackendError(
                f"CREST ensemble line {cursor + 3 + offset} is not an atom line: "
                f"{line!r}",
                code="crest_parse_failed",
            )
        try:
            values = [float(value) for value in parts[1:4]]
        except ValueError as exc:
            raise ConformerBackendError(
                f"CREST ensemble line {cursor + 3 + offset} holds unreadable "
                f"coordinates: {line!r}",
                code="crest_parse_failed",
            ) from exc
        elements.append(parts[0])
        coordinates.append((values[0], values[1], values[2]))

    geometry = ConformerGeometry(
        elements=tuple(elements),
        coordinates=tuple(coordinates),
    )
    conformer = Conformer(
        index=index,
        geometry=geometry,
        energy_kcal_mol=_parse_energy(lines[cursor + 1], factor),
    )
    return conformer, end


def _parse_energy(comment: str, factor: float) -> float | None:
    """Read the energy CREST wrote on a comment line, if it wrote one."""
    tokens = comment.split()
    if not tokens:
        return None
    try:
        value = float(tokens[0])
    except ValueError:
        return None
    return value * factor


__all__ = [
    "HARTREE_TO_KCAL_MOL",
    "CrestConformerSearch",
    "CrestRun",
    "CrestRunner",
    "parse_crest_ensemble",
]
