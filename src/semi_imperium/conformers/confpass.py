"""EXPERIMENTAL CONFPASS prioritization and its XYZ-to-SDF adaptation.

CONFPASS reads SDF, CREST writes XYZ, and the translation between them
is where a conformer ensemble is easiest to silently corrupt. The
adaptation here therefore preserves four things and is tested on all of
them: atom order, coordinates, connectivity, and the search provenance
that says where the structure came from.

Two boundaries are deliberate:

* the strategy is marked experimental, because CONFPASS prioritization
  has not been validated against this project's reference data;
* the PAS completeness classifier is recorded as advisory metadata and
  never as evidence, so it cannot justify a scientific selection.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

from semi_imperium.conformers.backends import (
    ConformerBackendError,
    ConfPassBackend,
    ConfPassCandidate,
    ConfPassRanking,
)
from semi_imperium.conformers.ensemble import (
    Conformer,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerSearchProvenance,
)
from semi_imperium.conformers.selection import (
    PAS_COMPLETENESS_LABEL_KEY,
    SelectionResult,
    require_strategy,
)
from semi_imperium.domain.configuration import ConformerSelectionSettings
from semi_imperium.domain.enums import ConformerSelectionStrategy

#: SD data field names written by :func:`to_sd_record`.
INDEX_FIELD = "CONFORMER_INDEX"
LABEL_FIELD = "CONFORMER_LABEL"
ENERGY_FIELD = "ENERGY_KCAL_MOL"
SOURCE_FIELD = "CONFORMER_SOURCE"
PROGRAM_FIELD = "SEARCH_PROGRAM"
VERSION_FIELD = "SEARCH_PROGRAM_VERSION"
RUN_FIELD = "RUN_ID"

_COUNTS_TAIL = "  0  0  0  0  0  0  0  0999 V2000"
_ATOM_TAIL = " 0  0  0  0  0  0  0  0  0  0  0  0"
_BOND_TAIL = "  0  0  0  0"
_END_TAG = "M  END"
_RECORD_TERMINATOR = "$$$$"


@dataclass(frozen=True)
class MoleculeTopology:
    """Connectivity shared by every conformer of one ensemble.

    Bonds address atoms by their position in the geometry, so a topology
    is only meaningful together with an ensemble that keeps that order.
    """

    atom_count: int
    bonds: tuple[tuple[int, int, int], ...] = ()
    """Zero-based ``(first_atom, second_atom, bond_order)`` triples."""

    def __post_init__(self) -> None:
        if self.atom_count < 1:
            raise ValueError(
                f"MoleculeTopology.atom_count must be >= 1, got {self.atom_count}"
            )
        seen: set[frozenset[int]] = set()
        for first, second, order in self.bonds:
            if first == second:
                raise ValueError(
                    f"MoleculeTopology bond {first}-{second} links atom {first} "
                    "to itself"
                )
            for atom in (first, second):
                if not 0 <= atom < self.atom_count:
                    raise ValueError(
                        f"MoleculeTopology bond atom {atom} is outside the "
                        f"molecule's {self.atom_count} atoms"
                    )
            if order < 1:
                raise ValueError(
                    f"MoleculeTopology bond {first}-{second} declares order "
                    f"{order}, which must be >= 1"
                )
            pair = frozenset({first, second})
            if pair in seen:
                raise ValueError(
                    f"MoleculeTopology declares bond {first}-{second} more " "than once"
                )
            seen.add(pair)


@dataclass(frozen=True)
class AdaptedStructure:
    """One SD record read back into geometry, connectivity and data."""

    title: str
    geometry: ConformerGeometry
    topology: MoleculeTopology
    data: dict[str, str] = field(default_factory=dict)


def to_sd_record(
    conformer: Conformer,
    topology: MoleculeTopology,
    provenance: ConformerSearchProvenance,
    *,
    molecule_id: str,
) -> str:
    """Render one conformer as an SD record CONFPASS can consume.

    Raises:
        ValueError: If ``topology`` describes a different atom count than
            the conformer's geometry, which would attach the wrong
            connectivity to the coordinates.
    """
    geometry = conformer.geometry
    if topology.atom_count != geometry.atom_count:
        raise ValueError(
            "Topology and geometry disagree on the molecule: "
            f"{topology.atom_count} atoms declared, {geometry.atom_count} given"
        )

    lines = [
        f"{molecule_id} {conformer.label}",
        f"  semi_imperium  {provenance.program}",
        f"conformer {conformer.index} of the {provenance.source.value} ensemble",
        f"{geometry.atom_count:>3}{len(topology.bonds):>3}{_COUNTS_TAIL}",
    ]
    for symbol, (x, y, z) in zip(geometry.elements, geometry.coordinates):
        lines.append(f"{x:>10.4f}{y:>10.4f}{z:>10.4f} {symbol:<3}{_ATOM_TAIL}")
    for first, second, order in topology.bonds:
        lines.append(f"{first + 1:>3}{second + 1:>3}{order:>3}{_BOND_TAIL}")
    lines.append(_END_TAG)

    data: dict[str, str] = {
        INDEX_FIELD: str(conformer.index),
        LABEL_FIELD: conformer.label,
        SOURCE_FIELD: provenance.source.value,
        PROGRAM_FIELD: provenance.program,
        VERSION_FIELD: provenance.program_version,
    }
    if conformer.energy_kcal_mol is not None:
        data[ENERGY_FIELD] = repr(conformer.energy_kcal_mol)
    if provenance.run_id is not None:
        data[RUN_FIELD] = provenance.run_id
    for key, value in data.items():
        lines.extend([f"> <{key}>", value, ""])
    lines.append(_RECORD_TERMINATOR)
    return "\n".join(lines)


def read_sd_record(text: str) -> AdaptedStructure:
    """Read back an SD record written by :func:`to_sd_record`.

    Reading is part of the contract on purpose: the adaptation claims to
    preserve order, coordinates, connectivity and provenance, and this
    is what makes that claim checkable instead of asserted.

    Raises:
        ConformerBackendError: If the record is not a readable molblock.
    """
    lines = text.splitlines()
    if len(lines) < 5:
        raise ConformerBackendError(
            "SD record is too short to contain a molblock",
            code="sdf_parse_failed",
        )

    counts = lines[3]
    try:
        atom_count = int(counts[0:3])
        bond_count = int(counts[3:6])
    except ValueError as exc:
        raise ConformerBackendError(
            f"SD record counts line is unreadable: {counts!r}",
            code="sdf_parse_failed",
        ) from exc

    atom_end = 4 + atom_count
    bond_end = atom_end + bond_count
    if bond_end >= len(lines):
        raise ConformerBackendError(
            f"SD record declares {atom_count} atoms and {bond_count} bonds "
            f"but holds only {len(lines)} lines",
            code="sdf_parse_failed",
        )

    elements: list[str] = []
    coordinates: list[tuple[float, float, float]] = []
    for line in lines[4:atom_end]:
        parts = line.split()
        if len(parts) < 4:
            raise ConformerBackendError(
                f"SD record atom line is unreadable: {line!r}",
                code="sdf_parse_failed",
            )
        try:
            values = [float(value) for value in parts[0:3]]
        except ValueError as exc:
            raise ConformerBackendError(
                f"SD record atom line holds unreadable coordinates: {line!r}",
                code="sdf_parse_failed",
            ) from exc
        coordinates.append((values[0], values[1], values[2]))
        elements.append(parts[3])

    bonds: list[tuple[int, int, int]] = []
    for line in lines[atom_end:bond_end]:
        try:
            first = int(line[0:3])
            second = int(line[3:6])
            order = int(line[6:9])
        except ValueError as exc:
            raise ConformerBackendError(
                f"SD record bond line is unreadable: {line!r}",
                code="sdf_parse_failed",
            ) from exc
        bonds.append((first - 1, second - 1, order))

    if lines[bond_end].strip() != _END_TAG:
        raise ConformerBackendError(
            f"SD record does not close its molblock with {_END_TAG!r}",
            code="sdf_parse_failed",
        )

    return AdaptedStructure(
        title=lines[0].strip(),
        geometry=ConformerGeometry(
            elements=tuple(elements),
            coordinates=tuple(coordinates),
        ),
        topology=MoleculeTopology(atom_count=atom_count, bonds=tuple(bonds)),
        data=_read_data_fields(lines[bond_end + 1 :]),
    )


def build_confpass_candidates(
    ensemble: ConformerEnsemble,
    topology: MoleculeTopology,
    *,
    molecule_id: str,
) -> tuple[ConfPassCandidate, ...]:
    """Adapt the whole ensemble to SDF, in the order the search wrote it.

    Every conformer is adapted, never a pre-filtered subset: CONFPASS is
    supposed to see the broad ensemble and decide, so any cut applied
    before this point would silently change what it is ranking.
    """
    return tuple(
        ConfPassCandidate(
            index=conformer.index,
            sd_record=to_sd_record(
                conformer,
                topology,
                ensemble.provenance,
                molecule_id=molecule_id,
            ),
            energy_kcal_mol=conformer.energy_kcal_mol,
        )
        for conformer in ensemble.conformers
    )


@dataclass(frozen=True)
class UnavailableConfPass:
    """Default CONFPASS backend: refuses instead of inventing a ranking."""

    reason: str = "no CONFPASS backend is configured"

    def prioritize(
        self,
        candidates: Sequence[ConfPassCandidate],
    ) -> Sequence[ConfPassRanking]:
        """Always fail, explaining why CONFPASS cannot run."""
        raise ConformerBackendError(
            f"CONFPASS prioritization was requested for {len(candidates)} "
            f"conformers, but {self.reason}",
            code="confpass_unavailable",
        )


@dataclass(frozen=True)
class ConfPassSelector:
    """EXPERIMENTAL: order the ensemble with CONFPASS, then keep N.

    The prioritizer receives the whole CREST ensemble before anything is
    discarded, so the ``top_n`` cut applies to CONFPASS's ordering rather
    than to a list something else already truncated.

    The PAS completeness class each ranking may carry is recorded in
    :attr:`SelectionResult.advisory_labels` only. It never enters the
    result's evidence, which :class:`SelectionResult` enforces.
    """

    backend: ConfPassBackend
    topology: MoleculeTopology
    molecule_id: str = "molecule"

    def select(
        self,
        ensemble: ConformerEnsemble,
        settings: ConformerSelectionSettings,
    ) -> SelectionResult:
        """Return the ``settings.top_n`` highest-priority conformers.

        Raises:
            ConformerBackendError: If CONFPASS is unavailable or returns
                a ranking that does not cover the ensemble exactly once.
        """
        strategy = ConformerSelectionStrategy.CONFPASS_PRIORITIZATION
        require_strategy(settings, strategy)

        candidates = build_confpass_candidates(
            ensemble,
            self.topology,
            molecule_id=self.molecule_id,
        )
        rankings = tuple(self.backend.prioritize(candidates))
        _validate_rankings(rankings, ensemble)

        by_index = {conformer.index: conformer for conformer in ensemble.conformers}
        ordered = tuple(
            by_index[ranking.index] for ranking in sorted(rankings, key=_priority_key)
        )
        advisory: dict[str, str] = {}
        for ranking in rankings:
            pas_class = ranking.pas_completeness_class
            if pas_class is None:
                continue
            label = by_index[ranking.index].label
            advisory[label] = f"{PAS_COMPLETENESS_LABEL_KEY}={pas_class}"

        return SelectionResult(
            strategy=strategy,
            selected=ordered[: settings.top_n],
            considered=ensemble.size,
            ranking_basis="confpass_priority",
            evidence=(
                "experimental_strategy",
                "confpass_priority_ordering",
                f"top_n={settings.top_n}",
            ),
            advisory_labels=advisory,
        )


def _priority_key(ranking: ConfPassRanking) -> tuple[int, int]:
    """Order by CONFPASS priority, breaking ties by ensemble order."""
    return (ranking.priority, ranking.index)


def _validate_rankings(
    rankings: tuple[ConfPassRanking, ...],
    ensemble: ConformerEnsemble,
) -> None:
    """Require a ranking for every conformer the ensemble contains."""
    ranked_indices = [ranking.index for ranking in rankings]
    expected = {conformer.index for conformer in ensemble.conformers}
    if len(set(ranked_indices)) != len(ranked_indices):
        raise ConformerBackendError(
            f"CONFPASS ranked a conformer twice: {sorted(ranked_indices)}",
            code="confpass_invalid_ranking",
        )
    if set(ranked_indices) != expected:
        raise ConformerBackendError(
            "CONFPASS must rank the whole ensemble before any selection: "
            f"ranked {sorted(ranked_indices)}, ensemble has {sorted(expected)}",
            code="confpass_invalid_ranking",
        )


def _read_data_fields(lines: Sequence[str]) -> dict[str, str]:
    """Read the ``> <NAME>`` / value pairs that follow the molblock."""
    data: dict[str, str] = {}
    cursor = 0
    while cursor < len(lines):
        header = lines[cursor].strip()
        if header in {"", _RECORD_TERMINATOR}:
            cursor += 1
            continue
        if not header.startswith("> <") or not header.endswith(">"):
            raise ConformerBackendError(
                f"SD record data header is unreadable: {header!r}",
                code="sdf_parse_failed",
            )
        cursor += 1
        if cursor >= len(lines):
            raise ConformerBackendError(
                f"SD record data field {header!r} carries no value",
                code="sdf_parse_failed",
            )
        name = header.removeprefix("> <").removesuffix(">")
        data[name] = lines[cursor].strip()
        cursor += 1
    return data


__all__ = [
    "AdaptedStructure",
    "ConfPassSelector",
    "MoleculeTopology",
    "UnavailableConfPass",
    "build_confpass_candidates",
    "read_sd_record",
    "to_sd_record",
]
