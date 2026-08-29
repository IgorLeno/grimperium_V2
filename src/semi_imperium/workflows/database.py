"""The Database area: an operational view over what was actually computed.

This is a report, never a second source of truth. Every field shown here
is read back from the persisted run manifests and calculation records,
and the Hamiltonian, CREST usage and selection strategy of a result come
from the effective configuration captured in *its own* run — not from
today's defaults. That is what makes the summary safe to read months
later, after the defaults have moved on.

Every cell is explicit. A molecule that was never computed under PM3
shows ``not calculated`` rather than a blank, and a finished result
always shows both its lifecycle state and its verification outcome.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

from semi_imperium.domain import (
    CalculationRecord,
    CalculationState,
    EffectiveConfiguration,
    MolecularIdentity,
    RunRecord,
    VerificationOutcome,
)
from semi_imperium.mopac.models import SUPPORTED_HAMILTONIANS
from semi_imperium.persistence import SemiImperiumStore

#: Shown wherever a molecule has no calculation for a Hamiltonian.
NOT_CALCULATED = "not calculated"

#: Shown when a record's run manifest cannot be read back.
UNKNOWN = "unknown"


@dataclass(frozen=True)
class HamiltonianCell:
    """One molecule's latest outcome under one Hamiltonian."""

    hamiltonian: str
    state: CalculationState
    verification: VerificationOutcome
    energy_hof_kcal_mol: float | None
    run_id: str
    calculation_id: str
    signature_digest: str
    crest_used: bool | None
    selection_strategy: str
    selection_experimental: bool
    updated_at: datetime

    @property
    def status_label(self) -> str:
        """Lifecycle state and verification outcome, always both."""
        return f"{self.state.value} / {self.verification.value}"

    @property
    def energy_label(self) -> str:
        """Heat of formation, or an explicit reason there is no number."""
        if self.energy_hof_kcal_mol is None:
            return "no value"
        return f"{self.energy_hof_kcal_mol:.2f}"

    @property
    def crest_label(self) -> str:
        """Whether the conformer search was used for this exact result."""
        if self.crest_used is None:
            return UNKNOWN
        return "yes" if self.crest_used else "no"


@dataclass(frozen=True)
class MoleculeSummary:
    """One row of the Database table: a molecule across AM1, PM3 and PM7."""

    molecule: MolecularIdentity
    cells: tuple[HamiltonianCell, ...]
    calculation_count: int
    last_run_id: str
    last_updated: datetime

    @property
    def label(self) -> str:
        """Best human label, falling back to the canonical SMILES."""
        return (
            self.molecule.display_name
            or self.molecule.resolved_name
            or self.molecule.canonical_smiles
        )

    def cell(self, hamiltonian: str) -> HamiltonianCell | None:
        """Return the latest cell for one Hamiltonian, if there is one."""
        name = hamiltonian.strip().upper()
        for cell in self.cells:
            if cell.hamiltonian == name:
                return cell
        return None

    def status_label(self, hamiltonian: str) -> str:
        """Explicit status for one Hamiltonian column."""
        cell = self.cell(hamiltonian)
        return NOT_CALCULATED if cell is None else cell.status_label

    @property
    def crest_label(self) -> str:
        """Whether CREST was used, across the results this molecule has."""
        used = {cell.crest_used for cell in self.cells}
        if not used:
            return NOT_CALCULATED
        if used == {True}:
            return "yes"
        if used == {False}:
            return "no"
        return "mixed"

    @property
    def selection_label(self) -> str:
        """The selection strategy behind this molecule's results."""
        strategies = sorted({cell.selection_strategy for cell in self.cells})
        if not strategies:
            return NOT_CALCULATED
        label = ", ".join(strategies)
        if any(cell.selection_experimental for cell in self.cells):
            label += " (experimental)"
        return label


@dataclass(frozen=True)
class CalculationDetail:
    """One stored calculation with the run context that explains it."""

    record: CalculationRecord
    run: RunRecord | None

    @property
    def configuration(self) -> EffectiveConfiguration | None:
        """The effective configuration captured for this calculation's run."""
        return None if self.run is None else self.run.configuration

    @property
    def hamiltonian(self) -> str:
        """The Hamiltonian this calculation used, per its own run."""
        configuration = self.configuration
        if configuration is None:
            return UNKNOWN
        return configuration.semiempirical.hamiltonian

    @property
    def crest_used(self) -> bool | None:
        """Whether the conformer search ran, per this calculation's run."""
        configuration = self.configuration
        if configuration is None:
            return None
        return configuration.conformer_search.enabled

    @property
    def selection_strategy(self) -> str:
        """The conformer selection strategy captured for this run."""
        configuration = self.configuration
        if configuration is None:
            return UNKNOWN
        return configuration.conformer_selection.strategy

    @property
    def selection_experimental(self) -> bool:
        """Whether that strategy is still marked experimental."""
        configuration = self.configuration
        return (
            False
            if configuration is None
            else configuration.conformer_selection.is_experimental
        )

    @property
    def updated_at(self) -> datetime:
        """When this calculation last changed state."""
        timestamps = self.record.timestamps
        return timestamps.completed_at or timestamps.created_at


@dataclass(frozen=True)
class DatabaseSummary:
    """The operational picture of everything Semi-Imperium has computed."""

    molecules: tuple[MoleculeSummary, ...] = ()
    runs: tuple[RunRecord, ...] = ()
    calculation_count: int = 0

    @property
    def molecule_count(self) -> int:
        """How many distinct molecular identities are stored."""
        return len(self.molecules)

    @property
    def run_count(self) -> int:
        """How many runs produced those calculations."""
        return len(self.runs)

    @property
    def is_empty(self) -> bool:
        """Whether nothing has been computed yet."""
        return not self.molecules

    def state_counts(self) -> dict[str, int]:
        """Count the latest cells by lifecycle state, states with zero included."""
        counts = {state.value: 0 for state in CalculationState}
        for molecule in self.molecules:
            for cell in molecule.cells:
                counts[cell.state.value] += 1
        return counts

    def find(self, molecule_id: str) -> MoleculeSummary | None:
        """Return one summary row by molecule id."""
        for molecule in self.molecules:
            if molecule.molecule.molecule_id == molecule_id:
                return molecule
        return None


def build_summary(
    store: SemiImperiumStore,
    *,
    hamiltonians: tuple[str, ...] = SUPPORTED_HAMILTONIANS,
) -> DatabaseSummary:
    """Read the store and build the molecule/calculation summary.

    The latest calculation wins per molecule and Hamiltonian, so the row
    describes the current state of that pair rather than an arbitrary
    historical attempt. Older attempts stay on disk and remain reachable
    through :func:`molecule_detail`.
    """
    runs = store.list_runs()
    by_run = {run.run_id: run for run in runs}
    requested = tuple(name.strip().upper() for name in hamiltonians)

    grouped: dict[str, list[CalculationDetail]] = {}
    identities: dict[str, MolecularIdentity] = {}
    total = 0
    for record in store.list_calculations():
        total += 1
        molecule_id = record.molecule.molecule_id
        identities.setdefault(molecule_id, record.molecule)
        grouped.setdefault(molecule_id, []).append(
            CalculationDetail(record=record, run=by_run.get(record.run_id))
        )

    summaries: list[MoleculeSummary] = []
    for molecule_id, details in grouped.items():
        latest: dict[str, CalculationDetail] = {}
        for detail in details:
            name = detail.hamiltonian
            if name not in requested:
                continue
            current = latest.get(name)
            if current is None or detail.updated_at > current.updated_at:
                latest[name] = detail

        cells = tuple(
            _to_cell(latest[name]) for name in requested if name in latest
        )
        newest = max(details, key=lambda item: item.updated_at)
        summaries.append(
            MoleculeSummary(
                molecule=identities[molecule_id],
                cells=cells,
                calculation_count=len(details),
                last_run_id=newest.record.run_id,
                last_updated=newest.updated_at,
            )
        )

    summaries.sort(key=lambda item: item.last_updated, reverse=True)
    return DatabaseSummary(
        molecules=tuple(summaries),
        runs=tuple(runs),
        calculation_count=total,
    )


def molecule_detail(
    store: SemiImperiumStore, molecule: MolecularIdentity | str
) -> tuple[CalculationDetail, ...]:
    """Return every stored calculation for one molecule, newest first.

    This is the drill-down behind a summary row: it keeps superseded
    attempts visible instead of collapsing them into the latest cell.
    """
    molecule_id = (
        molecule if isinstance(molecule, str) else molecule.molecule_id
    )
    by_run = {run.run_id: run for run in store.list_runs()}
    details = [
        CalculationDetail(record=record, run=by_run.get(record.run_id))
        for record in store.list_calculations()
        if record.molecule.molecule_id == molecule_id
    ]
    details.sort(key=lambda item: item.updated_at, reverse=True)
    return tuple(details)


def _to_cell(detail: CalculationDetail) -> HamiltonianCell:
    """Project one calculation and its run onto a summary cell."""
    record = detail.record
    result = record.result
    return HamiltonianCell(
        hamiltonian=detail.hamiltonian,
        state=record.state,
        verification=record.verification,
        energy_hof_kcal_mol=None if result is None else result.energy_hof_kcal_mol,
        run_id=record.run_id,
        calculation_id=record.calculation_id,
        signature_digest=record.signature.digest,
        crest_used=detail.crest_used,
        selection_strategy=detail.selection_strategy,
        selection_experimental=detail.selection_experimental,
        updated_at=detail.updated_at,
    )


__all__ = [
    "NOT_CALCULATED",
    "UNKNOWN",
    "CalculationDetail",
    "DatabaseSummary",
    "HamiltonianCell",
    "MoleculeSummary",
    "build_summary",
    "molecule_detail",
]
