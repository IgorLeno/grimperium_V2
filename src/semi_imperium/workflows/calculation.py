"""The Calculate area: a tabular session, reuse review and execution.

The area is one table of molecules. Each row carries its own identity,
its own charge and multiplicity, its own CREST switch and its own set of
independently requested Hamiltonians, because AM1, PM3 and PM7 are three
separate requests rather than three settings of one request.

Nothing here runs a program. Execution is reached through the
:class:`CalculationExecutor` boundary, so the whole flow — resolution,
reuse review, preparation, persistence and state transitions — is driven
in tests with an in-memory double.

Reuse is the part that has to be exact. After identity resolution the
review asks the store, per molecule *and per Hamiltonian*, whether a
compatible local result already exists under the same signature. What
happens next is always an explicit choice:

* :attr:`ReuseAction.REUSE_EXISTING` computes nothing;
* :attr:`ReuseAction.CALCULATE_MISSING` computes only the Hamiltonians
  that have no compatible result;
* :attr:`ReuseAction.RECALCULATE_ALL` computes everything again under a
  new run id, which means the reusable records stay on disk untouched.
"""

from __future__ import annotations

import uuid
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Protocol

from semi_imperium import __version__ as SEMI_IMPERIUM_VERSION
from semi_imperium.domain import (
    DEFAULT_REUSABLE_STATES,
    VERIFIED_ONLY_REUSABLE_STATES,
    CalculationRecord,
    CalculationResultData,
    CalculationSignature,
    CalculationState,
    EffectiveConfiguration,
    MolecularIdentity,
    MoleculeInputType,
    RunRecord,
    RunState,
    ScientificProvenance,
    VerificationOutcome,
    utc_now,
)
from semi_imperium.molecules.service import (
    MoleculeResolutionService,
    ResolutionOutcome,
    ResolutionStatus,
)
from semi_imperium.molecules.validation import ComputationalPath
from semi_imperium.mopac.models import SUPPORTED_HAMILTONIANS
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.settings import SemiImperiumSettings

#: What a freshly added row asks for until the user says otherwise.
DEFAULT_HAMILTONIANS: tuple[str, ...] = ("PM7",)


class CalculateSessionError(ValueError):
    """Raised when an edit does not apply to the row it names."""


def normalize_hamiltonians(values: Iterable[str]) -> tuple[str, ...]:
    """Return the requested Hamiltonians, deduplicated in canonical order."""
    requested = {value.strip().upper() for value in values}
    unsupported = sorted(requested - set(SUPPORTED_HAMILTONIANS))
    if unsupported:
        raise CalculateSessionError(
            f"Unsupported Hamiltonian(s): {', '.join(unsupported)}. "
            "Choose AM1, PM3 or PM7."
        )
    return tuple(name for name in SUPPORTED_HAMILTONIANS if name in requested)


@dataclass
class MoleculeEntry:
    """One editable row of the Calculate table."""

    entry_id: str
    raw_input: str
    input_type: MoleculeInputType = MoleculeInputType.CHEMICAL_NAME
    charge: int | None = None
    multiplicity: int = 1
    crest_enabled: bool = True
    hamiltonians: tuple[str, ...] = DEFAULT_HAMILTONIANS
    selected: bool = True
    outcome: ResolutionOutcome | None = None

    def __post_init__(self) -> None:
        self.hamiltonians = normalize_hamiltonians(self.hamiltonians)
        if self.multiplicity < 1:
            raise CalculateSessionError(
                f"Multiplicity must be >= 1, got {self.multiplicity}"
            )

    @property
    def identity(self) -> MolecularIdentity | None:
        """The resolved identity, or ``None`` while the row is not resolved."""
        return None if self.outcome is None else self.outcome.identity

    @property
    def status(self) -> str:
        """Explicit row state; ``"draft"`` means resolution never ran."""
        return "draft" if self.outcome is None else self.outcome.status.value

    @property
    def is_resolved(self) -> bool:
        """Whether this row currently carries a usable molecular identity."""
        return self.outcome is not None and (
            self.outcome.status is ResolutionStatus.RESOLVED
        )

    @property
    def label(self) -> str:
        """Best available human label, falling back to the raw input."""
        identity = self.identity
        if identity is not None:
            return identity.display_name or identity.resolved_name or self.raw_input
        return self.raw_input

    @property
    def message(self) -> str | None:
        """Why this row is not calculable yet, when that applies."""
        return None if self.outcome is None else self.outcome.message

    @property
    def computational_path(self) -> ComputationalPath:
        """The route this row asks for, used to validate its structure early."""
        return ComputationalPath(
            hamiltonians=self.hamiltonians or SUPPORTED_HAMILTONIANS,
            crest_enabled=self.crest_enabled,
        )


class CalculateSession:
    """The mutable Calculate table: one or many molecules, edited in place."""

    def __init__(
        self,
        *,
        resolution_service: MoleculeResolutionService | None = None,
        default_hamiltonians: Sequence[str] = DEFAULT_HAMILTONIANS,
        default_crest_enabled: bool = True,
    ) -> None:
        self._resolution_service = resolution_service
        self._entries: list[MoleculeEntry] = []
        self._counter = 0
        self.default_hamiltonians = normalize_hamiltonians(default_hamiltonians)
        self.default_crest_enabled = default_crest_enabled

    # ------------------------------------------------------------------
    # Resolution service
    # ------------------------------------------------------------------

    @property
    def resolution_service(self) -> MoleculeResolutionService:
        """The resolver, built against PubChem only when first needed."""
        if self._resolution_service is None:
            self._resolution_service = MoleculeResolutionService.with_pubchem()
        return self._resolution_service

    # ------------------------------------------------------------------
    # Table contents
    # ------------------------------------------------------------------

    @property
    def entries(self) -> tuple[MoleculeEntry, ...]:
        """Every row still in the table, in insertion order."""
        return tuple(self._entries)

    @property
    def selected_entries(self) -> tuple[MoleculeEntry, ...]:
        """The rows bulk operations and the review currently apply to."""
        return tuple(entry for entry in self._entries if entry.selected)

    @property
    def unresolved_entries(self) -> tuple[MoleculeEntry, ...]:
        """Rows that still need attention before they can be calculated."""
        return tuple(entry for entry in self._entries if not entry.is_resolved)

    def __len__(self) -> int:
        return len(self._entries)

    def get(self, entry_id: str) -> MoleculeEntry:
        """Return one row.

        Raises:
            CalculateSessionError: If no row carries ``entry_id``.
        """
        for entry in self._entries:
            if entry.entry_id == entry_id:
                return entry
        raise CalculateSessionError(f"No molecule row with id {entry_id!r}")

    def add(
        self,
        raw_input: str,
        input_type: MoleculeInputType | str = MoleculeInputType.CHEMICAL_NAME,
        *,
        charge: int | None = None,
        multiplicity: int = 1,
        crest_enabled: bool | None = None,
        hamiltonians: Sequence[str] | None = None,
    ) -> MoleculeEntry:
        """Append one molecule, by chemical name or by SMILES."""
        if not raw_input.strip():
            raise CalculateSessionError("A molecule needs a chemical name or SMILES")
        self._counter += 1
        entry = MoleculeEntry(
            entry_id=f"row-{self._counter}",
            raw_input=raw_input.strip(),
            input_type=MoleculeInputType(input_type),
            charge=charge,
            multiplicity=multiplicity,
            crest_enabled=(
                self.default_crest_enabled if crest_enabled is None else crest_enabled
            ),
            hamiltonians=(
                self.default_hamiltonians
                if hamiltonians is None
                else normalize_hamiltonians(hamiltonians)
            ),
        )
        self._entries.append(entry)
        return entry

    def add_many(
        self,
        raw_inputs: Iterable[str],
        input_type: MoleculeInputType | str = MoleculeInputType.CHEMICAL_NAME,
    ) -> tuple[MoleculeEntry, ...]:
        """Append several molecules at once, one per non-empty value."""
        return tuple(
            self.add(value, input_type)
            for value in raw_inputs
            if str(value).strip()
        )

    def remove(self, entry_id: str) -> MoleculeEntry:
        """Drop one row from the table."""
        entry = self.get(entry_id)
        self._entries.remove(entry)
        return entry

    def clear(self) -> None:
        """Empty the table without touching anything already persisted."""
        self._entries.clear()

    # ------------------------------------------------------------------
    # Per-row editing
    # ------------------------------------------------------------------

    def set_input(
        self,
        entry_id: str,
        raw_input: str,
        input_type: MoleculeInputType | str | None = None,
    ) -> MoleculeEntry:
        """Replace a row's identity input, discarding its stale resolution."""
        if not raw_input.strip():
            raise CalculateSessionError("A molecule needs a chemical name or SMILES")
        entry = self.get(entry_id)
        entry.raw_input = raw_input.strip()
        if input_type is not None:
            entry.input_type = MoleculeInputType(input_type)
        entry.outcome = None
        return entry

    def set_charge(self, entry_id: str, charge: int | None) -> MoleculeEntry:
        """Set the row's charge; the identity has to be resolved again."""
        entry = self.get(entry_id)
        entry.charge = charge
        entry.outcome = None
        return entry

    def set_multiplicity(self, entry_id: str, multiplicity: int) -> MoleculeEntry:
        """Set the row's spin multiplicity; forces a fresh resolution."""
        if multiplicity < 1:
            raise CalculateSessionError(
                f"Multiplicity must be >= 1, got {multiplicity}"
            )
        entry = self.get(entry_id)
        entry.multiplicity = multiplicity
        entry.outcome = None
        return entry

    def set_crest(self, entry_id: str, enabled: bool) -> MoleculeEntry:
        """Enable or disable the conformer search for one molecule."""
        entry = self.get(entry_id)
        entry.crest_enabled = enabled
        return entry

    def set_hamiltonians(
        self, entry_id: str, hamiltonians: Sequence[str]
    ) -> MoleculeEntry:
        """Replace the independent Hamiltonian requests of one molecule."""
        entry = self.get(entry_id)
        entry.hamiltonians = normalize_hamiltonians(hamiltonians)
        return entry

    def toggle_hamiltonian(self, entry_id: str, hamiltonian: str) -> MoleculeEntry:
        """Add or drop one Hamiltonian request for one molecule."""
        entry = self.get(entry_id)
        name = normalize_hamiltonians([hamiltonian])[0]
        requested = set(entry.hamiltonians)
        requested.symmetric_difference_update({name})
        entry.hamiltonians = normalize_hamiltonians(requested)
        return entry

    # ------------------------------------------------------------------
    # Selection and bulk operations
    # ------------------------------------------------------------------

    def set_selected(self, entry_id: str, selected: bool) -> MoleculeEntry:
        """Include or exclude one row from bulk operations and the review."""
        entry = self.get(entry_id)
        entry.selected = selected
        return entry

    def toggle_selection(self, entry_id: str) -> MoleculeEntry:
        """Flip one row's selection mark."""
        entry = self.get(entry_id)
        entry.selected = not entry.selected
        return entry

    def select_all(self) -> int:
        """Select every row; returns how many rows are now selected."""
        for entry in self._entries:
            entry.selected = True
        return len(self._entries)

    def clear_selection(self) -> int:
        """Deselect every row; returns how many rows were deselected."""
        cleared = sum(1 for entry in self._entries if entry.selected)
        for entry in self._entries:
            entry.selected = False
        return cleared

    def bulk_set_crest(self, enabled: bool) -> int:
        """Apply one CREST decision to every selected row."""
        targets = self.selected_entries
        for entry in targets:
            entry.crest_enabled = enabled
        return len(targets)

    def bulk_set_hamiltonian(self, hamiltonian: str, requested: bool) -> int:
        """Add or drop one Hamiltonian across every selected row."""
        name = normalize_hamiltonians([hamiltonian])[0]
        targets = self.selected_entries
        for entry in targets:
            current = set(entry.hamiltonians)
            if requested:
                current.add(name)
            else:
                current.discard(name)
            entry.hamiltonians = normalize_hamiltonians(current)
        return len(targets)

    def bulk_set_hamiltonians(self, hamiltonians: Sequence[str]) -> int:
        """Replace the Hamiltonian requests of every selected row."""
        normalized = normalize_hamiltonians(hamiltonians)
        targets = self.selected_entries
        for entry in targets:
            entry.hamiltonians = normalized
        return len(targets)

    # ------------------------------------------------------------------
    # Identity resolution
    # ------------------------------------------------------------------

    def resolve(self, entry_id: str) -> MoleculeEntry:
        """Resolve one row through the resolver-agnostic molecule service."""
        entry = self.get(entry_id)
        entry.outcome = self.resolution_service.resolve(
            entry.raw_input,
            entry.input_type,
            path=entry.computational_path,
            charge=entry.charge,
            multiplicity=entry.multiplicity,
        )
        return entry

    def resolve_all(self, *, only_selected: bool = True) -> tuple[MoleculeEntry, ...]:
        """Resolve every row that has no current outcome yet."""
        targets = self.selected_entries if only_selected else self.entries
        for entry in targets:
            if entry.outcome is None:
                self.resolve(entry.entry_id)
        return targets

    def select_candidate(
        self, entry_id: str, resolver_identifier: str
    ) -> MoleculeEntry:
        """Resolve an ambiguous row by naming one resolver candidate."""
        entry = self.get(entry_id)
        if entry.outcome is None:
            raise CalculateSessionError(
                f"Row {entry_id!r} has not been resolved yet"
            )
        entry.outcome = self.resolution_service.select_candidate(
            entry.outcome,
            resolver_identifier,
            path=entry.computational_path,
            charge=entry.charge,
            multiplicity=entry.multiplicity,
        )
        return entry

    def enter_manual_smiles(self, entry_id: str, smiles: str) -> MoleculeEntry:
        """Recover an unresolved name with an explicit structure."""
        entry = self.get(entry_id)
        if entry.outcome is None:
            raise CalculateSessionError(
                f"Row {entry_id!r} has not been resolved yet"
            )
        entry.outcome = self.resolution_service.enter_manual_smiles(
            entry.outcome,
            smiles,
            path=entry.computational_path,
            charge=entry.charge,
            multiplicity=entry.multiplicity,
        )
        return entry


# ---------------------------------------------------------------------------
# Reuse review
# ---------------------------------------------------------------------------


class ReuseAction(str, Enum):
    """What the user decided to do about compatible local results."""

    REUSE_EXISTING = "reuse_existing"
    """Adopt every compatible stored result and compute nothing."""

    CALCULATE_MISSING = "calculate_missing"
    """Compute only the Hamiltonians with no compatible stored result."""

    RECALCULATE_ALL = "recalculate_all"
    """Compute everything again; stored results are kept, never overwritten."""


@dataclass(frozen=True)
class HamiltonianPlan:
    """One molecule under one Hamiltonian, with its reuse verdict."""

    hamiltonian: str
    configuration: EffectiveConfiguration
    signature: CalculationSignature
    existing: CalculationRecord | None
    reason: str

    @property
    def is_reusable(self) -> bool:
        """Whether a compatible local result already answers this request."""
        return self.existing is not None

    @property
    def state_label(self) -> str:
        """Explicit reuse state; never a blank cell."""
        if self.existing is None:
            return "missing"
        return f"reusable ({self.existing.state.value})"


@dataclass(frozen=True)
class MoleculePlan:
    """The reuse verdict for one resolved row across its Hamiltonians."""

    entry_id: str
    label: str
    identity: MolecularIdentity
    crest_enabled: bool
    hamiltonian_plans: tuple[HamiltonianPlan, ...]

    @property
    def reusable(self) -> tuple[HamiltonianPlan, ...]:
        """Requests already answered by a compatible stored result."""
        return tuple(plan for plan in self.hamiltonian_plans if plan.is_reusable)

    @property
    def missing(self) -> tuple[HamiltonianPlan, ...]:
        """Requests with no compatible stored result."""
        return tuple(plan for plan in self.hamiltonian_plans if not plan.is_reusable)


@dataclass(frozen=True)
class BlockedEntry:
    """A row that cannot enter the plan, with the reason spelled out."""

    entry_id: str
    label: str
    status: str
    message: str


@dataclass(frozen=True)
class CalculatePlan:
    """What the Calculate area proposes to do, before anything is written."""

    molecules: tuple[MoleculePlan, ...] = ()
    blocked: tuple[BlockedEntry, ...] = ()

    @property
    def requested_count(self) -> int:
        """Total molecule/Hamiltonian pairs the selection asks for."""
        return sum(len(plan.hamiltonian_plans) for plan in self.molecules)

    @property
    def reusable_count(self) -> int:
        """Pairs a compatible local result can answer without computing."""
        return sum(len(plan.reusable) for plan in self.molecules)

    @property
    def missing_count(self) -> int:
        """Pairs that have no compatible local result."""
        return sum(len(plan.missing) for plan in self.molecules)

    @property
    def available_actions(self) -> tuple[ReuseAction, ...]:
        """The choices this plan actually offers, in decision order."""
        if not self.requested_count:
            return ()
        actions: list[ReuseAction] = []
        if self.reusable_count:
            actions.append(ReuseAction.REUSE_EXISTING)
        if self.missing_count:
            actions.append(ReuseAction.CALCULATE_MISSING)
        actions.append(ReuseAction.RECALCULATE_ALL)
        return tuple(actions)

    def plans_for(self, action: ReuseAction) -> tuple[HamiltonianPlan, ...]:
        """The pairs ``action`` would actually execute."""
        if action is ReuseAction.REUSE_EXISTING:
            return ()
        every = tuple(
            plan for molecule in self.molecules for plan in molecule.hamiltonian_plans
        )
        if action is ReuseAction.RECALCULATE_ALL:
            return every
        return tuple(plan for plan in every if not plan.is_reusable)

    def action_label(self, action: ReuseAction) -> str:
        """One-line description of what choosing ``action`` would do."""
        if action is ReuseAction.REUSE_EXISTING:
            return f"Reuse {self.reusable_count} compatible local result(s)"
        if action is ReuseAction.CALCULATE_MISSING:
            return (
                f"Calculate {self.missing_count} missing Hamiltonian(s), "
                f"reusing {self.reusable_count}"
            )
        return (
            f"Recalculate all {self.requested_count} (stored results are kept)"
        )


def review(
    session: CalculateSession,
    *,
    store: SemiImperiumStore,
    settings: SemiImperiumSettings,
    only_selected: bool = True,
    accepted_states: frozenset[CalculationState] | None = None,
) -> CalculatePlan:
    """Resolve the selected rows, then ask the store what already exists.

    Reuse is decided per molecule *and per Hamiltonian*: asking for PM7
    when only AM1 was computed before is a missing Hamiltonian, not a
    reusable result. Unless the caller supplies ``accepted_states``, a
    configuration that requires a minimum accepts verified records only;
    looser policies retain the general result-state reuse contract.
    """
    session.resolve_all(only_selected=only_selected)
    targets = session.selected_entries if only_selected else session.entries

    molecules: list[MoleculePlan] = []
    blocked: list[BlockedEntry] = []
    for entry in targets:
        identity = entry.identity
        if identity is None:
            blocked.append(
                BlockedEntry(
                    entry_id=entry.entry_id,
                    label=entry.label,
                    status=entry.status,
                    message=entry.message or "molecule identity is not resolved",
                )
            )
            continue
        if not entry.hamiltonians:
            blocked.append(
                BlockedEntry(
                    entry_id=entry.entry_id,
                    label=entry.label,
                    status="no_hamiltonian",
                    message="Request at least one of AM1, PM3 or PM7",
                )
            )
            continue

        hamiltonian_plans: list[HamiltonianPlan] = []
        for hamiltonian in entry.hamiltonians:
            configuration = settings.configuration_for(
                hamiltonian, crest_enabled=entry.crest_enabled
            )
            signature = configuration.signature()
            reusable_states = accepted_states
            if reusable_states is None:
                reusable_states = (
                    VERIFIED_ONLY_REUSABLE_STATES
                    if configuration.verification.requires_minimum
                    else DEFAULT_REUSABLE_STATES
                )
            decision = store.find_reusable(
                identity, signature, accepted_states=reusable_states
            )
            hamiltonian_plans.append(
                HamiltonianPlan(
                    hamiltonian=hamiltonian,
                    configuration=configuration,
                    signature=signature,
                    existing=decision.record,
                    reason=decision.reason,
                )
            )
        molecules.append(
            MoleculePlan(
                entry_id=entry.entry_id,
                label=entry.label,
                identity=identity,
                crest_enabled=entry.crest_enabled,
                hamiltonian_plans=tuple(hamiltonian_plans),
            )
        )

    return CalculatePlan(molecules=tuple(molecules), blocked=tuple(blocked))


# ---------------------------------------------------------------------------
# Preparation and persistence
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PreparedCalculation:
    """One pending calculation, tied back to the row that asked for it."""

    entry_id: str
    hamiltonian: str
    configuration: EffectiveConfiguration
    record: CalculationRecord


@dataclass(frozen=True)
class PreparedRun:
    """One run manifest plus the pending calculations it owns.

    A run holds exactly one effective configuration, so requests that
    differ in Hamiltonian or CREST usage land in separate runs; that is
    what keeps a run's results from being a mix of incomparable settings.
    """

    run: RunRecord
    calculations: tuple[PreparedCalculation, ...]

    @property
    def configuration(self) -> EffectiveConfiguration:
        """The effective configuration captured for this run."""
        return self.run.configuration


@dataclass(frozen=True)
class PreparedWork:
    """Everything one explicit reuse decision turned into."""

    action: ReuseAction
    runs: tuple[PreparedRun, ...] = ()
    reused: tuple[CalculationRecord, ...] = ()

    @property
    def pending_count(self) -> int:
        """How many calculations this decision will actually execute."""
        return sum(len(prepared.calculations) for prepared in self.runs)


def build_provenance(
    settings: SemiImperiumSettings,
    *,
    hostname: str | None = None,
    notes: str | None = None,
    program_versions: dict[str, str] | None = None,
) -> ScientificProvenance:
    """Capture who is producing the numbers of the run about to start."""
    return ScientificProvenance(
        method_id=settings.method_id,
        method_version=settings.method_version,
        property_id=settings.property_id,
        semi_imperium_version=SEMI_IMPERIUM_VERSION,
        grimperium_version=_grimperium_version(),
        program_versions=dict(program_versions or {}),
        hostname=hostname,
        source="semi_imperium",
        notes=notes,
    )


def prepare_runs(
    plan: CalculatePlan,
    action: ReuseAction,
    *,
    provenance: ScientificProvenance,
    run_id_factory: Callable[[], str] | None = None,
    label: str | None = None,
) -> PreparedWork:
    """Turn one explicit decision into runs and pending calculation records.

    ``RECALCULATE_ALL`` deliberately uses a fresh run id, so the new
    records get new calculation ids and land beside the stored ones under
    the same reuse key. A reusable result is therefore never overwritten,
    only superseded by something the user asked for explicitly.
    """
    executable = {id(item) for item in plan.plans_for(action)}
    reused = tuple(
        hamiltonian_plan.existing
        for molecule in plan.molecules
        for hamiltonian_plan in molecule.hamiltonian_plans
        if hamiltonian_plan.existing is not None
        and action is not ReuseAction.RECALCULATE_ALL
    )
    if not executable:
        return PreparedWork(action=action, runs=(), reused=reused)

    base_run_id = (run_id_factory or _default_run_id)()

    # Group by signature: one run manifest describes exactly one
    # effective configuration, so incomparable requests cannot share it.
    grouped: dict[str, list[tuple[MoleculePlan, HamiltonianPlan]]] = {}
    order: list[str] = []
    for molecule in plan.molecules:
        for hamiltonian_plan in molecule.hamiltonian_plans:
            if id(hamiltonian_plan) not in executable:
                continue
            digest = hamiltonian_plan.signature.digest
            if digest not in grouped:
                grouped[digest] = []
                order.append(digest)
            grouped[digest].append((molecule, hamiltonian_plan))

    runs: list[PreparedRun] = []
    for index, digest in enumerate(order, start=1):
        members = grouped[digest]
        configuration = members[0][1].configuration
        run_id = f"{base_run_id}-{index:02d}"
        run = RunRecord(
            run_id=run_id,
            configuration=configuration,
            provenance=provenance,
            state=RunState.PENDING,
            label=label,
            molecule_ids=tuple(
                dict.fromkeys(
                    molecule.identity.molecule_id for molecule, _ in members
                )
            ),
        )
        calculations = tuple(
            PreparedCalculation(
                entry_id=molecule.entry_id,
                hamiltonian=hamiltonian_plan.hamiltonian,
                configuration=hamiltonian_plan.configuration,
                record=CalculationRecord(
                    run_id=run_id,
                    molecule=molecule.identity,
                    signature=hamiltonian_plan.signature,
                    provenance=provenance,
                    state=CalculationState.PENDING,
                    verification=VerificationOutcome.NOT_REQUESTED,
                ),
            )
            for molecule, hamiltonian_plan in members
        )
        runs.append(PreparedRun(run=run, calculations=calculations))

    return PreparedWork(action=action, runs=tuple(runs), reused=reused)


def save_prepared_runs(
    work: PreparedWork, store: SemiImperiumStore
) -> tuple[Path, ...]:
    """Persist run manifests and their pending calculations before executing.

    Writing the pending records first is what makes an interrupted run
    legible afterwards: the work that was registered is on disk in state
    ``pending`` rather than missing.
    """
    paths: list[Path] = []
    for prepared in work.runs:
        paths.append(store.save_run(prepared.run))
        for calculation in prepared.calculations:
            paths.append(store.save_calculation(calculation.record))
    return tuple(paths)


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ExecutionRequest:
    """Everything an executor needs to compute one pending calculation."""

    run_id: str
    calculation_id: str
    identity: MolecularIdentity
    configuration: EffectiveConfiguration
    hamiltonian: str


@dataclass(frozen=True)
class ExecutionOutcome:
    """The terminal state one executed calculation reached."""

    state: CalculationState
    verification: VerificationOutcome = VerificationOutcome.NOT_REQUESTED
    result: CalculationResultData | None = None
    error_message: str | None = None


class CalculationExecutor(Protocol):
    """The boundary between the Calculate area and the scientific pipeline."""

    def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
        """Compute one molecule under one Hamiltonian."""
        ...


@dataclass(frozen=True)
class ExecutionSummary:
    """What actually happened, kept separate from what was reused."""

    action: ReuseAction
    runs: tuple[RunRecord, ...] = ()
    completed: tuple[CalculationRecord, ...] = ()
    failed: tuple[CalculationRecord, ...] = ()
    reused: tuple[CalculationRecord, ...] = ()

    @property
    def executed_count(self) -> int:
        """How many calculations were run, successfully or not."""
        return len(self.completed) + len(self.failed)

    def describe(self) -> str:
        """One line the Calculate area can show without interpreting it."""
        return (
            f"{self.action.value}: {len(self.completed)} completed, "
            f"{len(self.failed)} failed, {len(self.reused)} reused"
        )


def execute_prepared(
    work: PreparedWork,
    *,
    store: SemiImperiumStore,
    executor: CalculationExecutor,
    on_progress: Callable[[PreparedCalculation, CalculationRecord], None] | None = None,
) -> ExecutionSummary:
    """Execute every pending calculation, persisting each state transition.

    Reused records are reported untouched: they are never re-saved, so a
    finished scientific result cannot be rewritten by a later session.
    """
    completed: list[CalculationRecord] = []
    failed: list[CalculationRecord] = []
    finished_runs: list[RunRecord] = []

    for prepared in work.runs:
        run = prepared.run.transition_to(RunState.RUNNING)
        store.save_run(run)

        run_failures = 0
        for calculation in prepared.calculations:
            record = calculation.record.transition_to(CalculationState.RUNNING)
            store.save_calculation(record)

            try:
                outcome = executor.execute(
                    ExecutionRequest(
                        run_id=record.run_id,
                        calculation_id=record.calculation_id,
                        identity=record.molecule,
                        configuration=calculation.configuration,
                        hamiltonian=calculation.hamiltonian,
                    )
                )
            except Exception as exc:  # a crashing executor is a failed result
                outcome = ExecutionOutcome(
                    state=CalculationState.FAILED,
                    verification=VerificationOutcome.FAILED,
                    error_message=f"{type(exc).__name__}: {exc}",
                )

            record = record.transition_to(
                outcome.state,
                verification=outcome.verification,
                result=outcome.result,
                error_message=outcome.error_message,
            )
            store.save_calculation(record)

            if record.state is CalculationState.FAILED:
                run_failures += 1
                failed.append(record)
            else:
                completed.append(record)
            if on_progress is not None:
                on_progress(calculation, record)

        total = len(prepared.calculations)
        if run_failures == 0:
            finished = run.transition_to(RunState.COMPLETED)
        elif run_failures == total:
            finished = run.transition_to(
                RunState.FAILED,
                error_message=f"All {total} calculation(s) in this run failed",
            )
        else:
            finished = run.transition_to(RunState.PARTIAL)
        store.save_run(finished)
        finished_runs.append(finished)

    return ExecutionSummary(
        action=work.action,
        runs=tuple(finished_runs),
        completed=tuple(completed),
        failed=tuple(failed),
        reused=work.reused,
    )


def _default_run_id() -> str:
    """Return a fresh, path-safe run id."""
    stamp = utc_now().strftime("%Y%m%dT%H%M%S")
    return f"run-{stamp}-{uuid.uuid4().hex[:6]}"


def _grimperium_version() -> str:
    """Report the host Grimperium version without importing its CLI."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("grimperium")
    except PackageNotFoundError:  # pragma: no cover - source checkout without install
        return "unknown"


__all__ = [
    "DEFAULT_HAMILTONIANS",
    "BlockedEntry",
    "CalculatePlan",
    "CalculateSession",
    "CalculateSessionError",
    "CalculationExecutor",
    "ExecutionOutcome",
    "ExecutionRequest",
    "ExecutionSummary",
    "HamiltonianPlan",
    "MoleculeEntry",
    "MoleculePlan",
    "PreparedCalculation",
    "PreparedRun",
    "PreparedWork",
    "ReuseAction",
    "build_provenance",
    "execute_prepared",
    "normalize_hamiltonians",
    "prepare_runs",
    "review",
    "save_prepared_runs",
]
