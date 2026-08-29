"""Contract for the Calculate area, from the table to executed results.

The tests are grouped by the three properties the area has to hold:

1. the table really is a table — one or many molecules, entered by name
   or SMILES, each with its own editable identity, CREST switch and
   independent AM1/PM3/PM7 requests, plus bulk select/clear;
2. after resolution, reuse is decided per molecule *and* per Hamiltonian,
   and the three choices it produces are reachable from the normal
   interface, not only from the internal API;
3. choosing one of them actually prepares runs, persists them and runs
   the execution boundary — and a reusable result is never overwritten.
"""

from __future__ import annotations

import io
import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
from rich.console import Console

from semi_imperium.domain import (
    CalculationRecord,
    CalculationResultData,
    CalculationState,
    EffectiveConfiguration,
    MolecularIdentity,
    MoleculeInputType,
    RunRecord,
    RunState,
    ScientificProvenance,
    VerificationOutcome,
)
from semi_imperium.molecules.service import ResolutionOutcome, ResolutionStatus
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.settings import RuntimeSettings, SemiImperiumSettings
from semi_imperium.views import CalculateView
from semi_imperium.workflows.calculation import (
    CalculateSession,
    ExecutionOutcome,
    ExecutionRequest,
    ReuseAction,
    build_provenance,
    execute_prepared,
    prepare_runs,
    review,
    save_prepared_runs,
)
from semi_imperium.workspace import SemiImperiumWorkspace

ETHANOL = MolecularIdentity(
    canonical_smiles="CCO",
    charge=0,
    multiplicity=1,
    display_name="ethanol",
    original_input="ethanol",
    input_type=MoleculeInputType.CHEMICAL_NAME,
    resolved_name="Ethanol",
    resolver="stub",
)
ACETIC_ACID = MolecularIdentity(
    canonical_smiles="CC(=O)O",
    charge=0,
    multiplicity=1,
    original_input="CC(=O)O",
    input_type=MoleculeInputType.SMILES,
    resolver="stub",
)

KNOWN = {
    "ethanol": ETHANOL,
    "CCO": ETHANOL,
    "CC(=O)O": ACETIC_ACID,
    "acetic acid": ACETIC_ACID,
}


# ---------------------------------------------------------------------------
# Doubles
# ---------------------------------------------------------------------------


@dataclass
class StubResolutionService:
    """Resolves a fixed table of inputs without a network or RDKit."""

    known: dict[str, MolecularIdentity] = field(default_factory=lambda: dict(KNOWN))
    calls: list[str] = field(default_factory=list)

    def resolve(
        self,
        value: str,
        input_type: Any,
        *,
        path: Any = None,
        charge: int | None = None,
        multiplicity: int = 1,
    ) -> ResolutionOutcome:
        self.calls.append(value)
        identity = self.known.get(value.strip())
        parsed_type = MoleculeInputType(input_type)
        if identity is None:
            return ResolutionOutcome(
                original_input=value,
                input_type=parsed_type,
                status=ResolutionStatus.UNRESOLVED,
                message=f"No structure was found for {value!r}",
            )
        return ResolutionOutcome(
            original_input=value,
            input_type=parsed_type,
            status=ResolutionStatus.RESOLVED,
            identity=identity,
        )


@dataclass
class RecordingExecutor:
    """Execution boundary double: records requests, returns fixed outcomes."""

    energy: float = -56.1
    raises: bool = False
    requests: list[ExecutionRequest] = field(default_factory=list)

    def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
        self.requests.append(request)
        if self.raises:
            raise RuntimeError("MOPAC is not installed")
        return ExecutionOutcome(
            state=CalculationState.VERIFIED,
            verification=VerificationOutcome.MINIMUM_CONFIRMED,
            result=CalculationResultData(
                energy_hof_kcal_mol=self.energy,
                conformer_index=0,
                conformers_evaluated=3,
            ),
        )

    @property
    def hamiltonians(self) -> list[str]:
        return [request.hamiltonian for request in self.requests]


@dataclass
class ScriptedPrompter:
    """Answers the views' questions in order; exhausted answers cancel."""

    texts: list[str] = field(default_factory=list)
    choices: list[str] = field(default_factory=list)
    confirms: list[bool] = field(default_factory=list)
    pauses: int = 0

    def text(self, message: str, *, default: str = "") -> str | None:
        return self.texts.pop(0) if self.texts else None

    def choice(self, message: str, options: Any) -> str | None:
        return self.choices.pop(0) if self.choices else None

    def confirm(self, message: str, *, default: bool = False) -> bool:
        return self.confirms.pop(0) if self.confirms else default

    def pause(self) -> None:
        self.pauses += 1


# ---------------------------------------------------------------------------
# Fixtures and builders
# ---------------------------------------------------------------------------


@pytest.fixture
def settings(tmp_path: Path) -> SemiImperiumSettings:
    return replace(
        SemiImperiumSettings(),
        runtime=RuntimeSettings(
            work_dir=tmp_path / "work", store_root=tmp_path / "store"
        ),
    )


@pytest.fixture
def store(tmp_path: Path) -> SemiImperiumStore:
    return SemiImperiumStore(tmp_path / "store")


@pytest.fixture
def resolver() -> StubResolutionService:
    return StubResolutionService()


@pytest.fixture
def workspace(
    settings: SemiImperiumSettings,
    store: SemiImperiumStore,
    resolver: StubResolutionService,
) -> SemiImperiumWorkspace:
    return SemiImperiumWorkspace(
        settings=settings,
        store=store,
        executor=RecordingExecutor(),
        prompter=ScriptedPrompter(),
        session=CalculateSession(resolution_service=resolver),  # type: ignore[arg-type]
    )


def make_view(workspace: SemiImperiumWorkspace) -> CalculateView:
    """Build the Calculate view over a non-interactive console."""
    console = Console(file=io.StringIO(), width=200, force_terminal=False)
    controller = SimpleNamespace(console=console)
    return CalculateView(controller, workspace)  # type: ignore[arg-type]


def menu_values(view: CalculateView) -> list[str]:
    return [option.value for option in view.get_menu_options()]


def seed_result(
    store: SemiImperiumStore,
    identity: MolecularIdentity,
    configuration: EffectiveConfiguration,
    *,
    run_id: str = "run-seed-01",
    energy: float = -50.0,
) -> CalculationRecord:
    """Persist one finished, reusable calculation the way a real run would."""
    provenance = ScientificProvenance(
        method_id=configuration.method_id,
        method_version=configuration.method_version,
        property_id=configuration.property_id,
        semi_imperium_version="test",
        grimperium_version="test",
    )
    run = RunRecord(
        run_id=run_id,
        configuration=configuration,
        provenance=provenance,
        molecule_ids=(identity.molecule_id,),
    )
    store.save_run(run.transition_to(RunState.RUNNING).transition_to(RunState.COMPLETED))

    record = CalculationRecord(
        run_id=run_id,
        molecule=identity,
        signature=configuration.signature(),
        provenance=provenance,
    )
    record = record.transition_to(CalculationState.RUNNING).transition_to(
        CalculationState.VERIFIED,
        verification=VerificationOutcome.MINIMUM_CONFIRMED,
        result=CalculationResultData(
            energy_hof_kcal_mol=energy, conformer_index=0, conformers_evaluated=4
        ),
    )
    store.save_calculation(record)
    return record


# ---------------------------------------------------------------------------
# 1. The table
# ---------------------------------------------------------------------------


def test_table_holds_many_molecules_by_name_or_smiles(
    resolver: StubResolutionService,
) -> None:
    session = CalculateSession(resolution_service=resolver)  # type: ignore[arg-type]

    first = session.add("ethanol", MoleculeInputType.CHEMICAL_NAME)
    second = session.add("CC(=O)O", MoleculeInputType.SMILES)

    assert [entry.entry_id for entry in session.entries] == ["row-1", "row-2"]
    assert first.input_type is MoleculeInputType.CHEMICAL_NAME
    assert second.input_type is MoleculeInputType.SMILES
    assert first.hamiltonians == ("PM7",)
    assert first.crest_enabled is True
    assert first.status == "draft"


def test_rows_carry_independent_hamiltonians_and_crest_switches(
    resolver: StubResolutionService,
) -> None:
    session = CalculateSession(resolution_service=resolver)  # type: ignore[arg-type]
    session.add("ethanol")
    session.add("acetic acid")

    session.set_hamiltonians("row-1", ["pm7", "am1"])
    session.set_crest("row-1", False)
    session.toggle_hamiltonian("row-2", "PM3")

    first = session.get("row-1")
    second = session.get("row-2")
    assert first.hamiltonians == ("AM1", "PM7")
    assert first.crest_enabled is False
    assert second.hamiltonians == ("PM3", "PM7")
    assert second.crest_enabled is True


def test_bulk_select_and_clear_scope_bulk_edits(
    resolver: StubResolutionService,
) -> None:
    session = CalculateSession(resolution_service=resolver)  # type: ignore[arg-type]
    session.add("ethanol")
    session.add("acetic acid")

    assert session.clear_selection() == 2
    assert session.selected_entries == ()
    assert session.bulk_set_crest(False) == 0

    session.set_selected("row-2", True)
    assert session.bulk_set_crest(False) == 1
    assert session.get("row-1").crest_enabled is True
    assert session.get("row-2").crest_enabled is False

    assert session.select_all() == 2
    assert session.bulk_set_hamiltonian("AM1", True) == 2
    assert session.get("row-1").hamiltonians == ("AM1", "PM7")
    assert session.get("row-2").hamiltonians == ("AM1", "PM7")


def test_editing_identity_discards_the_stale_resolution(
    resolver: StubResolutionService,
) -> None:
    session = CalculateSession(resolution_service=resolver)  # type: ignore[arg-type]
    session.add("ethanol")
    session.resolve("row-1")
    assert session.get("row-1").is_resolved is True

    session.set_input("row-1", "CC(=O)O", MoleculeInputType.SMILES)

    assert session.get("row-1").outcome is None
    assert session.get("row-1").status == "draft"
    session.resolve("row-1")
    assert session.get("row-1").identity == ACETIC_ACID


def test_unresolved_rows_stay_visible_with_their_reason(
    resolver: StubResolutionService,
) -> None:
    session = CalculateSession(resolution_service=resolver)  # type: ignore[arg-type]
    session.add("unobtainium")
    session.resolve_all()

    entry = session.get("row-1")
    assert entry.is_resolved is False
    assert entry.status == ResolutionStatus.UNRESOLVED.value
    assert entry.message is not None
    assert session.unresolved_entries == (entry,)


# ---------------------------------------------------------------------------
# 2. Reuse review
# ---------------------------------------------------------------------------


def test_review_separates_reusable_from_missing_per_hamiltonian(
    workspace: SemiImperiumWorkspace, settings: SemiImperiumSettings
) -> None:
    seed_result(
        workspace.store,
        ETHANOL,
        settings.configuration_for("PM7", crest_enabled=True),
    )
    workspace.session.add("ethanol")
    workspace.session.set_hamiltonians("row-1", ["AM1", "PM7"])

    plan = review(
        workspace.session, store=workspace.store, settings=workspace.settings
    )

    molecule = plan.molecules[0]
    assert [item.hamiltonian for item in molecule.reusable] == ["PM7"]
    assert [item.hamiltonian for item in molecule.missing] == ["AM1"]
    assert plan.reusable_count == 1
    assert plan.missing_count == 1
    assert plan.available_actions == (
        ReuseAction.REUSE_EXISTING,
        ReuseAction.CALCULATE_MISSING,
        ReuseAction.RECALCULATE_ALL,
    )


def test_review_does_not_reuse_across_a_different_crest_decision(
    workspace: SemiImperiumWorkspace, settings: SemiImperiumSettings
) -> None:
    seed_result(
        workspace.store,
        ETHANOL,
        settings.configuration_for("PM7", crest_enabled=True),
    )
    workspace.session.add("ethanol")
    workspace.session.set_crest("row-1", False)

    plan = review(
        workspace.session, store=workspace.store, settings=workspace.settings
    )

    assert plan.reusable_count == 0
    assert plan.missing_count == 1
    assert ReuseAction.REUSE_EXISTING not in plan.available_actions


def test_review_blocks_rows_that_cannot_be_calculated(
    workspace: SemiImperiumWorkspace,
) -> None:
    workspace.session.add("ethanol")
    workspace.session.add("unobtainium")

    plan = review(
        workspace.session, store=workspace.store, settings=workspace.settings
    )

    assert [molecule.entry_id for molecule in plan.molecules] == ["row-1"]
    assert [blocked.entry_id for blocked in plan.blocked] == ["row-2"]
    assert plan.blocked[0].status == ResolutionStatus.UNRESOLVED.value


# ---------------------------------------------------------------------------
# 3. The interface offers the choices, and executing them writes results
# ---------------------------------------------------------------------------


def test_reuse_choices_are_absent_until_the_review_has_run(
    workspace: SemiImperiumWorkspace,
) -> None:
    view = make_view(workspace)
    workspace.session.add("ethanol")

    before = menu_values(view)
    assert "review" in before
    for action in ReuseAction:
        assert action.value not in before

    assert view.handle_action("review") is None

    after = menu_values(view)
    assert ReuseAction.CALCULATE_MISSING.value in after
    assert ReuseAction.RECALCULATE_ALL.value in after
    assert ReuseAction.REUSE_EXISTING.value not in after


def test_review_then_reuse_is_offered_and_recomputes_nothing(
    workspace: SemiImperiumWorkspace, settings: SemiImperiumSettings
) -> None:
    stored = seed_result(
        workspace.store,
        ETHANOL,
        settings.configuration_for("PM7", crest_enabled=True),
    )
    before = workspace.store.calculation_path(stored).read_text(encoding="utf-8")

    view = make_view(workspace)
    workspace.session.add("ethanol")
    view.handle_action("review")

    assert ReuseAction.REUSE_EXISTING.value in menu_values(view)
    assert view.handle_action(ReuseAction.REUSE_EXISTING.value) is None

    executor = workspace.executor
    assert isinstance(executor, RecordingExecutor)
    assert executor.requests == []
    assert workspace.store.calculation_path(stored).read_text(encoding="utf-8") == before
    assert workspace.store.list_runs() == [workspace.store.load_run("run-seed-01")]

    summary = view.last_summary
    assert summary is not None
    assert summary.action is ReuseAction.REUSE_EXISTING
    assert len(summary.reused) == 1
    assert summary.completed == ()


def test_calculate_missing_runs_only_the_missing_hamiltonian(
    workspace: SemiImperiumWorkspace, settings: SemiImperiumSettings
) -> None:
    stored = seed_result(
        workspace.store,
        ETHANOL,
        settings.configuration_for("PM7", crest_enabled=True),
    )
    before = workspace.store.calculation_path(stored).read_text(encoding="utf-8")

    view = make_view(workspace)
    workspace.session.add("ethanol")
    workspace.session.set_hamiltonians("row-1", ["AM1", "PM7"])
    view.handle_action("review")

    assert view.handle_action(ReuseAction.CALCULATE_MISSING.value) is None

    executor = workspace.executor
    assert isinstance(executor, RecordingExecutor)
    assert executor.hamiltonians == ["AM1"]

    # The reusable PM7 result was left exactly as it was found.
    assert workspace.store.calculation_path(stored).read_text(encoding="utf-8") == before

    am1_configuration = workspace.settings.configuration_for(
        "AM1", crest_enabled=True
    )
    decision = workspace.store.find_reusable(
        ETHANOL, am1_configuration.signature()
    )
    assert decision.can_reuse is True
    assert decision.record is not None
    assert decision.record.state is CalculationState.VERIFIED
    assert decision.record.result is not None
    assert decision.record.result.energy_hof_kcal_mol == pytest.approx(-56.1)

    summary = view.last_summary
    assert summary is not None
    assert len(summary.completed) == 1
    assert len(summary.reused) == 1
    assert view.plan is not None
    assert view.plan.missing_count == 0


def test_prepared_run_captures_the_effective_configuration_it_used(
    workspace: SemiImperiumWorkspace,
) -> None:
    view = make_view(workspace)
    workspace.session.add("ethanol")
    workspace.session.set_hamiltonians("row-1", ["AM1"])
    workspace.session.set_crest("row-1", False)
    view.handle_action("review")

    view.handle_action(ReuseAction.CALCULATE_MISSING.value)

    runs = workspace.store.list_runs()
    assert len(runs) == 1
    configuration = runs[0].configuration
    assert configuration.semiempirical.hamiltonian == "AM1"
    assert configuration.conformer_search.enabled is False
    assert runs[0].state is RunState.COMPLETED
    assert runs[0].molecule_ids == (ETHANOL.molecule_id,)

    # Changing the defaults afterwards must not re-describe the stored run.
    workspace.settings = replace(
        workspace.settings,
        semiempirical=replace(
            workspace.settings.semiempirical, hamiltonian="PM3", solvent="water"
        ),
    )
    reloaded = workspace.store.load_run(runs[0].run_id)
    assert reloaded.configuration.semiempirical.hamiltonian == "AM1"
    assert reloaded.configuration.semiempirical.solvent == "gas"


def test_recalculate_all_keeps_the_stored_result_and_adds_a_new_one(
    workspace: SemiImperiumWorkspace, settings: SemiImperiumSettings
) -> None:
    stored = seed_result(
        workspace.store,
        ETHANOL,
        settings.configuration_for("PM7", crest_enabled=True),
        energy=-49.5,
    )
    stored_path = workspace.store.calculation_path(stored)
    before = stored_path.read_text(encoding="utf-8")

    prompter = workspace.prompter
    assert isinstance(prompter, ScriptedPrompter)
    prompter.confirms.append(True)

    view = make_view(workspace)
    workspace.session.add("ethanol")
    view.handle_action("review")
    assert view.handle_action(ReuseAction.RECALCULATE_ALL.value) is None

    executor = workspace.executor
    assert isinstance(executor, RecordingExecutor)
    assert executor.hamiltonians == ["PM7"]

    assert stored_path.exists()
    assert stored_path.read_text(encoding="utf-8") == before

    signature = settings.configuration_for("PM7", crest_enabled=True).signature()
    records = workspace.store.list_calculations(molecule=ETHANOL, signature=signature)
    assert len(records) == 2
    assert {record.run_id for record in records} != {"run-seed-01"}
    energies = {
        record.result.energy_hof_kcal_mol
        for record in records
        if record.result is not None
    }
    assert energies == {-49.5, -56.1}


def test_recalculate_all_writes_nothing_when_it_is_not_confirmed(
    workspace: SemiImperiumWorkspace, settings: SemiImperiumSettings
) -> None:
    seed_result(
        workspace.store,
        ETHANOL,
        settings.configuration_for("PM7", crest_enabled=True),
    )
    prompter = workspace.prompter
    assert isinstance(prompter, ScriptedPrompter)
    prompter.confirms.append(False)

    view = make_view(workspace)
    workspace.session.add("ethanol")
    view.handle_action("review")
    view.handle_action(ReuseAction.RECALCULATE_ALL.value)

    executor = workspace.executor
    assert isinstance(executor, RecordingExecutor)
    assert executor.requests == []
    assert len(workspace.store.list_calculations()) == 1
    assert view.last_summary is None


def test_a_failing_executor_is_recorded_rather_than_lost(
    workspace: SemiImperiumWorkspace,
) -> None:
    workspace.executor = RecordingExecutor(raises=True)
    view = make_view(workspace)
    workspace.session.add("ethanol")
    view.handle_action("review")

    view.handle_action(ReuseAction.CALCULATE_MISSING.value)

    records = workspace.store.list_calculations()
    assert len(records) == 1
    assert records[0].state is CalculationState.FAILED
    assert records[0].verification is VerificationOutcome.FAILED
    assert "MOPAC is not installed" in (records[0].error_message or "")

    runs = workspace.store.list_runs()
    assert runs[0].state is RunState.FAILED
    assert runs[0].error_message is not None

    summary = view.last_summary
    assert summary is not None
    assert len(summary.failed) == 1
    assert summary.completed == ()


def test_editing_the_table_invalidates_the_previous_review(
    workspace: SemiImperiumWorkspace,
) -> None:
    view = make_view(workspace)
    workspace.session.add("ethanol")
    view.handle_action("review")
    assert view.plan is not None

    prompter = workspace.prompter
    assert isinstance(prompter, ScriptedPrompter)
    prompter.choices.append("on")
    view.handle_action("crest")

    assert view.plan is None
    assert ReuseAction.CALCULATE_MISSING.value not in menu_values(view)
    assert view.handle_action(ReuseAction.CALCULATE_MISSING.value) is None
    executor = workspace.executor
    assert isinstance(executor, RecordingExecutor)
    assert executor.requests == []


def test_the_interface_adds_and_removes_molecules(
    workspace: SemiImperiumWorkspace,
) -> None:
    view = make_view(workspace)
    prompter = workspace.prompter
    assert isinstance(prompter, ScriptedPrompter)

    prompter.texts.append("ethanol, acetic acid")
    view.handle_action("add_name")
    assert len(workspace.session) == 2

    prompter.texts.append("CC(=O)O")
    view.handle_action("add_smiles")
    assert len(workspace.session) == 3
    assert workspace.session.get("row-3").input_type is MoleculeInputType.SMILES

    prompter.choices.append("row-2")
    view.handle_action("remove")
    assert [entry.entry_id for entry in workspace.session.entries] == [
        "row-1",
        "row-3",
    ]

    view.handle_action("clear_selection")
    assert workspace.session.selected_entries == ()
    view.handle_action("select_all")
    assert len(workspace.session.selected_entries) == 2

    # The whole table renders without a terminal.
    view.render()


def test_the_interface_edits_one_row_in_place(
    workspace: SemiImperiumWorkspace,
) -> None:
    view = make_view(workspace)
    workspace.session.add("ethanol")
    prompter = workspace.prompter
    assert isinstance(prompter, ScriptedPrompter)

    prompter.choices.extend(["row-1", "hamiltonian", "AM1"])
    view.handle_action("edit")
    assert workspace.session.get("row-1").hamiltonians == ("AM1", "PM7")

    prompter.choices.extend(["row-1", "multiplicity"])
    prompter.texts.append("3")
    view.handle_action("edit")
    assert workspace.session.get("row-1").multiplicity == 3


# ---------------------------------------------------------------------------
# The same decisions, taken through the API the interface calls
# ---------------------------------------------------------------------------


def test_prepare_and_save_register_pending_work_before_executing(
    workspace: SemiImperiumWorkspace,
) -> None:
    workspace.session.add("ethanol")
    workspace.session.set_hamiltonians("row-1", ["AM1", "PM3"])
    plan = review(
        workspace.session, store=workspace.store, settings=workspace.settings
    )

    work = prepare_runs(
        plan,
        ReuseAction.CALCULATE_MISSING,
        provenance=build_provenance(workspace.settings),
        run_id_factory=lambda: "run-test",
    )
    save_prepared_runs(work, workspace.store)

    # One run per effective configuration: AM1 and PM3 are not comparable.
    assert [prepared.run.run_id for prepared in work.runs] == [
        "run-test-01",
        "run-test-02",
    ]
    pending = workspace.store.list_calculations()
    assert len(pending) == 2
    assert {record.state for record in pending} == {CalculationState.PENDING}

    summary = execute_prepared(
        work, store=workspace.store, executor=workspace.executor
    )
    assert len(summary.completed) == 2
    assert {
        record.state for record in workspace.store.list_calculations()
    } == {CalculationState.VERIFIED}


def test_reuse_action_prepares_no_run_at_all(
    workspace: SemiImperiumWorkspace, settings: SemiImperiumSettings
) -> None:
    seed_result(
        workspace.store,
        ETHANOL,
        settings.configuration_for("PM7", crest_enabled=True),
    )
    workspace.session.add("ethanol")
    plan = review(
        workspace.session, store=workspace.store, settings=workspace.settings
    )

    work = prepare_runs(
        plan,
        ReuseAction.REUSE_EXISTING,
        provenance=build_provenance(workspace.settings),
    )

    assert work.runs == ()
    assert work.pending_count == 0
    assert len(work.reused) == 1
    assert save_prepared_runs(work, workspace.store) == ()


def test_stored_records_round_trip_through_the_store(
    workspace: SemiImperiumWorkspace,
) -> None:
    view = make_view(workspace)
    workspace.session.add("ethanol")
    view.handle_action("review")
    view.handle_action(ReuseAction.CALCULATE_MISSING.value)

    records = workspace.store.list_calculations()
    assert len(records) == 1
    path = workspace.store.calculation_path(records[0])
    payload = json.loads(path.read_text(encoding="utf-8"))

    assert payload["state"] == "verified"
    assert payload["verification"] == "minimum_confirmed"
    assert payload["molecule"]["canonical_smiles"] == "CCO"
    assert payload["provenance"]["source"] == "semi_imperium"
    assert CalculationRecord.from_dict(payload) == records[0]
