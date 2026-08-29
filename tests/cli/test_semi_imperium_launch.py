"""Contract for the launched Semi-Imperium shell and its three areas.

What is asserted here is what a user actually sees after launching:

1. the shell registers Semi-Imperium's own areas and nothing borrowed
   from Grimperium's model, analytics or worker screens;
2. the Database area reports each Hamiltonian, its verification, CREST
   usage, the selection strategy and the run behind every number, and
   drills down to the individual calculations;
3. the Settings area exposes CREST, MOPAC and runtime/readiness defaults
   and says plainly that they apply to the next run, not to results that
   already exist.
"""

from __future__ import annotations

import io
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
from rich.console import Console

from semi_imperium import menu
from semi_imperium.app import SemiImperiumCLI
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
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.settings import RuntimeSettings, SemiImperiumSettings
from semi_imperium.views import (
    CalculateView,
    DatabaseView,
    HamiltonianView,
    SettingsView,
)
from semi_imperium.workflows.calculation import CalculateSession
from semi_imperium.workflows.database import build_summary
from semi_imperium.workspace import SemiImperiumWorkspace

#: Vocabulary from Grimperium that this focused application must not show.
FOREIGN_VOCABULARY = (
    "delta",
    "xgboost",
    "kernel ridge",
    "machine learning",
    "r2 score",
    "worker",
    "placeholder",
    "mock",
    "analytics",
    "prediction",
)

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


# ---------------------------------------------------------------------------
# Doubles and helpers
# ---------------------------------------------------------------------------


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


def make_workspace(tmp_path: Path) -> SemiImperiumWorkspace:
    settings = replace(
        SemiImperiumSettings(),
        runtime=RuntimeSettings(
            work_dir=tmp_path / "work", store_root=tmp_path / "store"
        ),
    )
    return SemiImperiumWorkspace(
        settings=settings,
        store=SemiImperiumStore(tmp_path / "store"),
        prompter=ScriptedPrompter(),
        session=CalculateSession(resolution_service=None),
    )


def make_view(view_class: Any, workspace: SemiImperiumWorkspace) -> Any:
    console = Console(file=io.StringIO(), width=400, force_terminal=False)
    return view_class(SimpleNamespace(console=console), workspace)


def render_text(view: Any) -> str:
    view.console.file.seek(0)
    view.console.file.truncate(0)
    view.render()
    return str(view.console.file.getvalue())


def store_calculation(
    store: SemiImperiumStore,
    configuration: EffectiveConfiguration,
    *,
    run_id: str,
    state: CalculationState,
    verification: VerificationOutcome,
    energy: float | None,
) -> CalculationRecord:
    """Persist one finished calculation together with its run manifest."""
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
        molecule_ids=(ETHANOL.molecule_id,),
    )
    store.save_run(run.transition_to(RunState.RUNNING).transition_to(RunState.COMPLETED))

    record = CalculationRecord(
        run_id=run_id,
        molecule=ETHANOL,
        signature=configuration.signature(),
        provenance=provenance,
    )
    record = record.transition_to(CalculationState.RUNNING).transition_to(
        state,
        verification=verification,
        result=CalculationResultData(
            energy_hof_kcal_mol=energy, conformer_index=0, conformers_evaluated=5
        ),
    )
    store.save_calculation(record)
    return record


@pytest.fixture
def seeded(tmp_path: Path) -> SemiImperiumWorkspace:
    """A workspace whose store holds a verified AM1 and a saddle PM7."""
    workspace = make_workspace(tmp_path)
    settings = workspace.settings
    store_calculation(
        workspace.store,
        settings.configuration_for("AM1", crest_enabled=True),
        run_id="run-a-01",
        state=CalculationState.VERIFIED,
        verification=VerificationOutcome.MINIMUM_CONFIRMED,
        energy=-56.12,
    )
    store_calculation(
        workspace.store,
        settings.configuration_for("PM7", crest_enabled=False),
        run_id="run-b-01",
        state=CalculationState.SADDLE,
        verification=VerificationOutcome.SADDLE_POINT,
        energy=-51.0,
    )
    return workspace


# ---------------------------------------------------------------------------
# 1. The launched shell
# ---------------------------------------------------------------------------


def test_shell_registers_only_the_focused_semi_imperium_areas() -> None:
    app = SemiImperiumCLI()

    assert set(app.controller._views) == {"calc", "methods", "databases", "settings"}
    assert isinstance(app.controller.get_view("calc"), CalculateView)
    assert isinstance(app.controller.get_view("databases"), DatabaseView)
    assert isinstance(app.controller.get_view("settings"), SettingsView)
    assert isinstance(app.controller.get_view("methods"), HamiltonianView)
    for absent in ("models", "results", "about", "batch"):
        assert app.controller.get_view(absent) is None


def test_main_menu_offers_exactly_three_areas(monkeypatch: Any) -> None:
    captured: dict[str, Any] = {}

    def fake_show_menu_with_separator(**kwargs: Any) -> str:
        captured.update(kwargs)
        return "calc"

    monkeypatch.setattr(
        menu, "show_menu_with_separator", fake_show_menu_with_separator
    )

    assert menu.show_main_menu() == "calc"
    assert [
        (option.label, option.value) for option in captured["options"]
    ] == [
        ("CALCULATE", "calc"),
        ("DATABASE", "databases"),
        ("SETTINGS", "settings"),
    ]


def test_main_launches_semi_imperium_without_preflight() -> None:
    with (
        patch("sys.argv", ["semi-imperium", "--skip-preflight"]),
        patch("grimperium.utils.logging.setup_logging"),
        patch("semi_imperium.app.SemiImperiumCLI") as cli_class,
    ):
        cli = MagicMock()
        cli.run.return_value = 0
        cli_class.return_value = cli

        from semi_imperium.app import main

        result = main()

    assert result == 0
    cli.run.assert_called_once_with(skip_preflight=True)


def test_status_line_describes_the_staged_work(tmp_path: Path) -> None:
    workspace = make_workspace(tmp_path)
    workspace.session.add("ethanol")
    workspace.session.add("acetic acid")
    workspace.session.set_selected("row-2", False)

    app = SemiImperiumCLI(workspace=workspace)
    line = app.status_line()

    assert "Molecules: 2 (1 selected)" in line
    assert "Hamiltonians: PM7" in line
    assert "CREST: on" in line
    assert "Store:" in line
    lowered = line.lower()
    for term in ("dataset", "model", "property"):
        assert term not in lowered


def test_normal_presentation_stays_focused(seeded: SemiImperiumWorkspace) -> None:
    seeded.session.add("ethanol")
    rendered: list[str] = []
    for view_class in (CalculateView, DatabaseView, SettingsView, HamiltonianView):
        view = make_view(view_class, seeded)
        rendered.append(render_text(view))
        rendered.extend(
            f"{option.label} {option.description}"
            for option in view.get_menu_options()
        )

    text = "\n".join(rendered).lower()
    for term in FOREIGN_VOCABULARY:
        assert term not in text, f"foreign concept leaked into the interface: {term}"
    assert "crest" in text
    assert "mopac" in text
    assert "am1" in text


# ---------------------------------------------------------------------------
# 2. The Database area
# ---------------------------------------------------------------------------


def test_database_summary_states_every_hamiltonian_explicitly(
    seeded: SemiImperiumWorkspace,
) -> None:
    summary = build_summary(seeded.store)

    assert summary.molecule_count == 1
    assert summary.calculation_count == 2
    assert summary.run_count == 2

    molecule = summary.molecules[0]
    assert molecule.status_label("AM1") == "verified / minimum_confirmed"
    assert molecule.status_label("PM7") == "saddle / saddle_point"
    assert molecule.status_label("PM3") == "not calculated"

    am1 = molecule.cell("AM1")
    assert am1 is not None
    assert am1.crest_used is True
    assert am1.run_id == "run-a-01"
    assert am1.energy_label == "-56.12"
    assert am1.selection_strategy == "crest_energy_top_n"
    assert am1.selection_experimental is False

    pm7 = molecule.cell("PM7")
    assert pm7 is not None
    assert pm7.crest_used is False
    assert molecule.crest_label == "mixed"


def test_database_area_renders_the_operational_table(
    seeded: SemiImperiumWorkspace,
) -> None:
    view = make_view(DatabaseView, seeded)

    text = render_text(view)

    for expected in (
        "AM1",
        "PM3",
        "PM7",
        "verified / minimum_confirmed",
        "saddle / saddle_point",
        "not calculated",
        "crest_energy_top_n",
        "Last run",
        "ethanol",
    ):
        assert expected in text
    assert "detail" in [option.value for option in view.get_menu_options()]


def test_database_area_drills_down_into_stored_calculations(
    seeded: SemiImperiumWorkspace,
) -> None:
    view = make_view(DatabaseView, seeded)
    render_text(view)

    prompter = seeded.prompter
    assert isinstance(prompter, ScriptedPrompter)
    prompter.choices.append(ETHANOL.molecule_id)

    view.console.file.seek(0)
    view.console.file.truncate(0)
    assert view.handle_action("detail") is None
    text = str(view.console.file.getvalue())

    assert "run-a-01" in text
    assert "run-b-01" in text
    assert "-56.12" in text
    assert "Provenance" in text
    assert prompter.pauses == 1


def test_database_area_says_so_when_nothing_is_stored(tmp_path: Path) -> None:
    workspace = make_workspace(tmp_path)
    view = make_view(DatabaseView, workspace)

    text = render_text(view)

    assert "No calculations stored yet" in text
    assert "detail" not in [option.value for option in view.get_menu_options()]


# ---------------------------------------------------------------------------
# 3. The Settings area
# ---------------------------------------------------------------------------


def test_settings_area_exposes_crest_mopac_and_runtime_readiness(
    tmp_path: Path,
) -> None:
    workspace = make_workspace(tmp_path)
    view = make_view(SettingsView, workspace)

    text = render_text(view)

    assert "CREST" in text
    assert "MOPAC" in text
    assert "Runtime and readiness" in text
    assert "Minimum verification" in text
    assert "crest_energy_top_n" in text
    assert "defaults for the next run" in text
    values = [option.value for option in view.get_menu_options()]
    assert values == [
        "crest",
        "selection",
        "mopac",
        "verification",
        "runtime",
        "back",
    ]


def test_settings_edits_apply_to_the_next_run_only(
    seeded: SemiImperiumWorkspace,
) -> None:
    before = seeded.store.load_run("run-a-01").configuration
    view = make_view(SettingsView, seeded)

    prompter = seeded.prompter
    assert isinstance(prompter, ScriptedPrompter)
    prompter.confirms.append(False)
    prompter.texts.extend(["gfnff", "4", "3.0"])

    assert view.handle_action("crest") is None

    search = seeded.settings.conformer_search
    assert search.enabled is False
    assert search.method == "gfnff"
    assert search.max_conformers == 4
    assert search.energy_window_kcal_mol == pytest.approx(3.0)

    # The stored run keeps the effective configuration it actually used.
    after = seeded.store.load_run("run-a-01").configuration
    assert after == before
    assert after.conformer_search.enabled is True
    assert after.conformer_search.method == "gfn2"


def test_settings_reports_readiness_without_guessing(tmp_path: Path) -> None:
    workspace = make_workspace(tmp_path)
    workspace.set_runtime(
        replace(
            workspace.settings.runtime,
            crest_executable=str(tmp_path / "no-such-crest"),
            mopac_executable=str(tmp_path / "no-such-mopac"),
        )
    )

    tools = workspace.settings.runtime.readiness()

    assert [tool.name for tool in tools] == ["CREST", "MOPAC"]
    assert all(tool.available is False for tool in tools)
    assert all(tool.label == "not found" for tool in tools)
    assert workspace.settings.runtime.is_ready is False


# ---------------------------------------------------------------------------
# 4. The Hamiltonian picker
# ---------------------------------------------------------------------------


def test_hamiltonian_area_presents_three_independent_requests(
    tmp_path: Path,
) -> None:
    workspace = make_workspace(tmp_path)
    workspace.session.add("ethanol")
    view = make_view(HamiltonianView, workspace)

    text = render_text(view)
    assert "AM1" in text
    assert "PM3" in text
    assert "PM7" in text
    assert "independent requests" in text

    prompter = workspace.prompter
    assert isinstance(prompter, ScriptedPrompter)
    prompter.choices.append("ALL")
    assert view.handle_action("apply") == "calc"
    assert workspace.session.get("row-1").hamiltonians == ("AM1", "PM3", "PM7")
