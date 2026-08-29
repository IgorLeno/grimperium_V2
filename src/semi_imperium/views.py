"""The three focused Semi-Imperium areas, plus the Hamiltonian picker.

These views are the whole application surface: a Calculate table, a
Database report and a Settings screen. Nothing else is reachable, and
none of them borrows Grimperium's model, dataset or analytics screens.

The Calculate area deliberately does not stop at a plan. Once the reuse
review has run, the choices it found — reuse, calculate only the missing
Hamiltonians, or recalculate everything — appear as ordinary menu
entries, and choosing one prepares the runs, persists them and executes
them through the workspace's execution boundary.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace
from typing import TYPE_CHECKING

from rich.markup import escape
from rich.panel import Panel
from rich.table import Table

from grimperium.cli.menu import MenuOption
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView
from semi_imperium.domain import (
    ConformerSelectionStrategy,
    MoleculeInputType,
    VerificationPolicy,
)
from semi_imperium.mopac.models import SUPPORTED_HAMILTONIANS
from semi_imperium.prompts import Prompter, QuestionaryPrompter
from semi_imperium.settings import SemiImperiumSettings
from semi_imperium.workflows.calculation import (
    CalculatePlan,
    CalculateSession,
    CalculateSessionError,
    ExecutionSummary,
    PreparedCalculation,
    ReuseAction,
    build_provenance,
    execute_prepared,
    prepare_runs,
    review,
    save_prepared_runs,
)
from semi_imperium.workflows.database import (
    CalculationDetail,
    DatabaseSummary,
    build_summary,
    molecule_detail,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from grimperium.cli.controller import CliController
    from semi_imperium.domain import CalculationRecord
    from semi_imperium.persistence import SemiImperiumStore
    from semi_imperium.workspace import SemiImperiumWorkspace

BACK = MenuOption(label="Back", value="back", icon=ICONS["back"])


class SemiImperiumView(BaseView):
    """Shared workspace access and the single input boundary."""

    def __init__(
        self, controller: CliController, workspace: SemiImperiumWorkspace
    ) -> None:
        super().__init__(controller)
        self.workspace = workspace
        self._prompter: Prompter | None = workspace.prompter

    @property
    def prompter(self) -> Prompter:
        """The prompt boundary, interactive unless the workspace set one."""
        if self._prompter is None:
            self._prompter = QuestionaryPrompter(self.console)
        return self._prompter

    @property
    def settings(self) -> SemiImperiumSettings:
        """The current defaults."""
        return self.workspace.settings

    @property
    def store(self) -> SemiImperiumStore:
        """The store the workspace is pointing at."""
        return self.workspace.store

    def handle_action(self, action: str) -> str | None:
        """Route one menu value to its handler."""
        if action == "back":
            return "main"
        handler = self.actions().get(action)
        if handler is None:
            self.show_error(f"Unknown action: {action}")
            return None
        return handler()

    def actions(self) -> dict[str, Callable[[], str | None]]:
        """Return this view's action table."""
        return {}

    def show_error(self, message: str) -> None:
        """Report a failure, escaping any molecular input it quotes."""
        super().show_error(escape(message))

    def show_success(self, message: str) -> None:
        """Confirm a state change, escaping any molecular input it quotes."""
        super().show_success(escape(message))

    def _notice(self, message: str) -> None:
        """Show a neutral, visible note about a state change.

        The message is escaped: molecular input such as ``C[N+](C)(C)C``
        would otherwise be read as console markup rather than chemistry.
        """
        self.console.print(
            f"[{COLORS['muted']}]{escape(message)}[/{COLORS['muted']}]"
        )


# ---------------------------------------------------------------------------
# Calculate
# ---------------------------------------------------------------------------


class CalculateView(SemiImperiumView):
    """One table of molecules, reviewed against local results, then executed."""

    name = "calc"
    title = "Calculate"
    icon = ICONS["calc"]
    color = COLORS["calc"]

    def __init__(
        self, controller: CliController, workspace: SemiImperiumWorkspace
    ) -> None:
        super().__init__(controller, workspace)
        self.plan: CalculatePlan | None = None
        self.last_summary: ExecutionSummary | None = None

    @property
    def session(self) -> CalculateSession:
        """The molecule table shared with the rest of the workspace."""
        return self.workspace.session

    # -- rendering ------------------------------------------------------

    def render(self) -> None:
        """Draw the molecule table, the reuse review and the last outcome."""
        self._render_molecules()
        if self.plan is not None:
            self._render_review(self.plan)
        if self.last_summary is not None:
            self._render_summary(self.last_summary)

    def _render_molecules(self) -> None:
        entries = self.session.entries
        if not entries:
            self.console.print(
                Panel(
                    "No molecules yet. Add one or many by chemical name or by "
                    "SMILES; each row keeps its own charge, multiplicity, CREST "
                    "switch and AM1/PM3/PM7 requests.",
                    border_style=COLORS["border"],
                    padding=(1, 2),
                )
            )
            return

        table = Table(title="Molecules", border_style=COLORS["border"])
        table.add_column("#")
        table.add_column("Sel")
        table.add_column("Input")
        table.add_column("Type")
        table.add_column("Status")
        table.add_column("Charge")
        table.add_column("Mult")
        table.add_column("CREST")
        for hamiltonian in SUPPORTED_HAMILTONIANS:
            table.add_column(hamiltonian)
        for entry in entries:
            row = [
                entry.entry_id,
                "x" if entry.selected else "-",
                escape(entry.raw_input),
                entry.input_type.value,
                entry.status,
                "auto" if entry.charge is None else str(entry.charge),
                str(entry.multiplicity),
                "on" if entry.crest_enabled else "off",
            ]
            row.extend(
                "yes" if name in entry.hamiltonians else "no"
                for name in SUPPORTED_HAMILTONIANS
            )
            table.add_row(*row)
        self.console.print(table)

        for entry in self.session.entries:
            if entry.message:
                self._notice(f"{entry.entry_id} ({entry.label}): {entry.message}")

    def _render_review(self, plan: CalculatePlan) -> None:
        table = Table(
            title="Compatible local results", border_style=COLORS["border"]
        )
        table.add_column("Molecule")
        table.add_column("Hamiltonian")
        table.add_column("CREST")
        table.add_column("Local result")
        table.add_column("Signature")
        for molecule in plan.molecules:
            for hamiltonian_plan in molecule.hamiltonian_plans:
                table.add_row(
                    escape(molecule.label),
                    hamiltonian_plan.hamiltonian,
                    "on" if molecule.crest_enabled else "off",
                    hamiltonian_plan.state_label,
                    hamiltonian_plan.signature.short,
                )
        self.console.print(table)
        self.console.print(
            f"[{COLORS['muted']}]{plan.reusable_count} reusable, "
            f"{plan.missing_count} missing, "
            f"{plan.requested_count} requested[/{COLORS['muted']}]"
        )
        for blocked in plan.blocked:
            self._notice(
                f"Not calculable — {blocked.label} ({blocked.status}): "
                f"{blocked.message}"
            )

    def _render_summary(self, summary: ExecutionSummary) -> None:
        self.console.print(
            Panel(
                summary.describe(),
                title="Last decision",
                border_style=COLORS["success"],
                padding=(0, 2),
            )
        )

    # -- menu -----------------------------------------------------------

    def get_menu_options(self) -> list[MenuOption]:
        """Offer table editing, the reuse review and the reuse decisions."""
        options = [
            MenuOption(
                label="Add molecule by name",
                value="add_name",
                description="One name, or several separated by commas",
            ),
            MenuOption(
                label="Add molecule by SMILES",
                value="add_smiles",
                description="One SMILES, or several separated by commas",
            ),
        ]
        if self.session.entries:
            options.extend(
                [
                    MenuOption(label="Edit molecule", value="edit"),
                    MenuOption(
                        label="CREST for selected molecules", value="crest"
                    ),
                    MenuOption(
                        label="AM1 / PM3 / PM7 for selected molecules",
                        value="hamiltonians",
                    ),
                    MenuOption(label="Select all molecules", value="select_all"),
                    MenuOption(label="Clear selection", value="clear_selection"),
                    MenuOption(label="Remove molecule", value="remove"),
                    MenuOption(
                        label="Review local results",
                        value="review",
                        description="Resolve identities and check the store",
                    ),
                ]
            )

        plan = self.plan
        if plan is not None:
            for action in plan.available_actions:
                options.append(
                    MenuOption(
                        label=plan.action_label(action),
                        value=action.value,
                        style="calc",
                    )
                )
        options.append(BACK)
        return options

    def actions(self) -> dict[str, Callable[[], str | None]]:
        """Map every Calculate menu value to its handler."""
        handlers: dict[str, Callable[[], str | None]] = {
            "add_name": lambda: self._add(MoleculeInputType.CHEMICAL_NAME),
            "add_smiles": lambda: self._add(MoleculeInputType.SMILES),
            "edit": self._edit,
            "crest": self._bulk_crest,
            "hamiltonians": self._bulk_hamiltonians,
            "select_all": self._select_all,
            "clear_selection": self._clear_selection,
            "remove": self._remove,
            "review": self._review,
        }
        for action in ReuseAction:
            handlers[action.value] = _bind(self._execute, action)
        return handlers

    # -- table editing --------------------------------------------------

    def _invalidate_plan(self) -> None:
        """Drop the review: an edited table is no longer what was reviewed."""
        self.plan = None

    def _add(self, input_type: MoleculeInputType) -> str | None:
        label = (
            "Chemical name(s)"
            if input_type is MoleculeInputType.CHEMICAL_NAME
            else "SMILES"
        )
        raw = self.prompter.text(f"{label}, separated by commas")
        if raw is None or not raw.strip():
            self._notice("No molecule added.")
            return None
        try:
            added = self.session.add_many(raw.split(","), input_type)
        except CalculateSessionError as exc:
            self.show_error(str(exc))
            return None
        self._invalidate_plan()
        self.show_success(f"Added {len(added)} molecule(s) to the table.")
        return None

    def _pick_entry(self, message: str) -> str | None:
        options = [
            (f"{entry.entry_id}  {entry.raw_input} [{entry.status}]", entry.entry_id)
            for entry in self.session.entries
        ]
        return self.prompter.choice(message, options)

    def _edit(self) -> str | None:
        entry_id = self._pick_entry("Which molecule?")
        if entry_id is None:
            return None
        field = self.prompter.choice(
            "What do you want to change?",
            [
                ("Identity input (name or SMILES)", "input"),
                ("Charge", "charge"),
                ("Multiplicity", "multiplicity"),
                ("CREST for this molecule", "crest"),
                ("AM1 / PM3 / PM7 for this molecule", "hamiltonian"),
            ],
        )
        if field is None:
            return None
        try:
            self._apply_edit(entry_id, field)
        except (CalculateSessionError, ValueError) as exc:
            self.show_error(str(exc))
            return None
        self._invalidate_plan()
        return None

    def _apply_edit(self, entry_id: str, field: str) -> None:
        entry = self.session.get(entry_id)
        if field == "input":
            input_type = self.prompter.choice(
                "How is it entered?",
                [("Chemical name", "chemical_name"), ("SMILES", "smiles")],
            )
            if input_type is None:
                return
            value = self.prompter.text("New value", default=entry.raw_input)
            if value is None or not value.strip():
                return
            self.session.set_input(entry_id, value, input_type)
            self.show_success(f"{entry_id} now resolves {value.strip()!r}.")
        elif field == "charge":
            value = self.prompter.text(
                "Charge (blank keeps the structure's own charge)",
                default="" if entry.charge is None else str(entry.charge),
            )
            if value is None:
                return
            charge = None if not value.strip() else int(value.strip())
            self.session.set_charge(entry_id, charge)
            self.show_success(
                f"{entry_id} charge set to "
                f"{'auto' if charge is None else charge}."
            )
        elif field == "multiplicity":
            value = self.prompter.text(
                "Multiplicity", default=str(entry.multiplicity)
            )
            if value is None or not value.strip():
                return
            self.session.set_multiplicity(entry_id, int(value.strip()))
            self.show_success(f"{entry_id} multiplicity set to {value.strip()}.")
        elif field == "crest":
            enabled = self.prompter.confirm(
                f"Run the CREST conformer search for {entry_id}?",
                default=entry.crest_enabled,
            )
            self.session.set_crest(entry_id, enabled)
            self.show_success(
                f"{entry_id} CREST {'enabled' if enabled else 'disabled'}."
            )
        else:
            hamiltonian = self.prompter.choice(
                "Which Hamiltonian?",
                [(name, name) for name in SUPPORTED_HAMILTONIANS],
            )
            if hamiltonian is None:
                return
            updated = self.session.toggle_hamiltonian(entry_id, hamiltonian)
            self.show_success(
                f"{entry_id} now requests: "
                f"{', '.join(updated.hamiltonians) or 'nothing'}."
            )

    def _bulk_crest(self) -> str | None:
        selected = self.session.selected_entries
        if not selected:
            self.show_error("Select at least one molecule first.")
            return None
        choice = self.prompter.choice(
            f"CREST for {len(selected)} selected molecule(s)",
            [("Enable CREST", "on"), ("Disable CREST", "off")],
        )
        if choice is None:
            return None
        count = self.session.bulk_set_crest(choice == "on")
        self._invalidate_plan()
        self.show_success(
            f"CREST {'enabled' if choice == 'on' else 'disabled'} "
            f"for {count} molecule(s)."
        )
        return None

    def _bulk_hamiltonians(self) -> str | None:
        selected = self.session.selected_entries
        if not selected:
            self.show_error("Select at least one molecule first.")
            return None
        options: list[tuple[str, str]] = []
        for name in SUPPORTED_HAMILTONIANS:
            options.append((f"Request {name}", f"{name}:on"))
            options.append((f"Drop {name}", f"{name}:off"))
        choice = self.prompter.choice(
            f"AM1 / PM3 / PM7 for {len(selected)} selected molecule(s)", options
        )
        if choice is None:
            return None
        hamiltonian, _, state = choice.partition(":")
        try:
            count = self.session.bulk_set_hamiltonian(hamiltonian, state == "on")
        except CalculateSessionError as exc:
            self.show_error(str(exc))
            return None
        self._invalidate_plan()
        self.show_success(
            f"{hamiltonian} {'requested for' if state == 'on' else 'dropped from'} "
            f"{count} molecule(s)."
        )
        return None

    def _select_all(self) -> str | None:
        count = self.session.select_all()
        self._invalidate_plan()
        self.show_success(f"Selected {count} molecule(s).")
        return None

    def _clear_selection(self) -> str | None:
        count = self.session.clear_selection()
        self._invalidate_plan()
        self.show_success(f"Cleared the selection of {count} molecule(s).")
        return None

    def _remove(self) -> str | None:
        entry_id = self._pick_entry("Remove which molecule?")
        if entry_id is None:
            return None
        try:
            removed = self.session.remove(entry_id)
        except CalculateSessionError as exc:
            self.show_error(str(exc))
            return None
        self._invalidate_plan()
        self.show_success(f"Removed {removed.raw_input!r} from the table.")
        return None

    # -- review and execution -------------------------------------------

    def _review(self) -> str | None:
        if not self.session.selected_entries:
            self.show_error(
                "Select at least one molecule before reviewing local results."
            )
            return None
        self._notice("Resolving identities and checking compatible local results...")
        plan = review(
            self.session, store=self.store, settings=self.settings
        )
        self.plan = plan
        self.last_summary = None
        if not plan.requested_count:
            self.show_error(
                "Nothing is calculable yet; fix the molecules listed above."
            )
            return None
        self.show_success(
            f"{plan.reusable_count} compatible local result(s), "
            f"{plan.missing_count} missing — choose what to do below."
        )
        return None

    def _execute(self, action: ReuseAction) -> str | None:
        plan = self.plan
        if plan is None:
            self.show_error("Review local results before choosing what to run.")
            return None
        if action not in plan.available_actions:
            self.show_error(f"{action.value} is not available for this table.")
            return None
        if action is ReuseAction.RECALCULATE_ALL and not self.prompter.confirm(
            f"Recalculate {plan.requested_count} molecule/Hamiltonian pair(s)? "
            "Stored results are kept, not overwritten.",
            default=False,
        ):
            self._notice("Recalculation cancelled; nothing was written.")
            return None

        work = prepare_runs(
            plan, action, provenance=build_provenance(self.settings)
        )
        save_prepared_runs(work, self.store)
        if work.pending_count:
            self._notice(
                f"Prepared {work.pending_count} calculation(s) in "
                f"{len(work.runs)} run(s); executing..."
            )
        summary = execute_prepared(
            work,
            store=self.store,
            executor=self.workspace.executor,
            on_progress=self._on_progress,
        )
        self.last_summary = summary
        # The store changed, so the previous review is stale by definition.
        self.plan = review(self.session, store=self.store, settings=self.settings)
        self._report(summary)
        self.prompter.pause()
        return None

    def _on_progress(
        self, calculation: PreparedCalculation, record: CalculationRecord
    ) -> None:
        """Show every state change as it is persisted."""
        self.console.print(
            f"[{COLORS['muted']}]{calculation.entry_id} "
            f"{calculation.hamiltonian}: {record.state.value} / "
            f"{record.verification.value}[/{COLORS['muted']}]"
        )

    def _report(self, summary: ExecutionSummary) -> None:
        if summary.failed:
            self.show_error(
                f"{len(summary.failed)} calculation(s) failed; "
                f"{len(summary.completed)} completed, {len(summary.reused)} reused."
            )
            for record in summary.failed:
                self._notice(
                    f"{record.calculation_id}: "
                    f"{record.error_message or 'no reason recorded'}"
                )
            return
        if summary.action is ReuseAction.REUSE_EXISTING:
            self.show_success(
                f"Reused {len(summary.reused)} compatible local result(s); "
                "nothing was recomputed."
            )
            return
        self.show_success(
            f"{len(summary.completed)} calculation(s) completed, "
            f"{len(summary.reused)} reused."
        )


# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------


class DatabaseView(SemiImperiumView):
    """An operational report over everything Semi-Imperium has computed."""

    name = "databases"
    title = "Database"
    icon = ICONS["databases"]
    color = COLORS["databases"]

    def __init__(
        self, controller: CliController, workspace: SemiImperiumWorkspace
    ) -> None:
        super().__init__(controller, workspace)
        self.summary: DatabaseSummary | None = None

    def render(self) -> None:
        """Draw one row per molecule with explicit per-Hamiltonian statuses."""
        summary = build_summary(self.store)
        self.summary = summary

        if summary.is_empty:
            self.console.print(
                Panel(
                    escape(
                        f"No calculations stored yet in {self.store.root}.\n"
                        "Molecules appear here as soon as Calculate runs them."
                    ),
                    border_style=COLORS["border"],
                    padding=(1, 2),
                )
            )
            return

        table = Table(title="Molecules", border_style=COLORS["border"])
        table.add_column("Molecule")
        table.add_column("SMILES")
        for hamiltonian in SUPPORTED_HAMILTONIANS:
            table.add_column(hamiltonian)
        table.add_column("CREST")
        table.add_column("Selection")
        table.add_column("Last run")
        table.add_column("Updated (UTC)")
        for molecule in summary.molecules:
            row = [
                escape(molecule.label),
                escape(molecule.molecule.canonical_smiles),
            ]
            row.extend(
                molecule.status_label(name) for name in SUPPORTED_HAMILTONIANS
            )
            row.extend(
                [
                    molecule.crest_label,
                    molecule.selection_label,
                    molecule.last_run_id,
                    molecule.last_updated.strftime("%Y-%m-%d %H:%M"),
                ]
            )
            table.add_row(*row)
        self.console.print(table)

        counts = summary.state_counts()
        spelled = ", ".join(f"{state}: {count}" for state, count in counts.items())
        self.console.print(
            f"[{COLORS['muted']}]{summary.molecule_count} molecule(s), "
            f"{summary.calculation_count} calculation(s), "
            f"{summary.run_count} run(s) — {spelled}[/{COLORS['muted']}]"
        )

    def get_menu_options(self) -> list[MenuOption]:
        """Offer the drill-down and a refresh."""
        options: list[MenuOption] = []
        if self.summary is not None and not self.summary.is_empty:
            options.append(
                MenuOption(
                    label="Open molecule detail",
                    value="detail",
                    description="Every stored calculation, newest first",
                )
            )
        options.append(MenuOption(label="Refresh", value="refresh"))
        options.append(BACK)
        return options

    def actions(self) -> dict[str, Callable[[], str | None]]:
        """Map Database menu values to handlers."""
        return {"detail": self._detail, "refresh": self._refresh}

    def _refresh(self) -> str | None:
        self._notice("Reloading from the store...")
        return None

    def _detail(self) -> str | None:
        summary = self.summary
        if summary is None or summary.is_empty:
            self.show_error("There is nothing stored to open yet.")
            return None
        molecule_id = self.prompter.choice(
            "Which molecule?",
            [
                (f"{molecule.label}  ({molecule.molecule.canonical_smiles})",
                 molecule.molecule.molecule_id)
                for molecule in summary.molecules
            ],
        )
        if molecule_id is None:
            return None
        details = molecule_detail(self.store, molecule_id)
        if not details:
            self.show_error("That molecule has no stored calculations.")
            return None
        self._render_detail(details)
        self.prompter.pause()
        return None

    def _render_detail(self, details: tuple[CalculationDetail, ...]) -> None:
        table = Table(
            title=escape(
                f"Calculations for {details[0].record.molecule.canonical_smiles}"
            ),
            border_style=COLORS["border"],
        )
        table.add_column("Calculation")
        table.add_column("Run")
        table.add_column("Hamiltonian")
        table.add_column("State / verification")
        table.add_column("H298 (kcal/mol)")
        table.add_column("CREST")
        table.add_column("Selection")
        table.add_column("Signature")
        table.add_column("Updated (UTC)")
        for detail in details:
            record = detail.record
            result = record.result
            crest = detail.crest_used
            table.add_row(
                record.calculation_id,
                record.run_id,
                detail.hamiltonian,
                f"{record.state.value} / {record.verification.value}",
                (
                    "no value"
                    if result is None or result.energy_hof_kcal_mol is None
                    else f"{result.energy_hof_kcal_mol:.2f}"
                ),
                "unknown" if crest is None else ("yes" if crest else "no"),
                detail.selection_strategy
                + (" (experimental)" if detail.selection_experimental else ""),
                record.signature.short,
                detail.updated_at.strftime("%Y-%m-%d %H:%M"),
            )
        self.console.print(table)

        newest = details[0]
        provenance = newest.record.provenance
        lines = [
            f"Method: {provenance.method_id} v{provenance.method_version}",
            f"Property: {provenance.property_id}",
            f"Produced by: {provenance.source} "
            f"(semi-imperium {provenance.semi_imperium_version})",
        ]
        if newest.record.error_message:
            lines.append(f"Last error: {newest.record.error_message}")
        self.console.print(
            Panel(
                escape("\n".join(lines)),
                title="Provenance",
                border_style=COLORS["border"],
                padding=(0, 2),
            )
        )


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class SettingsView(SemiImperiumView):
    """CREST, MOPAC and runtime defaults for the next prepared run."""

    name = "settings"
    title = "Settings"
    icon = ICONS["settings"]
    color = COLORS["settings"]

    def render(self) -> None:
        """Show the three groups of defaults and the live tool readiness."""
        settings = self.settings
        search = settings.conformer_search
        selection = settings.conformer_selection
        semiempirical = settings.semiempirical
        verification = settings.verification
        runtime = settings.runtime

        self.console.print(
            Panel(
                f"Conformer search: {'on' if search.enabled else 'off'}\n"
                f"Program / method: {search.program} / {search.method}\n"
                f"Quick mode: {search.quick_mode}   "
                f"Optimization level: {search.opt_level}\n"
                f"RMSD threshold: {search.rmsd_threshold} A   "
                f"Energy window: {search.energy_window_kcal_mol} kcal/mol\n"
                f"Max conformers: {search.max_conformers}   "
                f"Pre-optimizer: {search.preoptimizer}\n"
                f"Selection: {selection.strategy} (top {selection.top_n})"
                f"{' — experimental' if selection.is_experimental else ''}",
                title="CREST",
                border_style=COLORS["settings"],
                padding=(0, 2),
            )
        )
        self.console.print(
            Panel(
                f"Program: {semiempirical.program}\n"
                f"Default Hamiltonian: {semiempirical.hamiltonian}\n"
                f"Keywords: {', '.join(semiempirical.keywords) or 'none'}\n"
                f"SCF convergence: {semiempirical.scf_convergence}\n"
                f"Solvent: {semiempirical.solvent}\n"
                f"Minimum verification: {verification.policy.value}",
                title="MOPAC",
                border_style=COLORS["settings"],
                padding=(0, 2),
            )
        )
        readiness = "\n".join(
            f"{tool.name}: {tool.label} ({tool.detail})"
            for tool in runtime.readiness()
        )
        self.console.print(
            Panel(
                escape(
                    f"CREST executable: {runtime.crest_executable}\n"
                    f"MOPAC executable: {runtime.mopac_executable}\n"
                    f"CREST threads: {runtime.crest_threads}   "
                    f"CREST timeout: {runtime.crest_timeout_seconds}s\n"
                    f"MOPAC timeout: {runtime.mopac_timeout_seconds}s\n"
                    f"Work directory: {runtime.work_dir}\n"
                    f"Store: {runtime.store_root}\n\n"
                    f"{readiness}"
                ),
                title="Runtime and readiness",
                border_style=COLORS["settings"],
                padding=(0, 2),
            )
        )
        self.console.print(
            f"[{COLORS['muted']}]These are defaults for the next run. Every run "
            "stores the effective values it used, so changing them here never "
            f"re-describes a result that already exists.[/{COLORS['muted']}]"
        )

    def get_menu_options(self) -> list[MenuOption]:
        """Offer the focused groups of defaults."""
        return [
            MenuOption(label="CREST conformer search", value="crest"),
            MenuOption(label="Conformer selection strategy", value="selection"),
            MenuOption(label="MOPAC", value="mopac"),
            MenuOption(label="Minimum verification", value="verification"),
            MenuOption(
                label="Runtime (executables, threads, timeouts)", value="runtime"
            ),
            BACK,
        ]

    def actions(self) -> dict[str, Callable[[], str | None]]:
        """Map Settings menu values to handlers."""
        return {
            "crest": self._edit_crest,
            "selection": self._edit_selection,
            "mopac": self._edit_mopac,
            "verification": self._edit_verification,
            "runtime": self._edit_runtime,
        }

    def _edit_crest(self) -> str | None:
        search = self.settings.conformer_search
        enabled = self.prompter.confirm(
            "Run the CREST conformer search by default?", default=search.enabled
        )
        method = self.prompter.text("CREST method", default=search.method)
        if method is None:
            return None
        max_conformers = self.prompter.text(
            "Maximum conformers", default=str(search.max_conformers)
        )
        if max_conformers is None:
            return None
        window = self.prompter.text(
            "Energy window (kcal/mol)",
            default=str(search.energy_window_kcal_mol),
        )
        if window is None:
            return None
        try:
            updated = replace(
                search,
                enabled=enabled,
                method=method.strip() or search.method,
                max_conformers=int(max_conformers),
                energy_window_kcal_mol=float(window),
            )
        except ValueError as exc:
            self.show_error(f"Invalid CREST setting: {exc}")
            return None
        self.workspace.settings = replace(self.settings, conformer_search=updated)
        self.show_success("CREST defaults updated for the next run.")
        return None

    def _edit_selection(self) -> str | None:
        selection = self.settings.conformer_selection
        strategy = self.prompter.choice(
            "Which conformer selection strategy?",
            [
                (
                    f"{item.value}"
                    + (" (experimental)" if item.is_experimental else " (default)"),
                    item.value,
                )
                for item in ConformerSelectionStrategy
            ],
        )
        if strategy is None:
            return None
        top_n = self.prompter.text(
            "How many conformers reach MOPAC?", default=str(selection.top_n)
        )
        if top_n is None:
            return None
        try:
            updated = replace(selection, strategy=strategy, top_n=int(top_n))
        except ValueError as exc:
            self.show_error(f"Invalid selection setting: {exc}")
            return None
        self.workspace.settings = replace(self.settings, conformer_selection=updated)
        if updated.is_experimental:
            self._notice(
                "This strategy is experimental; results carry that mark."
            )
        self.show_success("Conformer selection updated for the next run.")
        return None

    def _edit_mopac(self) -> str | None:
        semiempirical = self.settings.semiempirical
        hamiltonian = self.prompter.choice(
            "Default Hamiltonian for new molecules",
            [(name, name) for name in SUPPORTED_HAMILTONIANS],
        )
        if hamiltonian is None:
            return None
        precise = self.prompter.confirm(
            "Use the PRECISE keyword?", default="PRECISE" in semiempirical.keywords
        )
        scf = self.prompter.text(
            "SCF convergence", default=str(semiempirical.scf_convergence)
        )
        if scf is None:
            return None
        try:
            updated = replace(
                semiempirical,
                hamiltonian=hamiltonian,
                keywords=("PRECISE",) if precise else (),
                scf_convergence=float(scf),
            )
        except ValueError as exc:
            self.show_error(f"Invalid MOPAC setting: {exc}")
            return None
        self.workspace.settings = replace(self.settings, semiempirical=updated)
        self.workspace.session.default_hamiltonians = (hamiltonian,)
        self.show_success(f"MOPAC defaults updated; new molecules ask for {hamiltonian}.")
        return None

    def _edit_verification(self) -> str | None:
        verification = self.settings.verification
        policy = self.prompter.choice(
            "How strongly must a geometry be proven to be a minimum?",
            [
                ("none — no frequencies are computed", VerificationPolicy.NONE.value),
                (
                    "advisory — frequencies recorded, never rejecting",
                    VerificationPolicy.ADVISORY.value,
                ),
                (
                    "require_minimum — only confirmed minima are accepted",
                    VerificationPolicy.REQUIRE_MINIMUM.value,
                ),
            ],
        )
        if policy is None:
            return None
        updated = replace(verification, policy=VerificationPolicy(policy))
        self.workspace.settings = replace(self.settings, verification=updated)
        self._notice(
            "Verification is part of the signature, so results computed under a "
            "looser policy are no longer reusable here."
        )
        self.show_success(f"Minimum verification set to {policy}.")
        return None

    def _edit_runtime(self) -> str | None:
        runtime = self.settings.runtime
        crest = self.prompter.text(
            "CREST executable", default=runtime.crest_executable
        )
        if crest is None:
            return None
        mopac = self.prompter.text(
            "MOPAC executable", default=runtime.mopac_executable
        )
        if mopac is None:
            return None
        threads = self.prompter.text(
            "CREST threads", default=str(runtime.crest_threads)
        )
        if threads is None:
            return None
        timeout = self.prompter.text(
            "MOPAC timeout (seconds)", default=str(runtime.mopac_timeout_seconds)
        )
        if timeout is None:
            return None
        try:
            self.workspace.set_runtime(
                replace(
                    runtime,
                    crest_executable=crest.strip() or runtime.crest_executable,
                    mopac_executable=mopac.strip() or runtime.mopac_executable,
                    crest_threads=int(threads),
                    mopac_timeout_seconds=float(timeout),
                )
            )
        except ValueError as exc:
            self.show_error(f"Invalid runtime setting: {exc}")
            return None
        for tool in self.workspace.settings.runtime.readiness():
            if not tool.available:
                self._notice(f"{tool.name} is not runnable yet: {tool.detail}")
        self.show_success("Runtime defaults updated for the next run.")
        return None


# ---------------------------------------------------------------------------
# Hamiltonians
# ---------------------------------------------------------------------------


class HamiltonianView(SemiImperiumView):
    """The AM1/PM3/PM7 picker that supports the Calculate table."""

    name = "methods"
    title = "Hamiltonians"
    icon = ICONS["methods"]
    color = COLORS["calc"]

    def render(self) -> None:
        """Explain that the three Hamiltonians are independent requests."""
        table = Table(title="Semiempirical Hamiltonians", border_style=COLORS["border"])
        table.add_column("Hamiltonian")
        table.add_column("Default for new molecules")
        table.add_column("Selected molecules requesting it")
        selected = self.workspace.session.selected_entries
        defaults = self.workspace.session.default_hamiltonians
        for name in SUPPORTED_HAMILTONIANS:
            table.add_row(
                name,
                "yes" if name in defaults else "no",
                str(sum(1 for entry in selected if name in entry.hamiltonians)),
            )
        self.console.print(table)
        self.console.print(
            f"[{COLORS['muted']}]AM1, PM3 and PM7 are independent requests: each "
            "one is optimized and verified on its own, and each one is stored "
            f"under its own signature.[/{COLORS['muted']}]"
        )

    def get_menu_options(self) -> list[MenuOption]:
        """Offer the default set and a bulk apply."""
        return [
            MenuOption(label="Set the default for new molecules", value="defaults"),
            MenuOption(label="Apply a set to the selected molecules", value="apply"),
            BACK,
        ]

    def actions(self) -> dict[str, Callable[[], str | None]]:
        """Map Hamiltonian menu values to handlers."""
        return {"defaults": self._set_defaults, "apply": self._apply_selected}

    def _pick_set(self, message: str) -> tuple[str, ...] | None:
        choice = self.prompter.choice(
            message,
            [(name, name) for name in SUPPORTED_HAMILTONIANS]
            + [("AM1 + PM3 + PM7", "ALL")],
        )
        if choice is None:
            return None
        return SUPPORTED_HAMILTONIANS if choice == "ALL" else (choice,)

    def _set_defaults(self) -> str | None:
        chosen = self._pick_set("Which Hamiltonians should new molecules request?")
        if chosen is None:
            return None
        self.workspace.session.default_hamiltonians = chosen
        self.show_success(f"New molecules will request {', '.join(chosen)}.")
        return None

    def _apply_selected(self) -> str | None:
        session = self.workspace.session
        if not session.selected_entries:
            self.show_error("Select at least one molecule in Calculate first.")
            return None
        chosen = self._pick_set("Which Hamiltonians should they request?")
        if chosen is None:
            return None
        count = session.bulk_set_hamiltonians(chosen)
        self.show_success(f"{', '.join(chosen)} requested for {count} molecule(s).")
        return "calc"


def _bind(
    handler: Callable[[ReuseAction], str | None], action: ReuseAction
) -> Callable[[], str | None]:
    """Bind one reuse action to its handler without a late-binding closure."""
    return lambda: handler(action)


__all__ = [
    "CalculateView",
    "DatabaseView",
    "HamiltonianView",
    "SemiImperiumView",
    "SettingsView",
]
