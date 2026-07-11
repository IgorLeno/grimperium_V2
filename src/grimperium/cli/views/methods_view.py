"""
Calculation Methods view for GRIMPERIUM CLI.

Lists available methods and sets the active method on the session controller.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from rich.panel import Panel
from rich.table import Table

from grimperium.calculation.methods import (
    CalculationMethodDefinition,
    get_calculation_method,
    list_calculation_methods,
)
from grimperium.cli.menu import MenuOption, show_back_menu, show_menu
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController

DEFAULT_PROPERTY_ID = "standard_enthalpy_of_formation"


class MethodsView(BaseView):
    """View for selecting and inspecting calculation methods."""

    name = "methods"
    title = "Calculation Methods"
    icon = ICONS["methods"]
    color = COLORS["calc"]

    def __init__(self, controller: CliController) -> None:
        super().__init__(controller)

    def _available_methods(self) -> list[CalculationMethodDefinition]:
        return list_calculation_methods(DEFAULT_PROPERTY_ID)

    def render(self) -> None:
        """Render the methods catalogue and active session method."""
        self.clear_screen()
        self.show_header()

        summary = self.controller.session_summary()
        active = (
            f"[{COLORS['muted']}]Active Property:[/{COLORS['muted']}] "
            f"[{COLORS['calc']}]{summary['property']}[/{COLORS['calc']}]\n"
            f"[{COLORS['muted']}]Active Method:[/{COLORS['muted']}] "
            f"[{COLORS['calc']}]{summary['method']}[/{COLORS['calc']}]\n"
            f"[{COLORS['muted']}]Model:[/{COLORS['muted']}] "
            f"[{COLORS['calc']}]{summary['model']}[/{COLORS['calc']}]\n"
            f"[{COLORS['muted']}]Status:[/{COLORS['muted']}] "
            f"[{COLORS['calc']}]{summary['status']}[/{COLORS['calc']}]"
        )
        self.console.print(
            Panel(
                active,
                title=f"[bold {COLORS['calc']}]Session[/bold {COLORS['calc']}]",
                border_style=COLORS["border"],
                padding=(1, 2),
            )
        )
        self.console.print()

        table = Table(
            title="Available Methods",
            show_header=True,
            header_style=f"bold {COLORS['calc']}",
            border_style=COLORS["border"],
        )
        table.add_column("Method ID", style="bold", min_width=28, no_wrap=True)
        table.add_column("Name")
        table.add_column("Model")
        table.add_column("Conformer Strategy")
        table.add_column("xTB")

        for method in self._available_methods():
            model = "Required" if method.model_requirement.model_required else "No"
            xtb = (
                "Default"
                if method.xtb.optional and method.xtb.enabled_by_default
                else "Disabled"
            )
            marker = ""
            if self.controller.current_method_id == method.method_id:
                marker = f" [{COLORS['success']}]●[/{COLORS['success']}]"
            table.add_row(
                f"{method.method_id}{marker}",
                method.display_name,
                model,
                method.conformer_selection.strategy,
                xtb,
            )

        self.console.print(table)
        self.console.print()

    def render_method_details(self, method: CalculationMethodDefinition) -> None:
        """Render a detailed panel for one method definition."""
        model_line = (
            "Required" if method.model_requirement.model_required else "Not required"
        )
        xtb_line = (
            "Optional (enabled by default)"
            if method.xtb.optional and method.xtb.enabled_by_default
            else "Disabled by default"
        )
        details = f"""
[bold]{method.display_name}[/bold]

[{COLORS['muted']}]Method ID:[/{COLORS['muted']}] {method.method_id}
[{COLORS['muted']}]Version:[/{COLORS['muted']}] {method.version}
[{COLORS['muted']}]Property:[/{COLORS['muted']}] {method.property_name}
[{COLORS['muted']}]Conformer strategy:[/{COLORS['muted']}] {method.conformer_selection.strategy}
[{COLORS['muted']}]Model:[/{COLORS['muted']}] {model_line}
[{COLORS['muted']}]xTB:[/{COLORS['muted']}] {xtb_line}
"""
        self.console.print(
            Panel(
                details.strip(),
                title=f"[bold {COLORS['calc']}]Method Details[/bold {COLORS['calc']}]",
                border_style=COLORS["calc"],
                padding=(1, 2),
            )
        )
        if method.model_requirement.model_required and self.controller.status in {
            "Model required",
            "Model incompatible",
        }:
            self.console.print(
                f"[{COLORS['warning']}]{ICONS['warning']} "
                f"{self.controller.status}. Select a compatible model in MODELS "
                f"before running Calculate.[/{COLORS['warning']}]"
            )
        self.console.print()

    def select_active_method(self) -> bool:
        """Prompt for a method and persist it on the controller.

        Returns True if a method was selected, False if cancelled.
        """
        previous = self.controller.current_method_definition
        methods = self._available_methods()
        selected = show_menu(
            [
                MenuOption(
                    label=method.display_name,
                    value=method.method_id,
                    icon=ICONS["methods"],
                )
                for method in methods
            ],
            title="Select Active Method",
        )
        if selected is None:
            # Preserve previous session method on cancel.
            if previous is not None:
                self.controller.set_method(previous)
            return False

        method = get_calculation_method(
            selected,
            property_id=DEFAULT_PROPERTY_ID,
        )
        self.controller.set_method(method)
        self.render_method_details(method)
        self.show_success(f"Active method: {method.display_name}")
        self.wait_for_enter()
        return True

    def show_active_details(self) -> None:
        """Show details for the currently active method."""
        method = self.controller.current_method_definition
        if method is None:
            self.show_error("No method selected. Choose an active method first.")
            self.wait_for_enter()
            return
        self.render_method_details(method)
        self.wait_for_enter()

    def get_menu_options(self) -> list[MenuOption]:
        return [
            MenuOption(
                label="Select Active Method",
                value="select",
                icon=ICONS["methods"],
            ),
            MenuOption(
                label="View Active Method Details",
                value="details",
                icon=ICONS["about"],
            ),
        ]

    def handle_action(self, action: str) -> str | None:
        if action == "back":
            return "main"
        if action == "select":
            self.select_active_method()
            return None
        if action == "details":
            self.show_active_details()
            return None
        return None

    def run(self) -> str | None:
        while True:
            self.render()
            result = show_back_menu(
                options=self.get_menu_options(),
                title="Actions",
            )
            if result is None or result == "back":
                return "main"
            next_view = self.handle_action(result)
            if next_view:
                return next_view
