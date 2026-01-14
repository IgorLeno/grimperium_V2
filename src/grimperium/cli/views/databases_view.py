"""
Databases view for GRIMPERIUM CLI.

Displays and manages molecular databases.
"""

import json
import os
from datetime import date, datetime
from typing import TYPE_CHECKING, Any

from rich.panel import Panel
from rich.table import Table

from grimperium.cli.constants import PHASE_A_RESULTS_FILE
from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.mock_data import DATABASES, Database
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController


class DatabasesView(BaseView):
    """View for managing molecular databases."""

    name = "databases"
    title = "Databases"
    icon = ICONS["databases"]
    color = COLORS["databases"]

    def __init__(self, controller: "CliController") -> None:
        """Initialize the databases view."""
        super().__init__(controller)
        self.selected_db: Database | None = None

    @staticmethod
    def load_real_phase_a_results() -> Database | None:
        """Load real Phase A CREST PM7 results from phase_a_results.json.

        Returns:
            Database object with real data, or None if file doesn't exist.
        """
        if not PHASE_A_RESULTS_FILE.exists():
            return None
        try:
            with open(PHASE_A_RESULTS_FILE, encoding="utf-8") as f:
                data: dict[str, Any] = json.load(f)
                data: dict[str, Any] = json.load(f)

            molecules_count = data.get("n_molecules", 0)
            results = data.get("results", [])

            last_updated_str = None
            if results:
                last_result = results[-1]
                timestamp = last_result.get("timestamp", "")
                if timestamp:
                    last_updated_str = timestamp[:10]

            last_updated: date | None
            if last_updated_str:
                last_updated = date.fromisoformat(last_updated_str)
            else:
                # Quando não há timestamp no JSON, inferimos via mtime do arquivo.
                # Se falhar, deixamos None para o chamador tratar como "desconhecido".
                try:
                    file_mtime = os.path.getmtime(PHASE_A_RESULTS_FILE)
                    last_updated = date.fromtimestamp(file_mtime)
                except OSError:
                    last_updated = None

            return Database(
                name="CREST PM7",
                description="CREST conformer search with PM7 optimization",
                molecules=molecules_count,
                last_updated=last_updated,
                status="ready" if molecules_count > 0 else "in_development",
                properties=["H298_pm7", "conformers", "smiles", "quality_grade"],
            )
        except (json.JSONDecodeError, OSError, ValueError):
            return None

    def get_databases(self) -> list[Database]:
        """Get list of databases, replacing CREST PM7 with real data if available.

        Returns:
            List of Database objects with CREST PM7 from real calculations.
        """
        real_pm7 = self.load_real_phase_a_results()
        databases = []

        for db in DATABASES:
            if db.name == "CREST PM7" and real_pm7 is not None:
                databases.append(real_pm7)
            else:
                databases.append(db)

        return databases

    def _format_last_updated(self, value: datetime | date | None) -> str:
        """Format database last_updated timestamp consistently.

        Args:
            value: A datetime, date, or None.

        Returns:
            Formatted string or "-" if value is None.

        Note:
            If the Database.last_updated is timezone-aware, the timezone info
            is preserved in the formatted output.
        """
        if value is None:
            return "-"
        if isinstance(value, datetime):
            return value.strftime("%Y-%m-%d %H:%M:%S")
        if isinstance(value, date):
            return value.strftime("%Y-%m-%d")
        return str(value)

    def render(self) -> None:
        """Render the databases overview."""
        self.clear_screen()
        self.show_header()

        # Databases table
        table = Table(
            title="Available Databases",
            show_header=True,
            header_style=f"bold {COLORS['databases']}",
            border_style=COLORS["border"],
        )
        table.add_column("Name", style="bold")
        table.add_column("Molecules", justify="right")
        table.add_column("Last Updated")
        table.add_column("Status")

        for db in self.get_databases():
            if db.status == "ready":
                status = f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
            else:
                status = (
                    f"[{COLORS['in_dev']}]{ICONS['in_dev']} In Dev[/{COLORS['in_dev']}]"
                )

            # Format last_updated consistently
            if isinstance(db.last_updated, datetime):
                last_updated_str = db.last_updated.strftime("%Y-%m-%d %H:%M:%S")
            elif isinstance(db.last_updated, date):
                last_updated_str = db.last_updated.strftime("%Y-%m-%d")
            elif db.last_updated is None:
                last_updated_str = "-"
            else:
                last_updated_str = str(db.last_updated)

            table.add_row(
                db.name,
                f"{db.molecules:,}" if db.molecules > 0 else "-",
                last_updated_str,
                status,
            )

        self.console.print(table)
        self.console.print()

    def render_database_detail(self, db: Database) -> None:
        """Render detailed view for a specific database."""
        self.clear_screen()
        self.show_header()

        # Database info panel
        status_text = (
            f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
            if db.status == "ready"
            else f"[{COLORS['in_dev']}]{ICONS['in_dev']} In Development[/{COLORS['in_dev']}]"
        )

        last_updated_str = self._format_last_updated(db.last_updated)

        info = f"""
[bold]Name:[/bold]         {db.name}
[bold]Description:[/bold]  {db.description}
[bold]Molecules:[/bold]    {db.molecules:,}
[bold]Last Updated:[/bold] {last_updated_str}
[bold]Status:[/bold]       {status_text}

[bold]Properties:[/bold]
"""
        for prop in db.properties:
            info += f"  • {prop}\n"

        self.console.print(
            Panel(
                info,
                title=f"[bold {COLORS['databases']}]{db.name}[/bold {COLORS['databases']}]",
                border_style=COLORS["databases"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the databases view."""
        options = []
        for db in self.get_databases():
            options.append(
                MenuOption(
                    label=f"View {db.name}",
                    value=f"view_{db.name}",
                    icon=ICONS["databases"],
                )
            )

        options.extend(
            [
                MenuOption(
                    label="Calculate PM7 Values",
                    value="calculate",
                    disabled=False,
                ),
                MenuOption(
                    label="Add New Database",
                    value="add",
                    disabled=True,
                    disabled_reason="In Development",
                ),
            ]
        )

        return options

    def get_detail_menu_options(self) -> list[MenuOption]:
        """Return menu options for database detail view."""
        return [
            MenuOption(
                label="Calculate PM7 Values",
                value="calculate",
                disabled=False,
            ),
            MenuOption(
                label="Edit Database",
                value="edit",
                disabled=True,
                disabled_reason="In Development",
            ),
            MenuOption(
                label="Delete Database",
                value="delete",
                disabled=True,
                disabled_reason="In Development",
            ),
        ]

    def handle_action(self, action: str) -> str | None:
        """Handle menu actions."""
        if action == "back":
            if self.selected_db:
                self.selected_db = None
                return None  # Stay in databases view
            return "main"

        if action and action.startswith("view_"):
            db_name = action.removeprefix("view_")
            found = False
            for db in self.get_databases():
                if db.name == db_name:
                    self.selected_db = db
                    found = True
                    break

            if not found:
                self.selected_db = None

            return None

        if action == "calculate":
            self.handle_calculate_pm7()
            return None

        if action in ["add", "edit", "delete"]:
            self.show_in_development(action.title())
            return None

        return None

    def handle_calculate_pm7(self) -> None:
        """Handle 'Calculate PM7 Values' action.

        Shows instructions for running CREST PM7 batch calculations.
        """
        self.console.print()
        self.console.print(
            Panel(
                "[bold]CREST PM7 Batch Calculations[/bold]\n\n"
                "To run CREST PM7 calculations on your dataset:\n\n"
                "1. Prepare CSV with molecules (SMILES column)\n"
                "2. Run: [cyan]python -m grimperium.crest_pm7.batch "
                "--input batch.csv[/cyan]\n"
                "3. Results saved to: [cyan]data/molecules_pm7/computed/[/cyan]\n"
                "4. Reload this view to see updated results\n\n"
                f"[dim]Current results: {PHASE_A_RESULTS_FILE}[/dim]",
                title="[bold green]Calculate PM7 Values[/bold green]",
                border_style="green",
            )
        )
        self.console.print()
        self.console.input("[dim]Press Enter to continue...[/dim]")

    def run(self) -> str | None:
        """Run the databases view interaction loop."""
        while True:
            if self.selected_db:
                self.render_database_detail(self.selected_db)
                result = show_back_menu(
                    options=self.get_detail_menu_options(),
                    title="Actions",
                )
            else:
                self.render()
                result = show_back_menu(
                    options=self.get_menu_options(),
                    title="Select Database",
                )

            if result is None or result == "back":  # Ctrl+C
                if self.selected_db:
                    self.selected_db = None
                else:
                    return "main"
            else:
                next_view = self.handle_action(result)
                if next_view:
                    return next_view
