"""
Databases view for GRIMPERIUM CLI.

Displays and manages molecular databases.
"""

from datetime import date, datetime
from typing import TYPE_CHECKING, Optional

from rich.panel import Panel
from rich.table import Table

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
        self.selected_db: Optional[Database] = None

    def _format_last_updated(self, value: Optional[datetime | date]) -> str:
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

        for db in DATABASES:
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
        for db in DATABASES:
            options.append(
                MenuOption(
                    label=f"View {db.name}",
                    value=f"view_{db.name}",
                    icon=ICONS["databases"],
                )
            )

        # Additional options (in development)
        options.extend(
            [
                MenuOption(
                    label="Calculate PM7 Values",
                    value="calculate",
                    disabled=True,
                    disabled_reason="In Development",
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
                disabled=True,
                disabled_reason="In Development",
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

    def handle_action(self, action: str) -> Optional[str]:
        """Handle menu actions."""
        if action == "back":
            if self.selected_db:
                self.selected_db = None
                return None  # Stay in databases view
            return "main"

        if action and action.startswith("view_"):
            db_name = action.removeprefix("view_")
            found = False
            for db in DATABASES:
                if db.name == db_name:
                    self.selected_db = db
                    found = True
                    break

            if not found:
                self.selected_db = None
                # Could log a warning here if logging was set up

            return None

        # Handle in-development features
        if action in ["calculate", "add", "edit", "delete"]:
            self.show_in_development(action.title())
            return None

        return None

    def run(self) -> Optional[str]:
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
