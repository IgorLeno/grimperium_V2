"""
About view for GRIMPERIUM CLI.

Displays application information, version, and system status.
"""

from typing import Optional

from rich.panel import Panel
from rich.table import Table

from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.mock_data import SYSTEM_INFO, get_ready_databases, get_ready_models
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView


class AboutView(BaseView):
    """View displaying application information and system status."""

    name = "about"
    title = "About"
    icon = ICONS["about"]
    color = COLORS["about"]

    def render(self) -> None:
        """Render the about information."""
        self.clear_screen()
        self.show_header()

        # Application info panel
        app_info = f"""
[bold]GRIMPERIUM[/bold]
ML-Enhanced PM7 Molecular Property Predictor

[{COLORS['muted']}]Version:[/{COLORS['muted']}]     {SYSTEM_INFO['version']}
[{COLORS['muted']}]Build Date:[/{COLORS['muted']}]  {SYSTEM_INFO['build_date']}
[{COLORS['muted']}]Python:[/{COLORS['muted']}]      {SYSTEM_INFO['python_version']}

[bold]Description[/bold]
Grimperium uses Delta-Learning to predict high-accuracy CBS
(Complete Basis Set) molecular properties by learning the
correction to semi-empirical PM7 calculations.

[bold]Core Hypothesis[/bold]
Learning Δ = (y_cbs - y_pm7) is easier than learning y_cbs directly.
"""
        self.console.print(
            Panel(
                app_info,
                title=f"[bold {COLORS['about']}]Application Info[/bold {COLORS['about']}]",
                border_style=COLORS["border"],
                padding=(1, 2),
            )
        )
        self.console.print()

        # System status table
        status_table = Table(
            title="System Status",
            show_header=True,
            header_style=f"bold {COLORS['about']}",
            border_style=COLORS["border"],
        )
        status_table.add_column("Component", style="bold")
        status_table.add_column("Status")
        status_table.add_column("Details")

        # Databases status
        db_ready = len(get_ready_databases())
        db_total = SYSTEM_INFO["databases_total"]
        
        # Determine database status: Ready / Not Ready / Partial
        if db_ready == db_total and db_total > 0:
            db_status = f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
        elif db_ready == 0:
            db_status = f"[{COLORS['error']}]{ICONS['error']} Not Ready[/{COLORS['error']}]"
        else:
            db_status = f"[{COLORS['warning']}]{ICONS['warning']} Partial[/{COLORS['warning']}]"
        
        status_table.add_row(
            f"{ICONS['databases']} Databases",
            db_status,
            f"{db_ready}/{db_total} available",
        )

        # Models status
        models_ready = len(get_ready_models())
        models_total = SYSTEM_INFO["models_total"]
        
        # Determine models status: Ready / Not Ready / Partial
        if models_ready == models_total and models_total > 0:
            models_status = f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
        elif models_ready == 0:
            models_status = f"[{COLORS['error']}]{ICONS['error']} Not Ready[/{COLORS['error']}]"
        else:
            models_status = f"[{COLORS['warning']}]{ICONS['warning']} Partial[/{COLORS['warning']}]"
        
        status_table.add_row(
            f"{ICONS['models']} Models",
            models_status,
            f"{models_ready}/{models_total} trained",
        )

        # Default model
        status_table.add_row(
            f"{ICONS['calc']} Default Model",
            f"[{COLORS['success']}]{ICONS['success']} Active[/{COLORS['success']}]",
            SYSTEM_INFO["default_model"],
        )

        self.console.print(status_table)
        self.console.print()

        # Links panel
        links = """
[bold]Documentation[/bold]
  https://grimperium.readthedocs.io

[bold]GitHub Repository[/bold]
  https://github.com/grimperium/grimperium

[bold]Report Issues[/bold]
  https://github.com/grimperium/grimperium/issues
"""
        self.console.print(
            Panel(
                links,
                title=f"[bold {COLORS['about']}]Links[/bold {COLORS['about']}]",
                border_style=COLORS["border"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the about view."""
        return []  # Only back option (added by show_back_menu)

    def handle_action(self, action: Optional[str]) -> Optional[str]:
        """Handle menu actions."""
        if action == "back" or action is None:
            return "main"
        return None

    def run(self) -> Optional[str]:
        """Run the about view interaction loop."""
        self.render()

        result = show_back_menu(
            options=[],
            title="",
        )

        return self.handle_action(result)
