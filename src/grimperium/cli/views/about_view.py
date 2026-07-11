"""
About view for GRIMPERIUM CLI.

Displays application information, version, and real session status.
"""

from __future__ import annotations

import sys

from rich.panel import Panel
from rich.table import Table

from grimperium import __version__
from grimperium.cli.menu import MenuOption, show_back_menu
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

        py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        app_info = f"""
[bold]GRIMPERIUM[/bold]
ML-Enhanced Molecular Property Calculation

[{COLORS['muted']}]Version:[/{COLORS['muted']}]     {__version__}
[{COLORS['muted']}]Python:[/{COLORS['muted']}]      {py_version}+

[bold]Description[/bold]
Grimperium runs declarative calculation methods (semiempirical and
delta-learning) to estimate molecular thermochemical properties.

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

        summary = self.controller.session_summary()
        status_table = Table(
            title="Session Status",
            show_header=True,
            header_style=f"bold {COLORS['about']}",
            border_style=COLORS["border"],
        )
        status_table.add_column("Component", style="bold")
        status_table.add_column("Value")

        status_table.add_row(f"{ICONS['calc']} Property", summary["property"])
        status_table.add_row(f"{ICONS['methods']} Method", summary["method"])
        status_table.add_row(f"{ICONS['models']} Model", summary["model"])
        status_table.add_row("Status", summary["status"])

        self.console.print(status_table)
        self.console.print()

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
        return []

    def handle_action(self, action: str | None) -> str | None:
        """Handle menu actions."""
        if action == "back" or action is None:
            return "main"
        return None

    def run(self) -> str | None:
        """Run the about view interaction loop."""
        self.render()

        result = show_back_menu(
            options=[],
            title="",
        )

        return self.handle_action(result)
