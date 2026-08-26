"""Semi-Imperium's independently launchable application shell."""

from __future__ import annotations

import sys

from rich.panel import Panel

from grimperium.cli.app import GrimperiumCLI
from grimperium.cli.styles import COLORS, FOOTER_NAVIGATION
from grimperium.cli.views import CalcView, DatabasesView, MethodsView, SettingsView
from semi_imperium import __version__
from semi_imperium.menu import show_main_menu


class SemiImperiumCLI(GrimperiumCLI):
    """Focused shell backed by Grimperium's calculation and data services."""

    def _register_views(self) -> None:
        """Register the three areas and CALCULATE's supporting method picker."""
        self.controller.register_view("calc", CalcView(self.controller))
        self.controller.register_view("methods", MethodsView(self.controller))
        self.controller.register_view("databases", DatabasesView(self.controller))
        self.controller.register_view("settings", SettingsView(self.controller))

    def show_welcome(self) -> None:
        """Display the Semi-Imperium identity without changing Grimperium's UI."""
        self.console.clear()
        self.console.print()
        self.console.print(
            Panel(
                "[bold cyan]SEMI-IMPERIUM[/bold cyan]\n"
                "[muted]Focused molecular calculation workspace[/muted]\n"
                f"[muted]v{__version__}[/muted]",
                border_style=COLORS["calc"],
                padding=(1, 4),
            )
        )
        self.console.print()

    def display_main_menu(self) -> str | None:
        """Display the focused three-area top-level menu."""
        self.console.print(
            Panel(
                f"[bold]MAIN MENU[/bold]\n\n{FOOTER_NAVIGATION}",
                border_style=COLORS["border"],
                padding=(0, 2),
            )
        )
        self.console.print()

        summary = self.controller.session_summary()
        return show_main_menu(
            property_label=summary["property"],
            method_label=summary["method"],
            dataset_label=summary["dataset"],
            model_label=summary["model"],
            status=summary["status"],
            width=self.console.size.width,
        )


def main() -> int:
    """Launch the Semi-Imperium interactive shell."""
    from grimperium.cli.constants import get_project_root
    from grimperium.utils.logging import setup_logging

    log_file = get_project_root() / "logs" / "semi-imperium.log"
    setup_logging(
        level="DEBUG",
        log_file=log_file,
        console_level="WARNING",
    )

    skip_preflight = "--skip-preflight" in sys.argv
    app = SemiImperiumCLI()
    return app.run(skip_preflight=skip_preflight)
