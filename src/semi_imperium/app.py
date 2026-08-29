"""Semi-Imperium's independently launchable application shell."""

from __future__ import annotations

import sys
from dataclasses import replace

from rich.panel import Panel

from grimperium.cli.app import GrimperiumCLI
from grimperium.cli.styles import COLORS, FOOTER_NAVIGATION
from semi_imperium import __version__
from semi_imperium.menu import format_status_line, show_main_menu
from semi_imperium.settings import SemiImperiumSettings
from semi_imperium.views import (
    CalculateView,
    DatabaseView,
    HamiltonianView,
    SettingsView,
)
from semi_imperium.workspace import SemiImperiumWorkspace


class SemiImperiumCLI(GrimperiumCLI):
    """Focused shell over Semi-Imperium's own calculation and data workflows."""

    def __init__(self, workspace: SemiImperiumWorkspace | None = None) -> None:
        self.workspace = workspace or SemiImperiumWorkspace(
            settings=_default_settings()
        )
        super().__init__()

    def _register_views(self) -> None:
        """Register the three areas and CALCULATE's supporting method picker."""
        self.controller.register_view(
            "calc", CalculateView(self.controller, self.workspace)
        )
        self.controller.register_view(
            "methods", HamiltonianView(self.controller, self.workspace)
        )
        self.controller.register_view(
            "databases", DatabaseView(self.controller, self.workspace)
        )
        self.controller.register_view(
            "settings", SettingsView(self.controller, self.workspace)
        )

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

    def status_line(self) -> str:
        """Summarize the staged molecules, their requests and the store."""
        session = self.workspace.session
        settings = self.workspace.settings
        return format_status_line(
            molecule_count=len(session),
            selected_count=len(session.selected_entries),
            hamiltonians=session.default_hamiltonians,
            crest_enabled=settings.conformer_search.enabled,
            store_root=str(settings.runtime.store_root),
            tools_ready=settings.runtime.is_ready,
        )

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
        return show_main_menu(self.status_line())


def _default_settings() -> SemiImperiumSettings:
    """Root the store inside the project rather than the current directory."""
    from grimperium.cli.constants import get_project_root

    settings = SemiImperiumSettings()
    runtime = replace(
        settings.runtime,
        store_root=get_project_root() / "data" / "semi_imperium",
    )
    return replace(settings, runtime=runtime)


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
