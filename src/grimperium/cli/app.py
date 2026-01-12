"""
Main application for GRIMPERIUM CLI.

This module provides the main entry point and application loop.
"""

import sys
from typing import Optional

from rich.panel import Panel

from grimperium.cli.controller import CliController
from grimperium.cli.menu import show_main_menu
from grimperium.cli.styles import (
    BANNER,
    COLORS,
    FOOTER_NAVIGATION,
    SUBTITLE,
    VERSION_LINE,
)
from grimperium.cli.views import (
    AboutView,
    CalcView,
    DatabasesView,
    ModelsView,
    ResultsView,
    SettingsView,
)


class GrimperiumCLI:
    """
    Main CLI application class.

    Manages the application lifecycle, view registration, and main loop.
    """

    def __init__(self) -> None:
        """Initialize the CLI application."""
        self.controller = CliController()
        self.console = self.controller.console
        self._register_views()

    def _register_views(self) -> None:
        """Register all views with the controller."""
        self.controller.register_view("calc", CalcView(self.controller))
        self.controller.register_view("databases", DatabasesView(self.controller))
        self.controller.register_view("models", ModelsView(self.controller))
        self.controller.register_view("results", ResultsView(self.controller))
        self.controller.register_view("settings", SettingsView(self.controller))
        self.controller.register_view("about", AboutView(self.controller))

    def show_welcome(self) -> None:
        """Display the welcome screen with ASCII banner."""
        self.console.clear()
        self.console.print()

        # Banner panel
        banner_content = f"{BANNER}\n{SUBTITLE}\n{VERSION_LINE}"
        self.console.print(
            Panel(
                banner_content,
                border_style=COLORS["calc"],
                padding=(1, 4),
            )
        )
        self.console.print()

    def display_main_menu(self) -> Optional[str]:
        """
        Display and handle the main menu.

        Returns:
            Selected menu option or None if cancelled
        """
        self.console.print(
            Panel(
                f"[bold]MAIN MENU[/bold]\n\n{FOOTER_NAVIGATION}",
                border_style=COLORS["border"],
                padding=(0, 2),
            )
        )
        self.console.print()

        return show_main_menu(
            current_model=self.controller.current_model,
            status=self.controller.status,
        )

    def run_view(self, view_name: str) -> Optional[str]:
        """
        Run a specific view.

        Args:
            view_name: Name of the view to run

        Returns:
            Next view to navigate to, or None to stay at main menu
        """
        view = self.controller.get_view(view_name)
        if view is None:
            self.console.print(f"[error]View '{view_name}' not found[/error]")
            return None

        try:
            return view.run()
        except KeyboardInterrupt:
            return "main"

    def run(self) -> int:
        """
        Run the main application loop.

        Returns:
            Exit code (0 for success, non-zero for errors)
        """
        self.controller.start()

        try:
            while self.controller.is_running():
                # Show welcome screen and main menu
                self.show_welcome()
                selection = self.display_main_menu()

                if selection is None:
                    # User pressed Ctrl+C at main menu - just break, stop() in finally
                    break

                # Navigate to selected view
                if selection == "main":
                    continue

                # Run the selected view
                next_view = self.run_view(selection)

                # Handle navigation from view
                while next_view and next_view != "main":
                    next_view = self.run_view(next_view)

        except KeyboardInterrupt:
            self.console.print()
            self.console.print("[muted]Interrupted by user[/muted]")

        except Exception as e:
            self.console.print()
            self.console.print(f"[error]Unexpected error: {e}[/error]")
            return 1

        finally:
            self.controller.stop()

        return 0


def main() -> int:
    """
    Main entry point for the GRIMPERIUM CLI.

    Returns:
        Exit code (0 for success, non-zero for errors)
    """
    app = GrimperiumCLI()
    return app.run()


if __name__ == "__main__":
    sys.exit(main())
