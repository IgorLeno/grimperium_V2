"""
Settings view for GRIMPERIUM CLI.

Currently marked as [IN DEVELOPMENT].
"""

import sys
from typing import Optional

from grimperium.cli.menu import MenuOption
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView


class SettingsView(BaseView):
    """View for application settings (currently in development)."""

    name = "settings"
    title = "Settings"
    icon = ICONS["settings"]
    color = COLORS["settings"]

    def render(self) -> None:
        """Render the settings view."""
        self.clear_screen()
        self.show_header()
        self.show_in_development("Settings")

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the settings view."""
        return [
            MenuOption(
                label="Select Default Model",
                value="default_model",
                disabled=True,
                disabled_reason="In Development",
            ),
            MenuOption(
                label="Theme (Light/Dark)",
                value="theme",
                disabled=True,
                disabled_reason="In Development",
            ),
            MenuOption(
                label="Display Options",
                value="display",
                disabled=True,
                disabled_reason="In Development",
            ),
            MenuOption(
                label="Data Directory",
                value="data_dir",
                disabled=True,
                disabled_reason="In Development",
            ),
        ]

    def handle_action(self, action: str) -> Optional[str]:
        """Handle menu actions."""
        if action == "back":
            return "main"
        # All actions are disabled
        self.show_in_development(action.replace("_", " ").title())
        return None

    def run(self) -> Optional[str]:
        """Run the settings view (IN DEVELOPMENT)."""
        self.render()
        
        # Only wait for user acknowledgement in interactive mode
        if sys.stdin.isatty():
            try:
                self.wait_for_enter()
            except (EOFError, KeyboardInterrupt):
                # In interactive mode, treat EOFError/Ctrl+C as continue
                pass
        # else: non-interactive (CI/piped), skip waiting
        
        return "main"
