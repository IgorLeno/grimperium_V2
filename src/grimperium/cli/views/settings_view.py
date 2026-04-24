"""
Settings view for GRIMPERIUM CLI.

Provides interactive configuration for CREST, MOPAC, and xTB settings.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from rich.panel import Panel

from grimperium.cli.menu import MenuOption
from grimperium.cli.settings_manager import SettingsManager
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController


class SettingsView(BaseView):
    """View for configuring CREST, MOPAC, and xTB settings.

    Provides interactive menus for adjusting computation parameters
    used in the molecular optimization pipeline.
    """

    name = "settings"
    title = "Settings"
    icon = ICONS["settings"]
    color = COLORS["settings"]

    def __init__(self, controller: CliController) -> None:
        """Initialize the settings view.

        Args:
            controller: The CLI controller managing navigation.
        """
        super().__init__(controller)
        self._settings_manager: SettingsManager | None = None

    @property
    def settings_manager(self) -> SettingsManager:
        """Get or create the settings manager instance.

        Returns:
            SettingsManager instance for this view.
        """
        if self._settings_manager is None:
            controller_manager = getattr(self.controller, "settings_manager", None)
            if controller_manager is not None:
                self._settings_manager = controller_manager
            else:
                self._settings_manager = SettingsManager(console=self.console)
        return self._settings_manager

    def render(self) -> None:
        """Render the settings view with current configuration summary."""
        self.console.print()
        self.console.print(
            Panel(
                "[dim]Configure computational parameters for the molecular "
                "optimization pipeline.[/dim]",
                title="[bold]Settings Overview[/bold]",
                border_style=self.color,
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the settings view.

        Returns:
            List of MenuOption objects for CREST, MOPAC configuration.
        """
        return [
            MenuOption(
                label="CREST Configuration",
                value="crest",
                icon="🔬",
                description="Conformer search settings (includes xTB pre-optimization)",
            ),
            MenuOption(
                label="MOPAC Configuration",
                value="mopac",
                icon="⚗️",
                description="PM7 optimization settings",
            ),
            MenuOption(
                label="View Current Settings",
                value="view_all",
                icon="📋",
                description="Display all current settings",
            ),
            MenuOption(
                label="Save as Default",
                value="save_as_default",
                icon="💾",
                description="Persist current settings to disk for future sessions",
            ),
            MenuOption(
                label="Reset All to Defaults",
                value="reset_all",
                icon="🔄",
                description="Reset all settings to default values",
            ),
            MenuOption(
                label="Distributed Settings",
                value="distributed",
                icon="🌐",
                description="Manage calculation profiles and default dispatch settings",
            ),
            MenuOption(
                label="Atualizar Documentação do Projeto",
                value="update_docs",
                icon="📊",
                description="Atualizar a seção dinâmica do README com dados reais",
            ),
            MenuOption(
                label="Back to Main Menu",
                value="back",
                icon="↩️",
                description="Return to main menu",
            ),
        ]

    def handle_action(self, action: str) -> str | None:
        """Handle menu actions for settings.

        Args:
            action: The action value from the selected menu option.

        Returns:
            Next view name to navigate to, or None to stay in current view.
        """
        if action == "back":
            return "main"
        if action == "crest":
            self.settings_manager.display_crest_menu()
            return None
        if action == "mopac":
            self.settings_manager.display_mopac_menu()
            return None
        if action == "view_all":
            self._display_all_settings()
            return None
        if action == "save_as_default":
            self._save_as_default()
            return None
        if action == "reset_all":
            self._reset_all_settings()
            return None
        if action == "update_docs":
            self.settings_manager.display_readme_updater_menu()
            return None
        if action == "distributed":
            self.settings_manager.display_distributed_menu()
            return None
        return None

    def _display_all_settings(self) -> None:
        """Display a summary of all current settings."""
        self.console.clear()
        self.console.print()

        self.console.print(
            Panel(
                self.settings_manager.show_crest_summary(),
                title="[bold]CREST Settings[/bold]",
                border_style=self.color,
            )
        )

        self.console.print(
            Panel(
                self.settings_manager.show_mopac_summary(),
                title="[bold]MOPAC Settings[/bold]",
                border_style=self.color,
            )
        )

        self.console.print(
            Panel(
                self.settings_manager.show_xtb_summary(),
                title="[bold]xTB Settings[/bold]",
                border_style=self.color,
            )
        )

        self.console.print()
        self.wait_for_enter()

    def _save_as_default(self) -> None:
        """Persist current settings to disk."""
        if self.settings_manager.save_to_file():
            self.show_success("Settings saved as default.")
        else:
            self.console.print("[red]Failed to save settings.[/red]")
        self.wait_for_enter()

    def _reset_all_settings(self) -> None:
        """Reset all settings to defaults with confirmation."""
        from grimperium.cli.menu import confirm

        if confirm("Reset all settings to defaults?", default=False):
            self.settings_manager.reset_all()
            self.settings_manager.delete_config_file()
            self.show_success("All settings reset to defaults.")
        else:
            self.console.print("[dim]Reset cancelled.[/dim]")
        self.wait_for_enter()
