"""
CLI Controller for GRIMPERIUM.

Manages navigation state and view transitions.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from rich.console import Console

from grimperium.cli.model_compatibility import (
    ModelCompatibilityError,
    validate_model_for_method,
)
from grimperium.cli.settings_manager import SettingsManager
from grimperium.cli.styles import CLI_THEME

if TYPE_CHECKING:
    from grimperium.calculation.methods import CalculationMethodDefinition
    from grimperium.cli.views.base_view import BaseView


class CliController:
    """
    Manages navigation state and view transitions.

    The controller is the central hub for:
    - Maintaining navigation history (breadcrumbs)
    - Managing the current view state
    - Providing a shared console for rendering
    - Tracking application-wide settings and scientific session context
    """

    def __init__(self) -> None:
        """Initialize the CLI controller."""
        self.history: list[str] = []
        self.current_view: str = "main"
        self.current_model: str | None = None
        self.current_model_path: Path | None = None
        self.current_csv_path: Path | None = None
        self.current_property_id: str | None = None
        self.current_method_id: str | None = None
        self.current_method_version: str | None = None
        self.current_method_definition: CalculationMethodDefinition | None = None
        self.status: str = "No method selected"
        self.console = Console(theme=CLI_THEME)
        self.settings_manager = SettingsManager(console=self.console)
        self.settings_manager.load_from_file()
        self._views: dict[str, BaseView] = {}
        self._running: bool = False

    def register_view(self, name: str, view: BaseView) -> None:
        """
        Register a view with the controller.

        Args:
            name: View identifier (e.g., "calc", "databases")
            view: View instance
        """
        self._views[name] = view

    def get_view(self, name: str) -> BaseView | None:
        """
        Get a registered view by name.

        Args:
            name: View identifier

        Returns:
            View instance or None if not found
        """
        return self._views.get(name)

    def navigate_to(self, view: str) -> None:
        """
        Navigate to a view.

        - If view is the current view, do nothing.
        - If view is "main", clear history and set current_view.
        - Otherwise, append current_view to history and navigate.
        """
        # Guard: same view, no-op
        if view == self.current_view:
            return

        # Special case: navigate to main
        if view == "main":
            self.history.clear()
        else:
            # Normal case: push current to history
            self.history.append(self.current_view)

        self.current_view = view

    def go_back(self) -> bool:
        """
        Navigate back to the previous view.

        Returns:
            True if navigated back, False if at root (main menu)
        """
        if self.history:
            self.current_view = self.history.pop()
            return True
        else:
            self.current_view = "main"
            return False

    def go_to_main(self) -> None:
        """Navigate directly to the main menu."""
        self.history.clear()
        self.current_view = "main"

    def get_breadcrumb(self) -> str:
        """
        Get the current navigation breadcrumb.

        Returns:
            Breadcrumb string (e.g., "Main > Databases > CBS Reference")
        """
        parts = ["Main"]
        for view_name in self.history:
            view = self.get_view(view_name)
            if view:
                parts.append(view.title)
            else:
                parts.append(view_name.title())

        current = self.get_view(self.current_view)
        if current and self.current_view != "main":
            parts.append(current.title)

        return " > ".join(parts)

    def set_model(self, model_name: str, model_path: Path | None = None) -> None:
        """
        Set the current active model and optional path to its file.

        Args:
            model_name: Name of the model to set as active
            model_path: Optional path to the model file
        """
        self.current_model = model_name
        self.current_model_path = model_path
        self._refresh_status()

    def set_method(self, method: CalculationMethodDefinition) -> None:
        """Set the active calculation method for this session."""
        self.current_method_definition = method
        self.current_method_id = method.method_id
        self.current_method_version = method.version
        self.current_property_id = method.property_id
        self._refresh_status()

    def clear_method(self) -> None:
        """Clear the active calculation method."""
        self.current_method_definition = None
        self.current_method_id = None
        self.current_method_version = None
        self.current_property_id = None
        self._refresh_status()

    def _refresh_status(self) -> None:
        """Recompute session status from method and model context."""
        method = self.current_method_definition
        if method is None:
            self.status = "No method selected"
            return

        if not method.model_requirement.model_required:
            self.status = "Ready"
            return

        model_path = self.current_model_path
        if model_path is None or not model_path.exists():
            self.status = "Model required"
            return

        try:
            validate_model_for_method(model_path, method)
        except ModelCompatibilityError:
            self.status = "Model incompatible"
            return

        self.status = "Ready"

    def session_summary(self) -> dict[str, str]:
        """Return display strings for the main-menu session header."""
        method = self.current_method_definition

        if method is not None:
            property_label = method.property_name
            method_label = method.display_name
        else:
            property_label = "Not selected"
            method_label = "Not selected"

        if method is not None and not method.model_requirement.model_required:
            model_label = "Not required"
        elif self.current_model:
            model_label = self.current_model
        else:
            model_label = "No model selected"

        return {
            "property": property_label,
            "method": method_label,
            "dataset": "Not selected",
            "model": model_label,
            "status": self.status,
        }

    def set_csv(self, csv_path: Path | None) -> None:
        """
        Set the current active CSV path for analysis.

        Args:
            csv_path: Optional path to the active CSV file
        """
        self.current_csv_path = csv_path

    def set_status(self, status: str) -> None:
        """
        Set the application status.

        Args:
            status: Status string (e.g., "Ready", "Processing")
        """
        self.status = status

    @property
    def running(self) -> bool:
        """Check if the application is running."""
        return self._running

    def is_running(self) -> bool:
        """Check if the application is running (compatibility wrapper)."""
        return self.running

    def start(self) -> None:
        """Start the application."""
        self._running = True

    def stop(self) -> None:
        """Stop the application."""
        self._running = False
        self.console.print()
        self.console.print("[muted]Goodbye![/muted]")
