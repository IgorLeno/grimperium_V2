"""
CLI Controller for GRIMPERIUM.

Manages navigation state and view transitions.
"""

from typing import TYPE_CHECKING, Any, Optional

from rich.console import Console

from grimperium.cli.mock_data import DEFAULT_MODEL
from grimperium.cli.styles import CLI_THEME

if TYPE_CHECKING:
    from grimperium.cli.views.base_view import BaseView


class CliController:
    """
    Manages navigation state and view transitions.

    The controller is the central hub for:
    - Maintaining navigation history (breadcrumbs)
    - Managing the current view state
    - Providing a shared console for rendering
    - Tracking application-wide settings
    """

    def __init__(self) -> None:
        """Initialize the CLI controller."""
        self.history: list[str] = []
        self.current_view: str = "main"
        self.current_model: str = DEFAULT_MODEL
        self.status: str = "Ready"
        self.console = Console(theme=CLI_THEME)
        self._views: dict[str, Any] = {}
        self._running: bool = False

    def register_view(self, name: str, view: "BaseView") -> None:  # type: ignore
        """
        Register a view with the controller.

        Args:
            name: View identifier (e.g., "calc", "databases")
            view: View instance
        """
        self._views[name] = view

    def get_view(self, name: str) -> Optional["BaseView"]:  # type: ignore
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
        Navigate to a different view.

        Args:
            view: Name of the view to navigate to
        """
        if self.current_view != "main":
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

    def set_model(self, model_name: str) -> None:
        """
        Set the current active model.

        Args:
            model_name: Name of the model to set as active
        """
        self.current_model = model_name

    def set_status(self, status: str) -> None:
        """
        Set the application status.

        Args:
            status: Status string (e.g., "Ready", "Processing")
        """
        self.status = status

    def is_running(self) -> bool:
        """Check if the application is running."""
        return self._running

    def start(self) -> None:
        """Start the application."""
        self._running = True

    def stop(self) -> None:
        """Stop the application."""
        self._running = False
        self.console.print()
        self.console.print("[muted]Goodbye![/muted]")
