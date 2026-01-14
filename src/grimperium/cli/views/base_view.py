"""
Base view class for GRIMPERIUM CLI.

All view modules inherit from this abstract base class.
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, ClassVar

from rich.console import Console
from rich.panel import Panel

from grimperium.cli.menu import MenuOption
from grimperium.cli.styles import COLORS, ICONS

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController


class BaseView(ABC):
    """
    Abstract base class for all view modules.

    Each view is responsible for:
    - Rendering its content using Rich
    - Providing menu options for navigation
    - Handling user interactions within the view
    """

    # View metadata (override in subclasses)
    name: ClassVar[str] = "base"
    title: ClassVar[str] = "Base View"
    icon: ClassVar[str] = ""
    color: ClassVar[str] = "#FFFFFF"

    def __init__(self, controller: "CliController") -> None:
        """
        Initialize the view.

        Args:
            controller: The CLI controller managing navigation
        """
        self.controller = controller
        self.console: Console = controller.console

    @abstractmethod
    def render(self) -> None:
        """
        Render the view content.

        This method should use self.console to display content.
        """
        pass

    @abstractmethod
    def get_menu_options(self) -> list[MenuOption]:
        """
        Return menu options available in this view.

        Returns:
            List of MenuOption objects for the view's menu
        """
        pass

    @abstractmethod
    def handle_action(self, action: str) -> str | None:
        """
        Handle a selected menu action.

        Args:
            action: The action value from the selected menu option

        Returns:
            Next view name to navigate to, or None to stay in current view
        """
        pass

    def run(self) -> str | None:
        """
        Run the view: render content, show menu, handle action.

        Returns:
            Next view name to navigate to, or None if no options/user cancelled
        """
        from grimperium.cli.menu import show_menu

        while True:
            self.clear_screen()
            self.show_header()
            self.render()

            options = self.get_menu_options()
            if not options:
                return None

            choice = show_menu(options)
            if choice is None:
                return None

            result = self.handle_action(choice)
            if result is not None:
                return result

    def show_header(self) -> None:
        """Display the view header with title and icon."""
        header_text = f"{self.icon}  {self.title}" if self.icon else self.title
        self.console.print()
        self.console.print(f"[bold {self.color}]{header_text}[/bold {self.color}]")

        # Adaptive separator width based on console width
        width = max(20, min(self.console.width, 80))
        self.console.print(f"[{COLORS['muted']}]{'─' * width}[/{COLORS['muted']}]")
        self.console.print()

    def show_in_development(self, feature: str) -> None:
        """
        Show an [IN DEVELOPMENT] message for incomplete features.

        Args:
            feature: Name of the feature that's in development
        """
        self.console.print()
        self.console.print(
            Panel(
                f"[{COLORS['in_dev']}]{ICONS['in_dev']} {feature} is under development.\n\n"
                f"This feature will be available in a future update.\n"
                f"Please check back soon![/{COLORS['in_dev']}]",
                title=f"[{COLORS['warning']}]In Development[/{COLORS['warning']}]",
                border_style=COLORS["in_dev"],
                padding=(1, 2),
            )
        )
        self.console.print()
        self.wait_for_enter()

    def show_error(self, message: str) -> None:
        """
        Display an error message.

        Args:
            message: Error message to display
        """
        self.console.print()
        self.console.print(
            Panel(
                f"[{COLORS['error']}]{ICONS['error']} {message}[/{COLORS['error']}]",
                title=f"[{COLORS['error']}]Error[/{COLORS['error']}]",
                border_style=COLORS["error"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def show_success(self, message: str) -> None:
        """
        Display a success message.

        Args:
            message: Success message to display
        """
        self.console.print()
        self.console.print(
            Panel(
                f"[{COLORS['success']}]{ICONS['success']} {message}[/{COLORS['success']}]",
                title=f"[{COLORS['success']}]Success[/{COLORS['success']}]",
                border_style=COLORS["success"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def clear_screen(self) -> None:
        """Clear the terminal screen."""
        self.console.clear()

    def wait_for_enter(self) -> None:
        """Wait for the user to press Enter."""
        self.console.input(
            f"[{COLORS['muted']}]Press Enter to continue...[/{COLORS['muted']}]"
        )
