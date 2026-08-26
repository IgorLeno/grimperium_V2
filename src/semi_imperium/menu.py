"""Focused top-level menu for Semi-Imperium."""

from grimperium.cli.menu import (
    MenuOption,
    format_session_header,
    show_menu_with_separator,
)
from grimperium.cli.styles import ICONS


def show_main_menu(
    *,
    property_label: str = "Not selected",
    method_label: str = "Not selected",
    dataset_label: str = "Not selected",
    model_label: str = "No model selected",
    status: str = "No method selected",
    width: int | None = None,
) -> str | None:
    """Show the three-area Semi-Imperium application menu."""
    options = [
        MenuOption(
            label="CALCULATE",
            value="calc",
            icon=ICONS["calc"],
            description="Run a molecular calculation",
            style="calc",
        ),
        MenuOption(
            label="DATABASE",
            value="databases",
            icon=ICONS["databases"],
            description="Manage molecular data",
            style="databases",
        ),
        MenuOption(
            label="SETTINGS",
            value="settings",
            icon=ICONS["settings"],
            description="Configure calculation tools",
            style="settings",
        ),
    ]
    return show_menu_with_separator(
        options=options,
        title=format_session_header(
            property_label=property_label,
            method_label=method_label,
            dataset_label=dataset_label,
            model_label=model_label,
            status=status,
            width=width,
        ),
    )
