"""
Menu system for GRIMPERIUM CLI.

This module provides menu rendering and selection using questionary.
"""

from dataclasses import dataclass
from typing import Callable, Optional

import questionary
from prompt_toolkit.styles import Style
from questionary import Choice, Separator
from rich.console import Console

from grimperium.cli.styles import COLORS, ICONS


@dataclass
class MenuOption:
    """Represents a menu option."""

    label: str
    value: str
    icon: str = ""
    description: str = ""
    disabled: bool = False
    disabled_reason: str = "In Development"
    style: str = ""  # Color name from COLORS


# Custom questionary style matching our theme
QUESTIONARY_STYLE = Style.from_dict(
    {
        "qmark": COLORS["calc"],
        "question": "bold",
        "answer": COLORS["success"],
        "pointer": f"{COLORS['calc']} bold",
        "highlighted": f"{COLORS['calc']} bold",
        "selected": COLORS["success"],
        "separator": COLORS["muted"],
        "instruction": COLORS["muted"],
        "text": "",
        "disabled": COLORS["in_dev"],
    }
)


def create_choice(option: MenuOption) -> Choice:
    """Create a questionary Choice from a MenuOption."""
    label = option.label
    if option.icon:
        label = f"{option.icon}  {label}"
    if option.disabled:
        label = f"{label} [{option.disabled_reason}]"

    return Choice(
        title=label,
        value=option.value,
        disabled=option.disabled_reason if option.disabled else None,
    )


def show_menu(
    options: list[MenuOption],
    title: str = "",
    instruction: str = "",
    pointer: str = ICONS["arrow"],
    console: Optional[Console] = None,
) -> Optional[str]:
    """
    Display an interactive menu and return the selected value.

    Args:
        options: List of MenuOption objects
        title: Menu title/question
        instruction: Navigation hint text
        pointer: Pointer symbol for current selection
        console: Rich console for rendering (optional)

    Returns:
        Selected option value, or None if cancelled (Ctrl+C)
    """
    if not options:
        return None

    choices = [create_choice(opt) for opt in options]

    if not instruction:
        instruction = "(↑↓ to move, Enter to select)"

    try:
        result = questionary.select(
            message=title,
            choices=choices,
            style=QUESTIONARY_STYLE,
            qmark="",
            pointer=pointer,
            instruction=instruction,
            use_arrow_keys=True,
            use_jk_keys=True,
            use_shortcuts=False,
        ).ask()
        return result
    except KeyboardInterrupt:
        return None


def show_menu_with_separator(
    options: list[MenuOption],
    title: str = "",
    separator_after: Optional[list[int]] = None,
) -> Optional[str]:
    """
    Display a menu with separators between option groups.

    Args:
        options: List of MenuOption objects
        title: Menu title/question
        separator_after: Indices after which to add separators

    Returns:
        Selected option value, or None if cancelled
    """
    if not options:
        return None

    separator_after = separator_after or []
    choices: list[Choice | Separator] = []

    for i, opt in enumerate(options):
        choices.append(create_choice(opt))
        if i in separator_after:
            choices.append(Separator())

    try:
        result = questionary.select(
            message=title,
            choices=choices,
            style=QUESTIONARY_STYLE,
            qmark="",
            pointer=ICONS["arrow"],
            instruction="(↑↓ to move, Enter to select)",
            use_arrow_keys=True,
            use_jk_keys=True,
        ).ask()
        return result
    except KeyboardInterrupt:
        return None


def confirm(
    message: str,
    default: bool = False,
) -> bool:
    """
    Ask for confirmation (yes/no).

    Args:
        message: Confirmation question
        default: Default answer if Enter is pressed

    Returns:
        True if confirmed, False otherwise
    """
    try:
        return questionary.confirm(
            message=message,
            default=default,
            style=QUESTIONARY_STYLE,
            qmark="",
        ).ask()
    except KeyboardInterrupt:
        return False


def text_input(
    message: str,
    default: str = "",
    validate: Optional[Callable[[str], bool | str]] = None,
) -> Optional[str]:
    """
    Get text input from the user.

    Args:
        message: Input prompt
        default: Default value
        validate: Validation function (returns True or error message)

    Returns:
        User input, or None if cancelled
    """
    try:
        return questionary.text(
            message=message,
            default=default,
            style=QUESTIONARY_STYLE,
            qmark="",
            validate=validate,
        ).ask()
    except KeyboardInterrupt:
        return None


def show_main_menu(
    current_model: str = "DeltaXGB_v1.0",
    status: str = "Ready",
) -> Optional[str]:
    """
    Display the main application menu.

    Args:
        current_model: Name of the currently selected model
        status: System status string

    Returns:
        Selected menu action, or None if cancelled
    """
    options = [
        MenuOption(
            label="CALC",
            value="calc",
            icon=ICONS["calc"],
            description="Predict molecular properties",
            style="calc",
        ),
        MenuOption(
            label="DATABASES",
            value="databases",
            icon=ICONS["databases"],
            description="Manage databases",
            style="databases",
        ),
        MenuOption(
            label="MODELS",
            value="models",
            icon=ICONS["models"],
            description="ML model management",
            style="models",
        ),
        MenuOption(
            label="RESULTS",
            value="results",
            icon=ICONS["results"],
            description="Performance analytics",
            style="results",
        ),
        MenuOption(
            label="SETTINGS",
            value="settings",
            icon=ICONS["settings"],
            description="Configuration",
            disabled=True,
            disabled_reason="In Development",
            style="settings",
        ),
        MenuOption(
            label="ABOUT",
            value="about",
            icon=ICONS["about"],
            description="Info & help",
            style="about",
        ),
    ]

    return show_menu_with_separator(
        options=options,
        title=f"[Model: {current_model} | Status: {status}]",
        separator_after=[3],  # Separator after RESULTS
    )


def show_back_menu(
    options: list[MenuOption],
    title: str = "",
) -> Optional[str]:
    """
    Display a menu with a "Back" option at the end.

    Args:
        options: List of MenuOption objects
        title: Menu title/question

    Returns:
        Selected option value, or None if cancelled
    """
    back_option = MenuOption(
        label="Back",
        value="back",
        icon=ICONS["back"],
    )

    all_options = options + [back_option]

    return show_menu_with_separator(
        options=all_options,
        title=title,
        separator_after=[len(options) - 1],  # Separator before Back
    )
