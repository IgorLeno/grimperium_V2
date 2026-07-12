"""
Menu system for GRIMPERIUM CLI.

This module provides menu rendering and selection using questionary.
"""

from collections.abc import Callable
from dataclasses import dataclass

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
    console: Console | None = None,
) -> str | None:
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
        return str(result) if result is not None else None
    except KeyboardInterrupt:
        return None


def show_menu_with_separator(
    options: list[MenuOption],
    title: str = "",
    separator_after: list[int] | None = None,
) -> str | None:
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
        return str(result) if result is not None else None
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
        result = questionary.confirm(
            message=message,
            default=default,
            style=QUESTIONARY_STYLE,
            qmark="",
        ).ask()
        return bool(result)  # Coerce None to False
    except KeyboardInterrupt:
        return False


def text_input(
    message: str,
    default: str = "",
    validate: Callable[[str], bool | str] | None = None,
) -> str | None:
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
        result = questionary.text(
            message=message,
            default=default,
            style=QUESTIONARY_STYLE,
            qmark="",
            validate=validate,
        ).ask()
        return str(result) if result is not None else None
    except KeyboardInterrupt:
        return None


def format_session_header(
    *,
    property_label: str = "Not selected",
    method_label: str = "Not selected",
    dataset_label: str = "Not selected",
    model_label: str = "No model selected",
    status: str = "No method selected",
    width: int | None = None,
) -> str:
    """Format the main-menu session context header."""
    if width is not None:
        available = max(8, (width - 72) // 5)
        property_label = _truncate_label(property_label, available)
        method_label = _truncate_label(method_label, available)
        dataset_label = _truncate_label(dataset_label, available)
        model_label = _truncate_label(model_label, available)
        status = _truncate_label(status, available)
    return (
        f"[Property: {property_label} | Method: {method_label} | "
        f"Dataset: {dataset_label} | Model: {model_label} | Status: {status}]"
    )


def _truncate_label(value: str, max_length: int) -> str:
    if len(value) <= max_length:
        return value
    if max_length <= 3:
        return "." * max_length
    return value[: max_length - 3] + "..."


def show_main_menu(
    *,
    property_label: str = "Not selected",
    method_label: str = "Not selected",
    dataset_label: str = "Not selected",
    model_label: str = "No model selected",
    status: str = "No method selected",
    current_model: str | None = None,
    width: int | None = None,
) -> str | None:
    """
    Display the main application menu.

    Args:
        property_label: Active scientific property display name
        method_label: Active calculation method display name
        model_label: Active model display name (or Not required)
        status: System status string
        current_model: Deprecated alias for model_label (tests/compat)

    Returns:
        Selected menu action, or None if cancelled
    """
    if current_model is not None:
        model_label = current_model

    options = [
        MenuOption(
            label="CALCULATE",
            value="calc",
            icon=ICONS["calc"],
            description="Run calculation with active method",
            style="calc",
        ),
        MenuOption(
            label="CALCULATION METHODS",
            value="methods",
            icon=ICONS["methods"],
            description="Select and configure calculation methods",
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
        title=format_session_header(
            property_label=property_label,
            method_label=method_label,
            dataset_label=dataset_label,
            model_label=model_label,
            status=status,
            width=width,
        ),
        separator_after=[4],  # Separator after RESULTS
    )


def show_back_menu(
    options: list[MenuOption],
    title: str = "",
) -> str | None:
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

    # Only add separator if there are options before "Back"
    if len(options) == 0:
        separator_after: list[int] = []
    else:
        separator_after = [len(options) - 1]

    return show_menu_with_separator(
        options=all_options,
        title=title,
        separator_after=separator_after,
    )
