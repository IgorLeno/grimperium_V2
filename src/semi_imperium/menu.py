"""Focused top-level menu for Semi-Imperium.

The header deliberately does not reuse Grimperium's session line: this
application has no model, no dataset and no generic property to report.
What matters at the top level is how many molecules are staged, what
they are asking for, and where their results are being written.
"""

from grimperium.cli.menu import MenuOption, show_menu_with_separator
from grimperium.cli.styles import ICONS


def show_main_menu(status_line: str = "") -> str | None:
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
    return show_menu_with_separator(options=options, title=status_line)


def format_status_line(
    *,
    molecule_count: int,
    selected_count: int,
    hamiltonians: tuple[str, ...],
    crest_enabled: bool,
    store_root: str,
    tools_ready: bool,
) -> str:
    """Format the focused status header shown above the three areas."""
    requested = ", ".join(hamiltonians) if hamiltonians else "none"
    return (
        f"[Molecules: {molecule_count} ({selected_count} selected) | "
        f"Hamiltonians: {requested} | "
        f"CREST: {'on' if crest_enabled else 'off'} | "
        f"Store: {store_root} | "
        f"Tools: {'ready' if tools_ready else 'incomplete'}]"
    )
