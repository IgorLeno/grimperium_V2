"""The one input boundary of the Semi-Imperium views.

Every question a view asks goes through :class:`Prompter`. Keeping the
prompts behind one small protocol is what lets the views be driven end
to end in tests — including the reuse decision — without a terminal, and
it keeps cancellation handling in a single place instead of scattered
through the views.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol

from rich.console import Console

from grimperium.cli.menu import MenuOption, confirm, show_menu, text_input
from grimperium.cli.styles import COLORS, ICONS

#: Value returned by :meth:`Prompter.choice` when the user backs out.
CANCEL = "__cancel__"


class Prompter(Protocol):
    """Asks the user for one value at a time, or reports a cancellation."""

    def text(self, message: str, *, default: str = "") -> str | None:
        """Ask for free text; ``None`` means the user cancelled."""
        ...

    def choice(
        self, message: str, options: Sequence[tuple[str, str]]
    ) -> str | None:
        """Ask the user to pick one ``(label, value)``; ``None`` cancels."""
        ...

    def confirm(self, message: str, *, default: bool = False) -> bool:
        """Ask a yes/no question."""
        ...

    def pause(self) -> None:
        """Hold the screen until the user acknowledges what is on it."""
        ...


class QuestionaryPrompter:
    """The interactive prompter used by the running application."""

    def __init__(self, console: Console) -> None:
        self.console = console

    def text(self, message: str, *, default: str = "") -> str | None:
        """Ask for free text through the shared questionary styling."""
        return text_input(message, default=default)

    def choice(
        self, message: str, options: Sequence[tuple[str, str]]
    ) -> str | None:
        """Show a list with an explicit cancel entry at the end."""
        if not options:
            return None
        menu = [MenuOption(label=label, value=value) for label, value in options]
        menu.append(MenuOption(label="Cancel", value=CANCEL, icon=ICONS["back"]))
        selected = show_menu(menu, title=message, console=self.console)
        if selected is None or selected == CANCEL:
            return None
        return selected

    def confirm(self, message: str, *, default: bool = False) -> bool:
        """Ask a yes/no question through the shared questionary styling."""
        return confirm(message, default=default)

    def pause(self) -> None:
        """Wait for Enter so a result stays readable before the next redraw."""
        self.console.input(
            f"[{COLORS['muted']}]Press Enter to continue...[/{COLORS['muted']}]"
        )


__all__ = ["CANCEL", "Prompter", "QuestionaryPrompter"]
