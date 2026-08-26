"""Load the optional questionary dependency without breaking CLI imports.

The interactive preflight reports a missing CLI extra before prompting. This
fallback keeps non-interactive commands and library imports available in lean
installations while preserving the real questionary API whenever it is
installed.
"""

from __future__ import annotations

from typing import Any, NoReturn


class _MissingQuestionaryPrompt:
    """Prompt returned when the optional interactive dependency is absent."""

    def ask(self) -> NoReturn:
        raise ModuleNotFoundError(
            "Interactive CLI prompts require the 'cli' extra: "
            "install grimperium[cli]."
        )


class _FallbackChoice:
    """Minimal choice value used while questionary is unavailable."""

    def __init__(
        self,
        title: str,
        value: Any = None,
        disabled: str | None = None,
    ) -> None:
        self.title = title
        self.value = value
        self.disabled = disabled


class _FallbackSeparator:
    """Minimal separator value used while questionary is unavailable."""

    def __init__(self, title: str = "────────") -> None:
        self.title = title


class _QuestionaryFallback:
    """Small import-safe subset of questionary's public module API."""

    Choice = _FallbackChoice
    Separator = _FallbackSeparator

    @staticmethod
    def select(*_args: Any, **_kwargs: Any) -> _MissingQuestionaryPrompt:
        return _MissingQuestionaryPrompt()

    @staticmethod
    def confirm(*_args: Any, **_kwargs: Any) -> _MissingQuestionaryPrompt:
        return _MissingQuestionaryPrompt()

    @staticmethod
    def text(*_args: Any, **_kwargs: Any) -> _MissingQuestionaryPrompt:
        return _MissingQuestionaryPrompt()


try:
    import questionary as _questionary
except ModuleNotFoundError as exc:
    if exc.name != "questionary":
        raise
    _questionary = _QuestionaryFallback()  # type: ignore[assignment]

questionary: Any = _questionary
Choice: Any = questionary.Choice
Separator: Any = questionary.Separator

__all__ = ["Choice", "Separator", "questionary"]
