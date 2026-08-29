"""Shared state behind Calculate, Database and Settings.

The three areas are views onto one workspace: the same defaults, the
same store and the same Calculate table. Holding them together here is
what lets Settings change a default and Calculate pick it up on the next
prepared run — while everything already persisted keeps the effective
configuration it was computed under.
"""

from __future__ import annotations

from dataclasses import replace

from semi_imperium.molecules.service import MoleculeResolutionService
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.prompts import Prompter
from semi_imperium.settings import RuntimeSettings, SemiImperiumSettings
from semi_imperium.workflows.calculation import CalculateSession, CalculationExecutor


class SemiImperiumWorkspace:
    """One session's defaults, store, molecule table and execution boundary."""

    def __init__(
        self,
        *,
        settings: SemiImperiumSettings | None = None,
        store: SemiImperiumStore | None = None,
        session: CalculateSession | None = None,
        executor: CalculationExecutor | None = None,
        prompter: Prompter | None = None,
        resolution_service: MoleculeResolutionService | None = None,
    ) -> None:
        self._settings = settings or SemiImperiumSettings()
        self._store = store
        self._store_pinned = store is not None
        self._executor = executor
        self._executor_pinned = executor is not None
        self.prompter = prompter
        self.session = session or CalculateSession(
            resolution_service=resolution_service
        )

    @property
    def settings(self) -> SemiImperiumSettings:
        """The defaults the next prepared run will be resolved from."""
        return self._settings

    @settings.setter
    def settings(self, value: SemiImperiumSettings) -> None:
        self._settings = value
        if not self._executor_pinned:
            # A default executor is built from the settings, so it must not
            # outlive them; the next execution rebuilds it.
            self._executor = None

    def set_runtime(self, runtime: RuntimeSettings) -> RuntimeSettings:
        """Adopt new runtime settings and return them."""
        self.settings = replace(self._settings, runtime=runtime)
        return runtime

    @property
    def store(self) -> SemiImperiumStore:
        """The store rooted at the configured location.

        An explicitly injected store is honoured as-is; otherwise the
        store follows ``runtime.store_root`` so changing it in Settings
        actually points the next run somewhere else.
        """
        if self._store_pinned and self._store is not None:
            return self._store
        root = self._settings.runtime.store_root
        if self._store is None or self._store.root != root:
            self._store = SemiImperiumStore(root)
        return self._store

    @property
    def executor(self) -> CalculationExecutor:
        """The execution boundary, built from the current settings on demand."""
        if self._executor is None:
            from semi_imperium.workflows.execution import (
                ScientificCalculationExecutor,
            )

            self._executor = ScientificCalculationExecutor(self._settings)
        return self._executor

    @executor.setter
    def executor(self, value: CalculationExecutor) -> None:
        self._executor = value
        self._executor_pinned = True


__all__ = ["SemiImperiumWorkspace"]
