"""Orchestration of the conformer stage: route first, then selection.

Two decisions live here, and nothing else does:

* whether the structures come from a CREST search or, when the search
  is disabled, from the initial-3D route — MOPAC always gets a geometry;
* which selection strategy narrows the resulting ensemble.

Both external programs are reached through the protocols in
:mod:`semi_imperium.conformers.backends`, so this orchestration is
exercised with in-memory doubles and never spawns anything.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from semi_imperium.conformers.backends import (
    ConformerBackendError,
    ConformerRequest,
    ConformerSearchBackend,
    ConfPassBackend,
    InitialStructureBackend,
)
from semi_imperium.conformers.confpass import ConfPassSelector, MoleculeTopology
from semi_imperium.conformers.ensemble import (
    Conformer,
    ConformerEnsemble,
    ConformerSearchProvenance,
)
from semi_imperium.conformers.selection import (
    ConformerSelector,
    EnergyTopNSelector,
    SelectionResult,
)
from semi_imperium.domain.configuration import (
    ConformerSearchSettings,
    ConformerSelectionSettings,
)
from semi_imperium.domain.enums import ConformerSelectionStrategy


@dataclass(frozen=True)
class ConformerPreparation:
    """What the conformer stage hands over to the MOPAC stage."""

    ensemble: ConformerEnsemble
    selection: SelectionResult

    @property
    def selected(self) -> tuple[Conformer, ...]:
        """The conformers that will actually be optimized."""
        return self.selection.selected

    @property
    def provenance(self) -> ConformerSearchProvenance:
        """How the structures were produced."""
        return self.ensemble.provenance

    @property
    def is_experimental(self) -> bool:
        """Whether an experimental strategy produced this selection."""
        return self.selection.is_experimental

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "provenance": self.provenance.to_dict(),
            "ensemble_size": self.ensemble.size,
            "selection": self.selection.to_dict(),
        }


class ConformerWorkflow:
    """Produces conformers for one molecule and narrows them to a subset."""

    def __init__(
        self,
        *,
        search_backend: ConformerSearchBackend,
        initial_structure_backend: InitialStructureBackend,
        confpass_backend: ConfPassBackend | None = None,
    ) -> None:
        self._search_backend = search_backend
        self._initial_structure_backend = initial_structure_backend
        self._confpass_backend = confpass_backend

    def prepare(
        self,
        request: ConformerRequest,
        *,
        search_settings: ConformerSearchSettings,
        selection_settings: ConformerSelectionSettings,
        topology: MoleculeTopology | None = None,
    ) -> ConformerPreparation:
        """Build the ensemble and apply the configured selection strategy.

        Args:
            request: The molecule to prepare structures for.
            search_settings: CREST settings; ``enabled=False`` routes the
                molecule through the initial-3D structure instead.
            selection_settings: Which strategy narrows the ensemble.
            topology: Connectivity matching the ensemble's atom order.
                Required only by CONFPASS, which needs SDF input.

        Raises:
            ConformerBackendError: If the chosen route fails.
            ValueError: If the configured strategy is missing something
                it needs, such as a CONFPASS backend or a topology.
        """
        ensemble = self.build_ensemble(request, search_settings)
        selector = self._selector_for(selection_settings, request, topology)
        selection = selector.select(ensemble, selection_settings)
        return ConformerPreparation(ensemble=ensemble, selection=selection)

    def build_ensemble(
        self,
        request: ConformerRequest,
        search_settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        """Return the ensemble for ``request`` from the configured route."""
        if search_settings.enabled:
            return self._search_backend.search(request, search_settings)
        return self._initial_structure_backend.build(request, search_settings)

    def _selector_for(
        self,
        settings: ConformerSelectionSettings,
        request: ConformerRequest,
        topology: MoleculeTopology | None,
    ) -> ConformerSelector:
        """Return the selector the configuration asked for."""
        strategy = settings.resolved_strategy
        if strategy is ConformerSelectionStrategy.CREST_ENERGY_TOP_N:
            return EnergyTopNSelector()
        if self._confpass_backend is None:
            raise ConformerBackendError(
                "CONFPASS prioritization was configured but no CONFPASS "
                "backend was provided to the workflow",
                code="confpass_unavailable",
            )
        if topology is None:
            raise ValueError(
                "CONFPASS prioritization needs the molecule topology to adapt "
                "the XYZ ensemble to SDF; pass one matching the atom order"
            )
        return ConfPassSelector(
            backend=self._confpass_backend,
            topology=topology,
            molecule_id=request.molecule_id,
        )


__all__ = ["ConformerPreparation", "ConformerWorkflow"]
