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

from collections.abc import Callable
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
from semi_imperium.conformers.folding import (
    FoldingFilterOutcome,
    apply_folding_filter,
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

#: Derives the topology of a request in the atom order of its ensemble;
#: raises :class:`ConformerBackendError` when it cannot vouch for that order.
TopologyProvider = Callable[[ConformerRequest, ConformerEnsemble], MoleculeTopology]


@dataclass(frozen=True)
class ConformerPreparation:
    """What the conformer stage hands over to the MOPAC stage."""

    ensemble: ConformerEnsemble
    selection: SelectionResult
    folding: FoldingFilterOutcome | None = None
    """Set only when the folding filter ran before the selection."""

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
        payload: dict[str, Any] = {
            "provenance": self.provenance.to_dict(),
            "ensemble_size": self.ensemble.size,
            "selection": self.selection.to_dict(),
        }
        if self.folding is not None:
            payload["folding_filter"] = self.folding.to_dict()
        return payload


class ConformerWorkflow:
    """Produces conformers for one molecule and narrows them to a subset."""

    def __init__(
        self,
        *,
        search_backend: ConformerSearchBackend,
        initial_structure_backend: InitialStructureBackend,
        confpass_backend: ConfPassBackend | None = None,
        topology_provider: TopologyProvider | None = None,
    ) -> None:
        self._search_backend = search_backend
        self._initial_structure_backend = initial_structure_backend
        self._confpass_backend = confpass_backend
        self._topology_provider = topology_provider

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
                Required by CONFPASS, which needs SDF input, and by the
                folding filter when it is enabled. When omitted, the
                workflow's topology provider derives it from the built
                ensemble, and only if one of those two needs it.

        Raises:
            ConformerBackendError: If the chosen route or the topology
                provider fails.
            ValueError: If the configured strategy or the folding filter is
                missing something it needs, such as a CONFPASS backend or a
                topology.
        """
        ensemble = self.build_ensemble(request, search_settings)
        if topology is None and self._needs_topology(selection_settings):
            topology = self._derive_topology(request, ensemble)
        folding_settings = selection_settings.folding_filter
        folding: FoldingFilterOutcome | None = None
        if folding_settings.enabled:
            if topology is None:
                raise ValueError(
                    "The folding filter needs the molecule topology to count "
                    "contacts; pass one matching the atom order"
                )
            folding = apply_folding_filter(ensemble, topology, folding_settings)
        selector = self._selector_for(selection_settings, request, topology)
        selection = selector.select(
            ensemble if folding is None else folding.kept, selection_settings
        )
        return ConformerPreparation(
            ensemble=ensemble, selection=selection, folding=folding
        )

    def build_ensemble(
        self,
        request: ConformerRequest,
        search_settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        """Return the ensemble for ``request`` from the configured route."""
        if search_settings.enabled:
            return self._search_backend.search(request, search_settings)
        return self._initial_structure_backend.build(request, search_settings)

    @staticmethod
    def _needs_topology(settings: ConformerSelectionSettings) -> bool:
        """Whether the folding filter or the strategy reads connectivity."""
        return (
            settings.folding_filter.enabled
            or settings.resolved_strategy
            is ConformerSelectionStrategy.CONFPASS_PRIORITIZATION
        )

    def _derive_topology(
        self,
        request: ConformerRequest,
        ensemble: ConformerEnsemble,
    ) -> MoleculeTopology | None:
        """Ask the provider, if any, for the topology of ``ensemble``."""
        if self._topology_provider is None:
            return None
        return self._topology_provider(request, ensemble)

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


__all__ = ["ConformerPreparation", "ConformerWorkflow", "TopologyProvider"]
