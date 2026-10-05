"""The three selection arms and the pieces that run them on one ensemble.

CREST runs once per molecule; every arm narrows that same ensemble, so
the search backend handed to :class:`ConformerWorkflow` here only ever
returns the ensemble it was built with.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from semi_imperium.conformers import (
    ConformerEnsemble,
    ConformerPreparation,
    ConformerRequest,
    ConformerWorkflow,
    CrestRun,
)
from semi_imperium.conformers.confpass_port import PortedConfPass
from semi_imperium.conformers.crest_runner import ENSEMBLE_FILENAME
from semi_imperium.conformers.initial_structure import RDKitInitialStructure
from semi_imperium.conformers.topology import SmilesTopology
from semi_imperium.domain import (
    ConformerSearchSettings,
    ConformerSelectionSettings,
    EffectiveConfiguration,
    FoldingFilterSettings,
    VerificationPolicy,
    VerificationSettings,
)

#: Also the order arms run in and the order cost ties are broken in.
ARM_NAMES = ("A", "B", "C")

ARM_SELECTIONS: dict[str, ConformerSelectionSettings] = {
    "A": ConformerSelectionSettings(),
    "B": ConformerSelectionSettings(strategy="confpass_prioritization"),
    "C": ConformerSelectionSettings(
        strategy="confpass_prioritization",
        folding_filter=FoldingFilterSettings(enabled=True, max_rg_ratio=None),
    ),
}

HAMILTONIAN = "PM7"


def configuration(arm: str) -> EffectiveConfiguration:
    """Effective configuration of ``arm``: default CREST, PM7, minimum required."""
    return EffectiveConfiguration(
        method_id="crest_pm7",
        method_version="1.0",
        property_id="standard_enthalpy_of_formation",
        conformer_search=ConformerSearchSettings(),
        conformer_selection=ARM_SELECTIONS[arm],
        verification=VerificationSettings(policy=VerificationPolicy.REQUIRE_MINIMUM),
    )


def parse_arms(text: str) -> tuple[str, ...]:
    """Parse ``A,B,C`` into a de-duplicated, ordered tuple of arm names."""
    names = [part.strip().upper() for part in text.split(",") if part.strip()]
    unknown = sorted(set(names) - set(ARM_NAMES))
    if unknown or not names:
        raise ValueError(f"Unknown arms {unknown}; choose from {ARM_NAMES}")
    return tuple(arm for arm in ARM_NAMES if arm in names)


@dataclass(frozen=True)
class FixedEnsembleSearch:
    """Search backend that hands every arm the one ensemble CREST produced."""

    ensemble: ConformerEnsemble

    def search(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        del request, settings
        return self.ensemble


@dataclass(frozen=True)
class StoredEnsembleRunner:
    """CREST runner for smoke runs: reads ``<root>/<mol_id>/crest_conformers.xyz``.

    The stored files carry no settings record, so the command says plainly
    that nothing was run and where the ensemble came from.
    """

    root: Path

    def run(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> CrestRun:
        del settings
        path = self.root / request.molecule_id / ENSEMBLE_FILENAME
        if not path.exists():
            raise FileNotFoundError(f"No stored ensemble at {path}")
        return CrestRun(
            ensemble_xyz=path.read_text(encoding="utf-8"),
            command=("stored", str(path)),
        )


def conformer_workflow(ensemble: ConformerEnsemble) -> ConformerWorkflow:
    """Production conformer workflow, fed the molecule's fixed ensemble."""
    return ConformerWorkflow(
        search_backend=FixedEnsembleSearch(ensemble),
        initial_structure_backend=RDKitInitialStructure(),
        confpass_backend=PortedConfPass(),
        topology_provider=SmilesTopology(),
    )


def prepare(
    ensemble: ConformerEnsemble,
    request: ConformerRequest,
    arm: str,
) -> ConformerPreparation:
    """Run ``arm``'s selection (and filter, for C) on ``ensemble``."""
    settings = configuration(arm)
    return conformer_workflow(ensemble).prepare(
        request,
        search_settings=settings.conformer_search,
        selection_settings=settings.conformer_selection,
    )


__all__ = [
    "ARM_NAMES",
    "ARM_SELECTIONS",
    "HAMILTONIAN",
    "FixedEnsembleSearch",
    "StoredEnsembleRunner",
    "configuration",
    "conformer_workflow",
    "parse_arms",
    "prepare",
]
