"""Protocol 2: the canonical Semi-Imperium protocol behind the training data.

Protocol 1 is the legacy ``crest_pm7`` pipeline that produced the first 8,232
rows of ``thermo_pm7.csv``. Protocol 2 selects conformers with CONFPASS plus
the folding filter (specification 2026-10-03) and demands a frequency check
that confirms a minimum. Both choices are part of the calculation signature,
so a number computed under another configuration is never mistaken for it.
"""

from __future__ import annotations

from dataclasses import replace

from semi_imperium.domain import (
    ConformerSelectionStrategy,
    EffectiveConfiguration,
    VerificationPolicy,
)
from semi_imperium.settings import SemiImperiumSettings

LEGACY_PROTOCOL = "1"
CANONICAL_PROTOCOL = "2"
HAMILTONIAN = "PM7"


def canonical_settings(
    base: SemiImperiumSettings | None = None,
) -> SemiImperiumSettings:
    """Return ``base`` (or the defaults) adjusted to protocol 2."""
    settings = base or SemiImperiumSettings()
    selection = settings.conformer_selection
    return replace(
        settings,
        conformer_search=replace(settings.conformer_search, enabled=True),
        conformer_selection=replace(
            selection,
            strategy=ConformerSelectionStrategy.CONFPASS_PRIORITIZATION.value,
            folding_filter=replace(selection.folding_filter, enabled=True),
        ),
        verification=replace(
            settings.verification, policy=VerificationPolicy.REQUIRE_MINIMUM
        ),
    )


def canonical_configuration(settings: SemiImperiumSettings) -> EffectiveConfiguration:
    """Resolve the one configuration every protocol-2 calculation runs under."""
    return settings.configuration_for(HAMILTONIAN, crest_enabled=True)


def canonical_signature(settings: SemiImperiumSettings) -> str:
    """Digest that identifies protocol 2 under ``settings``."""
    return canonical_configuration(settings).signature().digest
