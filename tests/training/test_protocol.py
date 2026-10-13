"""Protocol 2 is CONFPASS + folding filter + PM7 with minimum required."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from semi_imperium.domain import ConformerSelectionStrategy, VerificationPolicy
from semi_imperium.settings import RuntimeSettings, SemiImperiumSettings
from semi_imperium.training.protocol import (
    CANONICAL_PROTOCOL,
    HAMILTONIAN,
    LEGACY_PROTOCOL,
    canonical_configuration,
    canonical_settings,
    canonical_signature,
)


def test_protocol_labels_are_strings() -> None:
    assert (LEGACY_PROTOCOL, CANONICAL_PROTOCOL, HAMILTONIAN) == ("1", "2", "PM7")


def test_canonical_settings_select_with_confpass_and_folding_filter() -> None:
    settings = canonical_settings()
    selection = settings.conformer_selection
    assert (
        selection.resolved_strategy
        is ConformerSelectionStrategy.CONFPASS_PRIORITIZATION
    )
    assert selection.folding_filter.enabled is True
    assert settings.conformer_search.enabled is True
    assert settings.verification.policy is VerificationPolicy.REQUIRE_MINIMUM


def test_canonical_configuration_is_pm7_with_crest() -> None:
    configuration = canonical_configuration(canonical_settings())
    assert configuration.semiempirical.hamiltonian == "PM7"
    assert configuration.conformer_search.enabled is True


def test_signature_ignores_runtime_knobs() -> None:
    base = canonical_settings()
    other_runtime = replace(
        base,
        runtime=RuntimeSettings(crest_threads=9, work_dir=Path("/tmp/elsewhere")),
    )
    assert canonical_signature(base) == canonical_signature(other_runtime)


def test_signature_changes_with_the_selection_strategy() -> None:
    energy_top_n = replace(
        canonical_settings(),
        conformer_selection=replace(
            canonical_settings().conformer_selection,
            strategy=ConformerSelectionStrategy.CREST_ENERGY_TOP_N.value,
        ),
    )
    assert canonical_signature(energy_top_n) != canonical_signature(
        canonical_settings()
    )


def test_base_settings_are_not_mutated() -> None:
    base = SemiImperiumSettings()
    canonical_settings(base)
    assert (
        base.conformer_selection.resolved_strategy
        is ConformerSelectionStrategy.CREST_ENERGY_TOP_N
    )
    assert base.verification.policy is VerificationPolicy.NONE
