"""Bounded conformer selection: which structures reach the MOPAC stage.

Selection is what keeps a calculation affordable. The search may return
dozens of conformers; a strategy decides which of them are optimized,
and every strategy reports the same :class:`SelectionResult` so the
reasoning behind a number stays inspectable.

One rule is enforced here rather than left to convention: the PAS
completeness classifier may be *recorded* as advisory metadata, but it
can never appear as scientific evidence for a selection.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol

from semi_imperium.conformers.ensemble import Conformer, ConformerEnsemble
from semi_imperium.domain.configuration import ConformerSelectionSettings
from semi_imperium.domain.enums import ConformerSelectionStrategy

#: Key under which a PAS completeness label is kept as advisory metadata.
PAS_COMPLETENESS_LABEL_KEY = "pas_completeness_class"

#: Terms that must never appear in :attr:`SelectionResult.evidence`.
FORBIDDEN_EVIDENCE_TERMS = ("pas_completeness", "pas completeness")


@dataclass(frozen=True)
class SelectionResult:
    """The conformers carried forward and the reasoning behind them."""

    strategy: ConformerSelectionStrategy
    selected: tuple[Conformer, ...]
    considered: int
    ranking_basis: str
    evidence: tuple[str, ...] = ()
    advisory_labels: dict[str, str] = field(default_factory=dict)
    """Observations recorded for later study, never used to justify a pick."""

    def __post_init__(self) -> None:
        if not self.selected:
            raise ValueError("SelectionResult.selected must not be empty")
        if not self.ranking_basis.strip():
            raise ValueError("SelectionResult.ranking_basis must not be empty")
        if self.considered < len(self.selected):
            raise ValueError(
                "SelectionResult.considered must cover every selected conformer: "
                f"{self.considered} considered, {len(self.selected)} selected"
            )
        indices = [conformer.index for conformer in self.selected]
        if len(set(indices)) != len(indices):
            raise ValueError(
                f"SelectionResult selected duplicate conformers: {sorted(indices)}"
            )
        _reject_forbidden_evidence(self.evidence)

    @property
    def is_experimental(self) -> bool:
        """Whether the strategy that produced this result is experimental."""
        return self.strategy.is_experimental

    @property
    def selected_indices(self) -> tuple[int, ...]:
        """Ensemble indices of the selected conformers, in selection order."""
        return tuple(conformer.index for conformer in self.selected)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible primitives."""
        return {
            "strategy": self.strategy.value,
            "experimental": self.is_experimental,
            "considered": self.considered,
            "ranking_basis": self.ranking_basis,
            "selected": [conformer.label for conformer in self.selected],
            "selected_indices": list(self.selected_indices),
            "evidence": list(self.evidence),
            "advisory_labels": dict(self.advisory_labels),
        }


class ConformerSelector(Protocol):
    """Decides which conformers of an ensemble reach the MOPAC stage."""

    def select(
        self,
        ensemble: ConformerEnsemble,
        settings: ConformerSelectionSettings,
    ) -> SelectionResult:
        """Return the selection, with the reasoning that produced it."""
        ...


@dataclass(frozen=True)
class EnergyTopNSelector:
    """Default strategy: keep the N lowest-energy conformers of the search.

    The ranking uses the energies the *search* reported, so it happens
    entirely at the CREST level: no MOPAC number is needed to decide
    which structures are worth optimizing, and the decision can be
    replayed from the stored ensemble alone.
    """

    def select(
        self,
        ensemble: ConformerEnsemble,
        settings: ConformerSelectionSettings,
    ) -> SelectionResult:
        """Return the ``settings.top_n`` lowest-energy conformers.

        Raises:
            ValueError: If ``settings`` asks for another strategy, or if
                a multi-conformer ensemble carries no energies to rank.
        """
        strategy = ConformerSelectionStrategy.CREST_ENERGY_TOP_N
        require_strategy(settings, strategy)

        if not ensemble.has_energies:
            if ensemble.size > 1:
                raise ValueError(
                    "Energy Top-N needs an energy for every conformer, but the "
                    f"{ensemble.source.value!r} ensemble of {ensemble.size} "
                    "conformers reported none"
                )
            return SelectionResult(
                strategy=strategy,
                selected=ensemble.conformers,
                considered=ensemble.size,
                ranking_basis="single_structure_without_search_energies",
                evidence=("single_conformer_ensemble",),
            )

        ranked = ensemble.ranked_by_energy()
        window = settings.energy_window_kcal_mol
        evidence = ["search_energy_ranking", f"top_n={settings.top_n}"]
        if window is not None:
            lowest = ranked[0].require_energy()
            ranked = tuple(
                item for item in ranked if item.require_energy() - lowest <= window
            )
            evidence.append(f"energy_window_kcal_mol={window}")

        return SelectionResult(
            strategy=strategy,
            selected=ranked[: settings.top_n],
            considered=ensemble.size,
            ranking_basis="crest_search_energy",
            evidence=tuple(evidence),
        )


def require_strategy(
    settings: ConformerSelectionSettings,
    expected: ConformerSelectionStrategy,
) -> None:
    """Refuse to apply a selector the configuration did not ask for."""
    if settings.resolved_strategy is not expected:
        raise ValueError(
            f"Selector for {expected.value!r} was given settings configured for "
            f"{settings.strategy!r}"
        )


def _reject_forbidden_evidence(evidence: tuple[str, ...]) -> None:
    """Keep the PAS completeness classifier out of scientific evidence."""
    for entry in evidence:
        lowered = entry.lower()
        if any(term in lowered for term in FORBIDDEN_EVIDENCE_TERMS):
            raise ValueError(
                "The PAS completeness classifier is advisory metadata and must "
                f"not be recorded as selection evidence, got {entry!r}"
            )


__all__ = [
    "FORBIDDEN_EVIDENCE_TERMS",
    "PAS_COMPLETENESS_LABEL_KEY",
    "ConformerSelector",
    "EnergyTopNSelector",
    "SelectionResult",
    "require_strategy",
]
