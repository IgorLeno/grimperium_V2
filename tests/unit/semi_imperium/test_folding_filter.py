"""Extreme-folding filter: settings, descriptor, invariant and workflow wiring.

Geometries here are schematic: only the distances the descriptor reads
matter, so bond lengths are not chemically realistic.
"""

from __future__ import annotations

from collections.abc import Sequence

import pytest

from semi_imperium.conformers import (
    BYPASSED_ALL_FOLDED,
    Conformer,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerRequest,
    ConformerSearchProvenance,
    ConformerWorkflow,
    MoleculeTopology,
    apply_folding_filter,
)
from semi_imperium.domain import (
    ConformerSearchSettings,
    ConformerSelectionSettings,
    ConformerSource,
    EffectiveConfiguration,
    FoldingFilterSettings,
)

Point = tuple[float, float, float]

ENABLED = FoldingFilterSettings(enabled=True)

#: Seven-carbon chain: atoms 0 and 6 are six bonds apart (k = 6).
CHAIN7 = MoleculeTopology(atom_count=7, bonds=tuple((i, i + 1, 1) for i in range(6)))


def extended_chain(count: int) -> tuple[Point, ...]:
    return tuple((1.5 * i, 0.0, 0.0) for i in range(count))


def folded_chain(count: int, gap: float = 3.0) -> tuple[Point, ...]:
    """Chain whose last atom sits ``gap`` Å from the first."""
    return (*extended_chain(count - 1)[:-1], (4.5, 2.0, 0.0), (0.0, gap, 0.0))


def conformer(
    index: int,
    coordinates: Sequence[Point],
    elements: Sequence[str] | None = None,
    energy: float | None = None,
) -> Conformer:
    symbols = tuple(elements) if elements is not None else ("C",) * len(coordinates)
    return Conformer(
        index=index,
        geometry=ConformerGeometry(elements=symbols, coordinates=tuple(coordinates)),
        energy_kcal_mol=energy if energy is not None else float(index),
    )


def provenance() -> ConformerSearchProvenance:
    return ConformerSearchProvenance(
        source=ConformerSource.CREST,
        program="crest",
        program_version="3.0.1",
        settings=ConformerSearchSettings(),
    )


def ensemble(*conformers: Conformer) -> ConformerEnsemble:
    return ConformerEnsemble(conformers=conformers, provenance=provenance())


def configuration(selection: ConformerSelectionSettings) -> EffectiveConfiguration:
    return EffectiveConfiguration(
        method_id="crest_pm7",
        method_version="1.0",
        property_id="standard_enthalpy_of_formation",
        conformer_selection=selection,
    )


# ---------------------------------------------------------------------------
# Settings and signature
# ---------------------------------------------------------------------------


def test_defaults_are_the_calibrated_values_and_disabled() -> None:
    settings = FoldingFilterSettings()

    assert settings.enabled is False
    assert settings.min_topological_distance == 6
    assert settings.vdw_scale == pytest.approx(1.0)
    assert settings.min_contacts == 1
    assert settings.max_rg_ratio is None
    assert settings.hbond_max_h_acceptor_angstrom == pytest.approx(2.5)
    assert ConformerSelectionSettings().folding_filter == settings


@pytest.mark.parametrize(
    "kwargs",
    [
        {"min_topological_distance": 1},
        {"vdw_scale": 0.0},
        {"min_contacts": 0},
        {"max_rg_ratio": 0.0},
        {"max_rg_ratio": 1.5},
        {"hbond_max_h_acceptor_angstrom": 0.0},
    ],
)
def test_settings_reject_meaningless_values(kwargs: dict[str, float]) -> None:
    with pytest.raises(ValueError, match="FoldingFilterSettings"):
        FoldingFilterSettings(**kwargs)  # type: ignore[arg-type]


def test_disabled_filter_keeps_existing_signatures() -> None:
    plain = configuration(ConformerSelectionSettings())
    disabled_tuned = configuration(
        ConformerSelectionSettings(
            folding_filter=FoldingFilterSettings(enabled=False, vdw_scale=0.9)
        )
    )

    assert "folding_filter" not in plain.to_dict()["conformer_selection"]
    assert plain.signature() == disabled_tuned.signature()


def test_enabled_filter_and_its_values_enter_the_signature() -> None:
    plain = configuration(ConformerSelectionSettings())
    enabled = configuration(ConformerSelectionSettings(folding_filter=ENABLED))
    tuned = configuration(
        ConformerSelectionSettings(
            folding_filter=FoldingFilterSettings(enabled=True, vdw_scale=0.95)
        )
    )

    assert len({plain.signature(), enabled.signature(), tuned.signature()}) == 3


def test_selection_settings_round_trip_with_and_without_the_filter() -> None:
    enabled = ConformerSelectionSettings(
        folding_filter=FoldingFilterSettings(enabled=True, max_rg_ratio=0.8)
    )
    legacy_payload = {
        "strategy": "crest_energy_top_n",
        "top_n": 10,
        "energy_window_kcal_mol": None,
    }

    assert ConformerSelectionSettings.from_dict(enabled.to_dict()) == enabled
    assert ConformerSelectionSettings.from_dict(legacy_payload) == (
        ConformerSelectionSettings()
    )


# ---------------------------------------------------------------------------
# Descriptor
# ---------------------------------------------------------------------------


def test_folded_chain_end_contact_is_discarded() -> None:
    outcome = apply_folding_filter(
        ensemble(conformer(0, extended_chain(7)), conformer(1, folded_chain(7))),
        CHAIN7,
        ENABLED,
    )

    contacts = {item.index: item.fold_contacts for item in outcome.measurements}
    assert contacts == {0: 0, 1: 1}
    assert outcome.discarded_indices == (1,)
    assert [c.index for c in outcome.kept.conformers] == [0]
    assert outcome.bypassed is False
    assert "folding_filter_folded=1/2" in outcome.evidence


def test_contact_must_be_closer_than_the_scaled_vdw_sum() -> None:
    # C···C vdW sum is 3.40 Å; 3.45 Å is proximity, not contact, at f = 1.
    near_miss = ensemble(
        conformer(0, extended_chain(7)), conformer(1, folded_chain(7, gap=3.45))
    )

    assert apply_folding_filter(near_miss, CHAIN7, ENABLED).discarded_indices == ()

    loose = FoldingFilterSettings(enabled=True, vdw_scale=1.05)
    assert apply_folding_filter(near_miss, CHAIN7, loose).discarded_indices == (1,)


def test_pairs_closer_than_k_bonds_are_not_folding() -> None:
    chain6 = MoleculeTopology(
        atom_count=6, bonds=tuple((i, i + 1, 1) for i in range(5))
    )
    outcome = apply_folding_filter(
        ensemble(conformer(0, extended_chain(6)), conformer(1, folded_chain(6))),
        chain6,
        ENABLED,
    )

    assert outcome.discarded_indices == ()


def test_classic_hydrogen_bond_is_not_a_folding_contact() -> None:
    # HO-C-C-C-C-C-O: donor O0 carries H7; O0···O6 at 2.80 Å (< 3.10 vdW sum).
    elements = ("O", "C", "C", "C", "C", "C", "O", "H")
    topology = MoleculeTopology(
        atom_count=8,
        bonds=(*((i, i + 1, 1) for i in range(6)), (0, 7, 1)),
    )
    folded = folded_chain(7, gap=2.8)
    h_toward_acceptor = (*folded, (0.0, 0.97, 0.0))  # H···O6 = 1.83 Å
    h_away = (*folded, (-0.97, 0.0, 0.0))  # H···O6 ≈ 2.96 Å

    bonded = apply_folding_filter(
        ensemble(conformer(0, h_toward_acceptor, elements)), topology, ENABLED
    )
    not_bonded = apply_folding_filter(
        ensemble(conformer(0, h_away, elements)), topology, ENABLED
    )

    assert bonded.measurements[0].fold_contacts == 0
    assert not_bonded.measurements[0].fold_contacts == 1


def test_fragments_close_together_are_not_folding() -> None:
    # Two three-atom fragments whose ends touch: no bond path between them.
    topology = MoleculeTopology(
        atom_count=6, bonds=((0, 1, 1), (1, 2, 1), (3, 4, 1), (4, 5, 1))
    )
    coordinates = (
        (0.0, 0.0, 0.0),
        (1.5, 0.0, 0.0),
        (3.0, 0.0, 0.0),
        (3.0, 3.0, 0.0),
        (1.5, 3.0, 0.0),
        (0.0, 3.0, 0.0),
    )
    outcome = apply_folding_filter(
        ensemble(conformer(0, coordinates)), topology, ENABLED
    )

    assert outcome.measurements[0].fold_contacts == 0


def test_rg_condition_applies_only_when_configured() -> None:
    pair = ensemble(conformer(0, extended_chain(7)), conformer(1, folded_chain(7)))
    folded = next(
        item
        for item in apply_folding_filter(pair, CHAIN7, ENABLED).measurements
        if item.index == 1
    )
    assert 0 < folded.rg_ratio < 1

    strict = FoldingFilterSettings(enabled=True, max_rg_ratio=folded.rg_ratio / 2)
    lenient = FoldingFilterSettings(enabled=True, max_rg_ratio=1.0)

    assert apply_folding_filter(pair, CHAIN7, strict).discarded_indices == ()
    assert apply_folding_filter(pair, CHAIN7, lenient).discarded_indices == (1,)


# ---------------------------------------------------------------------------
# Invariant: the ensemble never becomes empty
# ---------------------------------------------------------------------------


def test_all_folded_ensemble_is_kept_and_the_bypass_is_recorded() -> None:
    original = ensemble(
        conformer(0, folded_chain(7)), conformer(1, folded_chain(7, gap=2.9))
    )
    outcome = apply_folding_filter(original, CHAIN7, ENABLED)

    assert outcome.bypassed is True
    assert outcome.kept == original
    assert outcome.discarded_indices == ()
    assert BYPASSED_ALL_FOLDED in outcome.evidence
    assert "folding_filter_folded=2/2" in outcome.evidence
    assert outcome.to_dict()["bypassed"] is True


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------


def test_filter_refuses_a_topology_of_another_molecule() -> None:
    with pytest.raises(ValueError, match="7 atoms but the ensemble has 6"):
        apply_folding_filter(ensemble(conformer(0, extended_chain(6))), CHAIN7, ENABLED)


def test_filter_refuses_an_element_without_a_radius() -> None:
    elements = ("C", "C", "C", "C", "C", "C", "Xe")
    with pytest.raises(ValueError, match="no van der Waals radius for 'Xe'"):
        apply_folding_filter(
            ensemble(conformer(0, extended_chain(7), elements)), CHAIN7, ENABLED
        )


def test_filter_refuses_to_run_disabled() -> None:
    with pytest.raises(ValueError, match="disabled"):
        apply_folding_filter(
            ensemble(conformer(0, extended_chain(7))),
            CHAIN7,
            FoldingFilterSettings(),
        )


# ---------------------------------------------------------------------------
# Workflow wiring
# ---------------------------------------------------------------------------


class FixedSearch:
    """Search double that returns a prepared ensemble."""

    def __init__(self, result: ConformerEnsemble) -> None:
        self.result = result

    def search(
        self, request: ConformerRequest, settings: ConformerSearchSettings
    ) -> ConformerEnsemble:
        return self.result


class NoInitialStructure:
    def build(
        self, request: ConformerRequest, settings: ConformerSearchSettings
    ) -> ConformerEnsemble:
        raise AssertionError("the initial-3D route ran while CREST was enabled")


def chain_workflow() -> ConformerWorkflow:
    # The folded conformer is the lowest in energy, so Top-N would pick it.
    return ConformerWorkflow(
        search_backend=FixedSearch(
            ensemble(
                conformer(0, folded_chain(7), energy=0.0),
                conformer(1, extended_chain(7), energy=1.0),
                conformer(2, extended_chain(7), energy=2.0),
            )
        ),
        initial_structure_backend=NoInitialStructure(),
    )


REQUEST = ConformerRequest(molecule_id="heptane", smiles="CCCCCCC")


def test_workflow_filters_before_the_selection_strategy() -> None:
    prepared = chain_workflow().prepare(
        REQUEST,
        search_settings=ConformerSearchSettings(),
        selection_settings=ConformerSelectionSettings(top_n=1, folding_filter=ENABLED),
        topology=CHAIN7,
    )
    payload = prepared.to_dict()

    assert prepared.selection.selected_indices == (1,)
    assert prepared.selection.considered == 2
    assert prepared.ensemble.size == 3
    assert payload["ensemble_size"] == 3
    assert payload["folding_filter"]["discarded_indices"] == [0]


def test_workflow_without_the_filter_is_unchanged() -> None:
    prepared = chain_workflow().prepare(
        REQUEST,
        search_settings=ConformerSearchSettings(),
        selection_settings=ConformerSelectionSettings(top_n=1),
    )

    assert prepared.folding is None
    assert prepared.selection.selected_indices == (0,)
    assert "folding_filter" not in prepared.to_dict()


def test_workflow_requires_a_topology_when_the_filter_is_enabled() -> None:
    with pytest.raises(ValueError, match="folding filter needs the molecule topology"):
        chain_workflow().prepare(
            REQUEST,
            search_settings=ConformerSearchSettings(),
            selection_settings=ConformerSelectionSettings(folding_filter=ENABLED),
        )
